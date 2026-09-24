# Method wrappers for the on-demand simulator.
#
# Ported from the DiNA pipeline (inspo_code/03_estimate.R, fit_method())
# with the same hyperparameters the grid was built with, so simulated cells
# are comparable to gridded ones. Each wrapper returns only what scoring
# needs:
#   score    : p x p differential score (higher = more likely differential)
#   pip_diff : p x p posterior P(edge differs) -- multiGGM only
#
# Method packages are Suggests (several are GitHub-only); they are checked
# once, up front, by .check_method_pkgs().

.METHOD_PKGS <- list(
  JGL_BIC     = "JGL",
  JGL_BIC_GGL = "JGL",
  SSJGL       = "spikeyglass",
  SSJGL_GGL   = "spikeyglass",
  multiGGM    = "multiGGM",
  JRF         = "iJRFNet",
  InvCov      = character(0)
)

.check_method_pkgs <- function(method) {
  pkgs <- .METHOD_PKGS[[method]]
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Simulating ", method, " requires package(s): ",
         paste(missing, collapse = ", "),
         ". See the DESCRIPTION 'Remotes' field for install sources.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Methods whose inference uses posterior probabilities, not a permutation null.
.BAYES_METHODS <- "multiGGM"

.fit_method <- function(method, dat) {
  switch(method,
    JGL_BIC     = .fit_jgl(dat, penalty = "fused"),
    JGL_BIC_GGL = .fit_jgl(dat, penalty = "group"),
    SSJGL       = .fit_ssjgl(dat, penalty = "fused"),
    SSJGL_GGL   = .fit_ssjgl(dat, penalty = "group"),
    multiGGM    = .fit_multiggm(dat),
    JRF         = .fit_jrf(dat),
    InvCov      = .fit_invcov(dat),
    stop("No simulation wrapper for method: ", method, call. = FALSE)
  )
}

.abs_diff <- function(theta) abs(theta[[1]] - theta[[2]])

# --- JGL, BIC-tuned (two-stage lambda1 then lambda2 grid search) ----------

.vbic <- function(dat, theta, thr = 1e-4) {
  sum(vapply(seq_along(dat), function(k) {
    n <- nrow(dat[[k]])
    S <- stats::cov(dat[[k]])
    Th <- theta[[k]]
    sum(diag(S %*% Th)) -
      determinant(Th, logarithm = TRUE)$modulus[1] +
      log(n) * sum(abs(Th) >= thr) / n
  }, numeric(1)))
}

.fit_jgl <- function(dat, penalty, n_lambda = 10L) {
  S_pool <- (stats::cov(dat[[1]]) + stats::cov(dat[[2]])) / 2
  signal <- stats::median(abs(S_pool[upper.tri(S_pool)]))
  lam_seq <- exp(seq(log(signal / 50), log(signal * 2), length.out = n_lambda))
  jgl <- function(l1, l2) JGL::JGL(Y = dat, penalty = penalty, lambda1 = l1,
                                   lambda2 = l2, return.whole.theta = TRUE)
  bic1 <- vapply(lam_seq, function(l) .vbic(dat, jgl(l, lam_seq[1])$theta),
                 numeric(1))
  l1 <- lam_seq[which.min(bic1)]
  bic2 <- vapply(lam_seq, function(l) .vbic(dat, jgl(l1, l)$theta), numeric(1))
  l2 <- lam_seq[which.min(bic2)]
  list(score = .abs_diff(jgl(l1, l2)$theta))
}

# --- SSJGL (spike-and-slab JGL, 3-step v0 ladder) --------------------------

.fit_ssjgl <- function(dat, penalty) {
  lambda1 <- 1
  v0s <- lambda1 / (lambda1 + (1:3) * 5)
  fit <- spikeyglass::ssjgl(Y = dat, penalty = penalty, lambda0 = 1,
                            lambda1 = lambda1, lambda2 = 1, v1 = 1, v0s = v0s,
                            a = 1, b = 1, doubly = TRUE, normalize = FALSE)
  list(score = .abs_diff(fit$thetalist[[length(v0s)]]))
}

# --- multiGGM (Bayesian SSVS; PIP-masked precision for ranking) ------------

.fit_multiggm <- function(dat, burnin = 5000L, nsave = 1000L, nchains = 2L,
                          pip_mask_thr = 0.05) {
  fit <- multiGGM::multiggm_mcmc(
    data_list = dat, method = "ssvs", burnin = burnin, nsave = nsave,
    nchains = nchains, parallel = FALSE, seed = 123,
    hyper = list(v0 = 0.0004, v1 = 1, lambda = 1, a = 1, b = 4,
                 alpha = 2, beta = 5, w = 0.8, alpha_prop = 1, beta_prop = 1))
  chains <- if (inherits(fit, "multiggm_fit_list")) seq_len(fit$nchains) else 1L
  avg <- function(f) Reduce(`+`, lapply(chains, function(ch) f(fit, chain = ch))) /
    length(chains)
  pip <- avg(multiGGM::pip_edges)
  pip_diff <- avg(multiGGM::pip_diff_edge)
  theta_ch <- lapply(chains, function(ch) multiGGM::posterior_precision(fit, chain = ch))
  theta <- lapply(1:2, function(k)
    Reduce(`+`, lapply(theta_ch, `[[`, k)) / length(chains))
  masked <- list(theta[[1]] * (pip[, , 1] >= pip_mask_thr),
                 theta[[2]] * (pip[, , 2] >= pip_mask_thr))
  list(score = .abs_diff(masked), pip_diff = pip_diff)
}

# --- JRF (joint random forests; |importance1 - importance2|) --------------

.fit_jrf <- function(dat, ntree = 1000L) {
  p <- ncol(dat[[1]])
  genes <- paste0("V", seq_len(p))
  X <- lapply(dat, function(d) {
    z <- scale(d)
    z[, !is.finite(colSums(z))] <- 0     # constant columns -> no signal
    t(z)                                 # iJRF wants p x n
  })
  fit <- iJRFNet::iJRF(X = X, ntree = ntree, mtry = round(sqrt(p - 1)),
                       genes.name = genes)
  idx <- stats::setNames(seq_len(p), genes)
  r <- idx[as.character(fit[, "gene1"])]
  c <- idx[as.character(fit[, "gene2"])]
  d <- abs(as.numeric(fit[, "importance1"]) - as.numeric(fit[, "importance2"]))
  score <- matrix(0, p, p)
  score[cbind(r, c)] <- d
  score[cbind(c, r)] <- d
  list(score = score)
}

# --- InvCov (unregularized inverse / pseudoinverse baseline) --------------

.fit_invcov <- function(dat) {
  theta <- lapply(dat, function(X) {
    S <- stats::cov(X)
    out <- tryCatch(solve(S), error = function(e) MASS::ginv(S))
    (out + t(out)) / 2
  })
  list(score = .abs_diff(theta))
}
