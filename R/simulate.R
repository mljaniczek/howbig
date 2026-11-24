#' Simulate power for a design
#'
#' Runs an on-demand simulation for a single design cell: generate networks
#' and data under the grid's data-generating family, fit the method, and
#' score its FDR-controlled differential calls. Use this directly to answer
#' designs the grid does not cover (other `p`, unequal group sizes, other FDR
#' targets), or via `power_for(..., simulate = TRUE)`.
#'
#' Inference matches how the grid was scored. Frequentist methods are called
#' against a pooled label-permutation null (one representative network,
#' `n_perm` permutations with fresh data, Benjamini-Hochberg); multiGGM uses
#' its posterior probabilities directly (direct-posterior FDR).
#'
#' Simulation is slow: each replicate and each permutation is a full method
#' fit. Start with small `reps` / `n_perm` to gauge runtime.
#'
#' @inheritParams power_for
#' @param n Per-group sample size; a single number, or `c(n1, n2)` for
#'   unequal groups.
#' @param fdr Target false-discovery rate for the differential calls.
#' @param reps Number of simulated replicates (fresh network + data each).
#' @param n_perm Number of label permutations for the null (ignored for
#'   multiGGM).
#' @param cores Parallel worker processes. Workers are fresh R sessions
#'   (PSOCK), not forks, so they are safe with multithreaded BLAS such as
#'   Apple Accelerate; `howbig` must be installed for `cores > 1`.
#' @param seed Base seed; results are reproducible for a given seed
#'   regardless of `cores`.
#' @param max_minutes Refuse to start if the estimated runtime exceeds this
#'   many minutes. The estimate comes from per-fit timings of the production
#'   grid runs; raise this to run long designs deliberately.
#' @param verbose Print progress messages.
#' @return A `howbig_power` object with `provenance = "simulated"`.
#' @examples
#' \dontrun{
#' simulate_power(p = 100, n = 500, method = "JGL_BIC", reps = 10)
#' }
#' @export
simulate_power <- function(p, n, method, perturbation = "structural",
                           fdr = 0.05, reps = 20L, n_perm = 30L,
                           cores = 1L, seed = 1L, max_minutes = 60,
                           verbose = TRUE) {
  stopifnot(length(p) == 1L, p >= 10, length(n) %in% 1:2, all(n >= 10),
            fdr > 0, fdr < 1, reps >= 2, n_perm >= 1)
  p <- as.integer(p)
  n <- .normalize_n(n)
  method <- .normalize_method(method)
  perturbation <- .normalize_perturbation(perturbation)
  .check_method_pkgs(method)
  cores <- max(1L, as.integer(cores))

  key <- .cache_key(p, n, method, perturbation, fdr, reps, n_perm, seed)
  if (!is.null(hit <- .howbig_cache[[key]])) {
    if (verbose) message("howbig: using cached simulation for this design.")
    return(hit)
  }

  bayes <- method %in% .BAYES_METHODS
  n_fits <- reps + if (bayes) 0L else n_perm
  eta <- .estimate_minutes(method, p, n, n_fits, cores)
  if (eta > max_minutes) {
    stop(sprintf(paste0(
      "Estimated runtime ~%s for %d %s fits at p = %d, n = %d/%d on %d core(s) ",
      "exceeds max_minutes = %g. Reduce reps/n_perm, add cores, or raise ",
      "max_minutes to run it anyway."),
      .fmt_minutes(eta), n_fits, method, p, n[1], n[2], cores, max_minutes),
      call. = FALSE)
  }
  if (verbose) message(sprintf("howbig: %d %s fits, estimated ~%s on %d core(s).",
                               n_fits, method, .fmt_minutes(eta), cores))

  t0 <- Sys.time()
  ut_n <- p * (p - 1L) / 2L
  cl <- .make_cluster(cores)
  if (!is.null(cl)) on.exit(parallel::stopCluster(cl), add = TRUE)

  # --- 1. Permutation null (frequentist methods) --------------------------
  null_pool <- NULL
  if (!bayes) {
    if (verbose) message(sprintf("howbig: building permutation null (%d refits)...",
                                 n_perm))
    set.seed(seed)
    null_gold <- .generate_gold(p, perturbation)
    perms <- .papply(seq_len(n_perm), cl, function(k) {
      set.seed(seed + k)
      dat <- .gen_data(null_gold, n)
      pooled <- rbind(dat[[1]], dat[[2]])
      idx <- sample(sum(n))
      dat_perm <- list(pooled[idx[seq_len(n[1])], , drop = FALSE],
                       pooled[idx[-seq_len(n[1])], , drop = FALSE])
      s <- .fit_method(method, dat_perm)$score
      s[upper.tri(s)]
    })
    ok <- vapply(perms, function(v) is.numeric(v) && length(v) == ut_n, logical(1))
    if (!any(ok)) stop("All permutation refits failed; cannot build a null.",
                       call. = FALSE)
    null_pool <- sort(unlist(perms[ok], use.names = FALSE))
    null_pool <- null_pool[is.finite(null_pool)]
  }

  # --- 2. Replicates: fresh network + data, fit, score --------------------
  if (verbose) message(sprintf("howbig: simulating %d replicates...", reps))
  rows <- .papply(seq_len(reps), cl, function(r) {
    set.seed(seed + 100000L + r)
    gold <- .generate_gold(p, perturbation)
    fit <- .fit_method(method, .gen_data(gold, n))
    .score_rep(fit, gold, null_pool, fdr)
  })
  ok <- vapply(rows, is.numeric, logical(1))
  if (!any(ok)) stop("All replicate fits failed.", call. = FALSE)
  M <- do.call(rbind, rows[ok])

  q <- function(x, pr) as.numeric(stats::quantile(x, pr, na.rm = TRUE))
  mn <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  sens <- M[, "sensitivity"]
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  res <- structure(list(
    method = method, p = p, n = n, perturbation = perturbation,
    scenario = "default", fdr_target = fdr, provenance = "simulated",
    sensitivity = mn(sens), sensitivity_lo = q(sens, .10),
    sensitivity_hi = q(sens, .90),
    sensitivity_mcse = stats::sd(sens, na.rm = TRUE) / sqrt(sum(!is.na(sens))),
    achieved_fdr = mn(M[, "achieved_fdr"]), aupr = mn(M[, "aupr"]),
    n_diff_edges = mn(M[, "n_diff"]), n_pairs = p * (p - 1) / 2,
    n_reps = sum(ok), n_perm = if (bayes) NA_integer_ else n_perm,
    runtime_sec = runtime, family = howbig_family,
    notes = if (sum(!ok)) sprintf("%d of %d replicate fits failed and were dropped.",
                                  sum(!ok), reps) else character(0)
  ), class = "howbig_power")
  if (.miscalibrated(res$achieved_fdr, fdr)) {
    res$notes <- c(res$notes, sprintf(
      "achieved FDR %.2f exceeds target %.2f \u2014 this method is not calibrated here.",
      res$achieved_fdr, fdr))
  }
  .howbig_cache[[key]] <- res
  if (verbose) message(sprintf("howbig: done in %.1f min.", runtime / 60))
  res
}

# Score one replicate: AUPR of the differential score, plus sensitivity and
# achieved FDR of the FDR-controlled calls.
.score_rep <- function(fit, gold, null_pool, fdr) {
  ut <- upper.tri(fit$score)
  obs <- fit$score[ut]
  truth <- gold$delta
  called <- if (!is.null(fit$pip_diff)) {
    .bayes_fdr_call(fit$pip_diff[ut], fdr)
  } else {
    # pooled empirical p-value: (1 + #{null >= obs}) / (1 + |null|)
    exceed <- length(null_pool) -
      findInterval(obs - .Machine$double.eps, null_pool)
    pv <- (1 + exceed) / (1 + length(null_pool))
    stats::p.adjust(pv, method = "BH") < fdr
  }
  TP <- sum(truth & called); FP <- sum(!truth & called)
  c(sensitivity  = if (any(truth)) TP / sum(truth) else NA_real_,
    achieved_fdr = if (TP + FP > 0) FP / (TP + FP) else NA_real_,
    aupr         = .aupr(truth, obs),
    n_diff       = sum(truth))
}

# Direct-posterior FDR: call the largest top-k set (by PIP) whose running
# mean of (1 - PIP) stays <= target.
.bayes_fdr_call <- function(pip, fdr) {
  called <- logical(length(pip))
  ord <- order(pip, decreasing = TRUE, na.last = NA)
  if (!length(ord)) return(called)
  running <- cumsum(1 - pip[ord]) / seq_along(ord)
  k <- sum(cummax(running) <= fdr)
  if (k > 0) called[ord[seq_len(k)]] <- TRUE
  called
}

# Trapezoid AUPR (as in the DiNA metrics), truth logical, score numeric.
.aupr <- function(truth, score) {
  ok <- !is.na(score)
  truth <- truth[ok]; score <- score[ok]
  n_pos <- sum(truth)
  if (n_pos == 0 || n_pos == length(truth)) return(NA_real_)
  t_sorted <- truth[order(score, decreasing = TRUE)]
  tp <- cumsum(t_sorted); fp <- cumsum(!t_sorted)
  prec <- c(1, tp / (tp + fp)); rec <- c(0, tp / n_pos)
  sum(diff(rec) * (prec[-length(prec)] + prec[-1]) / 2)
}

# Worker pool for cores > 1. PSOCK (fresh R processes) rather than fork:
# forking after a multithreaded BLAS (e.g. Apple Accelerate) has initialized
# segfaults the children. Each worker's BLAS is pinned to one thread so
# `cores` workers don't oversubscribe the machine.
.make_cluster <- function(cores) {
  if (cores <= 1L) return(NULL)
  cl <- parallel::makePSOCKcluster(cores)
  ok <- tryCatch({
    parallel::clusterEvalQ(cl, {
      Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
                 OMP_NUM_THREADS = "1")
      loadNamespace("howbig")
      NULL
    })
    TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    parallel::stopCluster(cl)
    stop("Parallel workers could not load howbig; install the package to use ",
         "cores > 1, or run with cores = 1.", call. = FALSE)
  }
  cl
}

# Apply `f` over `x`, serially or on the worker pool `cl`. Task errors become
# NULL (dropped by the caller) so one failed fit doesn't sink a whole run.
.papply <- function(x, cl, f) {
  g <- function(i) tryCatch(f(i), error = function(e) NULL)
  # static chunking: the closure (which may carry the null pool) is sent
  # once per worker rather than once per task
  if (is.null(cl)) lapply(x, g) else parallel::parLapply(cl, x, g)
}

# Pre-flight runtime estimate (minutes) from the per-method log-linear model
# of production fit times; uses the mean group size.
.estimate_minutes <- function(method, p, n, n_fits, cores) {
  r <- howbig_runtime[howbig_runtime$method == method, ]
  if (!nrow(r)) return(0)
  per_fit <- exp(r$b0 + r$b_logp * log(p) + r$b_logn * log(mean(n))) / 60
  per_fit * ceiling(n_fits / max(1L, cores))
}

.fmt_minutes <- function(m) {
  if (m < 1) "<1 min" else if (m < 120) sprintf("%.0f min", m)
  else sprintf("%.1f hours", m / 60)
}

# In-session cache: repeat queries for the same design are instant.
.howbig_cache <- new.env(parent = emptyenv())

.cache_key <- function(...) {
  paste(c(..., howbig_family$m_pa), collapse = "|")
}
