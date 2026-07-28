#------------------------------------------------------
# Estimation Wrappers
# Modular functions for fitting JGL, SSJGL, and
# computing tuning parameters
#------------------------------------------------------

library(JGL)
library(spikeyglass)
library(multiGGM)
library(iDINGO)
library(iJRFNet)

#--- Tune JGL via BIC grid search ----

tune_jgl <- function(dat_list,
                     penalty = "fused",
                     lambda1_seq = NULL,
                     lambda2_seq = NULL,
                     n_lambda = 10) {

  # Data-adaptive default: set lambda range centered around the actual
  # signal level in the pooled sample covariance. Using the median of
  # |off-diagonal entries| as the anchor ensures the range brackets
  # the signal, preventing BIC from selecting overly sparse solutions.
  if (is.null(lambda1_seq)) {
    S_pool <- (cov(dat_list[[1]]) + cov(dat_list[[2]])) / 2
    off_diag <- abs(S_pool[upper.tri(S_pool)])
    signal_level <- median(off_diag)
    lam_max <- signal_level * 2
    lam_min <- signal_level / 50
    lambda1_seq <- exp(seq(log(lam_min), log(lam_max), length.out = n_lambda))
  }

  # Stage 1: search lambda1, holding lambda2 at minimum
  if (is.null(lambda2_seq)) lambda2_seq <- lambda1_seq

  bic_l1 <- numeric(length(lambda1_seq))
  for (i in seq_along(lambda1_seq)) {
    fit <- JGL(Y = dat_list,
               penalty = penalty,
               lambda1 = lambda1_seq[i],
               lambda2 = lambda2_seq[1],
               return.whole.theta = TRUE)
    bic_l1[i] <- vBIC(dat_list, fit, thr = 1e-4)
  }

  best_l1 <- lambda1_seq[which.min(bic_l1)]

  # Stage 2: search lambda2, using best lambda1
  bic_l2 <- numeric(length(lambda2_seq))
  for (j in seq_along(lambda2_seq)) {
    fit <- JGL(Y = dat_list,
               penalty = penalty,
               lambda1 = best_l1,
               lambda2 = lambda2_seq[j],
               return.whole.theta = TRUE)
    bic_l2[j] <- vBIC(dat_list, fit, thr = 1e-4)
  }

  best_l2 <- lambda2_seq[which.min(bic_l2)]

  return(list(lambda1 = best_l1, lambda2 = best_l2,
              bic_l1 = bic_l1, bic_l2 = bic_l2,
              lambda1_seq = lambda1_seq, lambda2_seq = lambda2_seq))
}

#-----Run JGL with given or tuned lambdas-----------------------------------------------

run_jgl <- function(dat_list,
                    penalty = "fused",
                    lambda1 = NULL,
                    lambda2 = NULL,
                    tune = TRUE,
                    lambda1_seq = NULL,
                    lambda2_seq = NULL) {

  if (tune && (is.null(lambda1) || is.null(lambda2))) {
    tuned <- tune_jgl(dat_list, penalty, lambda1_seq, lambda2_seq)
    lambda1 <- tuned$lambda1
    lambda2 <- tuned$lambda2
  }

  fit <- JGL(Y = dat_list,
             penalty = penalty,
             lambda1 = lambda1,
             lambda2 = lambda2,
             return.whole.theta = TRUE)

  return(list(
    fit = fit,
    theta = fit$theta,
    lambda1 = lambda1,
    lambda2 = lambda2
  ))
}

#------- Tune JGL via cross-validation------------------------------------------------

cv_jgl <- function(dat_list, penalty = "fused", nfolds = 5,
                   lambda1_seq = NULL, n_lambda = 10) {

  n1 <- nrow(dat_list[[1]])
  n2 <- nrow(dat_list[[2]])

  # Data-adaptive lambda range
  if (is.null(lambda1_seq)) {
    S_pool <- (cov(dat_list[[1]]) + cov(dat_list[[2]])) / 2
    off_diag <- abs(S_pool[upper.tri(S_pool)])
    signal_level <- median(off_diag)
    lambda1_seq <- exp(seq(log(signal_level / 50), log(signal_level * 2),
                           length.out = n_lambda))
  }

  # Create folds
  folds1 <- sample(rep(1:nfolds, length.out = n1))
  folds2 <- sample(rep(1:nfolds, length.out = n2))

  # For each lambda, compute held-out negative log-likelihood
  cv_scores <- numeric(length(lambda1_seq))

  for (li in seq_along(lambda1_seq)) {
    lam <- lambda1_seq[li]
    fold_scores <- numeric(nfolds)

    for (k in 1:nfolds) {
      train1 <- dat_list[[1]][folds1 != k, , drop = FALSE]
      train2 <- dat_list[[2]][folds2 != k, , drop = FALSE]
      test1 <- dat_list[[1]][folds1 == k, , drop = FALSE]
      test2 <- dat_list[[2]][folds2 == k, , drop = FALSE]

      fit <- tryCatch(
        JGL(Y = list(train1, train2),
            penalty = penalty,
            lambda1 = lam,
            lambda2 = lam / 5,
            return.whole.theta = TRUE),
        error = function(e) NULL
      )

      if (is.null(fit)) {
        fold_scores[k] <- Inf
        next
      }

      # Held-out negative log-likelihood: -log|Theta| + tr(S_test * Theta)
      for (g in 1:2) {
        S_test <- cov(if (g == 1) test1 else test2)
        Theta_hat <- fit$theta[[g]]
        log_det <- determinant(Theta_hat, logarithm = TRUE)$modulus[1]
        fold_scores[k] <- fold_scores[k] + sum(diag(S_test %*% Theta_hat)) - log_det
      }
    }

    cv_scores[li] <- mean(fold_scores[is.finite(fold_scores)])
  }

  best_lam <- lambda1_seq[which.min(cv_scores)]

  # Refit on full data
  fit_final <- JGL(Y = dat_list, penalty = penalty,
                   lambda1 = best_lam, lambda2 = best_lam / 5,
                   return.whole.theta = TRUE)

  return(list(
    theta = fit_final$theta,
    lambda1 = best_lam,
    cv_scores = cv_scores,
    lambda1_seq = lambda1_seq
  ))
}

#---- Run SSJGL (spike-and-slab joint graphical lasso) ----
#
# SSJGL places a spike-and-slab prior on each off-diagonal
# element of the precision matrix. The "spike" (v0) shrinks
# excluded edges toward zero; the "slab" (v1) allows included
# edges to be estimated freely.
#
# v0s LADDER JUSTIFICATION:
#   The v0s parameter controls the spike variance in the
#   spike-and-slab prior. Smaller v0 = more aggressive
#   shrinkage of excluded edges. ssjgl() uses an annealing
#   ladder: starting from large v0 (loose) and decreasing
#   to small v0 (strict), warm-starting each step from the
#   previous solution. We use 3 steps:
#     v0s = lambda1 / (lambda1 + c(5, 10, 15))
#   which gives v0 ~ [0.167, 0.091, 0.0625] for lambda1=1.
#   Going beyond 3 steps (v0 < 0.05) tends to over-sparsify
#   the graph on p=100, n=100 data, removing true edges.
#   The 3-step ladder balances sparsity with edge retention.
#
# problist1 OUTPUT:
#   fit$problist1[[v]] is a SHARED p x p matrix of edge
#   inclusion probabilities at each v0 step. This is NOT
#   per-group — it represents P(gamma_ij = 1), where gamma
#   is the shared binary indicator for whether edge (i,j)
#   is present in the joint graph structure. Under the
#   fused/group penalty, SSJGL assumes a common sparsity
#   pattern across groups (shared gamma), with group-specific
#   edge weights. Thus problist1 can be used for AUC
#   computation as a continuous edge score, but cannot
#   directly identify differential edges (since it is shared).
#   Differential detection relies on |theta_g1 - theta_g2|.
#------------------------------------------------------

run_ssjgl <- function(dat_list,
                      penalty = "fused",
                      lambda0 = 1,
                      lambda1 = 1,
                      lambda2 = 1,
                      v1 = 1,
                      v0s = NULL,
                      a = 1, b = 1,
                      doubly = TRUE,
                      normalize = FALSE) {

  if (is.null(v0s)) {
    lam_eff <- lambda1 + (1:3) * 5
    v0s <- lambda1 / lam_eff
  }

  fit <- spikeyglass::ssjgl(Y = dat_list,
               penalty = penalty,
               lambda0 = lambda0,
               lambda1 = lambda1,
               lambda2 = lambda2,
               v1 = v1,
               v0s = v0s,
               a = a, b = b,
               doubly = doubly,
               normalize = normalize)

  return(list(
    fit = fit,
    thetalist = fit$thetalist,
    v0s = v0s
  ))
}

#------------------------------------------------------
# Method dispatcher — split into fit + score
#
# The pipeline separates two concerns:
#   fit_method(method_name, dat_list, params)
#     Runs the actual method. Returns raw fit output (theta
#     matrices, PIPs, partial correlations, importance, ...)
#     plus metadata. NO scoring / metrics. No mutation of the
#     raw output (e.g., no PIP masking applied here).
#     -> one file per rep saved to results/raw_fits/
#
#   score_method(method_name, fit_output, gold, <scoring opts>)
#     Takes a raw fit output and a gold standard; applies the
#     scoring rules (thresholds, PIP mask, inference at a given
#     FDR) to produce the one-row metrics data.frame that 05c
#     merges. Cheap; rerun whenever a scoring option changes.
#     -> metrics table saved to results/method_results/
#
# run_method() is retained as a thin orchestrator that calls
# fit_method() then score_method(). Useful for ad-hoc / notebook
# work where you just want the metrics, not the raw fit.
#
# Used by: 05b_run_methods.R (fit phase) and
#          05d_score_results.R (score phase).
#------------------------------------------------------

# --- Helper: build the NA-row schema used when a rep fails ---
.init_metrics_na_row <- function(method_name) {
  data.frame(
    method = method_name,
    # --- Group (aggregated) ---
    group_TPR = NA_real_, group_FPR = NA_real_,
    group_FDR = NA_real_, group_MCC = NA_real_,
    group_AUC = NA_real_, group_AUPR = NA_real_,
    group_frob_rel = NA_real_,
    group_TP = NA_integer_, group_FP = NA_integer_, group_FN = NA_integer_,
    group_edges_true = NA_integer_, group_edges_est = NA_integer_,
    group_KL = NA_real_,
    # --- Per-group diagnostics ---
    g1_TPR = NA_real_, g1_FPR = NA_real_,
    g1_TP = NA_integer_, g1_FP = NA_integer_, g1_FN = NA_integer_,
    g1_edges_true = NA_integer_, g1_edges_est = NA_integer_,
    g1_FDR = NA_real_, g1_MCC = NA_real_,
    g1_AUC = NA_real_, g1_AUPR = NA_real_,
    g1_frob_rel = NA_real_, g1_KL = NA_real_,
    g2_TPR = NA_real_, g2_FPR = NA_real_,
    g2_TP = NA_integer_, g2_FP = NA_integer_, g2_FN = NA_integer_,
    g2_edges_true = NA_integer_, g2_edges_est = NA_integer_,
    g2_FDR = NA_real_, g2_MCC = NA_real_,
    g2_AUC = NA_real_, g2_AUPR = NA_real_,
    g2_frob_rel = NA_real_, g2_KL = NA_real_,
    # --- Differential (primary) ---
    diff_TPR = NA_real_, diff_FPR = NA_real_,
    diff_TP = NA_integer_, diff_FP = NA_integer_, diff_FN = NA_integer_,
    diff_TN = NA_integer_,
    diff_edges_true = NA_integer_, diff_edges_est = NA_integer_,
    diff_FDR = NA_real_, diff_MCC = NA_real_,
    diff_AUC = NA_real_, diff_AUPR = NA_real_,
    diff_frob_rel = NA_real_,
    # --- Differential (secondary XOR diagnostic) ---
    diff_XOR_TPR = NA_real_, diff_XOR_FPR = NA_real_,
    # --- Inference columns ---
    # Source depends on method:
    #   Frequentist (JGL/SSJGL/JRF): cached permutation null
    #     via apply_cached_inference(). Per-group inference
    #     columns left NA.
    #   Bayesian (multiGGM): direct PIP-based FDR via
    #     apply_bayesian_fdr_pip(). ALL inference columns
    #     (diff_*, g1_*, g2_*, group_*) populated.
    diff_TPR_inference = NA_real_,
    diff_FPR_inference = NA_real_,
    diff_FDR_inference = NA_real_,
    diff_MCC_inference = NA_real_,
    diff_Power_inference = NA_real_,
    g1_TPR_inference = NA_real_, g1_FPR_inference = NA_real_,
    g1_FDR_inference = NA_real_, g1_MCC_inference = NA_real_,
    g2_TPR_inference = NA_real_, g2_FPR_inference = NA_real_,
    g2_FDR_inference = NA_real_, g2_MCC_inference = NA_real_,
    group_TPR_inference = NA_real_, group_FPR_inference = NA_real_,
    group_FDR_inference = NA_real_, group_MCC_inference = NA_real_,
    # --- Run metadata ---
    runtime_sec = NA_real_,
    lambda1 = NA_real_, lambda2 = NA_real_,
    status = "pending",
    stringsAsFactors = FALSE
  )
}

# --- Methods that produce estimated precision matrices ---
# (Used by score_method to route to calc_metrics.)
.precision_methods <- c("JGL_BIC", "JGL_CV", "JGL_BIC_GGL", "JGL_CV_GGL",
                         "SSJGL", "SSJGL_GGL",
                         "multiGGM", "multiGGM_long","multiGGM_short",
                         "InvCov")

#------------------------------------------------------
# fit_method: run a method and return its RAW output.
#
# Returns a list with at least:
#   status       : "success" or "error: <msg>"
#   runtime_sec  : wall-clock time of the fit call
#   method       : method_name (for reference)
#   params       : the hyperparameter list actually used
#   (plus method-specific fields: theta, pip, pcor, importance, ...)
#
# Never throws — errors are captured in $status.
#------------------------------------------------------

fit_method <- function(method_name, dat_list, params = NULL) {
  if (is.null(params)) params <- method_params(method_name)

  start <- proc.time()
  raw <- tryCatch({
    if (method_name == "JGL_BIC") {
      fit <- run_jgl(dat_list, penalty = "fused", tune = TRUE)
      list(theta = fit$theta,
           lambda1 = fit$lambda1, lambda2 = fit$lambda2)

    } else if (method_name == "JGL_CV") {
      fit <- cv_jgl(dat_list, penalty = "fused", nfolds = params$nfolds %||% 5)
      list(theta = fit$theta,
           lambda1 = fit$lambda1, lambda2 = fit$lambda1 / 5)

    } else if (method_name == "JGL_BIC_GGL") {
      fit <- run_jgl(dat_list, penalty = "group", tune = TRUE)
      list(theta = fit$theta,
           lambda1 = fit$lambda1, lambda2 = fit$lambda2)

    } else if (method_name == "JGL_CV_GGL") {
      fit <- cv_jgl(dat_list, penalty = "group", nfolds = params$nfolds %||% 5)
      list(theta = fit$theta,
           lambda1 = fit$lambda1, lambda2 = fit$lambda1 / 5)

    } else if (method_name %in% c("SSJGL", "SSJGL_GGL")) {
      penalty <- if (method_name == "SSJGL_GGL") "group" else "fused"
      fit <- run_ssjgl(dat_list, penalty = penalty)
      n_v0 <- length(fit$v0s)
      list(theta       = fit$thetalist[[n_v0]],
           problist1   = fit$fit$problist1[[n_v0]],
           problist2   = fit$fit$problist2[[n_v0]],
           v0s         = fit$v0s,
           lambda1     = NA_real_, lambda2 = NA_real_)

    } else if (method_name %in% c("multiGGM", "multiGGM_long","multiGGM_short")) {
      burnin  <- params$burnin  %||% 5000
      nsave   <- params$nsave   %||% 1000
      nchains <- params$nchains %||% 2
      fit <- run_multiggm(dat_list, burnin = burnin, nsave = nsave,
                           nchains = nchains)
      list(theta    = fit$theta,    # unmasked posterior mean
           pcor     = fit$pcor,
           pip      = fit$pip,
           pip_diff = fit$pip_diff)

    } else if (method_name == "DiNGO") {
      fit <- run_dingo(dat_list, B = params$B %||% 0)
      list(theta             = fit$theta,
           pcor              = fit$pcor,
           diff_score_matrix = fit$diff_score_matrix,
           diff_type         = fit$diff_type,
           genepair          = fit$genepair,
           p_val             = fit$p_val)

    } else if (method_name == "JRF") {
      fit <- run_jrf(dat_list, ntree = params$ntree %||% 1000)
      list(importance = fit$importance)

    } else if (method_name == "InvCov") {
      fit <- run_invcov(dat_list, regularize = params$regularize %||% "ginv")
      list(theta   = fit$theta,
           lambda1 = NA_real_, lambda2 = NA_real_,
           singular_g1 = fit$singular_g1,
           singular_g2 = fit$singular_g2)

    } else {
      stop(sprintf("Unknown method: %s", method_name))
    }
  }, error = function(e) list(.err = conditionMessage(e)))

  elapsed <- as.numeric((proc.time() - start)["elapsed"])

  if (!is.null(raw$.err)) {
    return(list(method      = method_name,
                params      = params,
                status      = paste0("error: ", raw$.err),
                runtime_sec = elapsed))
  }

  c(list(method      = method_name,
         params      = params,
         status      = "success",
         runtime_sec = elapsed),
    raw)
}

#------------------------------------------------------
# score_method: apply scoring rules to a raw fit.
#
# Takes a raw fit output (from fit_method or loaded from disk)
# plus a gold standard and returns the one-row metrics
# data.frame used by 05c / the report scripts.
#
# Scoring options are SEPARATE from fit-time hyperparameters:
# you can rerun score_method with different thresholds without
# re-fitting anything.
#
# Args:
#   method_name          : method label
#   fit                  : output of fit_method() (or loaded from RDS)
#   gold                 : gold-standard list (with $adj, $omega, $truth, ...)
#   setting              : setting label (used for cached-null lookup)
#   threshold            : |theta| > threshold => per-group edge.
#                           Was 1e-10; bumped default to 0.01 so the
#                           same rule applies to JGL / SSJGL / masked
#                           multiGGM (all post-shrinkage precisions).
#   pcor_threshold       : DiNGO |pcor| > this => per-group edge.
#   diff_threshold       : |theta^(1) - theta^(2)| > this => differential.
#   diff_score_threshold : DiNGO bootstrap |z| > this => differential
#                           (only when diff_score_matrix is z-scores).
#   pip_threshold        : retained for auxiliary PIP scoring.
#   pip_mask_thr         : multiGGM PIP support mask applied at score
#                           time. 0 disables; default 0.05 matches the
#                           value determined by diagnostics (max
#                           non-edge PIP ~ 0.053).
#   fdr_target           : target FDR for diff_*_inference columns.
#   null_correction      : "pooled" (default; new), "bh", or "none".
#------------------------------------------------------

score_method <- function(method_name, fit, gold,
                          setting             = NULL,
                          threshold           = 0.01,
                          pcor_threshold      = 0.01,
                          diff_threshold      = 1e-5,   # nominal-zero floor;
                                                          # matches 05d default
                          diff_score_threshold = 2.0,
                          pip_threshold       = 0.5,
                          pip_mask_thr        = 0.05,
                          fdr_target          = 0.05,
                          null_correction     = "pooled") {

  na_row <- .init_metrics_na_row(method_name)

  # --- Handle fit failures up-front ---
  if (isTRUE(startsWith(fit$status %||% "", "error"))) {
    na_row$runtime_sec <- fit$runtime_sec %||% NA_real_
    na_row$status      <- fit$status
    return(na_row)
  }

  # --- Apply method-specific post-hoc mutations at score time ---
  # multiGGM: PIP-based support mask (does NOT change the stored raw
  # fit; only affects what gets passed to calc_metrics).
  theta_for_scoring <- fit$theta
  if (method_name %in% c("multiGGM", "multiGGM_long", "multiGGM_short") &&
      pip_mask_thr > 0 && !is.null(fit$pip)) {
    theta_for_scoring <- list(
      fit$theta[[1]] * (fit$pip[,,1] >= pip_mask_thr),
      fit$theta[[2]] * (fit$pip[,,2] >= pip_mask_thr)
    )
  }

  # --- Dispatch to the right calc_metrics flavor ---
  metrics <- tryCatch({
    if (method_name %in% .precision_methods) {
      if (is.null(theta_for_scoring)) stop("theta is NULL")
      # SSJGL: surface shared PIP as an auxiliary continuous score
      edge_scores <- NULL
      if (!is.null(fit$problist1)) {
        edge_scores <- list(
          g1   = fit$problist1,
          g2   = fit$problist1,
          diff = abs(theta_for_scoring[[1]] - theta_for_scoring[[2]])
        )
      }
      calc_metrics(theta_for_scoring, gold,
                   threshold      = threshold,
                   diff_threshold = diff_threshold,
                   edge_scores    = edge_scores)

    } else if (method_name == "DiNGO") {
      if (is.null(fit$pcor)) stop("pcor is NULL")
      ds_thr <- if (identical(fit$diff_type, "raw_pcor_diff")) diff_threshold
                else diff_score_threshold
      calc_metrics_from_pcor(fit$pcor, gold,
                              pcor_threshold     = pcor_threshold,
                              diff_threshold     = diff_threshold,
                              diff_score_matrix  = fit$diff_score_matrix,
                              diff_score_threshold = ds_thr,
                              theta_list         = fit$theta)

    } else if (method_name == "JRF") {
      if (is.null(fit$importance)) stop("importance is NULL")
      calc_metrics_from_importance(fit$importance, gold)

    } else {
      stop(sprintf("No scoring pathway for method: %s", method_name))
    }
  }, error = function(e) NULL)

  if (is.null(metrics)) {
    na_row$runtime_sec <- fit$runtime_sec %||% NA_real_
    na_row$status      <- "metrics_error"
    return(na_row)
  }

  # --- Inference columns ---
  # Route depends on method. For multiGGM (Bayesian SSVS) we use
  # the posterior inclusion probabilities directly — no permutation
  # null. For frequentist methods we consult the cached permutation
  # null. See NOTES_bayesian_fdr.md for the theory.
  inf_result <- NULL

  if (method_name %in% c("multiGGM", "multiGGM_long", "multiGGM_short")) {
    # Bayesian FDR via direct posterior probability procedure.
    # Populates BOTH diff_*_inference AND per-group inference
    # columns (g1_*, g2_*, group_*).
    inf_result <- tryCatch(
      apply_bayesian_fdr_pip(pip_per_group = fit$pip,
                              pip_diff      = fit$pip_diff,
                              gold          = gold,
                              fdr_target    = fdr_target),
      error = function(e) NULL
    )
  } else {
    # Frequentist methods: two inference sources.
    #
    # (A) DIFFERENTIAL inference -- cached permutation null
    #     (apply_cached_inference). Populates the 5 diff_*_inference
    #     columns.
    # (B) PER-GROUP inference -- cached stability-selection
    #     frequencies (apply_stability_fdr). Populates the 12
    #     g1_/g2_/group_*_inference columns.
    #
    # Either cache (or both) can be missing; the corresponding
    # inference columns simply stay NA. See 15_compute_null.R for
    # the permutation null and 19_stability_selection.R for
    # stability frequencies, and NOTES_bayesian_fdr.md for the
    # underlying theory (same direct-posterior FDR procedure
    # drives both paths; only the quantity being ranked differs).

    # (A) Differential via permutation null
    observed_score <- NULL
    if (method_name %in% .precision_methods && !is.null(theta_for_scoring)) {
      observed_score <- abs(theta_for_scoring[[1]] - theta_for_scoring[[2]])
    } else if (method_name == "DiNGO" && !is.null(fit$pcor)) {
      observed_score <- abs(fit$pcor[[1]] - fit$pcor[[2]])
    } else if (method_name == "JRF" && !is.null(fit$importance)) {
      observed_score <- fit$importance
    }

    inf_diff <- tryCatch(
      apply_cached_inference(method_name     = method_name,
                              setting         = setting,
                              observed_score  = observed_score,
                              gold            = gold,
                              fdr_target      = fdr_target,
                              correction      = null_correction),
      error = function(e) NULL
    )

    # (B) Per-group via pooled permutation null (JGL / SSJGL only;
    #     other frequentist methods skip silently). Replaces the
    #     previous stability-selection approach so the inference
    #     works under fresh-GS mode (where stability frequencies
    #     would be network-specific and not transferable across
    #     reps). See 19_stability_selection.R header for details
    #     and apply_pooled_pergroup_inference() in 04_calc_metrics.R
    #     for the inference math.
    inf_grp <- NULL
    if (method_name %in% c("JGL_BIC", "JGL_CV",
                            "JGL_BIC_GGL", "JGL_CV_GGL",
                            "SSJGL", "SSJGL_GGL") &&
        !is.null(setting) && !is.null(theta_for_scoring)) {
      inf_grp <- tryCatch(
        apply_pooled_pergroup_inference(
          method_name       = method_name,
          setting           = setting,
          observed_score_g1 = abs(theta_for_scoring[[1]]),
          observed_score_g2 = abs(theta_for_scoring[[2]]),
          gold              = gold,
          fdr_target        = fdr_target,
          correction        = null_correction
        ),
        error = function(e) NULL
      )
    }

    # Merge the two sources
    inf_result <- c(inf_diff %||% list(), inf_grp %||% list())
    if (length(inf_result) == 0) inf_result <- NULL
  }

  if (!is.null(inf_result)) {
    # Splice in only the columns the inference function produced;
    # leave the rest at the NA default from .build_result_row().
    for (nm in names(inf_result)) {
      if (nm %in% names(metrics)) metrics[[nm]] <- inf_result[[nm]]
    }
  }

  # --- Populate run metadata ---
  metrics$method      <- method_name
  metrics$runtime_sec <- fit$runtime_sec %||% NA_real_
  metrics$lambda1     <- fit$lambda1 %||% NA_real_
  metrics$lambda2     <- fit$lambda2 %||% NA_real_
  metrics$status      <- "success"
  metrics
}

#------------------------------------------------------
# run_method: orchestrator (fit + score in one shot).
#
# Used by interactive / ad-hoc code paths. The production
# pipeline (05b fit + 05d score) calls fit_method and
# score_method directly.
#
# The signature below intentionally mirrors score_method's
# scoring options, so existing callers continue to work.
#------------------------------------------------------

run_method <- function(method_name, dat_list, gold,
                       setting = NULL,
                       threshold = 0.01,
                       pcor_threshold = 0.01,
                       diff_threshold = 1e-5,         # nominal-zero floor
                       pip_threshold = 0.5,
                       pip_mask_thr = 0.05,
                       diff_score_threshold = 2.0,
                       fdr_target = 0.05,
                       null_correction = "pooled") {

  fit <- fit_method(method_name, dat_list)
  score_method(method_name, fit, gold,
                setting              = setting,
                threshold            = threshold,
                pcor_threshold       = pcor_threshold,
                diff_threshold       = diff_threshold,
                diff_score_threshold = diff_score_threshold,
                pip_threshold        = pip_threshold,
                pip_mask_thr         = pip_mask_thr,
                fdr_target           = fdr_target,
                null_correction      = null_correction)
}

#-------DiNGO wrapper -----
# Package: iDINGO (CRAN; published version)
# Returns: group-specific partial correlations (p×p) and
#          differential scores (bootstrap z-scores if B > 0,
#          else raw |R1 - R2|).
#
# NOTE on bootstrapping:
#   Setting B > 0 enables bootstrap-standardized differential
#   scores (~5-10 min extra for p=100). With B = 0 (default),
#   we use raw |R1 - R2| as the differential score. The
#   published iDINGO package does not support the `pval`
#   argument; we retain the parameter name for forward
#   compatibility but only pass it if the installed version
#   accepts it (see dynamic argument construction below).
#------------------------------------------------------

run_dingo <- function(dat_list, B = 0, pval = FALSE,
                       rhoarray = NULL, cores = 1, ...) {

  n1 <- nrow(dat_list[[1]])
  n2 <- nrow(dat_list[[2]])
  p  <- ncol(dat_list[[1]])

  # DiNGO expects a single data frame with column names + binary group vector
  col_names <- colnames(dat_list[[1]])
  if (is.null(col_names)) col_names <- paste0("V", 1:p)
  dat_combined <- rbind(dat_list[[1]], dat_list[[2]])
  colnames(dat_combined) <- col_names
  dat_mat <- data.frame(dat_combined)

  group_vec <- c(rep(-1, n1), rep(1, n2))

  # Build dingo() argument list dynamically — different iDINGO
  # versions accept different arguments (e.g., `pval` was removed
  # from the published version). We pass only arguments the
  # installed version supports.
  #
  # IMPORTANT: When B = 0, we must pass diff.score = FALSE. The
  # published iDINGO has a bug where diff.score = TRUE with B = 0
  # produces NaN in diff.score (because boot.diff is an empty matrix)
  # and then getEfronp() -> quantile() throws "missing values and
  # NaN's not allowed if 'na.rm' is FALSE". When diff.score = FALSE,
  # dingo() still returns R1 and R2 (which is all we need — we compute
  # the differential score ourselves as |R1 - R2| below).
  dingo_formals <- names(formals(dingo))
  use_diff_score <- B > 0
  dingo_args <- list(
    dat = dat_mat, x = group_vec,
    rhoarray = rhoarray,
    diff.score = use_diff_score,
    B = B,
    verbose = FALSE, cores = cores
  )
  if ("pval" %in% dingo_formals) dingo_args$pval <- pval

  fit <- do.call(dingo, dingo_args)

  # Reconstruct p×p partial correlation matrices from R1/R2 vectors
  # genepair is a matrix with columns gene1, gene2 (character names)
  pcor_g1 <- matrix(0, p, p)
  pcor_g2 <- matrix(0, p, p)
  rownames(pcor_g1) <- colnames(pcor_g1) <- col_names
  rownames(pcor_g2) <- colnames(pcor_g2) <- col_names

  # Map gene names to numeric indices
  name_to_idx <- setNames(seq_len(p), col_names)

  for (i in seq_len(nrow(fit$genepair))) {
    r <- name_to_idx[fit$genepair[i, 1]]
    c <- name_to_idx[fit$genepair[i, 2]]
    pcor_g1[r, c] <- pcor_g1[c, r] <- fit$R1[i]
    pcor_g2[r, c] <- pcor_g2[c, r] <- fit$R2[i]
  }
  diag(pcor_g1) <- 1
  diag(pcor_g2) <- 1

  # Build p×p differential score matrix from R1, R2
  # Use |R1 - R2| as the differential metric. This works with or without
  # bootstrap (diff.score requires B > 0 for standardization).
  # When B > 0 and diff.score is available, use that instead (z-score-like).
  diff_score_mat <- matrix(0, p, p)
  rownames(diff_score_mat) <- colnames(diff_score_mat) <- col_names

  use_raw_diff <- all(is.na(fit$diff.score))
  for (i in seq_len(nrow(fit$genepair))) {
    r <- name_to_idx[fit$genepair[i, 1]]
    c <- name_to_idx[fit$genepair[i, 2]]
    if (use_raw_diff) {
      # B=0: use raw partial correlation difference
      diff_score_mat[r, c] <- diff_score_mat[c, r] <- abs(fit$R1[i] - fit$R2[i])
    } else {
      # B>0: use bootstrap-standardized diff.score
      diff_score_mat[r, c] <- diff_score_mat[c, r] <- fit$diff.score[i]
    }
  }

  # Back-calculate precision matrices from partial correlations
  # Pcor = -cov2cor(Theta), so standardized Theta_std = -Pcor + 2*I
  # i.e., Theta_std[i,j] = -Pcor[i,j] for i!=j, Theta_std[i,i] = 1
  # To recover actual Theta, use sample precision for diagonal scaling:
  #   Theta = D %*% Theta_std %*% D where D = diag(sqrt(diag(S^{-1})))
  theta_list <- tryCatch({
    dat_full <- rbind(dat_list[[1]], dat_list[[2]])
    S_pool <- cov(dat_full)
    S_inv <- solve(S_pool)
    diag_vals <- diag(S_inv)
    # Guard: if n ~= p, numerical round-off can make diag slightly negative
    # even though S_inv should be PD. Treat negatives as 0 -> sqrt(0) = 0,
    # giving a zero-scaled precision matrix rather than NaN. Return NULL
    # if too many diag entries are invalid (back-calc not meaningful).
    if (any(!is.finite(diag_vals)) ||
        sum(diag_vals < 0) > nrow(S_inv) * 0.1) {
      NULL
    } else {
      d_scale_vec <- sqrt(pmax(diag_vals, 0))

      lapply(list(pcor_g1, pcor_g2), function(pcor_mat) {
        Theta_std <- -pcor_mat
        diag(Theta_std) <- 1
        D <- diag(d_scale_vec)
        Theta <- D %*% Theta_std %*% D
        Theta <- (Theta + t(Theta)) / 2
        Theta
      })
    }
  }, error = function(e) NULL)

  return(list(
    theta             = theta_list,
    pcor              = list(pcor_g1, pcor_g2),
    diff_scores       = fit$diff.score,
    diff_score_matrix = diff_score_mat,
    diff_type         = if (use_raw_diff) "raw_pcor_diff" else "bootstrap_zscore",
    genepair          = fit$genepair,
    p_val             = if (!is.null(fit$p.val)) fit$p.val else NULL,
    method            = "DiNGO"
  ))
}

#------- multiGGM wrapper (Bayesian SSVS joint GGM) ------
# Package: multiGGM (mljaniczek/multiGGMr)
# Returns: precision matrices, partial correlations,
#          PIPs, and differential PIPs
#
# MCMC PARAMETERS (comparison to Shaddox et al. 2020):
#   Shaddox paper used 10K burnin + 30K saved iterations
#   (total 40K MCMC samples). Our defaults of 5K + 1K are
#   much shorter for computational feasibility across many
#   reps. For production runs, consider burnin=10000,
#   nsave=5000 or higher. The paper used a single long chain;
#   we use nchains=2 with averaging for better mixing.
#
# HYPERPARAMETERS:
#   v0: spike variance. Shaddox used 0.02^2 = 0.0004.
#       Controls sparsity: smaller v0 = more aggressive
#       shrinkage of excluded edges. We use v0=0.01 as a
#       compromise between sparsity and edge detection.
#   w: MRF coupling weight. Shaddox used w=0.5 in their
#       simulations. Higher w = more confidence that graphs
#       are similar across groups. Default 0.8 may over-
#       couple; consider w=0.5 for heterogeneous groups.
#   a=1, b=4: Beta(1,4) prior on edge inclusion, giving
#       E[p_edge] = 0.20 (expects ~20% of edges included).
#------------------------------------------------------

run_multiggm <- function(dat_list,
                          method = "ssvs",
                          burnin = 2000,
                          nsave = 1000,
                          nchains = 2,
                          parallel = FALSE,
                          hyper = list(
                            v0 = 0.0004, v1 = 1,
                            lambda = 1,
                            a = 1, b = 4,
                            alpha = 2, beta = 5,
                            w = 0.8,
                            alpha_prop = 1, beta_prop = 1
                          ),
                          seed = 123,
                          ...) {

  fit <- multiggm_mcmc(
    data_list = dat_list,
    method    = method,
    burnin    = burnin,
    nsave     = nsave,
    nchains   = nchains,
    parallel  = parallel,
    hyper     = hyper,
    seed      = seed,
    ...
  )

  # For multi-chain fits, average PIPs and precision across chains
  if (inherits(fit, "multiggm_fit_list") && fit$nchains > 1) {
    nc <- fit$nchains

    pip_list <- lapply(1:nc, function(ch) pip_edges(fit, chain = ch))
    pip <- Reduce("+", pip_list) / nc

    pip_diff_list <- lapply(1:nc, function(ch) pip_diff_edge(fit, chain = ch))
    pip_diff <- Reduce("+", pip_diff_list) / nc

    theta_per_chain <- lapply(1:nc, function(ch) posterior_precision(fit, chain = ch))
    K <- length(theta_per_chain[[1]])
    theta <- lapply(1:K, function(k) {
      Reduce("+", lapply(theta_per_chain, `[[`, k)) / nc
    })
    names(theta) <- names(theta_per_chain[[1]])

    pcor_per_chain <- lapply(1:nc, function(ch) {
      lapply(posterior_precision(fit, chain = ch), function(omega) {
        d <- sqrt(diag(omega))
        pc <- -omega / outer(d, d)
        diag(pc) <- 1
        pc
      })
    })
    pcor <- lapply(1:K, function(k) {
      Reduce("+", lapply(pcor_per_chain, `[[`, k)) / nc
    })
  } else {
    pip <- pip_edges(fit)
    pip_diff <- pip_diff_edge(fit)
    theta <- posterior_precision(fit)

    pcor <- lapply(theta, function(omega) {
      d <- sqrt(diag(omega))
      pc <- -omega / outer(d, d)
      diag(pc) <- 1
      pc
    })
  }

  return(list(
    theta    = theta,
    pcor     = pcor,
    pip      = pip,
    pip_diff = pip_diff,
    method   = "multiGGM"
  ))
}

#------------------------------------------------------
# JRF wrapper (Joint Random Forest)
# Package: iJRFNet
# Returns: p×p importance matrix for differential
#          edges (|importance1 - importance2|)
#
# NOTE: iJRF expects p×n matrices (genes as rows,
#       samples as columns), standardized to mean=0
#       and variance=1. This is method-specific
#       preprocessing done inside the wrapper only.
#------------------------------------------------------

run_jrf <- function(dat_list, ntree = 1000, mtry = NULL, ...) {

  p <- ncol(dat_list[[1]])
  genes_name <- colnames(dat_list[[1]])
  if (is.null(genes_name)) genes_name <- paste0("V", 1:p)

  if (is.null(mtry)) mtry <- round(sqrt(p - 1))

  # iJRF expects p×n matrices (genes as rows), standardized.
  # Guard against constant columns (sd = 0) which would produce NaN
  # and cause downstream failures (e.g., quantile() errors in iJRF).
  X <- lapply(dat_list, function(d) {
    d_t <- t(d)  # p × n
    t(apply(d_t, 1, function(x) {
      s <- sd(x)
      if (!is.finite(s) || s == 0) {
        # constant feature: return zero vector (no signal)
        rep(0, length(x))
      } else {
        (x - mean(x)) / s
      }
    }))
  })

  fit <- iJRF(X = X, ntree = ntree, mtry = mtry,
              genes.name = genes_name)

  # fit is a data.frame with columns: gene1, gene2, importance1, importance2
  # Build p×p importance matrix: |importance1 - importance2| captures
  # differential signal between groups
  importance <- matrix(0, p, p)
  rownames(importance) <- colnames(importance) <- genes_name
  name_to_idx <- setNames(seq_len(p), genes_name)

  for (i in seq_len(nrow(fit))) {
    r <- name_to_idx[fit[i, "gene1"]]
    c <- name_to_idx[fit[i, "gene2"]]
    diff_imp <- abs(as.numeric(fit[i, "importance1"]) -
                    as.numeric(fit[i, "importance2"]))
    importance[r, c] <- importance[c, r] <- diff_imp
  }

  return(list(
    theta      = NULL,
    importance = importance,
    raw_fit    = fit,
    method     = "JRF"
  ))
}

#------- InvCov wrapper (unregularized MLE / pseudoinverse) -----
# Per-group precision = inverse of sample covariance.
#
# When n > p the sample covariance is invertible and this gives the
# Gaussian MLE precision. When n <= p the sample covariance is rank-
# deficient (singular) and we fall back to the Moore-Penrose
# pseudoinverse (MASS::ginv) so the method always returns *something*
# — the resulting matrix is uninformative-by-construction in the
# rank-deficient regime, which is exactly the story this baseline
# tells. The $singular_g1 / $singular_g2 flags record which groups
# fell back to ginv() so downstream consumers can flag/audit.
#
# Useful as a "no-regularization floor" comparator: regularized
# methods (JGL, SSJGL, multiGGM) should beat this at small n/p,
# while at large n/p the gap should close.
#
# Inputs:
#   dat_list   : list of n x p data matrices (samples in rows)
#   regularize : "ginv" (default) -> pseudoinverse on singular S
#                "drop"            -> return NA matrix on singular S
#
# Returns: list(theta = list of K p x p precision matrices,
#               singular_g1, singular_g2 = logical flags,
#               method = "InvCov")
#------------------------------------------------------

run_invcov <- function(dat_list, regularize = c("ginv", "drop"), ...) {
  regularize <- match.arg(regularize)
  K <- length(dat_list)
  p <- ncol(dat_list[[1]])

  theta <- vector("list", K)
  singular_flags <- logical(K)

  for (k in seq_len(K)) {
    X  <- dat_list[[k]]
    Xc <- scale(X, center = TRUE, scale = FALSE)
    nk <- nrow(Xc)
    S  <- crossprod(Xc) / max(nk - 1L, 1L)

    out <- tryCatch(solve(S), error = function(e) NULL)

    if (is.null(out)) {
      singular_flags[k] <- TRUE
      if (regularize == "ginv") {
        if (!requireNamespace("MASS", quietly = TRUE)) {
          stop("InvCov: MASS package needed for ginv() fallback. ",
               "Install MASS or use regularize = 'drop'.")
        }
        out <- MASS::ginv(S)
      } else {
        out <- matrix(NA_real_, p, p)
      }
    }

    # Symmetrize and add genenames if available
    out <- (out + t(out)) / 2
    if (!is.null(colnames(X))) {
      rownames(out) <- colnames(out) <- colnames(X)
    }
    theta[[k]] <- out
  }

  if (!is.null(names(dat_list))) names(theta) <- names(dat_list)

  list(theta       = theta,
       singular_g1 = singular_flags[1],
       singular_g2 = if (K >= 2) singular_flags[2] else NA,
       method      = "InvCov")
}
