#------------------------------------------------------
# Performance Metrics
# Functions to compare estimated networks against gold
# standard truths.
#
# SCHEMA (as of Phase 0 refactor):
#
# Per-group diagnostics (kept in the result data frame but not
# typically shown in the main report tables/figures):
#   g1_*, g2_* — per-group versions of each metric
#
# Main reported columns (what the report uses):
#   group_*  — averaged across groups for rates, summed across
#              groups for counts (see convention below)
#   diff_*   — differential (cross-group) metric
#
# Convention for group_* aggregation:
#   - Rates (TPR, FPR, FDR, MCC, AUC, AUPR, frob_rel):  mean
#   - Counts (TP, FP, FN, edges_true, edges_est):       sum
#   - KL:  group_KL = g1_KL + g2_KL (summed; matches historical sum_KL)
#
# Differential edge primary rule (used for diff_TPR/FPR/FDR/MCC):
#   |theta_hat^(1) - theta_hat^(2)| > diff_threshold
# (Previously XOR of group adjacencies. The old XOR rule is
#  retained as a secondary diagnostic in diff_XOR_TPR/FPR.)
#------------------------------------------------------

# --- Helper: threshold precision matrix to adjacency ---
threshold_to_adjacency <- function(theta, threshold = 0.01) {
  adj <- (abs(theta) > threshold) * 1L
  diag(adj) <- 0L
  adj
}

# --- Helper: pick the XOR-based truth for the diff_XOR_* diagnostic ---
# Since the Delta_A redefinition, gold$truth$Delta_A is the primary
# magnitude-threshold truth, while gold$truth$Delta_A_XOR is the
# structural XOR truth. Older cached gold standards (pre-refactor)
# don't carry Delta_A_XOR; in that case we fall back to Delta_A so
# the diagnostic still runs (semantically equivalent on old caches
# because Delta_A was defined as XOR there).
.xor_truth <- function(gold) {
  if (!is.null(gold$truth$Delta_A_XOR)) gold$truth$Delta_A_XOR
  else gold$truth$Delta_A
}

# --- Helper: partial correlations from precision matrix ---
pcor_from_precision <- function(Theta) {
  d <- sqrt(diag(Theta))
  Pcor <- -Theta / outer(d, d)
  diag(Pcor) <- 1
  Pcor
}

#------ AUC via trapezoid rule ------

compute_auc <- function(truth_vec, score_vec) {
  valid <- !is.na(truth_vec) & !is.na(score_vec)
  truth_vec <- truth_vec[valid]
  score_vec <- score_vec[valid]
  n_pos <- sum(truth_vec == 1)
  n_neg <- sum(truth_vec == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  ord <- order(score_vec, decreasing = TRUE)
  truth_sorted <- truth_vec[ord]
  TPR <- cumsum(truth_sorted == 1) / n_pos
  FPR <- cumsum(truth_sorted == 0) / n_neg
  TPR <- c(0, TPR); FPR <- c(0, FPR)
  sum(diff(FPR) * (TPR[-length(TPR)] + TPR[-1]) / 2)
}

#------ AUPR via trapezoid rule ------

compute_aupr <- function(truth_vec, score_vec) {
  valid <- !is.na(truth_vec) & !is.na(score_vec)
  truth_vec <- truth_vec[valid]
  score_vec <- score_vec[valid]
  n_pos <- sum(truth_vec == 1)
  if (n_pos == 0 || n_pos == length(truth_vec)) return(NA_real_)
  ord <- order(score_vec, decreasing = TRUE)
  truth_sorted <- truth_vec[ord]
  TP_cum <- cumsum(truth_sorted == 1)
  FP_cum <- cumsum(truth_sorted == 0)
  precision <- TP_cum / (TP_cum + FP_cum)
  recall <- TP_cum / n_pos
  recall <- c(0, recall); precision <- c(1, precision)
  sum(diff(recall) * (precision[-length(precision)] + precision[-1]) / 2)
}

#------ Relative Frobenius norm ------

compute_relative_frobenius <- function(est_mat, truth_mat) {
  truth_norm <- sqrt(sum(truth_mat^2))
  if (truth_norm == 0) return(NA_real_)
  sqrt(sum((est_mat - truth_mat)^2)) / truth_norm
}

#------ FDR ------

compute_fdr <- function(TP, FP) {
  if (is.na(TP) || is.na(FP)) return(NA_real_)
  denom <- TP + FP
  if (denom == 0) return(NA_real_)
  FP / denom
}

#------ MCC ------

compute_mcc <- function(TP, FP, TN, FN) {
  if (is.na(TP) || is.na(FP) || is.na(TN) || is.na(FN)) return(NA_real_)
  num <- as.numeric(TP) * TN - as.numeric(FP) * FN
  denom <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  if (denom == 0) return(NA_real_)
  num / denom
}

#------ KL divergence for one group ------
# D_KL(est || true) = tr(Sigma_true * Theta_hat) - log|Sigma_true * Theta_hat| - p

compute_kl_divergence <- function(theta_hat, sigma_true) {
  p <- nrow(theta_hat)
  tryCatch({
    prod_mat <- sigma_true %*% theta_hat
    tr_val <- sum(diag(prod_mat))
    log_det <- determinant(prod_mat, logarithm = TRUE)
    if (log_det$sign <= 0) return(NA_real_)
    tr_val - log_det$modulus[1] - p
  }, error = function(e) NA_real_)
}

#------ TPR, FPR, TP, FP, TN, FN from binary matrices ------

compute_tpr_fpr <- function(truth_adj, estimated_adj) {
  truth_vec <- truth_adj[upper.tri(truth_adj)]
  est_vec <- estimated_adj[upper.tri(estimated_adj)]
  TP <- sum(truth_vec == 1 & est_vec == 1)
  FP <- sum(truth_vec == 0 & est_vec == 1)
  TN <- sum(truth_vec == 0 & est_vec == 0)
  FN <- sum(truth_vec == 1 & est_vec == 0)
  TPR <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
  FPR <- ifelse((FP + TN) > 0, FP / (FP + TN), NA)
  list(TPR = TPR, FPR = FPR, TP = TP, FP = FP, TN = TN, FN = FN)
}

#------ Helper: build a result row in the new schema ------
# Takes perf objects (from compute_tpr_fpr) for g1, g2, diff, and
# optional diff_xor; plus AUC/AUPR/Frob/KL scalars. Returns a
# one-row data.frame with the unified schema used by all methods.

.build_result_row <- function(
  perf_g1, perf_g2, perf_diff, perf_diff_xor,
  g1_AUC, g2_AUC, diff_AUC,
  g1_AUPR, g2_AUPR, diff_AUPR,
  g1_frob_rel, g2_frob_rel, diff_frob_rel,
  g1_KL, g2_KL,
  g1_edges_true, g2_edges_true, diff_edges_true,
  g1_edges_est, g2_edges_est, diff_edges_est
) {
  # Per-group rates
  g1_FDR <- compute_fdr(perf_g1$TP, perf_g1$FP)
  g2_FDR <- compute_fdr(perf_g2$TP, perf_g2$FP)
  g1_MCC <- compute_mcc(perf_g1$TP, perf_g1$FP, perf_g1$TN, perf_g1$FN)
  g2_MCC <- compute_mcc(perf_g2$TP, perf_g2$FP, perf_g2$TN, perf_g2$FN)

  # Differential rates (from primary rule)
  diff_FDR <- compute_fdr(perf_diff$TP, perf_diff$FP)
  diff_MCC <- compute_mcc(perf_diff$TP, perf_diff$FP, perf_diff$TN, perf_diff$FN)

  # --- Group-level aggregates ---
  # Rates: mean across groups (NA-safe)
  safe_mean <- function(a, b) {
    if (is.na(a) && is.na(b)) NA_real_
    else if (is.na(a)) b
    else if (is.na(b)) a
    else (a + b) / 2
  }
  safe_sum <- function(a, b) {
    if (is.na(a) && is.na(b)) NA_real_
    else if (is.na(a)) b
    else if (is.na(b)) a
    else a + b
  }

  group_TPR  <- safe_mean(perf_g1$TPR, perf_g2$TPR)
  group_FPR  <- safe_mean(perf_g1$FPR, perf_g2$FPR)
  group_FDR  <- safe_mean(g1_FDR, g2_FDR)
  group_MCC  <- safe_mean(g1_MCC, g2_MCC)
  group_AUC  <- safe_mean(g1_AUC, g2_AUC)
  group_AUPR <- safe_mean(g1_AUPR, g2_AUPR)
  group_frob_rel <- safe_mean(g1_frob_rel, g2_frob_rel)

  # Counts: sum across groups
  group_TP <- as.integer(safe_sum(perf_g1$TP, perf_g2$TP))
  group_FP <- as.integer(safe_sum(perf_g1$FP, perf_g2$FP))
  group_FN <- as.integer(safe_sum(perf_g1$FN, perf_g2$FN))
  group_edges_true <- as.integer(safe_sum(g1_edges_true, g2_edges_true))
  group_edges_est  <- as.integer(safe_sum(g1_edges_est, g2_edges_est))

  # KL: sum across groups
  group_KL <- if (!is.na(g1_KL) && !is.na(g2_KL)) g1_KL + g2_KL else NA_real_

  # Build data.frame (fixed column order)
  data.frame(
    # --- Group (aggregated, main display) ---
    group_TPR = group_TPR, group_FPR = group_FPR,
    group_FDR = group_FDR, group_MCC = group_MCC,
    group_AUC = group_AUC, group_AUPR = group_AUPR,
    group_frob_rel = group_frob_rel,
    group_TP = group_TP, group_FP = group_FP, group_FN = group_FN,
    group_edges_true = group_edges_true, group_edges_est = group_edges_est,
    group_KL = group_KL,
    # --- Per-group diagnostics ---
    g1_TPR = perf_g1$TPR, g1_FPR = perf_g1$FPR,
    g1_TP = perf_g1$TP, g1_FP = perf_g1$FP, g1_FN = perf_g1$FN,
    g1_edges_true = g1_edges_true, g1_edges_est = g1_edges_est,
    g1_FDR = g1_FDR, g1_MCC = g1_MCC,
    g1_AUC = g1_AUC, g1_AUPR = g1_AUPR,
    g1_frob_rel = g1_frob_rel, g1_KL = g1_KL,
    g2_TPR = perf_g2$TPR, g2_FPR = perf_g2$FPR,
    g2_TP = perf_g2$TP, g2_FP = perf_g2$FP, g2_FN = perf_g2$FN,
    g2_edges_true = g2_edges_true, g2_edges_est = g2_edges_est,
    g2_FDR = g2_FDR, g2_MCC = g2_MCC,
    g2_AUC = g2_AUC, g2_AUPR = g2_AUPR,
    g2_frob_rel = g2_frob_rel, g2_KL = g2_KL,
    # --- Differential (primary = magnitude threshold) ---
    diff_TPR = perf_diff$TPR, diff_FPR = perf_diff$FPR,
    diff_TP = perf_diff$TP, diff_FP = perf_diff$FP, diff_FN = perf_diff$FN,
    diff_TN = perf_diff$TN,
    diff_edges_true = diff_edges_true, diff_edges_est = diff_edges_est,
    diff_FDR = diff_FDR, diff_MCC = diff_MCC,
    diff_AUC = diff_AUC, diff_AUPR = diff_AUPR,
    diff_frob_rel = diff_frob_rel,
    # --- Differential (secondary = XOR diagnostic) ---
    diff_XOR_TPR = if (!is.null(perf_diff_xor)) perf_diff_xor$TPR else NA_real_,
    diff_XOR_FPR = if (!is.null(perf_diff_xor)) perf_diff_xor$FPR else NA_real_,
    # --- Inference columns (populated by run_method / score_method) ---
    # Sources:
    #   - Frequentist methods (JGL / SSJGL / JRF): cached permutation
    #     null via apply_cached_inference() in 04_calc_metrics.R.
    #   - Bayesian multiGGM: direct PIP-based FDR via
    #     apply_bayesian_fdr_pip() (no permutation needed; see
    #     NOTES_bayesian_fdr.md for the theory).
    # Per-group inference columns are only populated when the method
    # provides a per-group inclusion score (currently: multiGGM's
    # pip[,,g]). Frequentist methods leave them NA.
    diff_TPR_inference = NA_real_,
    diff_FPR_inference = NA_real_,
    diff_FDR_inference = NA_real_,
    diff_MCC_inference = NA_real_,
    diff_Power_inference = NA_real_,
    g1_TPR_inference = NA_real_,
    g1_FPR_inference = NA_real_,
    g1_FDR_inference = NA_real_,
    g1_MCC_inference = NA_real_,
    g2_TPR_inference = NA_real_,
    g2_FPR_inference = NA_real_,
    g2_FDR_inference = NA_real_,
    g2_MCC_inference = NA_real_,
    group_TPR_inference = NA_real_,
    group_FPR_inference = NA_real_,
    group_FDR_inference = NA_real_,
    group_MCC_inference = NA_real_,
    stringsAsFactors = FALSE
  )
}

#======================================================
# calc_metrics: JGL, SSJGL, GGL variants, multiGGM
# Uses estimated precision matrices for both per-group
# edge detection AND differential detection (via magnitude
# threshold on |theta^(1) - theta^(2)|).
#
# Primary differential rule: |theta^(1) - theta^(2)| > diff_threshold
# Continuous score for AUC/AUPR: |theta^(1) - theta^(2)|
# Secondary XOR rule: XOR of per-group adjacencies -> diff_XOR_*
#======================================================

calc_metrics <- function(theta_list, gold,
                         threshold = 1e-10,          # edge = nonzero precision
                         diff_threshold = 1e-5,      # primary differential rule:
                                                      # nominal-zero floor; matches
                                                      # 05d --diff-threshold default
                         edge_scores = NULL) {
  ut <- upper.tri(theta_list[[1]])

  # --- Per-group adjacency detection ---
  adj_hat_g1 <- threshold_to_adjacency(theta_list[[1]], threshold)
  adj_hat_g2 <- threshold_to_adjacency(theta_list[[2]], threshold)
  perf_g1 <- compute_tpr_fpr(gold$adj$g1, adj_hat_g1)
  perf_g2 <- compute_tpr_fpr(gold$adj$g2, adj_hat_g2)

  # --- Differential detection (primary: magnitude threshold) ---
  diff_theta <- abs(theta_list[[1]] - theta_list[[2]])
  delta_hat_primary <- (diff_theta > diff_threshold) * 1L
  diag(delta_hat_primary) <- 0L
  perf_diff <- compute_tpr_fpr(gold$truth$Delta_A, delta_hat_primary)

  # --- Differential detection (secondary: XOR diagnostic) ---
  delta_hat_xor <- (adj_hat_g1 != adj_hat_g2) * 1L
  diag(delta_hat_xor) <- 0L
  perf_diff_xor <- compute_tpr_fpr(.xor_truth(gold), delta_hat_xor)

  # --- Frobenius norms ---
  g1_frob_rel <- compute_relative_frobenius(theta_list[[1]], gold$omega$g1)
  g2_frob_rel <- compute_relative_frobenius(theta_list[[2]], gold$omega$g2)
  diff_est  <- theta_list[[1]] - theta_list[[2]]
  diff_true <- gold$omega$g1 - gold$omega$g2
  diff_frob_rel <- compute_relative_frobenius(diff_est, diff_true)

  # --- AUC / AUPR ---
  g1_score   <- abs(theta_list[[1]])[ut]
  g2_score   <- abs(theta_list[[2]])[ut]
  diff_score <- diff_theta[ut]

  # Optional alternative continuous scores (e.g., SSJGL inclusion probs)
  if (!is.null(edge_scores)) {
    if (!is.null(edge_scores$g1))   g1_score   <- edge_scores$g1[ut]
    if (!is.null(edge_scores$g2))   g2_score   <- edge_scores$g2[ut]
    if (!is.null(edge_scores$diff)) diff_score <- edge_scores$diff[ut]
  }

  g1_AUC   <- compute_auc(gold$adj$g1[ut], g1_score)
  g2_AUC   <- compute_auc(gold$adj$g2[ut], g2_score)
  diff_AUC <- compute_auc(gold$truth$Delta_A[ut], diff_score)
  g1_AUPR   <- compute_aupr(gold$adj$g1[ut], g1_score)
  g2_AUPR   <- compute_aupr(gold$adj$g2[ut], g2_score)
  diff_AUPR <- compute_aupr(gold$truth$Delta_A[ut], diff_score)

  # --- KL divergence ---
  g1_KL <- tryCatch(compute_kl_divergence(theta_list[[1]], gold$sigma$g1),
                    error = function(e) NA_real_)
  g2_KL <- tryCatch(compute_kl_divergence(theta_list[[2]], gold$sigma$g2),
                    error = function(e) NA_real_)

  .build_result_row(
    perf_g1, perf_g2, perf_diff, perf_diff_xor,
    g1_AUC, g2_AUC, diff_AUC,
    g1_AUPR, g2_AUPR, diff_AUPR,
    g1_frob_rel, g2_frob_rel, diff_frob_rel,
    g1_KL, g2_KL,
    g1_edges_true = sum(gold$adj$g1[ut]),
    g2_edges_true = sum(gold$adj$g2[ut]),
    diff_edges_true = sum(gold$truth$Delta_A[ut]),
    g1_edges_est = sum(adj_hat_g1[ut]),
    g2_edges_est = sum(adj_hat_g2[ut]),
    diff_edges_est = sum(delta_hat_primary[ut])
  )
}

#======================================================
# calc_metrics_from_pcor: DiNGO
# Partial correlation matrices + optional back-calculated
# precision matrices for Frobenius/KL.
#======================================================

calc_metrics_from_pcor <- function(pcor_list, gold,
                                    pcor_threshold = 0.01,
                                    diff_threshold = 1e-5,    # nominal-zero floor
                                    diff_score_matrix = NULL,
                                    diff_score_threshold = 2.0,
                                    theta_list = NULL) {
  ut <- upper.tri(pcor_list[[1]])

  # --- Per-group adjacency from pcor magnitude ---
  adj_hat_g1 <- (abs(pcor_list[[1]]) > pcor_threshold) * 1L; diag(adj_hat_g1) <- 0L
  adj_hat_g2 <- (abs(pcor_list[[2]]) > pcor_threshold) * 1L; diag(adj_hat_g2) <- 0L
  perf_g1 <- compute_tpr_fpr(gold$adj$g1, adj_hat_g1)
  perf_g2 <- compute_tpr_fpr(gold$adj$g2, adj_hat_g2)

  # --- Differential detection (primary) ---
  if (!is.null(diff_score_matrix)) {
    # Bootstrap z-score route (B > 0)
    delta_hat_primary <- (abs(diff_score_matrix) > diff_score_threshold) * 1L
  } else {
    # Raw partial-correlation difference route (B = 0)
    diff_pcor <- abs(pcor_list[[1]] - pcor_list[[2]])
    delta_hat_primary <- (diff_pcor > diff_threshold) * 1L
  }
  diag(delta_hat_primary) <- 0L
  perf_diff <- compute_tpr_fpr(gold$truth$Delta_A, delta_hat_primary)

  # --- Differential XOR diagnostic ---
  delta_hat_xor <- (adj_hat_g1 != adj_hat_g2) * 1L; diag(delta_hat_xor) <- 0L
  perf_diff_xor <- compute_tpr_fpr(.xor_truth(gold), delta_hat_xor)

  # --- Frobenius norms ---
  if (!is.null(theta_list)) {
    g1_frob_rel <- compute_relative_frobenius(theta_list[[1]], gold$omega$g1)
    g2_frob_rel <- compute_relative_frobenius(theta_list[[2]], gold$omega$g2)
    diff_est  <- theta_list[[1]] - theta_list[[2]]
    diff_true <- gold$omega$g1 - gold$omega$g2
    diff_frob_rel <- compute_relative_frobenius(diff_est, diff_true)
  } else {
    g1_frob_rel <- compute_relative_frobenius(pcor_list[[1]], gold$truth$Pcor_g1)
    g2_frob_rel <- compute_relative_frobenius(pcor_list[[2]], gold$truth$Pcor_g2)
    diff_est  <- pcor_list[[1]] - pcor_list[[2]]
    diff_true <- gold$truth$Pcor_g1 - gold$truth$Pcor_g2
    diff_frob_rel <- compute_relative_frobenius(diff_est, diff_true)
  }

  # --- AUC / AUPR ---
  g1_AUC <- compute_auc(gold$adj$g1[ut], abs(pcor_list[[1]])[ut])
  g2_AUC <- compute_auc(gold$adj$g2[ut], abs(pcor_list[[2]])[ut])
  diff_continuous <- if (!is.null(diff_score_matrix)) {
    abs(diff_score_matrix)[ut]
  } else {
    abs(pcor_list[[1]] - pcor_list[[2]])[ut]
  }
  diff_AUC <- compute_auc(gold$truth$Delta_A[ut], diff_continuous)
  g1_AUPR <- compute_aupr(gold$adj$g1[ut], abs(pcor_list[[1]])[ut])
  g2_AUPR <- compute_aupr(gold$adj$g2[ut], abs(pcor_list[[2]])[ut])
  diff_AUPR <- compute_aupr(gold$truth$Delta_A[ut], diff_continuous)

  # --- KL divergence (only if back-calculated theta available) ---
  if (!is.null(theta_list)) {
    g1_KL <- tryCatch(compute_kl_divergence(theta_list[[1]], gold$sigma$g1),
                      error = function(e) NA_real_)
    g2_KL <- tryCatch(compute_kl_divergence(theta_list[[2]], gold$sigma$g2),
                      error = function(e) NA_real_)
  } else {
    g1_KL <- NA_real_
    g2_KL <- NA_real_
  }

  .build_result_row(
    perf_g1, perf_g2, perf_diff, perf_diff_xor,
    g1_AUC, g2_AUC, diff_AUC,
    g1_AUPR, g2_AUPR, diff_AUPR,
    g1_frob_rel, g2_frob_rel, diff_frob_rel,
    g1_KL, g2_KL,
    g1_edges_true = sum(gold$adj$g1[ut]),
    g2_edges_true = sum(gold$adj$g2[ut]),
    diff_edges_true = sum(gold$truth$Delta_A[ut]),
    g1_edges_est = sum(adj_hat_g1[ut]),
    g2_edges_est = sum(adj_hat_g2[ut]),
    diff_edges_est = sum(delta_hat_primary[ut])
  )
}

#======================================================
# calc_metrics_from_importance: JRF
# Only differential metrics are defined — per-group are NA.
# AUC / AUPR use raw importance scores. Threshold-dependent
# diff metrics are NA in the main analysis (they populate
# only when cached permutation inference is available, via
# run_method's inference path).
#======================================================

calc_metrics_from_importance <- function(importance_matrix, gold) {
  p <- nrow(importance_matrix)
  ut <- upper.tri(importance_matrix)
  scores <- importance_matrix[ut]

  # Per-group metrics: undefined for JRF
  na_perf <- list(TPR = NA_real_, FPR = NA_real_,
                  TP = NA_integer_, FP = NA_integer_,
                  TN = NA_integer_, FN = NA_integer_)
  perf_g1 <- na_perf; perf_g2 <- na_perf

  # Differential: threshold-free AUC/AUPR only in main analysis
  diff_AUC  <- compute_auc(gold$truth$Delta_A[ut], scores)
  diff_AUPR <- compute_aupr(gold$truth$Delta_A[ut], scores)

  # No threshold-based diff metrics in main analysis for JRF —
  # they populate only when cached permutation inference is used.
  perf_diff     <- na_perf
  perf_diff_xor <- na_perf

  .build_result_row(
    perf_g1, perf_g2, perf_diff, perf_diff_xor,
    g1_AUC = NA_real_, g2_AUC = NA_real_, diff_AUC = diff_AUC,
    g1_AUPR = NA_real_, g2_AUPR = NA_real_, diff_AUPR = diff_AUPR,
    g1_frob_rel = NA_real_, g2_frob_rel = NA_real_, diff_frob_rel = NA_real_,
    g1_KL = NA_real_, g2_KL = NA_real_,
    g1_edges_true = sum(gold$adj$g1[ut]),
    g2_edges_true = sum(gold$adj$g2[ut]),
    diff_edges_true = sum(gold$truth$Delta_A[ut]),
    g1_edges_est = NA_integer_,
    g2_edges_est = NA_integer_,
    diff_edges_est = NA_integer_
  )
}

#======================================================
# apply_cached_inference: per-edge inference via cached null
#
# Loads the pre-computed null distribution from
# results/null_distributions/<setting>_<method>.rds (produced
# by 15_compute_null.R) and uses it to call differential
# edges at a target FDR.
#
# Three correction modes for permutation nulls:
#
# (A) correction = "pooled" (DEFAULT, recommended)
#     Concatenate all n_perm x n_pairs null values into one pooled
#     null vector. Each edge's empirical p-value is computed
#     against that pooled distribution:
#       p_ij = (1 + #{null_pool >= obs_ij}) / (1 + |null_pool|)
#     This gives p-value resolution of 1 / (1 + n_perm * n_pairs),
#     which at n_perm = 30, p = 250 is ~1e-6 — fine enough to
#     survive BH correction on ~31k tests at FDR 0.05.
#     Assumes exchangeability of null distributions across edges,
#     which is approximately true for label-permutation nulls.
#     BH-adjusted; call edges with p_adj < fdr_target.
#
# (B) correction = "bh" (per-edge null, like the original)
#     p_ij = (1 + #{null_scores[[k]][i,j] >= obs_ij}) / (1 + n_perm)
#     Then BH-adjust. This is very conservative in high-dimensional
#     settings: at n_perm = 30, min p = 1/31 = 0.032, and BH at
#     FDR 0.05 on ~31k tests requires raw p < 0.05/31k = 1.6e-6,
#     so nothing passes regardless of the true signal. Retained
#     for comparison / diagnostics.
#
# (C) correction = "none"
#     Raw empirical p < fdr_target. No multiple-testing correction.
#     Liberal; matches what many method papers actually report.
#     Mode is named "none" but fdr_target is still the alpha cutoff.
#
# For bootstrap nulls (DiNGO) the cache stores per-edge SE. Then:
#   z_ij = obs_ij / SE_ij (Wald-type); two-sided p-value;
#   correction applied as above ("pooled" treats all |z| values
#   as one pool; "bh" is per-edge; "none" uses raw p).
#
# Returns a named list with 5 components (diff_TPR_inference,
# diff_FPR_inference, diff_FDR_inference, diff_MCC_inference,
# diff_Power_inference). All NA if no cache file is available
# or if inference fails.
#======================================================

apply_cached_inference <- function(method_name, setting, observed_score,
                                    gold, fdr_target = 0.05,
                                    null_dir = "results/null_distributions",
                                    correction = c("pooled", "bh", "none")) {
  correction <- match.arg(correction)
  na_result <- list(
    diff_TPR_inference   = NA_real_,
    diff_FPR_inference   = NA_real_,
    diff_FDR_inference   = NA_real_,
    diff_MCC_inference   = NA_real_,
    diff_Power_inference = NA_real_
  )

  if (is.null(setting) || is.null(observed_score)) return(na_result)

  null_file <- file.path(null_dir, sprintf("%s_%s.rds", setting, method_name))
  if (!file.exists(null_file)) return(na_result)

  null_rec <- tryCatch(readRDS(null_file), error = function(e) NULL)
  if (is.null(null_rec) || is.null(null_rec$data)) return(na_result)

  p <- nrow(observed_score)
  ut <- upper.tri(observed_score)
  obs_vec <- observed_score[ut]

  pvals <- tryCatch({
    if (identical(null_rec$data$type, "permutation")) {
      null_list <- null_rec$data$null_scores
      # Drop permutations where the refit failed (all NA or NULL)
      valid <- vapply(null_list,
                       function(m) !is.null(m) && !all(is.na(m)),
                       logical(1))
      null_list <- null_list[valid]
      if (length(null_list) == 0) return(NULL)

      if (correction == "pooled") {
        # Pool all null values into one long vector and rank obs
        # against it. Fine p-value resolution: 1/(1 + n_perm * n_edges).
        null_pool <- unlist(lapply(null_list, function(m) m[ut]),
                             use.names = FALSE)
        null_pool <- null_pool[is.finite(null_pool)]
        # ecdf-style empirical p: count pool exceedances for each obs
        sorted_pool <- sort(null_pool, decreasing = TRUE)
        # For each obs_ij, (# pool >= obs) = number of sorted_pool
        # entries that are >= obs, found via binary search.
        # findInterval on a non-increasing vector doesn't exist, so
        # invert: number >= obs = length(pool) - findInterval(obs-eps, sort(pool, ascending))
        asc_pool <- sort(null_pool)
        # # of entries in asc_pool that are < obs  =>  findInterval(obs, asc_pool, left.open = FALSE)
        # Therefore # >= obs = length(asc_pool) - # < obs
        n_lt <- findInterval(obs_vec - .Machine$double.eps, asc_pool)
        exceed <- length(asc_pool) - n_lt
        (1 + exceed) / (1 + length(asc_pool))

      } else {  # correction %in% c("bh", "none") — per-edge null
        # Stack upper triangles into an (n_perm x n_edge) matrix
        null_mat <- do.call(rbind, lapply(null_list, function(m) m[ut]))
        n_perm_eff <- nrow(null_mat)
        obs_row_rep <- matrix(obs_vec, n_perm_eff, length(obs_vec),
                               byrow = TRUE)
        exceed <- colSums(null_mat >= obs_row_rep, na.rm = TRUE)
        (1 + exceed) / (1 + n_perm_eff)
      }

    } else if (identical(null_rec$data$type, "bootstrap_se")) {
      se_mat <- null_rec$data$se_matrix
      if (is.null(se_mat) || !all(dim(se_mat) == dim(observed_score))) {
        return(NULL)
      }
      se_vec <- se_mat[ut]
      z <- ifelse(is.finite(se_vec) & se_vec > 1e-12,
                  obs_vec / se_vec, NA_real_)
      2 * (1 - pnorm(abs(z)))
    } else {
      return(NULL)
    }
  }, error = function(e) NULL)

  if (is.null(pvals)) return(na_result)

  # Apply correction and threshold
  called <- if (correction == "none") {
    as.integer(!is.na(pvals) & pvals < fdr_target)
  } else {
    # Both "pooled" and "bh" use BH on the vector of per-edge p's
    # (for "pooled", the p's are already drawn from a finer null
    # so BH still controls FDR at target but more edges pass it).
    p_adj <- p.adjust(pvals, method = "BH")
    as.integer(!is.na(p_adj) & p_adj < fdr_target)
  }

  # Reconstruct symmetric adjacency for compute_tpr_fpr
  delta_hat <- matrix(0L, p, p)
  delta_hat[ut] <- called
  delta_hat <- delta_hat + t(delta_hat)
  perf <- compute_tpr_fpr(gold$truth$Delta_A, delta_hat)

  list(
    diff_TPR_inference   = perf$TPR,
    diff_FPR_inference   = perf$FPR,
    diff_FDR_inference   = compute_fdr(perf$TP, perf$FP),
    diff_MCC_inference   = compute_mcc(perf$TP, perf$FP, perf$TN, perf$FN),
    diff_Power_inference = perf$TPR  # Power = TPR; alias for clarity
  )
}

#======================================================
# apply_bayesian_fdr_pip: direct-posterior FDR control
# from posterior inclusion probabilities (PIPs).
#
# Uses the "direct posterior probability" approach (Newton 2004,
# Müller 2004; Wakefield 2007 BFDP is the same thing with
# BFDP = 1 - PIP). For a Bayesian method that produces PIPs:
#   1. Rank edges by PIP descending.
#   2. Running mean of (1 - PIP) across the top-k edges is the
#      posterior expected FDR for that selection.
#   3. Largest k with running mean <= fdr_target gives the
#      call set; posterior E[FDR] <= target by construction.
# No permutation null, no refitting. See NOTES_bayesian_fdr.md.
#
# .bayes_fdr_call() does the selection; returns a logical vector
# of length n_edges (upper-tri) indicating which edges are
# declared. Used by apply_bayesian_fdr_pip() and reusable
# elsewhere if needed.
#======================================================

.bayes_fdr_call <- function(pip_vec, fdr_target) {
  ok <- !is.na(pip_vec)
  if (!any(ok)) return(rep(FALSE, length(pip_vec)))

  ord <- order(pip_vec, decreasing = TRUE, na.last = NA)
  running_fdr <- cumsum(1 - pip_vec[ord]) / seq_along(ord)

  # Largest k with running_fdr <= target. If no such k, call nothing.
  pass <- running_fdr <= fdr_target
  # Only call the contiguous "top block" that satisfies the rule.
  # i.e., find the largest k where running_fdr[1..k] all <= target.
  # (Equivalent to: stop as soon as running mean exceeds target.)
  if (!pass[1]) {
    k_star <- 0L
  } else {
    # running_fdr is monotonic under tied PIPs; we still compute
    # safely using cummax.
    k_star <- sum(cummax(running_fdr) <= fdr_target)
  }

  called <- logical(length(pip_vec))
  if (k_star > 0) called[ord[seq_len(k_star)]] <- TRUE
  called
}

# Helper: score a single call set (logical edge-vec) against truth
# and return TPR/FPR/FDR/MCC/Power.
.score_bayes_call <- function(truth_vec, called_vec, return_power = FALSE) {
  TP <- sum(truth_vec == 1 & called_vec)
  FP <- sum(truth_vec == 0 & called_vec)
  TN <- sum(truth_vec == 0 & !called_vec)
  FN <- sum(truth_vec == 1 & !called_vec)
  TPR <- if (TP + FN > 0) TP / (TP + FN) else NA_real_
  FPR <- if (FP + TN > 0) FP / (FP + TN) else NA_real_
  FDR <- compute_fdr(TP, FP)
  MCC <- compute_mcc(TP, FP, TN, FN)
  out <- list(TPR = TPR, FPR = FPR, FDR = FDR, MCC = MCC)
  if (return_power) out$Power <- TPR
  out
}

#------------------------------------------------------
# apply_bayesian_fdr_pip: compute the 17 inference columns
# (diff_* + g1_* + g2_* + group_*) from multiGGM PIPs directly,
# with NO permutation null required.
#
# Inputs:
#   pip_per_group   : p x p x 2 array (group-specific PIPs).
#                     pip_per_group[,,g] = P(edge in group g | data).
#                     Used for g1_* and g2_* inference columns.
#   pip_diff        : p x p matrix; pip_diff[i,j] = P(edge (i,j)
#                     differs between groups | data). Used for
#                     diff_* inference columns.
#   gold            : gold-standard list (with $adj$g1, $adj$g2,
#                     $truth$Delta_A).
#   fdr_target      : posterior expected FDR cap (default 0.05).
#
# Returns a named list of 17 numeric values ready to be spliced
# into the result row.
#------------------------------------------------------

apply_bayesian_fdr_pip <- function(pip_per_group, pip_diff, gold,
                                    fdr_target = 0.05) {
  na_result <- list(
    diff_TPR_inference   = NA_real_,
    diff_FPR_inference   = NA_real_,
    diff_FDR_inference   = NA_real_,
    diff_MCC_inference   = NA_real_,
    diff_Power_inference = NA_real_,
    g1_TPR_inference   = NA_real_,
    g1_FPR_inference   = NA_real_,
    g1_FDR_inference   = NA_real_,
    g1_MCC_inference   = NA_real_,
    g2_TPR_inference   = NA_real_,
    g2_FPR_inference   = NA_real_,
    g2_FDR_inference   = NA_real_,
    g2_MCC_inference   = NA_real_,
    group_TPR_inference = NA_real_,
    group_FPR_inference = NA_real_,
    group_FDR_inference = NA_real_,
    group_MCC_inference = NA_real_
  )

  if (is.null(gold) || is.null(gold$adj) || is.null(gold$truth)) {
    return(na_result)
  }

  # Safe NA-mean helper
  safe_mean <- function(a, b) {
    if (is.na(a) && is.na(b)) NA_real_
    else if (is.na(a)) b
    else if (is.na(b)) a
    else (a + b) / 2
  }

  out <- na_result

  # --- Differential inference ---
  if (!is.null(pip_diff) &&
      is.matrix(pip_diff) &&
      identical(dim(pip_diff), dim(gold$truth$Delta_A))) {
    ut_d <- upper.tri(pip_diff)
    called <- .bayes_fdr_call(pip_diff[ut_d], fdr_target)
    perf   <- .score_bayes_call(gold$truth$Delta_A[ut_d], called,
                                 return_power = TRUE)
    out$diff_TPR_inference   <- perf$TPR
    out$diff_FPR_inference   <- perf$FPR
    out$diff_FDR_inference   <- perf$FDR
    out$diff_MCC_inference   <- perf$MCC
    out$diff_Power_inference <- perf$Power
  }

  # --- Per-group inference (requires pip_per_group array) ---
  if (!is.null(pip_per_group) && length(dim(pip_per_group)) == 3 &&
      dim(pip_per_group)[3] >= 2) {
    # Group 1
    pip_g1 <- pip_per_group[,,1]
    if (identical(dim(pip_g1), dim(gold$adj$g1))) {
      ut1 <- upper.tri(pip_g1)
      called_g1 <- .bayes_fdr_call(pip_g1[ut1], fdr_target)
      perf_g1   <- .score_bayes_call(gold$adj$g1[ut1], called_g1)
      out$g1_TPR_inference <- perf_g1$TPR
      out$g1_FPR_inference <- perf_g1$FPR
      out$g1_FDR_inference <- perf_g1$FDR
      out$g1_MCC_inference <- perf_g1$MCC
    }
    # Group 2
    pip_g2 <- pip_per_group[,,2]
    if (identical(dim(pip_g2), dim(gold$adj$g2))) {
      ut2 <- upper.tri(pip_g2)
      called_g2 <- .bayes_fdr_call(pip_g2[ut2], fdr_target)
      perf_g2   <- .score_bayes_call(gold$adj$g2[ut2], called_g2)
      out$g2_TPR_inference <- perf_g2$TPR
      out$g2_FPR_inference <- perf_g2$FPR
      out$g2_FDR_inference <- perf_g2$FDR
      out$g2_MCC_inference <- perf_g2$MCC
    }
    # Group-level = safe mean of g1 and g2 (matches the convention
    # used for the primary group_* metrics, averaging rates).
    out$group_TPR_inference <- safe_mean(out$g1_TPR_inference, out$g2_TPR_inference)
    out$group_FPR_inference <- safe_mean(out$g1_FPR_inference, out$g2_FPR_inference)
    out$group_FDR_inference <- safe_mean(out$g1_FDR_inference, out$g2_FDR_inference)
    out$group_MCC_inference <- safe_mean(out$g1_MCC_inference, out$g2_MCC_inference)
  }

  out
}

#======================================================
# apply_stability_fdr: per-group FDR control via
# stability-selection frequencies.
#
# For frequentist methods without a native uncertainty output
# (JGL, SSJGL), we estimate per-edge selection frequencies by
# subsampling + refitting (see 19_stability_selection.R). Those
# frequencies are interpretable as selection probabilities in
# [0, 1] — calibrated to the method's sampling behavior — and
# plug into the SAME direct-posterior procedure we use for
# multiGGM PIPs.
#
# This is the frequentist analog of apply_bayesian_fdr_pip():
# same .bayes_fdr_call() math, different source of the
# "inclusion" quantity.
#
# NOTE on approximation: the cached stability frequencies are
# computed ONCE per (method, setting, config) on a representative
# fixed-GS replicate, analogous to how the permutation null is
# cached. Each rep's per-group inference reuses those frequencies.
# This assumes the method's stability pattern is roughly constant
# across reps from the same DGP — an approximation we accept
# for compute feasibility. See NOTES_bayesian_fdr.md for details.
#
# Inputs:
#   freq_per_group : p x p x 2 array. freq_per_group[,,g] is the
#                    per-edge selection frequency for group g:
#                    P(edge detected by the method on a random
#                    subsample of the data).
#   gold           : gold-standard list (uses $adj$g1, $adj$g2).
#   fdr_target     : target posterior expected FDR (default 0.05).
#
# Returns a named list with 12 components (g1_*_inference,
# g2_*_inference, group_*_inference — 4 metrics each). The
# differential inference columns are NOT populated here;
# differential inference for frequentist methods still comes
# from the cached permutation null via apply_cached_inference().
#======================================================

#======================================================
# apply_pooled_pergroup_inference: pooled per-group permutation
# null + Benjamini-Hochberg for FDR-controlled per-group edge
# inference. Used by JGL family + SSJGL family (replaces the
# stability-selection path that only worked under fixed-GS mode).
#
# Inputs:
#   method_name        : string label (e.g. "JGL_BIC")
#   setting            : setting label (used to find the cache)
#   observed_score_g1  : p x p matrix of |theta_g1| at the rep being scored
#   observed_score_g2  : p x p matrix of |theta_g2|
#   gold               : gold-standard list (uses $adj$g1, $adj$g2)
#   fdr_target         : BH FDR target (default 0.05)
#   null_dir           : directory holding the per-group null caches
#                         (default "results/null_distributions_pergroup")
#   correction         : "pooled" (default), "bh", or "none"; mirrors
#                         apply_cached_inference's options.
#
# Returns a named list with the 12 per-group inference columns
# (g1_*_inference, g2_*_inference, group_*_inference x 4 metrics).
# Returns NA values when the cache file does not exist or fails
# to load. Cache file produced by 19_stability_selection.R.
#======================================================

apply_pooled_pergroup_inference <- function(method_name, setting,
                                             observed_score_g1,
                                             observed_score_g2,
                                             gold, fdr_target = 0.05,
                                             null_dir = "results/null_distributions_pergroup",
                                             correction = c("pooled", "bh", "none")) {
  correction <- match.arg(correction)
  na_result <- list(
    g1_TPR_inference = NA_real_, g1_FPR_inference = NA_real_,
    g1_FDR_inference = NA_real_, g1_MCC_inference = NA_real_,
    g2_TPR_inference = NA_real_, g2_FPR_inference = NA_real_,
    g2_FDR_inference = NA_real_, g2_MCC_inference = NA_real_,
    group_TPR_inference = NA_real_, group_FPR_inference = NA_real_,
    group_FDR_inference = NA_real_, group_MCC_inference = NA_real_
  )
  if (is.null(setting) || is.null(observed_score_g1) || is.null(observed_score_g2)) {
    return(na_result)
  }
  null_file <- file.path(null_dir, sprintf("%s_%s.rds", setting, method_name))
  if (!file.exists(null_file)) return(na_result)
  null_rec <- tryCatch(readRDS(null_file), error = function(e) NULL)
  if (is.null(null_rec) ||
      is.null(null_rec$null_scores_g1) ||
      is.null(null_rec$null_scores_g2)) {
    return(na_result)
  }

  safe_mean <- function(a, b) {
    if (is.na(a) && is.na(b)) NA_real_
    else if (is.na(a)) b
    else if (is.na(b)) a
    else (a + b) / 2
  }

  # ----- Inner helper: pvals + BH on one group -----
  compute_called_one_group <- function(observed, null_list) {
    p <- nrow(observed)
    ut <- upper.tri(observed)
    obs_vec <- observed[ut]
    valid <- vapply(null_list,
                     function(m) !is.null(m) && !all(is.na(m)),
                     logical(1))
    null_list <- null_list[valid]
    if (length(null_list) == 0) {
      return(matrix(0L, p, p))
    }
    if (correction == "pooled") {
      null_pool <- unlist(lapply(null_list, function(m) m[ut]),
                           use.names = FALSE)
      null_pool <- null_pool[is.finite(null_pool)]
      asc_pool  <- sort(null_pool)
      n_lt   <- findInterval(obs_vec - .Machine$double.eps, asc_pool)
      exceed <- length(asc_pool) - n_lt
      pvals  <- (1 + exceed) / (1 + length(asc_pool))
    } else {
      null_mat   <- do.call(rbind, lapply(null_list, function(m) m[ut]))
      n_perm_eff <- nrow(null_mat)
      obs_row_rep <- matrix(obs_vec, n_perm_eff, length(obs_vec),
                             byrow = TRUE)
      exceed <- colSums(null_mat >= obs_row_rep, na.rm = TRUE)
      pvals  <- (1 + exceed) / (1 + n_perm_eff)
    }
    if (correction == "none") {
      called <- as.integer(!is.na(pvals) & pvals < fdr_target)
    } else {
      p_adj  <- p.adjust(pvals, method = "BH")
      called <- as.integer(!is.na(p_adj) & p_adj < fdr_target)
    }
    delta_hat <- matrix(0L, p, p)
    delta_hat[ut] <- called
    delta_hat <- delta_hat + t(delta_hat)
    delta_hat
  }

  delta_hat_g1 <- tryCatch(
    compute_called_one_group(observed_score_g1, null_rec$null_scores_g1),
    error = function(e) NULL
  )
  delta_hat_g2 <- tryCatch(
    compute_called_one_group(observed_score_g2, null_rec$null_scores_g2),
    error = function(e) NULL
  )
  if (is.null(delta_hat_g1) || is.null(delta_hat_g2)) return(na_result)

  perf_g1 <- compute_tpr_fpr(gold$adj$g1, delta_hat_g1)
  perf_g2 <- compute_tpr_fpr(gold$adj$g2, delta_hat_g2)
  fdr_g1  <- compute_fdr(perf_g1$TP, perf_g1$FP)
  fdr_g2  <- compute_fdr(perf_g2$TP, perf_g2$FP)
  mcc_g1  <- compute_mcc(perf_g1$TP, perf_g1$FP, perf_g1$TN, perf_g1$FN)
  mcc_g2  <- compute_mcc(perf_g2$TP, perf_g2$FP, perf_g2$TN, perf_g2$FN)

  list(
    g1_TPR_inference = perf_g1$TPR, g1_FPR_inference = perf_g1$FPR,
    g1_FDR_inference = fdr_g1,      g1_MCC_inference = mcc_g1,
    g2_TPR_inference = perf_g2$TPR, g2_FPR_inference = perf_g2$FPR,
    g2_FDR_inference = fdr_g2,      g2_MCC_inference = mcc_g2,
    group_TPR_inference = safe_mean(perf_g1$TPR, perf_g2$TPR),
    group_FPR_inference = safe_mean(perf_g1$FPR, perf_g2$FPR),
    group_FDR_inference = safe_mean(fdr_g1, fdr_g2),
    group_MCC_inference = safe_mean(mcc_g1, mcc_g2)
  )
}

#======================================================
# apply_stability_fdr (DEPRECATED).
#
# Bayesian FDR on subsample-stability frequencies for per-group
# inference. Replaced by apply_pooled_pergroup_inference() because
# stability frequencies are network-specific and only valid under
# fixed-GS mode. Kept here for backward compatibility with any old
# results/stability_frequencies/ caches; new pipelines use the
# pooled-permutation null path.
#======================================================

apply_stability_fdr <- function(freq_per_group, gold, fdr_target = 0.05) {
  na_result <- list(
    g1_TPR_inference = NA_real_, g1_FPR_inference = NA_real_,
    g1_FDR_inference = NA_real_, g1_MCC_inference = NA_real_,
    g2_TPR_inference = NA_real_, g2_FPR_inference = NA_real_,
    g2_FDR_inference = NA_real_, g2_MCC_inference = NA_real_,
    group_TPR_inference = NA_real_, group_FPR_inference = NA_real_,
    group_FDR_inference = NA_real_, group_MCC_inference = NA_real_
  )

  if (is.null(gold) || is.null(gold$adj) ||
      is.null(freq_per_group) ||
      length(dim(freq_per_group)) != 3 ||
      dim(freq_per_group)[3] < 2) {
    return(na_result)
  }

  safe_mean <- function(a, b) {
    if (is.na(a) && is.na(b)) NA_real_
    else if (is.na(a)) b
    else if (is.na(b)) a
    else (a + b) / 2
  }

  out <- na_result

  # Group 1
  freq_g1 <- freq_per_group[,,1]
  if (identical(dim(freq_g1), dim(gold$adj$g1))) {
    ut1 <- upper.tri(freq_g1)
    called_g1 <- .bayes_fdr_call(freq_g1[ut1], fdr_target)
    perf_g1   <- .score_bayes_call(gold$adj$g1[ut1], called_g1)
    out$g1_TPR_inference <- perf_g1$TPR
    out$g1_FPR_inference <- perf_g1$FPR
    out$g1_FDR_inference <- perf_g1$FDR
    out$g1_MCC_inference <- perf_g1$MCC
  }
  # Group 2
  freq_g2 <- freq_per_group[,,2]
  if (identical(dim(freq_g2), dim(gold$adj$g2))) {
    ut2 <- upper.tri(freq_g2)
    called_g2 <- .bayes_fdr_call(freq_g2[ut2], fdr_target)
    perf_g2   <- .score_bayes_call(gold$adj$g2[ut2], called_g2)
    out$g2_TPR_inference <- perf_g2$TPR
    out$g2_FPR_inference <- perf_g2$FPR
    out$g2_FDR_inference <- perf_g2$FDR
    out$g2_MCC_inference <- perf_g2$MCC
  }
  # Group-level = safe mean (same convention as the primary group_*)
  out$group_TPR_inference <- safe_mean(out$g1_TPR_inference, out$g2_TPR_inference)
  out$group_FPR_inference <- safe_mean(out$g1_FPR_inference, out$g2_FPR_inference)
  out$group_FDR_inference <- safe_mean(out$g1_FDR_inference, out$g2_FDR_inference)
  out$group_MCC_inference <- safe_mean(out$g1_MCC_inference, out$g2_MCC_inference)

  out
}

#======================================================
# calc_metrics_from_pip: multiGGM PIP-based metrics
# (Retained as an auxiliary function for diagnostic use.
#  multiGGM's main analysis now uses calc_metrics on the
#  posterior-mean precision matrices; call this function
#  explicitly if you want PIP-based metrics.)
#======================================================

calc_metrics_from_pip <- function(pip_array, pip_diff_matrix, gold,
                                   pip_threshold = 0.5,
                                   theta_list = NULL) {
  ut <- upper.tri(pip_array[,,1])

  adj_hat_g1 <- (pip_array[,,1] > pip_threshold) * 1L; diag(adj_hat_g1) <- 0L
  adj_hat_g2 <- (pip_array[,,2] > pip_threshold) * 1L; diag(adj_hat_g2) <- 0L
  perf_g1 <- compute_tpr_fpr(gold$adj$g1, adj_hat_g1)
  perf_g2 <- compute_tpr_fpr(gold$adj$g2, adj_hat_g2)

  # Primary differential = pip_diff threshold
  delta_hat_primary <- (pip_diff_matrix > pip_threshold) * 1L
  diag(delta_hat_primary) <- 0L
  perf_diff <- compute_tpr_fpr(gold$truth$Delta_A, delta_hat_primary)
  # Secondary XOR
  delta_hat_xor <- (adj_hat_g1 != adj_hat_g2) * 1L; diag(delta_hat_xor) <- 0L
  perf_diff_xor <- compute_tpr_fpr(.xor_truth(gold), delta_hat_xor)

  if (!is.null(theta_list)) {
    g1_frob_rel <- compute_relative_frobenius(theta_list[[1]], gold$omega$g1)
    g2_frob_rel <- compute_relative_frobenius(theta_list[[2]], gold$omega$g2)
    diff_est  <- theta_list[[1]] - theta_list[[2]]
    diff_true <- gold$omega$g1 - gold$omega$g2
    diff_frob_rel <- compute_relative_frobenius(diff_est, diff_true)
    g1_KL <- tryCatch(compute_kl_divergence(theta_list[[1]], gold$sigma$g1),
                      error = function(e) NA_real_)
    g2_KL <- tryCatch(compute_kl_divergence(theta_list[[2]], gold$sigma$g2),
                      error = function(e) NA_real_)
  } else {
    g1_frob_rel <- NA_real_; g2_frob_rel <- NA_real_; diff_frob_rel <- NA_real_
    g1_KL <- NA_real_; g2_KL <- NA_real_
  }

  g1_AUC   <- compute_auc(gold$adj$g1[ut], pip_array[,,1][ut])
  g2_AUC   <- compute_auc(gold$adj$g2[ut], pip_array[,,2][ut])
  diff_AUC <- compute_auc(gold$truth$Delta_A[ut], pip_diff_matrix[ut])
  g1_AUPR   <- compute_aupr(gold$adj$g1[ut], pip_array[,,1][ut])
  g2_AUPR   <- compute_aupr(gold$adj$g2[ut], pip_array[,,2][ut])
  diff_AUPR <- compute_aupr(gold$truth$Delta_A[ut], pip_diff_matrix[ut])

  .build_result_row(
    perf_g1, perf_g2, perf_diff, perf_diff_xor,
    g1_AUC, g2_AUC, diff_AUC,
    g1_AUPR, g2_AUPR, diff_AUPR,
    g1_frob_rel, g2_frob_rel, diff_frob_rel,
    g1_KL, g2_KL,
    g1_edges_true = sum(gold$adj$g1[ut]),
    g2_edges_true = sum(gold$adj$g2[ut]),
    diff_edges_true = sum(gold$truth$Delta_A[ut]),
    g1_edges_est = sum(adj_hat_g1[ut]),
    g2_edges_est = sum(adj_hat_g2[ut]),
    diff_edges_est = sum(delta_hat_primary[ut])
  )
}
