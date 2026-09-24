#' Power (sensitivity) for a given design
#'
#' Estimates the power of a differential network analysis at a given sample
#' size and feature count, where **power is the sensitivity** -- the fraction
#' of truly differential edges recovered -- at a false-discovery-rate target.
#' The achieved FDR is reported alongside, because several methods are not
#' calibrated in every regime.
#'
#' Answers come from a precomputed simulation grid where possible. Queries
#' that hit a grid cell are exact; queries on a swept arm are interpolated.
#' Anything else (other designs, unequal group sizes, FDR targets other than
#' the grid's 0.05) is flagged `simulate_required`, and answered by an
#' on-demand simulation when `simulate = TRUE` (see [simulate_power()]).
#'
#' @param p Number of features.
#' @param n Per-group sample size; a single number, or `c(n1, n2)` for
#'   unequal groups (always simulated).
#' @param method Estimation method label. See [howbig_coverage()] for the
#'   available methods.
#' @param perturbation Type of differential signal: `"structural"` (edges
#'   appear/disappear), `"magnitude"` (edge weights change), or `"both"`.
#' @param fdr Target false-discovery rate. The grid is scored at `0.05`;
#'   other targets require simulation.
#' @param simulate If `TRUE`, answer off-grid designs by on-demand
#'   simulation (slow). If `FALSE` (default), they are flagged instead.
#' @param ... Passed to [simulate_power()] when simulating (e.g. `reps`,
#'   `n_perm`, `cores`, `seed`).
#'
#' @return An object of class `howbig_power`: a list with the sensitivity
#'   estimate and 80% band, achieved FDR, AUPR, provenance, and any notes.
#' @examples
#' power_for(p = 250, n = 1000, method = "JGL_BIC")
#' power_for(p = 250, n = 1500, method = "JGL_BIC")   # interpolated
#' power_for(p = 150, n = 1500, method = "JGL_BIC")   # off-grid: flagged
#' @export
power_for <- function(p, n, method,
                      perturbation = "structural",
                      fdr = .HOWBIG_FDR_TARGET,
                      simulate = FALSE, ...) {
  stopifnot(length(p) == 1L, length(n) %in% 1:2, is.numeric(p), is.numeric(n),
            p > 0, all(n > 0))
  method <- .normalize_method(method)
  perturbation <- .normalize_perturbation(perturbation)

  # Reasons the grid cannot answer, independent of where (p, n) falls.
  off_reason <- c(
    if (length(n) == 2L && n[1] != n[2]) "unequal group sizes",
    if (!isTRUE(all.equal(fdr, .HOWBIG_FDR_TARGET)))
      sprintf("FDR target %.2f (grid is scored at %.2f)", fdr, .HOWBIG_FDR_TARGET)
  )
  n_groups <- n                      # keep c(n1, n2) for the simulator
  n <- n[1]
  cls <- if (length(off_reason)) list(kind = "simulate_required")
         else .classify(method, p, n, perturbation)
  if (cls$kind == "simulate_required" && !length(off_reason)) {
    off_reason <- sprintf("(p = %g, n = %g) is off the grid for the %s perturbation",
                          p, n, perturbation)
  }

  if (cls$kind == "simulate_required" && isTRUE(simulate)) {
    return(simulate_power(p = p, n = n_groups, method = method,
                          perturbation = perturbation, fdr = fdr, ...))
  }
  notes <- character(0)

  vals <- if (cls$kind == "in_grid") {
    as.list(.grid_cell(method, p, n, perturbation))
  } else if (cls$kind == "interpolated") {
    .interp_cell(method, perturbation, cls$axis, at = if (cls$axis == "n") n else p)
  } else {
    NULL
  }

  res <- list(
    method = method, p = p, n = n_groups, perturbation = perturbation,
    fdr_target = fdr, provenance = cls$kind,
    sensitivity = NA_real_, sensitivity_lo = NA_real_, sensitivity_hi = NA_real_,
    achieved_fdr = NA_real_, aupr = NA_real_,
    n_reps = NA_integer_, notes = notes
  )

  if (is.null(vals)) {
    res$notes <- c(res$notes,
      paste0("not answerable from the grid: ", paste(off_reason, collapse = "; "), "."),
      "rerun with simulate = TRUE to simulate this design (slow; see ?simulate_power).")
    return(structure(res, class = "howbig_power"))
  }

  res$sensitivity    <- vals$sensitivity_mean
  res$sensitivity_lo <- vals$sensitivity_q10
  res$sensitivity_hi <- vals$sensitivity_q90
  res$achieved_fdr   <- vals$achieved_fdr_mean
  res$aupr           <- vals$aupr_mean
  res$n_reps         <- vals$n_reps

  # Metric hole: some cells (e.g. p = 500) have AUPR but no scored inference.
  if (!is.finite(res$sensitivity) && is.finite(res$aupr)) {
    res$notes <- c(res$notes, paste0(
      "sensitivity was not scored at this cell (no inference columns); ",
      "AUPR is shown as a threshold-free fallback. Use simulate_power() ",
      "to estimate sensitivity here."))
  }
  if (.miscalibrated(res$achieved_fdr, res$fdr_target)) {
    res$notes <- c(res$notes, sprintf(
      "achieved FDR %.2f exceeds target %.2f \u2014 this method is not calibrated here.",
      res$achieved_fdr, res$fdr_target))
  }
  structure(res, class = "howbig_power")
}

#' @export
print.howbig_power <- function(x, ...) {
  tag <- switch(x$provenance,
                in_grid = "in-grid", interpolated = "interpolated",
                simulated = "simulated", simulate_required = "OFF-GRID")
  cat(sprintf("howbig power estimate  [%s]\n", tag))
  n_txt <- if (length(x$n) == 2L && x$n[1] != x$n[2])
    sprintf("n = %g / %g", x$n[1], x$n[2]) else sprintf("n = %g", x$n[1])
  cat(sprintf("  %s | %s perturbation | p = %g, %s\n",
              x$method, x$perturbation, x$p, n_txt))
  if (x$provenance == "simulate_required") {
    cat("  Not available from the grid.\n")
  } else if (is.finite(x$sensitivity)) {
    cat(sprintf("  Sensitivity @ FDR %.2f:  %.2f   (80%% band %.2f\u2013%.2f)\n",
                x$fdr_target, x$sensitivity, x$sensitivity_lo, x$sensitivity_hi))
    if (is.finite(x$sensitivity_mcse %||% NA_real_))
      cat(sprintf("  Monte Carlo SE of mean:  %.3f\n", x$sensitivity_mcse))
    fdr_flag <- if (.miscalibrated(x$achieved_fdr, x$fdr_target)) " \u26a0" else " \u2713"
    cat(sprintf("  Achieved FDR:            %.3f%s\n", x$achieved_fdr, fdr_flag))
    cat(sprintf("  AUPR (threshold-free):   %.2f\n", x$aupr))
  } else if (is.finite(x$aupr)) {
    cat(sprintf("  AUPR (threshold-free):   %.2f   (sensitivity not scored here)\n",
                x$aupr))
  }
  if (is.finite(x$n_reps)) {
    extra <- if (x$provenance == "simulated")
      sprintf(", %s; %.1f min", if (is.na(x$n_perm)) "posterior FDR"
              else sprintf("%d-permutation null", x$n_perm), x$runtime_sec / 60)
      else ""
    cat(sprintf("  based on %d reps%s\n", x$n_reps, extra))
  }
  for (nt in x$notes) cat(strwrap(paste("  note:", nt), exdent = 8), sep = "\n")
  invisible(x)
}
