#' Power (sensitivity) for a given design
#'
#' Estimates the power of a differential network analysis at a given sample
#' size and feature count, where **power is the sensitivity** — the fraction
#' of truly differential edges recovered — at a false-discovery-rate target.
#' The achieved FDR is reported alongside, because several methods are not
#' calibrated in every regime.
#'
#' Answers come from a precomputed simulation grid. Queries that hit a grid
#' cell are exact; queries that fall on a swept arm are interpolated; queries
#' off the grid are flagged `simulate_required` (the simulation fallback is
#' not yet wired in this version).
#'
#' @param p Number of features.
#' @param n Per-group sample size.
#' @param method Estimation method label. See [howbig_coverage()] for the
#'   available methods.
#' @param perturbation Type of differential signal: `"structural"` (edges
#'   appear/disappear), `"magnitude"` (edge weights change), or `"both"`.
#' @param fdr Target false-discovery rate. Only `0.05` is available from the
#'   grid in this version; other values warn and use `0.05`.
#'
#' @return An object of class `howbig_power`: a list with the sensitivity
#'   estimate and 80% band, achieved FDR, AUPR, provenance, and any notes.
#' @examples
#' power_for(p = 250, n = 1000, method = "JGL_BIC")
#' power_for(p = 250, n = 1500, method = "JGL_BIC")   # interpolated
#' @export
power_for <- function(p, n, method,
                      perturbation = "structural",
                      fdr = .HOWBIG_FDR_TARGET) {
  stopifnot(length(p) == 1L, length(n) == 1L, is.numeric(p), is.numeric(n),
            p > 0, n > 0)
  method <- .normalize_method(method)
  perturbation <- .normalize_perturbation(perturbation)

  notes <- character(0)
  if (!isTRUE(all.equal(fdr, .HOWBIG_FDR_TARGET))) {
    notes <- c(notes, sprintf(
      "grid is scored at FDR = %.2f; requested %.2f ignored (needs simulation).",
      .HOWBIG_FDR_TARGET, fdr))
  }

  cls <- .classify(method, p, n, perturbation)

  vals <- if (cls$kind == "in_grid") {
    as.list(.grid_cell(method, p, n, perturbation))
  } else if (cls$kind == "interpolated") {
    .interp_cell(method, perturbation, cls$axis, at = if (cls$axis == "n") n else p)
  } else {
    NULL
  }

  res <- list(
    method = method, p = p, n = n, perturbation = perturbation,
    fdr_target = .HOWBIG_FDR_TARGET, provenance = cls$kind,
    sensitivity = NA_real_, sensitivity_lo = NA_real_, sensitivity_hi = NA_real_,
    achieved_fdr = NA_real_, aupr = NA_real_,
    n_reps = NA_integer_, notes = notes
  )

  if (is.null(vals)) {
    res$notes <- c(res$notes, paste0(
      "off-grid for perturbation '", perturbation,
      "': the simulation fallback is not available in this version."))
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
      "AUPR is shown as a threshold-free fallback."))
  }
  if (.miscalibrated(res$achieved_fdr, res$fdr_target)) {
    res$notes <- c(res$notes, sprintf(
      "achieved FDR %.2f exceeds target %.2f — this method is not calibrated here.",
      res$achieved_fdr, res$fdr_target))
  }
  structure(res, class = "howbig_power")
}

#' @export
print.howbig_power <- function(x, ...) {
  tag <- switch(x$provenance,
                in_grid = "in-grid", interpolated = "interpolated",
                simulate_required = "OFF-GRID")
  cat(sprintf("howbig power estimate  [%s]\n", tag))
  cat(sprintf("  %s | %s perturbation | p = %g, n = %g\n",
              x$method, x$perturbation, x$p, x$n))
  if (x$provenance == "simulate_required") {
    cat("  Not available from the grid.\n")
  } else if (is.finite(x$sensitivity)) {
    cat(sprintf("  Sensitivity @ FDR %.2f:  %.2f   (80%% band %.2f–%.2f)\n",
                x$fdr_target, x$sensitivity, x$sensitivity_lo, x$sensitivity_hi))
    fdr_flag <- if (.miscalibrated(x$achieved_fdr, x$fdr_target)) " ⚠" else " ✓"
    cat(sprintf("  Achieved FDR:            %.3f%s\n", x$achieved_fdr, fdr_flag))
    cat(sprintf("  AUPR (threshold-free):   %.2f\n", x$aupr))
  } else if (is.finite(x$aupr)) {
    cat(sprintf("  AUPR (threshold-free):   %.2f   (sensitivity not scored here)\n",
                x$aupr))
  }
  if (is.finite(x$n_reps)) cat(sprintf("  based on %d reps\n", x$n_reps))
  for (nt in x$notes) cat(strwrap(paste("  note:", nt), exdent = 8), sep = "\n")
  invisible(x)
}
