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
#' Anything else (other designs, other unequal group sizes, FDR targets
#' other than the grid's 0.05) is flagged `simulate_required`, and answered
#' by an on-demand simulation when `simulate = TRUE` (see [simulate_power()]).
#'
#' @param p Number of features.
#' @param n Per-group sample size; a single number, or `c(n1, n2)` for
#'   unequal groups. Group 1 is the reference group (its network is intact;
#'   group 2 carries the change), so order matters.
#' @param method Estimation method label. See [howbig_methods()].
#' @param perturbation Type of differential signal: `"structural"` (edges
#'   appear/disappear), `"magnitude"` (edge strengths change), or `"both"`.
#' @param fdr Target false-discovery rate. The grid is scored at `0.05`;
#'   other targets require simulation.
#' @param scenario Assumption about the networks and data, for checking how
#'   sensitive an answer is to them. `"default"` (block-structured
#'   Barabasi-Albert tree networks, Gaussian data, empirical edge weights) is
#'   the main grid. The others are in the grid only at p = 100, n = 1000 with
#'   a structural change: `"scale_free"` (denser Barabasi-Albert blocks),
#'   `"erdos_renyi"` (random blocks), `"ar1_banded"` (non-modular banded
#'   network), `"censored_5"`, `"censored_25"`, `"censored_50"` (that share
#'   of each feature's values set to zero, as at a detection limit), and
#'   `"uniform_weights"`.
#' @param simulate If `TRUE`, answer off-grid designs by on-demand
#'   simulation (slow; default scenario and simulatable methods only). If
#'   `FALSE` (default), they are flagged instead.
#' @param ... Passed to [simulate_power()] when simulating (e.g. `reps`,
#'   `n_perm`, `cores`, `seed`).
#'
#' @return An object of class `howbig_power`: a list with the sensitivity
#'   estimate and 80% band, achieved FDR, AUPR, the number of truly
#'   differential edges, provenance, and any notes.
#' @examples
#' power_for(p = 250, n = 1000, method = "JGL_BIC")
#' power_for(p = 250, n = 1500, method = "JGL_BIC")          # interpolated
#' power_for(p = 100, n = c(500, 1500), method = "JGL_BIC")  # unequal groups
#' power_for(p = 100, n = 1000, method = "JGL_BIC", scenario = "censored_25")
#' power_for(p = 150, n = 1500, method = "JGL_BIC")          # off-grid: flagged
#' @export
power_for <- function(p, n, method,
                      perturbation = "structural",
                      fdr = .HOWBIG_FDR_TARGET,
                      scenario = "default",
                      simulate = FALSE, ...) {
  stopifnot(length(p) == 1L, is.numeric(p), p > 0)
  n <- .normalize_n(n)
  method <- .normalize_method(method)
  perturbation <- .normalize_perturbation(perturbation)
  scenario <- .normalize_scenario(scenario)

  cls <- if (!isTRUE(all.equal(fdr, .HOWBIG_FDR_TARGET))) {
    list(kind = "simulate_required",
         reason = sprintf("FDR target %.2f (grid is scored at %.2f)", fdr,
                          .HOWBIG_FDR_TARGET))
  } else {
    .classify(method, p, n, perturbation, scenario)
  }

  if (cls$kind == "simulate_required" && isTRUE(simulate)) {
    if (scenario != "default") {
      stop("Simulation supports only scenario = \"default\"; scenario '",
           scenario, "' is available from the grid at p = 100, n = 1000.",
           call. = FALSE)
    }
    return(simulate_power(p = p, n = n, method = method,
                          perturbation = perturbation, fdr = fdr, ...))
  }

  vals <- switch(cls$kind,
    in_grid = as.list(.grid_cell(method, p, n, perturbation, scenario)),
    interpolated = .interp_cell(method, perturbation, cls$axis,
                                at = if (cls$axis == "n") n[1] else p),
    NULL)

  res <- list(
    method = method, p = p, n = n, perturbation = perturbation,
    scenario = scenario, fdr_target = fdr, provenance = cls$kind,
    sensitivity = NA_real_, sensitivity_lo = NA_real_, sensitivity_hi = NA_real_,
    achieved_fdr = NA_real_, aupr = NA_real_,
    n_diff_edges = NA_real_, n_pairs = p * (p - 1) / 2,
    n_reps = NA_integer_, notes = character(0)
  )

  if (is.null(vals)) {
    can_sim <- scenario == "default" &&
      howbig_method_info$simulate[howbig_method_info$method == method]
    res$notes <- c(paste0("not answerable from the grid: ", cls$reason, "."),
      if (can_sim) "rerun with simulate = TRUE to simulate this design (slow; see ?simulate_power)."
      else if (scenario != "default") "simulation supports only the default scenario."
      else sprintf("%s is available from the grid only; it cannot be simulated.", method))
    return(structure(res, class = "howbig_power"))
  }

  res$sensitivity    <- vals$sensitivity_mean
  res$sensitivity_lo <- vals$sensitivity_q10
  res$sensitivity_hi <- vals$sensitivity_q90
  res$achieved_fdr   <- vals$achieved_fdr_mean
  res$aupr           <- vals$aupr_mean
  res$n_diff_edges   <- vals$n_diff_edges
  res$n_reps         <- vals$n_reps

  if (is.finite(res$n_reps) && res$n_reps < 20) {
    res$notes <- c(res$notes, sprintf(
      "this grid cell rests on only %d replicates; treat it as rough.",
      res$n_reps))
  }
  # Metric hole: a cell can have AUPR but no scored inference.
  if (!is.finite(res$sensitivity) && is.finite(res$aupr)) {
    res$notes <- c(res$notes, paste0(
      "sensitivity was not scored at this cell; AUPR is shown as a ",
      "threshold-free fallback."))
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
  scen <- if (!is.null(x$scenario) && x$scenario != "default")
    sprintf(" | scenario: %s", x$scenario) else ""
  cat(sprintf("  %s | %s perturbation | p = %g, %s%s\n",
              x$method, x$perturbation, x$p, n_txt, scen))
  if (x$provenance == "simulate_required") {
    cat("  Not available from the grid.\n")
  } else {
    if (is.finite(x$n_diff_edges %||% NA_real_)) {
      cat(sprintf("  Truly differential edges: ~%s of %s pairs (%.2f%%)\n",
                  format(round(x$n_diff_edges), big.mark = ","),
                  format(x$n_pairs, big.mark = ","),
                  100 * x$n_diff_edges / x$n_pairs))
    }
    if (is.finite(x$sensitivity)) {
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
