#' Required sample size for a target power
#'
#' Inverts [power_for()]: finds the smallest per-group sample size `n` that
#' reaches a target power (sensitivity) at a feature count `p`.
#'
#' From the grid, the inversion uses the precomputed sample-size arm, which
#' exists at `p = 250` for the `structural` perturbation at FDR 0.05. Other
#' designs are flagged `simulate_required`; with `simulate = TRUE` they are
#' answered by simulating power at `n_points` log-spaced sample sizes across
#' `n_range`, smoothing the curve to be monotone, and inverting it.
#'
#' @param p Number of features.
#' @param power Target sensitivity in (0, 1).
#' @param method Estimation method label.
#' @param perturbation Type of differential signal; see [power_for()].
#' @param fdr Target false-discovery rate. The grid is scored at `0.05`;
#'   other targets require simulation.
#' @param simulate If `TRUE`, answer off-grid designs by simulation (slow:
#'   `n_points` full [simulate_power()] runs).
#' @param n_range,n_points Per-group sample sizes searched when simulating.
#' @param ... Passed to [simulate_power()] (e.g. `reps`, `n_perm`, `cores`).
#'
#' @return An object of class `howbig_n`: the required `n`, provenance, the
#'   `n` range searched, the power curve used, and any notes.
#' @examples
#' n_for(p = 250, power = 0.8, method = "JGL_BIC")
#' @export
n_for <- function(p, power = 0.8, method,
                  perturbation = "structural",
                  fdr = .HOWBIG_FDR_TARGET,
                  simulate = FALSE, n_range = c(100, 5000), n_points = 5L,
                  ...) {
  stopifnot(length(p) == 1L, is.numeric(p), p > 0,
            length(power) == 1L, power > 0, power < 1)
  method <- .normalize_method(method)
  perturbation <- .normalize_perturbation(perturbation)

  res <- list(method = method, p = p, power = power,
              perturbation = perturbation, fdr_target = fdr,
              provenance = "simulate_required", n = NA_real_,
              n_range = c(NA_real_, NA_real_), curve = NULL,
              notes = character(0))

  on_arm <- p == 250 && perturbation == "structural" &&
    isTRUE(all.equal(fdr, .HOWBIG_FDR_TARGET))
  if (!on_arm) {
    if (isTRUE(simulate)) {
      return(.n_for_simulated(res, n_range, n_points, ...))
    }
    res$notes <- c(res$notes, paste0(
      "the grid's sample-size arm covers only p = 250, structural, FDR 0.05; ",
      "rerun with simulate = TRUE to search n by simulation (slow)."))
    return(structure(res, class = "howbig_n"))
  }

  arm <- .arm_table(method, perturbation, "n")
  ok <- is.finite(arm$sensitivity_mean)
  arm <- arm[ok, ]
  if (nrow(arm) < 2L) {
    res$notes <- c(res$notes, "too few scored points on the arm to invert.")
    return(structure(res, class = "howbig_n"))
  }
  arm <- arm[order(arm$n), ]
  res$n_range <- range(arm$n)
  s <- arm$sensitivity_mean

  if (power <= min(s)) {
    res$n <- min(arm$n); res$provenance <- "interpolated"
    res$notes <- c(res$notes, sprintf(
      "target power is reached at or below the smallest gridded n (%g).",
      min(arm$n)))
    return(structure(res, class = "howbig_n"))
  }
  if (power > max(s)) {
    res$notes <- c(res$notes, sprintf(
      "target power %.2f exceeds the best gridded sensitivity (%.2f at n = %g); ",
      power, max(s), max(arm$n)))
    res$notes <- c(res$notes, "a larger n is needed than the grid covers (simulation required).")
    return(structure(res, class = "howbig_n"))
  }

  f <- stats::splinefun(log(arm$n), s, method = "monoH.FC")
  root <- stats::uniroot(function(ln) f(ln) - power,
                         interval = log(range(arm$n)))$root
  res$n <- ceiling(exp(root))
  res$provenance <- if (res$n %in% arm$n) "in_grid" else "interpolated"
  res$curve <- data.frame(n = arm$n, sensitivity = s)
  structure(res, class = "howbig_n")
}

# Simulated search: power at log-spaced n, isotonic smoothing (simulated
# points are noisy and power must be non-decreasing in n), linear inversion
# on log(n).
.n_for_simulated <- function(res, n_range, n_points, ...) {
  stopifnot(length(n_range) == 2L, n_range[1] < n_range[2], n_points >= 3)
  n_grid <- unique(round(exp(seq(log(n_range[1]), log(n_range[2]),
                                 length.out = n_points))))
  sens <- vapply(n_grid, function(nn) {
    message(sprintf("howbig: n_for simulating n = %g ...", nn))
    simulate_power(p = res$p, n = nn, method = res$method,
                   perturbation = res$perturbation, fdr = res$fdr_target,
                   ...)$sensitivity
  }, numeric(1))
  ok <- is.finite(sens)
  mono <- rep(NA_real_, length(sens))
  if (sum(ok) >= 2L) mono[ok] <- stats::isoreg(log(n_grid[ok]), sens[ok])$yf
  res$curve <- data.frame(n = n_grid, sensitivity = sens, monotone = mono)
  res$n_range <- range(n_grid)
  res$provenance <- "simulated"

  if (sum(ok) < 2L) {
    res$notes <- c(res$notes, "too few successful simulations to invert.")
  } else if (res$power > max(mono, na.rm = TRUE)) {
    res$notes <- c(res$notes, sprintf(
      "target power %.2f not reached by n = %g (best %.2f); widen n_range.",
      res$power, max(n_grid), max(mono, na.rm = TRUE)))
  } else if (res$power <= min(mono, na.rm = TRUE)) {
    res$n <- min(n_grid[ok])
    res$notes <- c(res$notes, sprintf(
      "target power is reached at or below the smallest n searched (%g).",
      min(n_grid[ok])))
  } else {
    # first crossing of the monotone curve, interpolated on log(n)
    ln <- log(n_grid[ok]); y <- mono[ok]
    i <- which(y >= res$power)[1]
    t <- (res$power - y[i - 1]) / (y[i] - y[i - 1])
    res$n <- ceiling(exp(ln[i - 1] + t * (ln[i] - ln[i - 1])))
  }
  structure(res, class = "howbig_n")
}

#' @export
print.howbig_n <- function(x, ...) {
  tag <- switch(x$provenance, in_grid = "in-grid", interpolated = "interpolated",
                simulated = "simulated", simulate_required = "OFF-GRID")
  cat(sprintf("howbig sample-size estimate  [%s]\n", tag))
  cat(sprintf("  %s | %s perturbation | p = %g | target power = %.2f @ FDR %.2f\n",
              x$method, x$perturbation, x$p, x$power, x$fdr_target))
  if (is.finite(x$n)) {
    cat(sprintf("  Required per-group n:  %g\n", x$n))
  } else {
    cat("  Required n: not available.\n")
  }
  if (identical(x$provenance, "simulated") && !is.null(x$curve)) {
    cat("  simulated power curve:\n")
    cv <- x$curve
    cat(sprintf("    n = %-6g sensitivity = %.2f\n", cv$n, cv$sensitivity), sep = "")
  }
  for (nt in x$notes) cat(strwrap(paste("  note:", nt), exdent = 8), sep = "\n")
  invisible(x)
}
