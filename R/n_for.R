#' Required sample size for a target power
#'
#' Inverts [power_for()]: finds the smallest per-group sample size `n` that
#' reaches a target power (sensitivity) at a feature count `p`.
#'
#' In this version the inversion uses the precomputed sample-size arm, which
#' exists at `p = 250` for the `structural` perturbation. Other `p` or
#' perturbations return `simulate_required` (the simulation fallback that
#' would sweep `n` at an arbitrary `p` is not yet wired).
#'
#' @param p Number of features.
#' @param power Target sensitivity in (0, 1).
#' @param method Estimation method label.
#' @param perturbation Type of differential signal; see [power_for()].
#' @param fdr Target false-discovery rate (only `0.05` is available).
#'
#' @return An object of class `howbig_n`: the required `n`, provenance, the
#'   `n` range searched, and any notes.
#' @examples
#' n_for(p = 250, power = 0.8, method = "JGL_BIC")
#' @export
n_for <- function(p, power = 0.8, method,
                  perturbation = "structural",
                  fdr = .HOWBIG_FDR_TARGET) {
  stopifnot(length(p) == 1L, is.numeric(p), p > 0,
            length(power) == 1L, power > 0, power < 1)
  method <- .normalize_method(method)
  perturbation <- .normalize_perturbation(perturbation)

  res <- list(method = method, p = p, power = power,
              perturbation = perturbation, fdr_target = .HOWBIG_FDR_TARGET,
              provenance = "simulate_required", n = NA_real_,
              n_range = c(NA_real_, NA_real_), notes = character(0))

  # v0: invertible only along the p = 250 sample-size arm, structural.
  if (!(p == 250 && perturbation == "structural")) {
    res$notes <- c(res$notes, paste0(
      "n_for currently inverts only the p = 250 structural sample-size arm; ",
      "other designs need the simulation fallback (not in this version)."))
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
  structure(res, class = "howbig_n")
}

#' @export
print.howbig_n <- function(x, ...) {
  tag <- switch(x$provenance, in_grid = "in-grid",
                interpolated = "interpolated", simulate_required = "OFF-GRID")
  cat(sprintf("howbig sample-size estimate  [%s]\n", tag))
  cat(sprintf("  %s | %s perturbation | p = %g | target power = %.2f @ FDR %.2f\n",
              x$method, x$perturbation, x$p, x$power, x$fdr_target))
  if (is.finite(x$n)) {
    cat(sprintf("  Required per-group n:  %g\n", x$n))
  } else {
    cat("  Required n: not available from the grid.\n")
  }
  for (nt in x$notes) cat(strwrap(paste("  note:", nt), exdent = 8), sep = "\n")
  invisible(x)
}
