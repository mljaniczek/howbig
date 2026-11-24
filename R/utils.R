# Internal utilities: input normalization and shared constants.

`%||%` <- function(x, y) if (is.null(x)) y else x

# The FDR target the shipped grid was scored at. Other targets require
# simulation.
.HOWBIG_FDR_TARGET <- 0.05

# Perturbation vocabulary: what changes between the two groups.
.HOWBIG_PERTS <- c("structural", "magnitude", "both")

# Scenario vocabulary: assumptions about the networks and data. "default" is
# the main grid (block-structured BA-tree networks, Gaussian data, empirical
# weights); the others exist only at the anchor design (p = 100, n = 1000,
# structural change).
.HOWBIG_SCENARIOS <- c("default", "scale_free", "erdos_renyi", "ar1_banded",
                       "censored_5", "censored_25", "censored_50",
                       "uniform_weights")

.normalize_perturbation <- function(perturbation) {
  match.arg(perturbation, .HOWBIG_PERTS)
}

.normalize_scenario <- function(scenario) {
  match.arg(scenario, .HOWBIG_SCENARIOS)
}

.normalize_method <- function(method) {
  methods <- howbig_method_info$method
  if (length(method) != 1L || !method %in% methods) {
    stop("`method` must be one of: ", paste(methods, collapse = ", "),
         call. = FALSE)
  }
  method
}

# Per-group sample sizes as an integer pair c(n1, n2).
.normalize_n <- function(n) {
  stopifnot(is.numeric(n), length(n) %in% 1:2, all(n > 0))
  as.integer(round(rep_len(n, 2L)))
}

# Monotone (or linear, for 2 points) interpolation of `y` at `at`, with NO
# extrapolation: returns NA when `at` falls outside the finite support.
# `x` is expected already on the interpolation scale (e.g. log(n)).
.interp <- function(x, y, at) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2L) return(NA_real_)
  o <- order(x); x <- x[o]; y <- y[o]
  if (at < min(x) - 1e-9 || at > max(x) + 1e-9) return(NA_real_)
  if (length(x) == 2L) return(stats::approx(x, y, xout = at)$y)
  stats::splinefun(x, y, method = "monoH.FC")(at)
}

# Is achieved FDR far enough from target to warn the user?
.miscalibrated <- function(achieved_fdr, target = .HOWBIG_FDR_TARGET) {
  is.finite(achieved_fdr) && achieved_fdr > 2 * target
}
