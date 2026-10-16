# Internal utilities: input normalization and shared constants.

`%||%` <- function(x, y) if (is.null(x)) y else x

# The FDR target the shipped grid was scored at. In v0 this is the only
# operating point available from the grid; other targets require the raw
# fits (deferred to the simulation path).
.HOWBIG_FDR_TARGET <- 0.05

# Perturbation vocabulary. Users speak in structural / magnitude / both;
# the grid stores Deletion / Activation / Deletion + Activation.
.HOWBIG_PERTS <- c("structural", "magnitude", "both")

.normalize_perturbation <- function(perturbation) {
  perturbation <- match.arg(perturbation, .HOWBIG_PERTS)
  perturbation
}

.normalize_method <- function(method) {
  methods <- sort(unique(howbig_grid$method))
  if (length(method) != 1L || !method %in% methods) {
    stop("`method` must be one of: ", paste(methods, collapse = ", "),
         call. = FALSE)
  }
  method
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
