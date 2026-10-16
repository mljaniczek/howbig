# Access layer over the internal `howbig_grid` lookup table.
#
# The grid is a "plus sign": a sample-size arm (varying n at p = 250) and a
# feature arm (varying p at n = 1000) that cross at a single anchor, plus a
# handful of extra cells. Only the `structural` perturbation has full arms;
# `magnitude` / `both` exist only at the p = 250, n = 1000 anchor.

# Metric column stems present in the grid (paired with _mean/_sd/_q10/_q90).
.HOWBIG_METRICS <- c("sensitivity", "achieved_fdr", "aupr")

# Return the single matching grid row, or NULL.
.grid_cell <- function(method, p, n, perturbation) {
  hit <- howbig_grid[howbig_grid$method == method &
                       howbig_grid$p == p &
                       howbig_grid$n == n &
                       howbig_grid$perturbation == perturbation, ]
  if (nrow(hit) == 0L) return(NULL)
  hit[1L, ]
}

# All grid rows on one arm for a (method, perturbation): axis "n" -> the
# p = 250 sweep; axis "p" -> the n = 1000 sweep. Ordered by the swept axis.
.arm_table <- function(method, perturbation, axis) {
  base <- howbig_grid[howbig_grid$method == method &
                        howbig_grid$perturbation == perturbation, ]
  arm <- if (axis == "n") base[base$p == 250, ] else base[base$n == 1000, ]
  arm[order(arm[[axis]]), ]
}

# Classify a query. Returns list(kind, axis):
#   kind "in_grid"            exact cell exists (axis NA)
#   kind "interpolated"       on an arm, between nodes (axis "n" or "p")
#   kind "simulate_required"  off-grid (axis NA)
.classify <- function(method, p, n, perturbation) {
  if (!is.null(.grid_cell(method, p, n, perturbation)))
    return(list(kind = "in_grid", axis = NA_character_))

  # arms only exist for the structural perturbation
  if (perturbation == "structural") {
    narm <- .arm_table(method, perturbation, "n")
    if (p == 250 && n >= min(narm$n) && n <= max(narm$n))
      return(list(kind = "interpolated", axis = "n"))
    parm <- .arm_table(method, perturbation, "p")
    if (n == 1000 && p >= min(parm$p) && p <= max(parm$p))
      return(list(kind = "interpolated", axis = "p"))
  }
  list(kind = "simulate_required", axis = NA_character_)
}

# Interpolate every metric summary along an arm at the query point.
.interp_cell <- function(method, perturbation, axis, at) {
  arm <- .arm_table(method, perturbation, axis)
  x <- log(arm[[axis]])
  lx <- log(at)
  out <- list()
  for (m in .HOWBIG_METRICS) {
    for (stat in c("mean", "q10", "q90")) {
      col <- paste0(m, "_", stat)
      out[[col]] <- .interp(x, arm[[col]], lx)
    }
  }
  out$n_reps <- NA_integer_   # interpolated, not a single cell
  out
}

#' Coverage of the built-in grid
#'
#' Returns the precomputed cells `howbig` can answer directly (before any
#' interpolation or simulation), one row per (method, p, n, perturbation).
#'
#' @param method Optional method label to filter to.
#' @return A data frame of grid cells with mean sensitivity, achieved FDR
#'   and AUPR, and the replicate count behind each cell.
#' @export
howbig_coverage <- function(method = NULL) {
  g <- howbig_grid[, c("method", "perturbation", "p", "n",
                        "sensitivity_mean", "achieved_fdr_mean",
                        "aupr_mean", "n_reps")]
  if (!is.null(method)) g <- g[g$method == .normalize_method(method), ]
  rownames(g) <- NULL
  g
}
