# Access layer over the internal `howbig_grid` lookup table.
#
# Cells are keyed by (method, scenario, perturbation, p, n1, n2). The
# default scenario has two swept arms that interpolation may use:
#   sample-size arm : n = 100 ... 5000 at p = 250 (balanced groups)
#   feature arm     : p = 50 ... 500 at n = 1000 (balanced groups)
# plus single cells off the arms (small n at p = 100, unequal groups, the
# perturbation types and the assumption scenarios at the p = 100, n = 1000
# anchor). Off-arm cells are answered exactly or not at all.

.HOWBIG_METRICS <- c("sensitivity", "achieved_fdr", "aupr")
.N_ARM_P <- 250L
.P_ARM_N <- 1000L

# Return the single matching grid row, or NULL. `n` is c(n1, n2).
.grid_cell <- function(method, p, n, perturbation, scenario = "default") {
  g <- howbig_grid
  hit <- g[g$method == method & g$scenario == scenario &
             g$perturbation == perturbation & g$p == p &
             g$n1 == n[1] & g$n2 == n[2], ]
  if (nrow(hit) == 0L) return(NULL)
  hit[1L, ]
}

# Rows on one arm of the default scenario for (method, perturbation):
# axis "n" -> the p = 250 sweep; axis "p" -> the n = 1000 sweep. Balanced
# groups only; ordered by the swept axis. `n` holds the per-group size.
.arm_table <- function(method, perturbation, axis) {
  g <- howbig_grid
  base <- g[g$method == method & g$scenario == "default" &
              g$perturbation == perturbation & g$n1 == g$n2, ]
  arm <- if (axis == "n") base[base$p == .N_ARM_P, ] else base[base$n1 == .P_ARM_N, ]
  arm$n <- arm$n1
  arm[order(arm[[axis]]), ]
}

# Classify a query. Returns list(kind, axis, reason):
#   kind "in_grid"            exact cell exists
#   kind "interpolated"       on an arm, between nodes (axis "n" or "p")
#   kind "simulate_required"  otherwise; `reason` says why
.classify <- function(method, p, n, perturbation, scenario) {
  if (!is.null(.grid_cell(method, p, n, perturbation, scenario)))
    return(list(kind = "in_grid", axis = NA_character_, reason = NULL))

  off <- function(reason) list(kind = "simulate_required", axis = NA_character_,
                               reason = reason)
  if (scenario != "default")
    return(off(sprintf(paste0("scenario '%s' is in the grid only at p = 100, ",
                              "n = 1000 with a structural change"), scenario)))
  if (n[1] != n[2])
    return(off(paste0("these unequal group sizes are not in the grid ",
                      "(see howbig_coverage() for the cells that are)")))
  if (perturbation == "structural") {
    narm <- .arm_table(method, perturbation, "n")
    if (p == .N_ARM_P && nrow(narm) && n[1] >= min(narm$n) && n[1] <= max(narm$n))
      return(list(kind = "interpolated", axis = "n", reason = NULL))
    parm <- .arm_table(method, perturbation, "p")
    if (n[1] == .P_ARM_N && nrow(parm) && p >= min(parm$p) && p <= max(parm$p))
      return(list(kind = "interpolated", axis = "p", reason = NULL))
  }
  off(sprintf("(p = %g, n = %g) is off the grid for the %s perturbation",
              p, n[1], perturbation))
}

# Interpolate every metric summary along an arm at the query point.
.interp_cell <- function(method, perturbation, axis, at) {
  arm <- .arm_table(method, perturbation, axis)
  x <- log(arm[[axis]])
  lx <- log(at)
  cols <- c(as.vector(outer(.HOWBIG_METRICS, c("mean", "q10", "q90"), paste,
                            sep = "_")), "n_diff_edges")
  out <- lapply(stats::setNames(cols, cols), function(col) .interp(x, arm[[col]], lx))
  out$n_reps <- NA_integer_   # interpolated, not a single cell
  out
}

#' Coverage of the built-in grid
#'
#' Returns the precomputed cells `howbig` can answer directly, before any
#' interpolation or simulation: one row per method, scenario, perturbation,
#' number of features, and pair of group sizes.
#'
#' @param method Optional method label to filter to.
#' @param scenario Optional scenario to filter to (see [power_for()]).
#' @return A data frame of grid cells with mean sensitivity, achieved FDR,
#'   AUPR, the mean number of truly differential edges, the number of
#'   feature pairs, and the replicate count behind each cell.
#' @examples
#' head(howbig_coverage("JGL_BIC"))
#' howbig_coverage("JGL_BIC", scenario = "erdos_renyi")
#' @export
howbig_coverage <- function(method = NULL, scenario = NULL) {
  g <- howbig_grid[, c("method", "scenario", "perturbation", "p", "n1", "n2",
                        "sensitivity_mean", "achieved_fdr_mean", "aupr_mean",
                        "n_diff_edges", "n_pairs", "n_reps")]
  if (!is.null(method)) g <- g[g$method == .normalize_method(method), ]
  if (!is.null(scenario)) g <- g[g$scenario == .normalize_scenario(scenario), ]
  rownames(g) <- NULL
  g
}

#' Methods available in howbig
#'
#' Lists the estimation methods in the grid, how each one decides which
#' edges are differential at the FDR target, and whether [simulate_power()]
#' can run it for designs the grid does not cover.
#'
#' @return A data frame with columns `method`, `description`, `inference`
#'   and `simulate`.
#' @examples
#' howbig_methods()
#' @export
howbig_methods <- function() {
  howbig_method_info
}
