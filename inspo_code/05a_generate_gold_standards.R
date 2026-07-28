#------------------------------------------------------
# Gold Standard Generation
# Generates and saves all gold standards to disk so
# that methods can be run incrementally on the same
# networks.
#
# Usage:
#   Rscript 05a_generate_gold_standards.R                     # full: all settings x 50 reps
#   Rscript 05a_generate_gold_standards.R --test              # test: 1 setting x 2 reps
#   Rscript 05a_generate_gold_standards.R --settings "SF-p125n100,SF-p250n1000" # specific settings only
#   Rscript 05a_generate_gold_standards.R --force             # overwrite existing files
#   Rscript 05a_generate_gold_standards.R --settings "Danaher,SSJGL-paper,PSV,Shaddox,DINGO-paper" # paper settings

#------------------------------------------------------

# TODO add ability to overwrite existing ("force" like in the other scripts)
library(igraph)
library(MASS)
library(future)
library(furrr)

source("00_helper_functions.R")
source("01_generate_gold_standard_networks.R")
source("01b_generate_paper_gold_standards.R")

# Load empirical weight distribution
empirical_dist <- read.csv("mxDist.csv")

#------------------------------------------------------
# Parse command-line arguments
#------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
test_mode <- "--test" %in% args
force_mode <- "--force" %in% args

# Parse --settings flag
settings_arg <- NULL
if ("--settings" %in% args) {
  idx <- which(args == "--settings")
  if (idx < length(args)) {
    settings_arg <- strsplit(args[idx + 1], ",")[[1]]
  }
}

# Parse --gold-mode flag: "fresh" (default) or "fixed"
#   fresh: each replicate draws a new gold standard (seed_gold varies
#          per rep). Matches 4 of 5 source papers.
#   fixed: all replicates within a setting share the SAME gold
#          standard (seed_gold is constant per-setting). Required
#          for per-edge permutation inference in Phase 2.
gold_mode <- "fresh"
if ("--gold-mode" %in% args) {
  idx <- which(args == "--gold-mode")
  if (idx < length(args)) {
    gold_mode <- args[idx + 1]
    if (!gold_mode %in% c("fresh", "fixed"))
      stop("--gold-mode must be 'fresh' or 'fixed'")
  }
}


all_settings <- c(
  # Primary tutorial settings (body of 09 report). SF-p50n1000
  # also serves as the small-p anchor of the feature-size sweep.
  "SF-p50n1000", "SF-p125n100", "SF-p250n1000", "SF+p250n1000",
  # Sample-size sweep (p=250 fixed, vary n)
  "SF-p250n250", "SF-p250n500", "SF-p250n2000", "SF-p250n5000",
  # Feature-size sweep (n=1000 fixed, vary p; p=50 is in primaries above)
  "SF-p100n1000", "SF-p200n1000", "SF-p300n1000", "SF-p500n1000",
  # Additive perturbation (Setting C: g2 = g1 + new Block-1 edges)
  "SFCp250n1000",
  # Extreme high-dim, JRF-paper-style (p >> n)
  "SF-p500n50",
  # Non-linear machinery (kept, not in report)
  "SF-nonlinear-saturation",
  # Paper-replicating settings
  "Danaher", "Danaher-full", "SSJGL-paper", "PSV", "Shaddox", "DINGO-paper"
)
# Note: "SF-p20n1000" is still defined in get_setting_params() for
# verification / spot-checks, but is not part of the production
# roster — superseded by SF-p50n1000 (commit on 2026-04-25).

if (test_mode) {
  n_reps <- 2
  settings <- if (!is.null(settings_arg)) settings_arg else c("SF-p50n1000")
} else {
  n_reps <- 100
  settings <- if (!is.null(settings_arg)) settings_arg else all_settings
}


#------------------------------------------------------
# Build parameter grid with two seeds per replicate
# seed_gold: for generate_gold_standard()
# seed_data: for gen_dat() at estimation time
#
# The grid is always built for ALL settings x ALL reps
# so that seeds are deterministic regardless of which
# subset is generated. We then filter to the requested
# settings/reps.
#------------------------------------------------------

set.seed(20240101)
full_grid <- expand.grid(
  setting_label = all_settings,
  rep_id = 1:1000,
  stringsAsFactors = FALSE
)
# Sort to match expand_grid order: setting varies slowest
full_grid <- full_grid[order(full_grid$setting_label, full_grid$rep_id), ]
rownames(full_grid) <- NULL

# seed_gold first (matches original 'seed' column from 05_run_simulation.R)
full_grid$seed_gold <- round(runif(nrow(full_grid)) * 2^30) + seq_len(nrow(full_grid))
# seed_data second (new, for reproducible data sampling)
full_grid$seed_data <- round(runif(nrow(full_grid)) * 2^30) + seq_len(nrow(full_grid))

# Fixed-GS mode: override seed_gold to a setting-specific constant
# so every replicate within a setting shares the same gold standard.
# seed_data still varies per rep, so data is re-drawn each time.
if (gold_mode == "fixed") {
  # Deterministic per-setting seed: hash the setting label to an integer.
  # Note: 2654435761 overflows int32, so we compute in double precision
  # and cast back to integer after modulo (which guarantees the result
  # fits in int32).
  setting_seed <- function(label) {
    # Parens matter: %% binds tighter than *, so without them we would
    # compute sum * (multiplier %% int.max) which overflows int32.
    as.integer((sum(utf8ToInt(label)) * 2654435761) %% .Machine$integer.max)
  }
  full_grid$seed_gold <- vapply(
    full_grid$setting_label, setting_seed, FUN.VALUE = integer(1)
  )
  cat("Gold-standard mode: FIXED (one GS per setting shared across reps)\n")
} else {
  cat("Gold-standard mode: FRESH (new GS per replicate)\n")
}

# Filter to requested settings and reps
params_grid <- full_grid[full_grid$setting_label %in% settings &
                         full_grid$rep_id <= n_reps, ]
rownames(params_grid) <- NULL

cat(sprintf("Settings: %s\n", paste(settings, collapse = ", ")))
cat(sprintf("Reps per setting: %d\n", n_reps))
cat(sprintf("Total gold standards to generate: %d\n", nrow(params_grid)))
if (force_mode) cat("Force mode: ON (overwriting existing gold standards)\n")

#------------------------------------------------------
# Save manifest (always save the full grid so that
# 05b_run_methods.R can look up seed_data for any rep)
#------------------------------------------------------

dir.create("gold_standards", showWarnings = FALSE)
saveRDS(full_grid, "gold_standards/manifest.rds")
cat("Manifest saved to gold_standards/manifest.rds\n\n")

#------------------------------------------------------
# Generate gold standards in parallel
#------------------------------------------------------

n_workers <- max(1, availableCores() - 1)
cat(sprintf("Using %d parallel workers\n\n", n_workers))
plan(multisession, workers = n_workers)

start_time <- Sys.time()

for (setting in settings) {
  setting_dir <- file.path("gold_standards", setting)
  dir.create(setting_dir, showWarnings = FALSE, recursive = TRUE)

  sub_grid <- params_grid[params_grid$setting_label == setting, ]
  params <- get_setting_params(setting)

  # Check how many already exist (for resumability)
  if (force_mode) {
    existing <- 0L
  } else {
    existing <- sum(file.exists(
      file.path(setting_dir, sprintf("rep_%04d.rds", sub_grid$rep_id))
    ))
  }
  to_generate <- nrow(sub_grid) - existing

  cat(sprintf("[%s] %d reps total, %d already exist, %d to generate%s\n",
              setting, nrow(sub_grid), existing, to_generate,
              if (force_mode) " (--force: overwriting)" else ""))

  if (to_generate == 0) {
    cat(sprintf("[%s] All done — skipping\n\n", setting))
    next
  }

  # Generate in parallel across reps
  future_pwalk(
    list(
      rep_id = sub_grid$rep_id,
      seed_gold = sub_grid$seed_gold
    ),
    function(rep_id, seed_gold) {
      outfile <- file.path(setting_dir, sprintf("rep_%04d.rds", rep_id))
      if (!force_mode && file.exists(outfile)) return(invisible(NULL))

      # Route to paper-specific generator or standard pipeline
      if (isTRUE(params$generator == "paper")) {
        gold <- generate_paper_gold_standard(
          paper = params$paper,
          seed = seed_gold,
          params = params
        )
      } else {
        gold <- generate_gold_standard(
          p = params$p,
          B = params$B %||% 1,
          p_block = params$p_block %||% params$p,
          topology = params$topology,
          # Renamed from `within_density`; the old field is silently
          # ignored. See NOTES_network_density.md.
          within_density_random = params$within_density_random %||%
                                   params$within_density %||% 0.15,
          between_prob = params$between_prob %||% 0,
          setting = params$setting,
          wt_method = params$wt_method,
          empirical_dist = empirical_dist,
          epsilon = params$epsilon,
          precision_method = params$precision_method %||% "row_sum_normalize",
          d_scale = params$d_scale,
          scale_factor = params$scale_factor,
          m_pa = params$m_pa,
          # AR(1) fields (ignored unless topology = "ar1_banded")
          ar1_rho       = params$ar1_rho       %||% 0.4,
          ar1_bandwidth = params$ar1_bandwidth %||% 1,
          seed = seed_gold
        )
      }

      saveRDS(gold, outfile)
    },
    .options = furrr_options(
      seed = TRUE,
      globals = c("generate_gold_standard", "generate_gold_standard_ar1",
                   "generate_block_graph",
                   "assemble_block_network", "assign_edge_weights",
                   "build_precision_matrix",
                   "build_precision_eigenvalue_boost",
                   "build_precision_row_sum",
                   "apply_perturbation",
                   "compute_gold_standard_truths",
                   "generate_paper_gold_standard",
                   "generate_gs_danaher", "generate_gs_ssjgl",
                   "generate_gs_psv", "generate_gs_shaddox",
                   "generate_gs_dingo",
                   "build_precision_danaher", "build_precision_shaddox",
                   "empirical_dist", "params", "setting_dir", "force_mode"),
      packages = c("igraph", "MASS", "BDgraph")
    ),
    .progress = TRUE
  )

  # Verify
  n_saved <- length(list.files(setting_dir, pattern = "\\.rds$"))
  cat(sprintf("[%s] Done — %d files saved\n\n", setting, n_saved))
}

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

cat(sprintf("Gold standard generation complete in %.1f minutes\n", as.numeric(elapsed)))
cat(sprintf("Output directory: gold_standards/\n"))

# Summary of what's on disk
for (s in settings) {
  n_files <- length(list.files(file.path("gold_standards", s), pattern = "\\.rds$"))
  cat(sprintf("  %s: %d files\n", s, n_files))
}
