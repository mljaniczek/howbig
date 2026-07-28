#------------------------------------------------------
# 05d_score_results.R — Phase 2: SCORE saved raw fits.
#
# Reads per-rep raw fit files from
#   results/raw_fits/<method>/<setting>/_config_<hash>/rep_XXXX.rds
# (produced by 05b_run_methods.R), applies the current scoring
# rules, and writes the aggregated metrics table to
#   results/method_results/<method>/<setting>.rds
# alongside a meta record
#   results/method_results/<method>/<setting>.meta.rds
# documenting exactly which config + scoring options produced it.
#
# Scoring is cheap (~seconds per setting), so iterate freely:
# change a threshold, rerun this script, regenerate metrics
# without re-fitting anything.
#
# Default config picking rule: "latest" = the _config_<hash>
# with the most reps on disk, tie-break by mtime. Override
# with --config <hash> to score a specific configuration.
#
# Usage:
#   # Score everything at default options (pip_mask=0.05,
#   # diff_threshold=1e-5 (nominal-zero floor; matches manuscript),
#   # fdr_target=0.05, pooled null)
#   Rscript 05d_score_results.R \
#     --methods "JGL_BIC,JGL_BIC_GGL,SSJGL,SSJGL_GGL,multiGGM,JRF" \
#     --settings "SF-p50n1000,SF-p125n100,SF-p250n1000,SF+p250n1000"
#
#   # Try a different PIP mask threshold without refitting
#   Rscript 05d_score_results.R --pip-mask-thr 0.03 \
#     --methods "multiGGM" --settings "SF-p250n1000"
#
#   # Use a looser FDR target for the inference columns
#   Rscript 05d_score_results.R --fdr-target 0.10 \
#     --methods "SSJGL" --settings "SF-p250n1000"
#
#   # Try raw empirical p-values instead of pooled BH
#   Rscript 05d_score_results.R --null-correction none \
#     --methods "SSJGL" --settings "SF-p250n1000"
#
#   # Score a specific config explicitly
#   Rscript 05d_score_results.R --config 3bbfee91 \
#     --methods "multiGGM" --settings "SF-p250n1000"
#------------------------------------------------------

suppressPackageStartupMessages({
  library(future)
  library(furrr)
  library(digest)
})

source("00_helper_functions.R")
source("01_generate_gold_standard_networks.R")
source("01b_generate_paper_gold_standards.R")
source("03_estimate.R")
source("04_calc_metrics.R")

empirical_dist <- read.csv("mxDist.csv")

#------------------------------------------------------
# CLI
#------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
parse_flag_value <- function(flag, default = NULL) {
  if (!(flag %in% args)) return(default)
  idx <- which(args == flag)[1]
  if (idx >= length(args)) return(default)
  args[idx + 1]
}

test_mode    <- "--test" %in% args
methods_arg  <- parse_flag_value("--methods")
settings_arg <- parse_flag_value("--settings")
config_arg   <- parse_flag_value("--config", "latest")
pip_mask_thr <- as.numeric(parse_flag_value("--pip-mask-thr", "0.05"))
threshold    <- as.numeric(parse_flag_value("--threshold", "0.00000001"))
diff_thr     <- as.numeric(parse_flag_value("--diff-threshold", "0.00001"))
pcor_thr     <- as.numeric(parse_flag_value("--pcor-threshold", "0.01"))
fdr_target   <- as.numeric(parse_flag_value("--fdr-target", "0.05"))
null_corr    <- parse_flag_value("--null-correction", "pooled")

method_names <- if (!is.null(methods_arg))
  strsplit(methods_arg, ",")[[1]] else NULL
settings_req <- if (!is.null(settings_arg))
  strsplit(settings_arg, ",")[[1]] else NULL

if (test_mode) {
  if (is.null(method_names))
    method_names <- c("JGL_BIC", "JGL_BIC_GGL", "SSJGL", "SSJGL_GGL",
                       "multiGGM", "JRF")
  if (is.null(settings_req)) settings_req <- c("SF-p50n1000")
}
if (is.null(method_names))
  stop("Must specify --methods (or pass --test).")
if (is.null(settings_req))
  stop("Must specify --settings (or pass --test).")

cat("=== Scoring raw fits ===\n")
cat(sprintf("  methods:     %s\n", paste(method_names, collapse = ", ")))
cat(sprintf("  settings:    %s\n", paste(settings_req, collapse = ", ")))
cat(sprintf("  config:      %s\n", config_arg))
cat(sprintf("  thresholds:  group=%.3g  diff=%.3g  pcor=%.3g  pip_mask=%.3g\n",
            threshold, diff_thr, pcor_thr, pip_mask_thr))
cat(sprintf("  inference:   fdr_target=%.3g  null_correction=%s\n",
            fdr_target, null_corr))

#------------------------------------------------------
# Utilities
#------------------------------------------------------

# Pick a config subdir for (method, setting): if user passed an
# explicit hash, use it; if "latest", pick the one with the most
# rep files (tie-break by mtime).
pick_config_dir <- function(method, setting, arg = "latest") {
  base <- file.path("results/raw_fits", method, setting)
  if (!dir.exists(base)) return(NULL)
  configs <- list.dirs(base, recursive = FALSE)
  if (length(configs) == 0) return(NULL)

  if (arg != "latest") {
    want <- file.path(base, paste0("_config_", arg))
    return(if (dir.exists(want)) want else NULL)
  }

  # Latest: most reps, break by mtime
  scores <- vapply(configs, function(d) {
    n_reps <- length(list.files(d, pattern = "^rep_\\d+\\.rds$"))
    n_reps + as.numeric(file.info(d)$mtime) / 1e12
  }, numeric(1))
  configs[which.max(scores)]
}

#------------------------------------------------------
# Main loop
#------------------------------------------------------

dir.create("results/method_results", showWarnings = FALSE, recursive = TRUE)

# Parallel scoring — cheap per rep, but worth parallelizing when
# there are lots of files. Keep workers modest.
default_workers <- max(1, min(availableCores() - 1, 6))
plan(multisession, workers = default_workers)
cat(sprintf("  workers:     %d\n\n", default_workers))

for (method_name in method_names) {
  for (setting in settings_req) {
    cfg_dir <- pick_config_dir(method_name, setting, config_arg)
    if (is.null(cfg_dir)) {
      cat(sprintf("[%s / %s] no raw_fits found — skipping\n",
                  method_name, setting))
      next
    }
    cfg_hash <- sub(".*_config_", "", basename(cfg_dir))

    # List all rep files in this config
    rep_files <- sort(list.files(cfg_dir, pattern = "^rep_\\d+\\.rds$",
                                  full.names = TRUE))
    if (length(rep_files) == 0) {
      cat(sprintf("[%s / %s] config %s has 0 reps — skipping\n",
                  method_name, setting, cfg_hash))
      next
    }

    # Warn if multiple configs exist
    n_configs <- length(list.dirs(file.path("results/raw_fits", method_name, setting),
                                   recursive = FALSE))
    extra_config_note <- if (n_configs > 1)
      sprintf(" (1 of %d configs; use --config to pick another)", n_configs) else ""

    cat(sprintf("[%s / %s] cfg=%s, %d reps%s — scoring...\n",
                method_name, setting, cfg_hash, length(rep_files),
                extra_config_note))

    # Score in parallel — each worker reads its own rep, loads the
    # matching gold standard, calls score_method(), returns one row.
    setting_params <- get_setting_params(setting)
    rows <- future_map(rep_files, function(rep_file) {
      rec <- readRDS(rep_file)
      gold_file <- file.path("gold_standards", setting,
                              sprintf("rep_%04d.rds", rec$rep_id))
      if (!file.exists(gold_file)) {
        na_row <- data.frame(method = method_name, setting = setting,
                              rep_id = rec$rep_id,
                              status = "gold_missing",
                              stringsAsFactors = FALSE)
        return(na_row)
      }
      gold <- readRDS(gold_file)

      row <- score_method(method_name, rec$fit, gold,
                          setting              = setting,
                          threshold            = threshold,
                          pcor_threshold       = pcor_thr,
                          diff_threshold       = diff_thr,
                          pip_mask_thr         = pip_mask_thr,
                          fdr_target           = fdr_target,
                          null_correction      = null_corr)

      # Augment with rep-level metadata carried from the raw record
      row$setting     <- setting
      row$rep_id      <- rec$rep_id
      row$seed        <- rec$seed_gold
      row$seed_data   <- rec$seed_data
      row$config_hash <- rec$config_hash
      row$n           <- setting_params$n
      row$topology    <- setting_params$topology    %||% NA_character_
      row$perturbation <- setting_params$setting    %||% NA_character_
      row$within_density <- setting_params$within_density %||% NA_real_
      row$between_prob   <- setting_params$between_prob   %||% NA_real_
      row$wt_method   <- setting_params$wt_method   %||% NA_character_
      row$generator   <- setting_params$generator   %||% "standard"
      row
    }, .options = furrr_options(
      seed = TRUE,
      globals = c("score_method", "method_name", "setting",
                  "threshold", "pcor_thr", "diff_thr", "pip_mask_thr",
                  "fdr_target", "null_corr", "setting_params",
                  "get_setting_params",
                  "calc_metrics", "calc_metrics_from_pcor",
                  "calc_metrics_from_importance", "calc_metrics_from_pip",
                  ".build_result_row", ".xor_truth", ".init_metrics_na_row",
                  ".precision_methods", "apply_cached_inference",
                  "apply_bayesian_fdr_pip", "apply_stability_fdr",
                  ".bayes_fdr_call", ".score_bayes_call",
                  "threshold_to_adjacency", "pcor_from_precision",
                  "compute_tpr_fpr", "compute_auc", "compute_aupr",
                  "compute_fdr", "compute_mcc", "compute_kl_divergence",
                  "compute_relative_frobenius"),
      packages = c("digest")
    ))

    results_df <- do.call(rbind, rows)
    outfile <- file.path("results/method_results", method_name,
                         paste0(setting, ".rds"))
    dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
    saveRDS(results_df, outfile)

    # Side-car meta file so downstream scripts know exactly what
    # scoring options produced this metrics table.
    meta <- list(
      method_name     = method_name,
      setting         = setting,
      config_hash     = cfg_hash,
      n_reps_scored   = nrow(results_df),
      scoring_options = list(
        threshold       = threshold,
        pcor_threshold  = pcor_thr,
        diff_threshold  = diff_thr,
        pip_mask_thr    = pip_mask_thr,
        fdr_target      = fdr_target,
        null_correction = null_corr
      ),
      timestamp = Sys.time()
    )
    meta_file <- file.path("results/method_results", method_name,
                           paste0(setting, ".meta.rds"))
    saveRDS(meta, meta_file)

    n_success <- sum(results_df$status == "success", na.rm = TRUE)
    cat(sprintf("[%s / %s] wrote %s (%d/%d success)\n\n",
                method_name, setting, outfile, n_success, nrow(results_df)))
  }
}

cat("Scoring complete.\n")
cat("Next: Rscript 05c_merge_results.R\n")
