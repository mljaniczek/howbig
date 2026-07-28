#------------------------------------------------------
# 05b_run_methods.R — Phase 1: FIT and save raw outputs.
#
# For each (method, setting, rep) this script:
#   1. Loads the gold standard and regenerates data deterministically
#      from the saved seeds.
#   2. Calls fit_method() to produce raw fit output (theta / pip /
#      pcor / importance / ...).
#   3. Saves that raw output to
#      results/raw_fits/<method>/<setting>/_config_<hash>/rep_XXXX.rds
#      ...with the rep's metadata embedded.
#
# Scoring is a SEPARATE phase; 05d_score_results.R reads these raw
# files and produces the metrics table. That way you can tweak
# thresholds / PIP masks / FDR targets without re-fitting.
#
# The script is resumable: each rep is an independent file, so
# rerunning after a crash / interruption picks up where it left
# off via the `--target-reps N` flag (default 100). You can start
# with a handful of reps, see how long it takes, and grow to the
# target count later.
#
# Configuration hashing: the full hyperparameter config for each
# method (from method_params() in 00_helper_functions.R) is
# hashed into the path via _config_<hash>. Different hyperparameter
# sets produce different subdirectories, so results never silently
# mix. To change hyperparameters, add a new method label in
# method_params() — never reuse an existing label.
#
# Usage:
#   # Fit 100 reps of 4 primary settings on the 6-method roster
#   # (DiNGO currently deferred — see CLAUDE.md "Open investigation:
#   # DiNGO behavior at p=250")
#   Rscript 05b_run_methods.R \
#     --methods "JGL_BIC,JGL_BIC_GGL,SSJGL,SSJGL_GGL,multiGGM,JRF" \
#     --settings "SF-p50n1000,SF-p125n100,SF-p250n1000,SF+p250n1000" \
#     --target-reps 100
#
#   # Start with 5 reps to verify pipeline; then run 95 more
#   Rscript 05b_run_methods.R --methods "multiGGM" \
#     --settings "SF-p250n1000" --target-reps 5
#   Rscript 05b_run_methods.R --methods "multiGGM" \
#     --settings "SF-p250n1000" --target-reps 100   # picks up from 6
#
#   # Quick smoke test
#   Rscript 05b_run_methods.R --test
#
#   # Force re-fit of existing reps
#   Rscript 05b_run_methods.R --methods "SSJGL" \
#     --settings "SF-p50n1000" --target-reps 10 --force
#
#   # List configs on disk for a (method, setting)
#   Rscript 05b_run_methods.R --list-configs \
#     --methods "multiGGM" --settings "SF-p250n1000"
#------------------------------------------------------

suppressPackageStartupMessages({
  library(igraph)
  library(MASS)
  library(JGL)
  library(future)
  library(furrr)
  library(spikeyglass)
  library(digest)
})

source("00_helper_functions.R")
source("01_generate_gold_standard_networks.R")
source("01b_generate_paper_gold_standards.R")
source("03_estimate.R")
source("04_calc_metrics.R")

empirical_dist <- read.csv("mxDist.csv")

#------------------------------------------------------
# CLI parsing
#------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

parse_flag_value <- function(flag, default = NULL) {
  if (!(flag %in% args)) return(default)
  idx <- which(args == flag)[1]
  if (idx >= length(args)) return(default)
  args[idx + 1]
}

test_mode    <- "--test" %in% args
force_mode   <- "--force" %in% args
list_configs <- "--list-configs" %in% args

methods_arg   <- parse_flag_value("--methods")
settings_arg  <- parse_flag_value("--settings")
target_reps   <- as.integer(parse_flag_value("--target-reps", "100"))
start_rep     <- as.integer(parse_flag_value("--start-rep",     "1"))
n_workers_arg <- parse_flag_value("--n-workers", NA)

method_names <- if (!is.null(methods_arg))
  strsplit(methods_arg, ",")[[1]] else NULL
settings_req <- if (!is.null(settings_arg))
  strsplit(settings_arg, ",")[[1]] else NULL

all_settings <- c(
  # Primary tutorial settings (body of 09 report). SF-p50n1000 also
  # serves as the small-p anchor of the feature-size sweep.
  "SF-p50n1000", "SF-p125n100", "SF-p250n1000", "SF+p250n1000",
  # Sample-size sweep (p=250 fixed, vary n)
  "SF-p250n250", "SF-p250n500", "SF-p250n2000", "SF-p250n5000",
  # Feature-size sweep (n=1000 fixed, vary p; p=50 is in primaries above)
  "SF-p100n1000", "SF-p200n1000", "SF-p300n1000", "SF-p500n1000",
  # Non-linear machinery (kept, not in report)
  "SF-nonlinear-saturation",
  # Paper-replicating settings
  "Danaher", "Danaher-full", "SSJGL-paper", "PSV", "Shaddox", "DINGO-paper"
)
all_methods <- c("JGL_BIC", "JGL_CV", "JGL_BIC_GGL", "JGL_CV_GGL",
                 "SSJGL", "SSJGL_GGL",
                 "DiNGO",
                 "multiGGM", "multiGGM_short", "multiGGM_long",
                 "JRF",
                 "InvCov")
# Note: DiNGO is currently deferred from the production roster
# while we resolve the poor diff_AUC at p=250 (see CLAUDE.md "Open
# investigation"). It is still in all_methods so it can be invoked
# explicitly via --methods DiNGO.

if (test_mode) {
  cat("=== TEST MODE ===\n\n")
  if (is.null(method_names))
    method_names <- c("JGL_BIC", "JGL_BIC_GGL", "SSJGL", "SSJGL_GGL",
                       "multiGGM", "JRF")
  if (is.null(settings_req)) settings_req <- c("SF-p50n1000")
  if (is.na(target_reps) || target_reps > 5) target_reps <- 2L
}

if (!list_configs) {
  if (is.null(method_names))
    stop("Must specify --methods (e.g., --methods \"JGL_BIC,SSJGL\"). Use --test for defaults.")
  if (is.null(settings_req)) settings_req <- all_settings

  bad <- setdiff(method_names, all_methods)
  if (length(bad) > 0) {
    stop(sprintf("Unknown method(s): %s\nValid: %s",
                 paste(bad, collapse = ", "),
                 paste(all_methods, collapse = ", ")))
  }
}

#------------------------------------------------------
# --list-configs mode: report what's on disk and exit
#------------------------------------------------------

raw_root <- "results/raw_fits"

if (list_configs) {
  for (m in method_names %||% all_methods) {
    for (s in settings_req %||% all_settings) {
      mdir <- file.path(raw_root, m, s)
      if (!dir.exists(mdir)) next
      configs <- list.dirs(mdir, recursive = FALSE)
      if (length(configs) == 0) next
      cat(sprintf("\n=== %s / %s ===\n", m, s))
      for (cfg in configs) {
        reps <- length(list.files(cfg, pattern = "^rep_\\d+\\.rds$"))
        cfg_file <- file.path(cfg, "config.rds")
        cfg_short <- sub(".*_config_", "", basename(cfg))
        note <- ""
        if (file.exists(cfg_file)) {
          params <- tryCatch(readRDS(cfg_file), error = function(e) NULL)
          current_hash <- method_config_hash(m)
          if (!is.null(params) && identical(cfg_short, current_hash)) {
            note <- " (current default)"
          }
        }
        cat(sprintf("  _config_%s  %3d reps%s\n", cfg_short, reps, note))
      }
    }
  }
  quit(save = "no")
}

#------------------------------------------------------
# Load the rep manifest (seeds) produced by 05a
#------------------------------------------------------

manifest_file <- "gold_standards/manifest.rds"
if (!file.exists(manifest_file)) {
  stop(sprintf("Manifest not found at %s. Run 05a_generate_gold_standards.R first.",
               manifest_file))
}
manifest <- readRDS(manifest_file)

cat(sprintf("Methods:      %s\n", paste(method_names, collapse = ", ")))
cat(sprintf("Settings:     %s\n", paste(settings_req, collapse = ", ")))
cat(sprintf("Target reps:  %d (start at %d)\n", target_reps, start_rep))
if (force_mode) cat("Force mode:   ON (overwriting existing rep files)\n")

#------------------------------------------------------
# Parallel plan — let user cap the worker count for memory-
# heavy methods like multiGGM on p=250. Default is min of
# (availableCores-1, 6) for safety.
#------------------------------------------------------

default_workers <- max(1, min(availableCores() - 1, 6))
n_workers <- if (!is.na(n_workers_arg)) as.integer(n_workers_arg) else default_workers
cat(sprintf("Workers:      %d\n\n", n_workers))
plan(multisession, workers = n_workers)

dir.create(raw_root, showWarnings = FALSE, recursive = TRUE)

#------------------------------------------------------
# Main loop: for each (method, setting), identify which reps
# need fitting, kick off future_pwalk so each worker writes its
# own RDS atomically.
#------------------------------------------------------

.pad <- function(x) sprintf("rep_%04d.rds", x)

start_time <- Sys.time()

for (method_name in method_names) {
  params_m <- method_params(method_name)
  cfg_hash <- method_config_hash(method_name, params_m)

  # ----------------------------------------------------------
  # Build ONE flat job list across all settings for this method.
  # The previous structure dispatched one future_pwalk per
  # (method, setting), which left workers idle while waiting for
  # the slowest rep in each batch — most acutely for multiGGM and
  # JRF where a single rep can be 5-15 minutes. Flattening the
  # dispatch lets a worker that finishes setting A's last rep
  # immediately pick up setting B's first rep, eliminating the
  # tail-end idle period at every (method, setting) boundary.
  #
  # Per-method (not cross-method) flattening keeps method-specific
  # --n-workers limits intact (multiGGM at 4 for memory, others
  # higher). Run different methods in separate invocations to
  # preserve those.
  # ----------------------------------------------------------
  jobs_rep_id     <- integer(0)
  jobs_seed_gold  <- numeric(0)
  jobs_seed_data  <- numeric(0)
  jobs_setting    <- character(0)
  jobs_gold_dir   <- character(0)
  jobs_out_dir    <- character(0)
  jobs_params_set <- list()
  out_dirs_seen   <- character(0)   # for the post-dispatch summary

  for (setting in settings_req) {
    gold_dir <- file.path("gold_standards", setting)
    if (!dir.exists(gold_dir)) {
      cat(sprintf("[%s / %s] No gold standards at %s — skipping\n",
                  method_name, setting, gold_dir))
      next
    }

    # Which reps have gold standards?
    gs_files <- list.files(gold_dir, pattern = "^rep_\\d+\\.rds$")
    available_reps <- as.integer(gsub("rep_0*(\\d+)\\.rds", "\\1", gs_files))
    eligible <- available_reps[available_reps >= start_rep &
                                 available_reps <= target_reps]
    if (length(eligible) == 0) {
      cat(sprintf("[%s / %s] 0 gold standards in range [%d, %d] — skipping (run 05a first)\n",
                  method_name, setting, start_rep, target_reps))
      next
    }

    # Where to save raw fits for this (method, setting, config)
    out_dir <- file.path(raw_root, method_name, setting,
                         paste0("_config_", cfg_hash))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    # Always (re)write the human-readable config.rds record
    saveRDS(params_m, file.path(out_dir, "config.rds"))

    # Resume: skip any reps whose output file already exists
    existing <- list.files(out_dir, pattern = "^rep_\\d+\\.rds$")
    existing_ids <- as.integer(
      gsub("rep_0*(\\d+)\\.rds", "\\1",
           setdiff(existing, "config.rds")))
    to_run <- if (force_mode) eligible else setdiff(eligible, existing_ids)
    to_run <- sort(to_run)

    cat(sprintf("[%s / %s] cfg=%s | %3d existing, %3d queued, target %d\n",
                method_name, setting, cfg_hash,
                length(existing_ids), length(to_run), target_reps))

    out_dirs_seen <- c(out_dirs_seen,
                       setNames(out_dir, setting))

    if (length(to_run) == 0) next

    # Look up seed_data per rep from the manifest
    man_sub <- manifest[manifest$setting_label == setting, ]
    man_sub <- man_sub[match(to_run, man_sub$rep_id), ]
    if (any(is.na(man_sub$seed_data))) {
      warning(sprintf("[%s / %s] %d reps missing from manifest; skipping those",
                      method_name, setting, sum(is.na(man_sub$seed_data))))
      man_sub <- man_sub[!is.na(man_sub$seed_data), ]
    }
    if (nrow(man_sub) == 0) next

    params_setting <- get_setting_params(setting)

    # Append this setting's reps to the flat job list
    n_new <- nrow(man_sub)
    jobs_rep_id     <- c(jobs_rep_id,    man_sub$rep_id)
    jobs_seed_gold  <- c(jobs_seed_gold, man_sub$seed_gold)
    jobs_seed_data  <- c(jobs_seed_data, man_sub$seed_data)
    jobs_setting    <- c(jobs_setting,   rep(setting,  n_new))
    jobs_gold_dir   <- c(jobs_gold_dir,  rep(gold_dir, n_new))
    jobs_out_dir    <- c(jobs_out_dir,   rep(out_dir,  n_new))
    jobs_params_set <- c(jobs_params_set, replicate(n_new, params_setting,
                                                    simplify = FALSE))
  }

  total_jobs <- length(jobs_rep_id)
  if (total_jobs == 0) {
    cat(sprintf("[%s] All settings up-to-date — nothing to dispatch\n\n",
                method_name))
    next
  }

  n_settings <- length(unique(jobs_setting))
  cat(sprintf("\n[%s] Dispatching %d total reps across %d settings as one queue (workers fill the tail across settings)\n\n",
              method_name, total_jobs, n_settings))

  # ----------------------------------------------------------
  # Single flat future_pwalk. Per-rep state (setting, gold_dir,
  # out_dir, params_setting) is passed as columns rather than
  # captured from the enclosing environment.
  # ----------------------------------------------------------
  future_pwalk(
    list(rep_id         = jobs_rep_id,
         seed_gold      = jobs_seed_gold,
         seed_data      = jobs_seed_data,
         setting        = jobs_setting,
         gold_dir       = jobs_gold_dir,
         out_dir        = jobs_out_dir,
         params_setting = jobs_params_set),
    function(rep_id, seed_gold, seed_data, setting, gold_dir,
             out_dir, params_setting) {
      outfile <- file.path(out_dir, .pad(rep_id))
      # Last-chance skip in case another worker got to it
      if (!force_mode && file.exists(outfile)) return(invisible(NULL))

      gold_file <- file.path(gold_dir, .pad(rep_id))
      gold <- readRDS(gold_file)

      set.seed(seed_data)
      dat <- gen_dat(gold, params_setting$n,
                     transform = params_setting$transform)
      dat_list <- list(dat$dat1, dat$dat2)

      fit <- fit_method(method_name, dat_list, params = params_m)

      record <- list(
        rep_id      = rep_id,
        method      = method_name,
        setting     = setting,
        config_hash = cfg_hash,
        seed_gold   = seed_gold,
        seed_data   = seed_data,
        params      = params_m,
        timestamp   = Sys.time(),
        r_version   = R.version.string,
        status      = fit$status,
        runtime_sec = fit$runtime_sec,
        fit         = fit
      )

      # Atomic write: write to a temp file, then rename
      tmpfile <- paste0(outfile, ".tmp")
      saveRDS(record, tmpfile)
      file.rename(tmpfile, outfile)

      # Free memory before the next rep runs on this worker
      rm(record, fit, gold, dat, dat_list); gc(verbose = FALSE)
      invisible(NULL)
    },
    .options = furrr_options(
      seed = TRUE,
      # scheduling = Inf: dispatch each job individually rather than
      # pre-chunking the queue. The default scheduling = 1.0 splits
      # the job list into n_workers consecutive chunks at dispatch
      # time, so workers that get assigned chunks of fast settings
      # (or fewer slow reps) finish early and sit idle for the rest
      # of the run while workers with heavier chunks grind through.
      # With scheduling = Inf, workers pull the next available job
      # from the queue as soon as they finish their current one --
      # the dynamic load balancing that the flat-dispatch redesign
      # was meant to enable in the first place. Per-job dispatch
      # overhead is in the millisecond range, negligible for fits
      # that take minutes to hours.
      scheduling = Inf,
      globals = c("fit_method", "score_method", "run_method",
                  "gen_dat", "get_setting_params", "method_params",
                  "method_config_hash",
                  "run_jgl", "tune_jgl", "cv_jgl", "run_ssjgl",
                  "run_dingo", "run_multiggm", "run_jrf",
                  "vBIC", "empirical_dist",
                  "method_name", "params_m", "cfg_hash",
                  "force_mode", ".pad"),
      packages = c("igraph", "MASS", "JGL", "spikeyglass",
                   "multiGGM", "iDINGO", "iJRFNet", "BDgraph",
                   "digest")
    ),
    .progress = TRUE
  )

  # Post-dispatch summary, one line per setting
  for (i in seq_along(out_dirs_seen)) {
    s <- names(out_dirs_seen)[i]
    od <- out_dirs_seen[[i]]
    n_done <- length(list.files(od, pattern = "^rep_\\d+\\.rds$"))
    cat(sprintf("[%s / %s] done — %d rep files now on disk (target %d)\n",
                method_name, s, n_done, target_reps))
  }
  cat("\n")
}

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("\nFit phase complete in %.1f min\n", as.numeric(elapsed)))
cat(sprintf("Next: Rscript 05d_score_results.R --methods \"...\" --settings \"...\"\n"))
