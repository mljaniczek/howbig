#------------------------------------------------------
# Merge Method Results
# Collects per-method result files into a single
# unified dataframe for use by 06_summarize_results.R.
#
# Usage:
#   Rscript 05c_merge_results.R
#------------------------------------------------------

#------------------------------------------------------
# Discover result files
#------------------------------------------------------

result_dir <- "results/method_results"

if (!dir.exists(result_dir)) {
  stop("No results directory found at ", result_dir,
       "\nRun 05b_run_methods.R first.")
}

# Pick up metrics files but SKIP the *.meta.rds side-cars written
# by 05d_score_results.R (those record which config + scoring
# options produced each metrics file; they are not part of the
# aggregated results table).
result_files <- list.files(result_dir, pattern = "\\.rds$",
                           recursive = TRUE, full.names = TRUE)
result_files <- result_files[!grepl("\\.meta\\.rds$", result_files)]

if (length(result_files) == 0) {
  stop("No result files found in ", result_dir,
       "\nHint: run 05b_run_methods.R to fit, then 05d_score_results.R to score.")
}

cat(sprintf("Found %d result files:\n", length(result_files)))
for (f in result_files) {
  cat(sprintf("  %s\n", sub("results/method_results/", "", f)))
}

#------------------------------------------------------
# Read and combine
#------------------------------------------------------

results_list <- lapply(result_files, readRDS)
# Use bind_rows for robustness: paper settings add a 'generator' column
# that standard settings lack; bind_rows fills missing columns with NA
results <- dplyr::bind_rows(results_list)

#------------------------------------------------------
# Translate perturbation codes to human-readable labels
#------------------------------------------------------
# The simulation pipeline writes single-letter codes for the
# three perturbation types. We translate them here so downstream
# consumers (06_results_tables.R, ad-hoc analysis, the manuscript)
# get the descriptive labels directly.
#
# Convention:
#   A = Deletion              (within-module-1 edges deleted)
#   B = Deletion + Activation (module 1 deleted, module 2 rescaled)
#   C = Activation            (within-module-1 edges rescaled)
#
# The mapping is idempotent: rows that already have full labels
# (e.g., from a re-run after the labels were already translated)
# pass through unchanged.

if ("perturbation" %in% colnames(results)) {
  pert_map <- c(
    "A" = "Deletion",
    "B" = "Deletion + Activation",
    "C" = "Activation"
  )
  unknown_codes <- setdiff(unique(results$perturbation),
                            c(names(pert_map), unname(pert_map), NA))
  if (length(unknown_codes) > 0) {
    warning(sprintf(
      "Unrecognized perturbation codes (left as-is): %s",
      paste(unknown_codes, collapse = ", ")
    ))
  }
  results$perturbation <- ifelse(
    results$perturbation %in% names(pert_map),
    pert_map[results$perturbation],
    results$perturbation
  )
}

cat(sprintf("\nTotal rows: %d\n", nrow(results)))
cat(sprintf("Methods: %s\n", paste(sort(unique(results$method)), collapse = ", ")))
cat(sprintf("Settings: %s\n", paste(sort(unique(results$setting)), collapse = ", ")))
cat(sprintf("Perturbations: %s\n", paste(sort(unique(results$perturbation)), collapse = ", ")))

# Status summary
status_tab <- table(results$status)
cat("\nStatus summary:\n")
for (s in names(status_tab)) {
  cat(sprintf("  %s: %d\n", s, status_tab[s]))
}

#------------------------------------------------------
# Save merged results
#------------------------------------------------------

dir.create("results", showWarnings = FALSE)

# Save in same format that 06_summarize_results.R expects
save(results, file = "results/simulation_results.RData")
cat(sprintf("\nSaved to results/simulation_results.RData\n"))

# Also save timestamped RDS backup
timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
rds_file <- sprintf("results/simulation_results_%s.rds", timestamp_str)
saveRDS(results, rds_file)
cat(sprintf("Saved to %s\n", rds_file))

cat("\nRun 06_summarize_results.R to generate manuscript tables.\n")
