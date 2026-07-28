#------------------------------------------------------
# Setup Script: Install and Load All Required Packages
#
# Run this once on any new machine before using the simulation
# pipeline (05a, 05b, 05c, 06, reports, etc.).
#
# Usage (from within the dina_simulation project directory):
#   Rscript 00_setup.R
#
# Or interactively:
#   source("00_setup.R")
#
# This script:
#   1. Installs 'remotes' if missing (needed for GitHub packages)
#   2. Installs all CRAN packages used by the pipeline
#   3. Installs GitHub-only packages (spikeyglass, multiGGM, iJRFNet)
#   4. Loads and verifies all packages
#   5. Prints a summary showing which packages are ready
#
# If any install fails, the script continues and reports the
# failures at the end. You may need to install system dependencies
# (e.g., gfortran, gmp, mpfr) for some packages.
#------------------------------------------------------

cat("\n===========================================================\n")
cat("  DiNA Simulation — Environment Setup\n")
cat("===========================================================\n\n")

cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("Platform: %s\n\n", R.version$platform))

#------------------------------------------------------
# Helper functions
#------------------------------------------------------

install_cran_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s from CRAN...\n", pkg))
    tryCatch({
      install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
      return(requireNamespace(pkg, quietly = TRUE))
    }, error = function(e) {
      cat(sprintf("    FAILED: %s\n", conditionMessage(e)))
      return(FALSE)
    })
  } else {
    return(TRUE)
  }
}

install_github_if_missing <- function(pkg, repo) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s from GitHub (%s)...\n", pkg, repo))
    tryCatch({
      remotes::install_github(repo, upgrade = "never", quiet = TRUE)
      return(requireNamespace(pkg, quietly = TRUE))
    }, error = function(e) {
      cat(sprintf("    FAILED: %s\n", conditionMessage(e)))
      return(FALSE)
    })
  } else {
    return(TRUE)
  }
}

#------------------------------------------------------
# Step 1: Install 'remotes' (needed for GitHub packages)
#------------------------------------------------------

cat("Step 1: Bootstrap 'remotes' package\n")
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
cat("  remotes: OK\n\n")

#------------------------------------------------------
# Step 2: Install CRAN packages
#------------------------------------------------------

cat("Step 2: Install CRAN packages\n")

cran_pkgs <- c(
  # Core simulation
  "igraph",      # network generation
  "MASS",        # multivariate normal sampling
  "future",      # parallel backend
  "furrr",       # purrr + future
  "BDgraph",     # G-Wishart sampling for SSJGL-paper setting
  # Method: DiNGO (available on CRAN)
  "iDINGO",      # differential network analysis in genomics
  # Data wrangling & analysis
  "tidyverse",   # dplyr, tidyr, ggplot2, readr, etc.
  "reshape2",    # melt/dcast for plotting
  "gridExtra",   # multi-panel plots
  # Reports
  "knitr",
  "kableExtra",
  "rmarkdown"
)

cran_status <- sapply(cran_pkgs, install_cran_if_missing)
cat(sprintf("\n  CRAN packages installed: %d/%d\n\n",
            sum(cran_status), length(cran_status)))

#------------------------------------------------------
# Step 3: Install GitHub-only packages
#------------------------------------------------------

cat("Step 3: Install GitHub-only packages\n")

github_pkgs <- list(
  # JGL: joint graphical lasso (install from GitHub for latest)
  JGL         = "patrickjdanaher/JGL-package",
  # SSJGL: spike-and-slab joint graphical lasso
  spikeyglass = "mljaniczek/spikeyglass",
  # multiGGM: Bayesian SSVS joint GGM
  multiGGM    = "mljaniczek/multiGGMr",
  # JRF: joint random forest (repo name is 'iJRF', package name is 'iJRFNet')
  iJRFNet     = "petraf01/iJRF"
)

gh_status <- mapply(install_github_if_missing,
                     names(github_pkgs), github_pkgs)
cat(sprintf("\n  GitHub packages installed: %d/%d\n\n",
            sum(gh_status), length(gh_status)))

#------------------------------------------------------
# Step 4: Load and verify all packages
#------------------------------------------------------

cat("Step 4: Load and verify packages\n")

all_pkgs <- c(cran_pkgs, names(github_pkgs))
load_status <- sapply(all_pkgs, function(pkg) {
  ok <- suppressPackageStartupMessages(
    tryCatch({
      requireNamespace(pkg, quietly = TRUE)
    }, error = function(e) FALSE)
  )
  status <- if (ok) "OK" else "FAILED"
  version <- if (ok) as.character(packageVersion(pkg)) else "--"
  cat(sprintf("  %-15s %-8s v%s\n", pkg, status, version))
  ok
})

#------------------------------------------------------
# Step 5: Summary
#------------------------------------------------------

cat("\n===========================================================\n")
cat("  Setup Summary\n")
cat("===========================================================\n")

n_ok <- sum(load_status)
n_total <- length(load_status)
cat(sprintf("  %d/%d packages ready\n", n_ok, n_total))

if (n_ok == n_total) {
  cat("\n  All packages installed and loaded successfully.\n")
  cat("  You can now run the simulation pipeline:\n")
  cat("    Rscript 05a_generate_gold_standards.R --test\n")
  cat("    Rscript 05b_run_methods.R --test --methods 'JGL_BIC,SSJGL,multiGGM,DiNGO,JRF'\n")
  cat("    Rscript 05c_merge_results.R\n")
} else {
  failed <- names(load_status)[!load_status]
  cat("\n  The following packages failed to install:\n")
  for (pkg in failed) cat(sprintf("    - %s\n", pkg))
  cat("\n  Troubleshooting tips:\n")
  cat("    * BDgraph needs gmp and mpfr system libraries\n")
  cat("      (macOS: brew install gmp mpfr; ubuntu: apt install libgmp-dev libmpfr-dev)\n")
  cat("    * iDINGO may need 'glasso' and 'scales' installed first\n")
  cat("    * For GitHub packages, ensure Rtools (Windows) or Xcode CLT (macOS)\n")
  cat("    * Set GITHUB_PAT env var if you hit rate limits\n")
}
cat("\n")
