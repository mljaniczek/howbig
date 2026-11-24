#!/usr/bin/env Rscript
# Build the internal lookup tables shipped in R/sysdata.rda from the DiNA
# production results. Run from the package root:
#   Rscript data-raw/build_sysdata.R
#
# Objects written:
#   howbig_grid    one row per (method, scenario, perturbation, p, n1, n2) cell
#   howbig_method_info method metadata (inference route, whether simulatable)
#   howbig_family  data-generating family the simulator uses (matches grid)
#   howbig_weights empirical edge-weight distribution (mxDist.csv)
#   howbig_runtime per-method log-linear model of fit time, for ETAs

SOURCE <- "results/simulation_results_20261117_112116.rds"
d <- readRDS(SOURCE)
d <- d[d$status == "success", ]

# v1 methods (agreed Nov 2026): six core methods plus three grid-only ones.
# Excluded: JRF (sensitivity scored in under half the cells, breaks under
# unequal groups), DTRACE/PNJGL (no FDR-controlled calls, so no power),
# JGL_CV*/multiGGM_long (run in only 2-3 settings).
howbig_method_info <- data.frame(
  method = c("JGL_BIC", "JGL_BIC_GGL", "SSJGL", "SSJGL_GGL", "multiGGM",
             "InvCov", "mGHS", "DiffNetFDR", "ONDSA"),
  description = c(
    "Joint graphical lasso, fused penalty, BIC-tuned",
    "Joint graphical lasso, group penalty, BIC-tuned",
    "Spike-and-slab joint graphical lasso, fused penalty",
    "Spike-and-slab joint graphical lasso, group penalty",
    "Bayesian joint GGM (stochastic search variable selection)",
    "Unregularized inverse covariance (baseline)",
    "Bayesian multiple graphical horseshoe",
    "Partial-correlation difference test (DiffNetFDR)",
    "Debiased (FastGGM-based) differential network test"),
  inference = c(rep("permutation null + BH", 4), "posterior FDR",
                "permutation null + BH", "posterior FDR",
                "built-in FDR procedure", "built-in FDR procedure"),
  simulate = c(rep(TRUE, 6), rep(FALSE, 3)),
  stringsAsFactors = FALSE
)
d <- d[d$method %in% howbig_method_info$method, ]

# --- Parse setting labels into design + scenario ---------------------------
# Labels look like <family><perturbation>p<P>n<N>[-<variant>], e.g.
# "SF-p100n1000", "SFCp100n1000", "SF2-p100n1000", "SF-p100n1000-zi25".
# NOTE: the results' `n` column is the AVERAGE group size. Unequal-group
# settings must take n1/n2 from the DiNA registry
# (dina_simulation/00_helper_functions.R, get_setting_params()), otherwise
# they silently merge into the balanced cell.
UNEQUAL <- list(unbal3 = c(500L, 1500L))   # "SF-p100n1000-unbal3"
VARIANT_SCENARIO <- c(zi05 = "censored_5", zi25 = "censored_25",
                      zi50 = "censored_50", uniform = "uniform_weights")
FAMILY_SCENARIO <- c(SF2 = "scale_free", ER = "erdos_renyi",
                     AR1 = "ar1_banded")
PERT <- c("-" = "structural", "+" = "both", "C" = "magnitude", "0" = "null")

rx <- "^(SF2|SF|ER|AR1)([-+C0])p([0-9]+)n([0-9]+)(-(.+))?$"
settings <- unique(d$setting)
bad <- settings[!grepl(rx, settings)]
if (length(bad)) stop("Unparsed setting labels: ", paste(bad, collapse = ", "))

parsed <- do.call(rbind, lapply(settings, function(s) {
  m <- regmatches(s, regexec(rx, s))[[1]]
  family <- m[2]; variant <- m[7]; n <- as.integer(m[5])
  scenario <- if (family != "SF") FAMILY_SCENARIO[[family]]
              else if (nzchar(variant) && variant %in% names(VARIANT_SCENARIO))
                VARIANT_SCENARIO[[variant]]
              else "default"
  if (nzchar(variant) && !variant %in% c(names(VARIANT_SCENARIO), names(UNEQUAL)))
    stop("Unknown setting variant: ", s)
  nn <- if (nzchar(variant) && variant %in% names(UNEQUAL)) UNEQUAL[[variant]]
        else c(n, n)
  data.frame(setting = s, scenario = scenario, perturbation = PERT[[m[3]]],
             p = as.integer(m[4]), n1 = nn[1], n2 = nn[2],
             stringsAsFactors = FALSE)
}))
# The global-null setting has no true differences: no power to report.
parsed <- parsed[parsed$perturbation != "null", ]
d <- merge(d[, setdiff(names(d), c("p", "perturbation"))], parsed, by = "setting")

metrics <- c(sensitivity = "diff_TPR_inference",
             achieved_fdr = "diff_FDR_inference",
             aupr = "diff_AUPR")
keys <- c("method", "scenario", "perturbation", "p", "n1", "n2")
split_idx <- split(seq_len(nrow(d)), interaction(d[keys], drop = TRUE, sep = "|"))

howbig_grid <- do.call(rbind, lapply(split_idx, function(ix) {
  sub <- d[ix, ]
  out <- sub[1, keys]
  for (mn in names(metrics)) {
    v <- sub[[metrics[mn]]]
    v <- v[is.finite(v)]
    out[[paste0(mn, "_mean")]] <- if (length(v)) mean(v) else NA_real_
    out[[paste0(mn, "_sd")]]   <- if (length(v) > 1) sd(v) else NA_real_
    out[[paste0(mn, "_q10")]]  <- if (length(v)) unname(quantile(v, .10)) else NA_real_
    out[[paste0(mn, "_q90")]]  <- if (length(v)) unname(quantile(v, .90)) else NA_real_
  }
  out$n_diff_edges <- mean(sub$diff_edges_true, na.rm = TRUE)
  out$n_pairs <- sub$p[1] * (sub$p[1] - 1) / 2
  out$n_reps <- nrow(sub)
  out
}))
rownames(howbig_grid) <- NULL
howbig_grid <- howbig_grid[order(howbig_grid$method, howbig_grid$scenario,
                                 howbig_grid$perturbation, howbig_grid$p,
                                 howbig_grid$n1, howbig_grid$n2), ]

attr(howbig_grid, "howbig_meta") <- list(
  source = basename(SOURCE),
  source_family = "SF1: Barabasi-Albert tree blocks (m_pa = 1)",
  fdr_target = 0.05,
  metric_cols = as.list(metrics),
  built = as.character(Sys.Date())
)

# Data-generating family for the on-demand simulator; must match the grid's
# "default" scenario so simulated and grid answers describe the same networks.
howbig_family <- list(
  topology       = "scale_free",
  m_pa           = 1L,        # SF1 = BA tree
  p_block        = 25L,
  wt_method      = "empirical",
  scale_factor   = 1.5,
  diff_threshold = 0.01       # |Theta_g1 - Theta_g2| > this => true differential
)

# Empirical partial-correlation weight distribution (histogram: mids,
# density) used to draw edge weights; from the DiNA simulation project.
howbig_weights <- read.csv("data-raw/mxDist.csv")

# Per-method runtime model for simulation pre-flight estimates:
# log(seconds per fit) ~ log(p) + log(n), on default-scenario timings.
rt <- d[d$scenario == "default" & is.finite(d$runtime_sec) & d$runtime_sec > 0, ]
rt$n <- (rt$n1 + rt$n2) / 2
howbig_runtime <- do.call(rbind, lapply(split(rt, rt$method), function(s) {
  cf <- if (length(unique(s$p)) > 1 && length(unique(s$n)) > 1) {
    coef(lm(log(runtime_sec) ~ log(p) + log(n), data = s))
  } else {
    c(log(stats::median(s$runtime_sec)), 0, 0)
  }
  data.frame(method = s$method[1], b0 = cf[[1]], b_logp = cf[[2]],
             b_logn = cf[[3]], row.names = NULL)
}))

dir.create("R", showWarnings = FALSE)
save(howbig_grid, howbig_method_info, howbig_family, howbig_weights, howbig_runtime,
     file = "R/sysdata.rda", compress = "xz")

cat("Wrote R/sysdata.rda:", nrow(howbig_grid), "cells,",
    length(unique(howbig_grid$method)), "methods,",
    file.info("R/sysdata.rda")$size, "bytes\n")
print(table(howbig_grid$scenario, howbig_grid$perturbation))
