#!/usr/bin/env Rscript
# Preprocess results/simulation_results.RData into a compact internal
# lookup table shipped as R/sysdata.rda (object: howbig_grid).
#
# One row per (method, p, n, perturbation) cell with mean / sd / q10 / q90
# and rep count for the three metrics howbig serves.

load("results/simulation_results.RData")   # -> results
d <- results

parse_p <- function(s) suppressWarnings(as.integer(sub(".*p(\\d+)n.*", "\\1", s)))
d$p <- parse_p(d$setting)

# keep only the SF- BA-family sweep/anchor cells (drop paper-replication
# generators Danaher/Shaddox which have no parseable p and a different DGP)
d <- d[!is.na(d$p), ]

# user-facing perturbation label <- grid perturbation
# structural = Deletion (A); magnitude = Activation (C); both = Deletion + Activation (B)
pert_map <- c("Deletion" = "structural",
              "Activation" = "magnitude",
              "Deletion + Activation" = "both")
d$pert <- pert_map[d$perturbation]

metrics <- c(sensitivity = "diff_TPR_inference",
             achieved_fdr = "diff_FDR_inference",
             aupr         = "diff_AUPR")

key <- interaction(d$method, d$p, d$n, d$pert, drop = TRUE, sep = "\u0001")
split_idx <- split(seq_len(nrow(d)), key)

rows <- lapply(split_idx, function(ix) {
  sub <- d[ix, ]
  out <- data.frame(method = sub$method[1], p = sub$p[1], n = sub$n[1],
                    perturbation = sub$pert[1], stringsAsFactors = FALSE)
  for (mn in names(metrics)) {
    v <- sub[[metrics[mn]]]
    v <- v[is.finite(v)]
    out[[paste0(mn, "_mean")]] <- if (length(v)) mean(v) else NA_real_
    out[[paste0(mn, "_sd")]]   <- if (length(v)) sd(v)   else NA_real_
    out[[paste0(mn, "_q10")]]  <- if (length(v)) as.numeric(quantile(v, .10)) else NA_real_
    out[[paste0(mn, "_q90")]]  <- if (length(v)) as.numeric(quantile(v, .90)) else NA_real_
  }
  out$n_reps <- nrow(sub)
  out
})
howbig_grid <- do.call(rbind, rows)
rownames(howbig_grid) <- NULL
howbig_grid <- howbig_grid[order(howbig_grid$method, howbig_grid$perturbation,
                                 howbig_grid$p, howbig_grid$n), ]

# provenance metadata travels with the object
attr(howbig_grid, "howbig_meta") <- list(
  source_family = "SF- (BA-tree, m_pa=1) — INTERIM; ship target is SF2-",
  fdr_target = 0.05,
  metric_cols = as.list(metrics),
  built = as.character(Sys.Date())
)

# Data-generating family the grid was simulated under. The on-demand
# simulator reads this so simulated answers describe the same networks as
# grid answers. When the SF2- grid replaces this one, set m_pa = 2 here.
howbig_family <- list(
  topology     = "scale_free",
  m_pa         = 1L,        # SF- = BA tree; SF2- = 2
  p_block      = 25L,
  wt_method    = "empirical",
  scale_factor = 1.5,
  diff_threshold = 0.01     # |Theta_g1 - Theta_g2| > this => true differential
)

# Empirical partial-correlation weight distribution (histogram: mids,
# density) used to draw edge weights; from the DiNA simulation project.
howbig_weights <- read.csv("data-raw/mxDist.csv")

# Per-method runtime model for simulation pre-flight estimates:
# log(seconds per fit) ~ log(p) + log(n), fit on the production timings.
# InvCov is effectively instant, so it gets a flat floor.
rt <- d[is.finite(d$runtime_sec) & d$runtime_sec > 0, ]
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
save(howbig_grid, howbig_family, howbig_weights, howbig_runtime,
     file = "R/sysdata.rda", compress = "xz")

cat("Wrote R/sysdata.rda:", nrow(howbig_grid), "cells,",
    length(unique(howbig_grid$method)), "methods\n")
cat("perturbations:", paste(unique(howbig_grid$perturbation), collapse=", "), "\n")
cat("file size:", file.info("R/sysdata.rda")$size, "bytes\n")
cat("\nsample rows (JGL_BIC, structural):\n")
print(head(howbig_grid[howbig_grid$method=="JGL_BIC" &
                       howbig_grid$perturbation=="structural",
                       c("p","n","sensitivity_mean","achieved_fdr_mean","aupr_mean","n_reps")], 12),
      row.names = FALSE, digits = 3)
