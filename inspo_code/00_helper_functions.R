# this function from Koyejo Lab Github https://github.com/koyejo-lab/JointGraphicalLasso
# input is list of matrices same as you would put into `JGL` function, as well as results of `JGL`
# output is AIC for particular estimated graph.
# use this function when tuning lambdas
vBIC <- function(X, est_graph, thr=0.001){
  num <- length(X)
  BIC_acc <- 0.
  for(i in 1:num){
    
    data_num <- dim(X[[i]])[1]
    sample_cov <- cov(X[[i]], X[[i]])
    tr_sum <- sum(diag(sample_cov %*% est_graph$theta[[i]]))
    
    log_det <- determinant(est_graph$theta[[i]], logarithm = TRUE)$modulus[1][1]
    
    E <- sum(sum(abs(est_graph$theta[[i]]) >= thr))
    BIC_acc <- BIC_acc + (tr_sum - log_det) + (log(data_num)*E/data_num)
  }
  return(BIC_acc)
}

# function to make the plots I want below. input data, penalty type and lambdas.
# made function to shorten code later.
lambda_plots <- function(data, penalty = "group", lambda1, lambda2, ...){
  res <- JGL(Y=data,
             penalty=penalty,
             lambda1=lambda1,
             lambda2=lambda2, 
             return.whole.theta=TRUE)
  test1 <- vBIC(data, res, thr=0.0001)
  resgraph <- graph_from_adjacency_matrix(
    -cov2cor(res$theta[[1]]),
    weighted = T,
    mode = "undirected",
    diag = FALSE
  )
  resplot <- plot_jgl(resgraph,
                      multiplier = 3,
                      sub = paste("AIC", round(test1, 2)),
                      ...
  )
}

# NOTE: assignRealDataEdgeWeights() and boost() have been superseded by
# assign_edge_weights() and build_precision_matrix() in 01_generate_gold_standard_networks.R

#------------------------------------------------------
# Simulation settings lookup
# Maps setting labels to generate_gold_standard() args.
# Shared by 05a_generate_gold_standards.R and
# 05b_run_methods.R.
#------------------------------------------------------

get_setting_params <- function(setting_label) {
  # Setting-naming convention (this revision):
  #   SF{±}p{P}n{N}  where
  #     {±}       A (delete Block 1 edges)  = "-"
  #               B (delete + rescale wts) = "+"
  #     p{P}      feature count
  #     n{N}      per-group sample size
  #   Topology is always scale-free in this family. Block structure
  #   holds p_block = 25 fixed; B = ceiling(p / p_block). Within-block
  #   density is 0.15. Between-block probability scales as
  #     between_prob = 0.01 * (5 / B)
  #   so the expected *number* of between-block edges stays comparable
  #   across p (≈31 at p=125). p=20 uses B=2, p_block=10 as a small-p
  #   special case that still supports the delete-a-block perturbation.
  #
  # Shared defaults across the legacy SF- family (Barabási–Albert
  # *tree* topology — m_pa = 1, ~8% within-block density). The new
  # SF2- family below overrides m_pa = 2 to produce proper scale-free
  # blocks (~15.7% within-block density, hubs and cycles); the AR1-
  # family overrides `topology` to an entirely different,
  # non-modular, banded structure.
  #
  # `within_density_random` is the renamed field formerly known as
  # `within_density`. It controls the Erdős–Rényi edge probability
  # for `topology = "random"` and is IGNORED for the scale-free and
  # AR(1) branches. Renamed so it's obvious it doesn't apply to all
  # topologies. See NOTES_network_density.md for the discovery.
  base <- list(
    wt_method = "empirical", epsilon = 0.1,
    precision_method = "row_sum_normalize", d_scale = NULL,
    scale_factor = 1.5, m_pa = 1,
    within_density_random = 0.15,
    topology = "scale_free",
    # Gold-standard mode:
    #   "fresh" (default): each replicate draws a new gold standard.
    #     Matches 4 of 5 source papers (JGL/SSJGL/multiGGM/DiNGO);
    #     gives the full variance decomposition (between-GS + within-GS).
    #   "fixed": all replicates share the same gold standard, with
    #     only the data being re-drawn. Matches JRF paper. Narrower
    #     variance estimates, but edge identities match across reps
    #     — required for per-edge permutation inference.
    gold_standard_mode = "fresh"
  )

  # SF2- family: same modular structure as the legacy SF- family,
  # but with m_pa = 2 -> proper Barabási–Albert scale-free graphs
  # (mean within-block degree ~3.8, density ~15.7%) instead of the
  # legacy BA-tree (mean degree ~1.9, density ~8%). All other
  # parameters are shared with `base`.
  base_sf2 <- modifyList(base, list(m_pa = 2))

  # AR1- family: AR(1) banded precision matrix (non-modular). The
  # topology branch in generate_gold_standard() routes to a
  # dedicated AR(1) generator; the `B`, `p_block`, `between_prob`,
  # and `m_pa` fields are not used and are set to inert values for
  # schema-consistency with the rest of the registry.
  #   ar1_rho       : geometric-decay coefficient on off-diagonals
  #                   (|i-j| -> rho^|i-j|, clamped to bandwidth)
  #   ar1_bandwidth : how many off-diagonals are non-zero on each
  #                   side of the diagonal. bandwidth = 1 gives a
  #                   true AR(1); larger gives "banded AR"
  base_ar1 <- modifyList(base, list(
    topology      = "ar1_banded",
    B             = 1, p_block = NA, between_prob = 0,
    m_pa          = NA,
    ar1_rho       = 0.4,
    ar1_bandwidth = 1
  ))

  switch(setting_label,
    # ========================================================
    # PRIMARY TUTORIAL SETTINGS (body of 09 report; 4 settings)
    # ========================================================
    # SF-p125n100: mirrors our applied TCGA RNA-seq setup
    # (p ~ 100, n ~ 100). High-dimensional p > n regime.
    "SF-p125n100" = modifyList(base, list(
      setting = "A",
      p = 125, B = 5, p_block = 25,
      between_prob = 0.01, n = 100
    )),
    # SF-p250n1000 / SF+p250n1000: metabolomics-scale setting
    # (p ~ 250, n ~ 1000). Low-dimensional p < n regime.
    "SF-p250n1000" = modifyList(base, list(
      setting = "A",
      p = 250, B = 10, p_block = 25,
      between_prob = 0.005, n = 1000
    )),
    "SF+p250n1000" = modifyList(base, list(
      setting = "B",
      p = 250, B = 10, p_block = 25,
      between_prob = 0.005, n = 1000
    )),
    # SF-p20n1000: small-network data-rich regime. 2 blocks of 10
    # (replacing the old 4-blocks-of-5 structure) so Setting A
    # deletes ~7 differential edges — enough signal for stable
    # per-edge metrics.
    "SF-p20n1000" = modifyList(base, list(
      setting = "A",
      p = 20, B = 2, p_block = 10,
      between_prob = 0.025, n = 1000
    )),

    # ========================================================
    # SAMPLE-SIZE SWEEP: p = 250 fixed, vary n
    # ========================================================
    "SF-p250n250"  = modifyList(base, list(setting = "A", p = 250, B = 10, p_block = 25, between_prob = 0.005, n = 250)),
    "SF-p250n500"  = modifyList(base, list(setting = "A", p = 250, B = 10, p_block = 25, between_prob = 0.005, n = 500)),
    "SF-p250n2000" = modifyList(base, list(setting = "A", p = 250, B = 10, p_block = 25, between_prob = 0.005, n = 2000)),
    "SF-p250n5000" = modifyList(base, list(setting = "A", p = 250, B = 10, p_block = 25, between_prob = 0.005, n = 5000)),
    # (SF-p250n1000 is defined above as a primary and doubles as a
    # sweep data point.)

    # ========================================================
    # FEATURE-SIZE SWEEP: n = 1000 fixed, vary p
    # between_prob = 0.01 * (5 / B) so E[#between-block edges] is
    # stable across p.
    # ========================================================
    "SF-p50n1000"  = modifyList(base, list(setting = "A", p = 50,  B = 2,  p_block = 25, between_prob = 0.025,    n = 1000)),
    "SF-p100n1000" = modifyList(base, list(setting = "A", p = 100, B = 4,  p_block = 25, between_prob = 0.0125,   n = 1000)),
    "SF-p200n1000" = modifyList(base, list(setting = "A", p = 200, B = 8,  p_block = 25, between_prob = 0.00625,  n = 1000)),
    "SF-p300n1000" = modifyList(base, list(setting = "A", p = 300, B = 12, p_block = 25, between_prob = 0.004167, n = 1000)),
    "SF-p500n1000" = modifyList(base, list(setting = "A", p = 500, B = 20, p_block = 25, between_prob = 0.0025,   n = 1000)),

    # ========================================================
    # ADDITIVE PERTURBATION (Setting C)
    # g2 = g1 + new within-Block-1 edges (no deletions, no
    # rescaling). Number of added edges = within-Block-1 edge
    # count of g1 so signal magnitude is comparable to Setting A.
    # See apply_perturbation() Setting C branch.
    # ========================================================
    "SFCp250n1000" = modifyList(base, list(
      setting = "C",
      p = 250, B = 10, p_block = 25,
      between_prob = 0.005, n = 1000
    )),

    # ========================================================
    # JRF-PAPER-STYLE EXTREME HIGH-DIM
    # p >> n regime to mirror the JRF paper's home-turf scale.
    # Setting A perturbation. Note: multiGGM at p=500 is far too
    # slow to be practical; this setting is intended for the
    # JGL family and JRF only.
    # ========================================================
    "SF-p500n50" = modifyList(base, list(
      setting = "A",
      p = 500, B = 20, p_block = 25,
      between_prob = 0.0025, n = 50
    )),

    # ========================================================
    # NON-GAUSSIAN MACHINERY (kept but not used in current report)
    # Same network structure as SF-p125n100; applies
    # f(x) = x^3 / (1 + x^2) coordinate-wise at data-sampling time.
    # ========================================================
    "SF-nonlinear-saturation" = modifyList(base, list(
      setting = "A",
      p = 125, B = 5, p_block = 25,
      between_prob = 0.01, n = 100,
      transform = "saturation"
    )),

    # ========================================================
    # SF2- FAMILY: Barabási–Albert SCALE-FREE topology (m_pa = 2)
    #
    # Mean within-block degree ~3.8, within-block density ~15.7%.
    # Matches the original JGL/Danaher (2014) scale-free topology
    # (their mean degree ~3.9) and the pathway-restricted GGM
    # literature on metabolomics (Krumsiek et al. 2011; 5-15%
    # within-pathway density). See NOTES_network_density.md.
    #
    # Identical naming convention to SF- but with "2" suffix:
    #   SF2{±|C}p{P}n{N}
    # so legacy and new gold standards / fits don't collide on disk.
    # ========================================================
    # Anchor + perturbation sweep
    "SF2-p250n1000" = modifyList(base_sf2, list(
      setting = "A",
      p = 250, B = 10, p_block = 25,
      between_prob = 0.005, n = 1000
    )),
    "SF2+p250n1000" = modifyList(base_sf2, list(
      setting = "B",
      p = 250, B = 10, p_block = 25,
      between_prob = 0.005, n = 1000
    )),
    "SF2Cp250n1000" = modifyList(base_sf2, list(
      setting = "C",
      p = 250, B = 10, p_block = 25,
      between_prob = 0.005, n = 1000
    )),

    # Sample-size sweep (p = 250 fixed, vary n; perturbation A)
    "SF2-p250n250"  = modifyList(base_sf2, list(setting = "A", p = 250, B = 10, p_block = 25, between_prob = 0.005, n = 250)),
    "SF2-p250n500"  = modifyList(base_sf2, list(setting = "A", p = 250, B = 10, p_block = 25, between_prob = 0.005, n = 500)),
    "SF2-p250n2000" = modifyList(base_sf2, list(setting = "A", p = 250, B = 10, p_block = 25, between_prob = 0.005, n = 2000)),
    "SF2-p250n5000" = modifyList(base_sf2, list(setting = "A", p = 250, B = 10, p_block = 25, between_prob = 0.005, n = 5000)),

    # Feature-size sweep (n = 1000 fixed, vary p; perturbation A)
    # between_prob = 0.01 * (5 / B) so E[between-block edges] stays
    # comparable across p.
    "SF2-p50n1000"  = modifyList(base_sf2, list(setting = "A", p = 50,  B = 2,  p_block = 25, between_prob = 0.025,    n = 1000)),
    "SF2-p100n1000" = modifyList(base_sf2, list(setting = "A", p = 100, B = 4,  p_block = 25, between_prob = 0.0125,   n = 1000)),
    "SF2-p200n1000" = modifyList(base_sf2, list(setting = "A", p = 200, B = 8,  p_block = 25, between_prob = 0.00625,  n = 1000)),
    "SF2-p300n1000" = modifyList(base_sf2, list(setting = "A", p = 300, B = 12, p_block = 25, between_prob = 0.004167, n = 1000)),
    "SF2-p500n1000" = modifyList(base_sf2, list(setting = "A", p = 500, B = 20, p_block = 25, between_prob = 0.0025,   n = 1000)),

    # SF2-p500n250: small-n x large-p stress test (p/n = 2.0).
    # Pushes the p > n regime further than SF-p125n100 (p/n = 1.25)
    # for the discussion's "regularization is essential" claim.
    "SF2-p500n250" = modifyList(base_sf2, list(
      setting = "A",
      p = 500, B = 20, p_block = 25,
      between_prob = 0.0025, n = 250
    )),

    # SF2-p250n1000-uniform: same anchor as SF2-p250n1000 but with
    # uniform edge weights (U[-0.4,-0.1] ∪ [0.1, 0.4]) instead of
    # the empirical metabolomics-style weights. Carries forward the
    # legacy "SF- Uniform" comparison from the original CLAUDE.md
    # design table into the new naming convention. Lets us claim the
    # calibration findings hold under uniform weights too.
    "SF2-p250n1000-uniform" = modifyList(base_sf2, list(
      setting = "A",
      p = 250, B = 10, p_block = 25,
      between_prob = 0.005, n = 1000,
      wt_method = "uniform"
    )),

    # ========================================================
    # AR1- FAMILY: AR(1) banded precision (NON-modular topology)
    #
    # Theta is banded with off-diagonal entries decaying
    # geometrically: Theta[i,j] = rho^|i-j| for |i-j| <= bandwidth,
    # 0 otherwise (then diagonal-boosted for PD).
    #
    # This is the canonical "non-modular" GGM benchmark in the
    # literature (Friedman et al. 2008; Cai 2017). Adds a genuinely
    # different topology to the comparison rather than another
    # point in the BA family.
    #
    # Perturbation semantics:
    #   A : zero out off-diagonals among the first p_block_ar1 = 25
    #       nodes (mimics "delete a block")
    #   B : delete first p_block_ar1 nodes + scale next p_block_ar1
    #   C : scale first p_block_ar1 off-diagonal magnitudes
    # ========================================================
    "AR1-p250n1000" = modifyList(base_ar1, list(
      setting = "A",
      p = 250, n = 1000
    )),
    "AR1+p250n1000" = modifyList(base_ar1, list(
      setting = "B",
      p = 250, n = 1000
    )),
    "AR1Cp250n1000" = modifyList(base_ar1, list(
      setting = "C",
      p = 250, n = 1000
    )),

    # ========================================================
    # PAPER-REPLICATING SETTINGS
    # Route to 01b_generate_paper_gold_standards.R via
    # generator = "paper".
    # ========================================================
    "Danaher" = list(
      generator = "paper", paper = "danaher",
      p = 100, B = 10, p_block = 10, n = 150,
      m_pa = 2, d_scale = NULL
    ),
    # Full p=500 version of Danaher matching the paper exactly.
    # Use only for methods that can handle this size (JGL, SSJGL,
    # JRF); multiGGM and DiNGO at p=500 are too slow to be practical.
    "Danaher-full" = list(
      generator = "paper", paper = "danaher",
      p = 500, B = 10, p_block = 50, n = 150,
      m_pa = 2, d_scale = NULL
    ),
    "SSJGL-paper" = list(
      generator = "paper", paper = "ssjgl",
      p = 100, B = 10, p_block = 10, n = 150,
      b_gwish = 3
    ),
    "PSV" = list(
      generator = "paper", paper = "psv",
      p = 20, B = 1, p_block = 20, n = 100,
      n_changes = 5
    ),
    "Shaddox" = list(
      generator = "paper", paper = "shaddox",
      p = 25, B = 1, p_block = 25, n = 100,
      pct_retain = 0.70
    ),
    "DINGO-paper" = list(
      generator = "paper", paper = "dingo",
      p = 100, B = 10, p_block = 10, n = 75,
      m_pa = 2
    ),

    stop(paste("Unknown setting:", setting_label))
  )
}

#------------------------------------------------------
# Sweep-membership registry
#
# Names of the canonical sweeps used by the manuscript and poster.
# Each entry is a character vector of setting labels that belong to
# that sweep, in display order. Downstream scripts
# (06_results_tables.R, 30_poster_figures.R) can call this to avoid
# hardcoding setting names.
#
# Sweeps (per manuscript Results section ordering):
#   1. sample_size      : vary n at p = 250, BA scale-free, Setting A
#   2. feature_size     : vary p at n = 1000, BA scale-free, Setting A
#   3. perturbation     : structural / magnitude / both at the anchor
#   4. topology (NEW)   : BA-tree vs BA scale-free vs AR(1) banded
#                         at the anchor, Setting A
# The legacy SF- (BA-tree) sweeps are retained as a parallel series
# so the original results stay valid and comparable.
#------------------------------------------------------

get_sweep_settings <- function(sweep_name, topology_family = "scale_free") {
  prefix <- switch(topology_family,
    "tree"        = "SF-",     # legacy BA-tree (m_pa = 1)
    "scale_free"  = "SF2-",    # new BA scale-free (m_pa = 2)
    stop(sprintf("Unknown topology_family: %s", topology_family))
  )

  switch(sweep_name,
    "sample_size" = paste0(prefix, c("p250n250", "p250n500",
                                       "p250n1000", "p250n2000",
                                       "p250n5000")),
    "feature_size" = paste0(prefix, c("p50n1000", "p100n1000",
                                        "p200n1000", "p250n1000",
                                        "p300n1000", "p500n1000")),
    "perturbation" = c(
      paste0(prefix,                "p250n1000"),  # Setting A
      paste0(sub("-$", "+", prefix), "p250n1000"), # Setting B
      paste0(sub("-$", "C", prefix), "p250n1000")  # Setting C
    ),
    "topology" = c(
      "SF-p250n1000",   # BA-tree (legacy)
      "SF2-p250n1000",  # BA scale-free (new)
      "AR1-p250n1000"   # AR(1) banded
    ),
    stop(sprintf("Unknown sweep: %s", sweep_name))
  )
}

#------------------------------------------------------
# Method hyperparameter registry
#
# Each method has a *fixed* hyperparameter configuration that's
# used by fit_method() when it's invoked under that label. If you
# want to try different hyperparameters (e.g., multiGGM with
# burnin = 20000), ADD A NEW METHOD LABEL here — do not reuse an
# existing label. That way results never silently mix configs.
#
# The return value is hashed to form the subdirectory where
# that method's raw fits live:
#   results/raw_fits/<method>/<setting>/_config_<hash>/rep_XXXX.rds
# so different configurations coexist cleanly on disk.
#------------------------------------------------------

method_params <- function(method_name) {
  p <- switch(method_name,
    # --- JGL family (frequentist joint graphical lasso) ---
    "JGL_BIC" = list(
      family = "JGL", penalty = "fused", tune = "BIC",
      n_lambda = 10
    ),
    "JGL_CV" = list(
      family = "JGL", penalty = "fused", tune = "CV",
      nfolds = 5, n_lambda = 10
    ),
    "JGL_BIC_GGL" = list(
      family = "JGL", penalty = "group", tune = "BIC",
      n_lambda = 10
    ),
    "JGL_CV_GGL" = list(
      family = "JGL", penalty = "group", tune = "CV",
      nfolds = 5, n_lambda = 10
    ),

    # --- SSJGL family (spike-and-slab JGL via spikeyglass) ---
    "SSJGL" = list(
      family = "SSJGL", penalty = "fused",
      lambda0 = 1, lambda1 = 1, lambda2 = 1,
      v1 = 1, v0s = "auto",    # v0s = NULL => auto from lambda1
      a = 1, b = 1,
      doubly = TRUE, normalize = FALSE
    ),
    "SSJGL_GGL" = list(
      family = "SSJGL", penalty = "group",
      lambda0 = 1, lambda1 = 1, lambda2 = 1,
      v1 = 1, v0s = "auto",
      a = 1, b = 1,
      doubly = TRUE, normalize = FALSE
    ),

    # --- multiGGM (Bayesian SSVS joint GGM) ---
    # NOTE: pip_mask_thr is NOT a fit-time hyperparameter; it's a
    # scoring option (controls the post-hoc support mask applied to
    # the stored posterior-mean theta). Swept at score time without
    # re-fitting. See method_params() docstring.
    "multiGGM" = list(
      family = "multiGGM", method = "ssvs",
      burnin = 5000, nsave = 1000, nchains = 2,
      v0 = 0.01, v1 = 1, lambda = 1,
      a = 1, b = 4, alpha = 2, beta = 5,
      w = 0.8, alpha_prop = 1, beta_prop = 1,
      seed = 123
    ),
    "multiGGM_short" = list(
      family = "multiGGM", method = "ssvs",
      burnin = 2000, nsave = 500, nchains = 2,
      v0 = 0.5, v1 = 1, lambda = 1,
      a = 1, b = 4, alpha = 2, beta = 5,
      w = 0.8, alpha_prop = 1, beta_prop = 1,
      seed = 123
    ),
    "multiGGM_long" = list(
      family = "multiGGM", method = "ssvs",
      burnin = 10000, nsave = 500, nchains = 2,
      v0 = 0.0004, v1 = 1, lambda = 1,
      a = 1, b = 4, alpha = 2, beta = 5,
      w = 0.8, alpha_prop = 1, beta_prop = 1,
      seed = 123
    ),

    # --- DiNGO (regression-based partial correlation decomposition) ---
    "DiNGO" = list(
      family = "DiNGO", B = 0, pval = FALSE, cores = 1,
      rho_strategy = "bic_tune"    # or "pilot_cache" (future work)
    ),

    # --- JRF (joint random forests) ---
    "JRF" = list(
      family = "JRF", ntree = 1000, mtry = "sqrt"    # mtry=NULL => sqrt(p-1)
    ),

    # --- InvCov (unregularized MLE baseline) ---
    # Per-group precision matrix as the inverse of the sample covariance.
    # This is the "no regularization" floor: when n > p it's the maximum-
    # likelihood estimator; when n <= p the sample covariance is singular
    # so we fall back to the Moore-Penrose pseudoinverse via MASS::ginv().
    # Useful as a baseline that quantifies how much regularized methods
    # actually buy you across (n, p) regimes.
    "InvCov" = list(
      family = "InvCov", regularize = "ginv"    # "ginv" or "drop"
    ),

    stop(sprintf("Unknown method label: '%s'", method_name))
  )
  p$method_name <- method_name
  p
}

#------------------------------------------------------
# Config hash: stable identifier for a (method, params) combo.
# Used to name the _config_<hash> subdir so different
# hyperparameter configurations coexist cleanly on disk.
#------------------------------------------------------

method_config_hash <- function(method_name, params = NULL) {
  if (is.null(params)) params <- method_params(method_name)
  # Strip non-fit fields that shouldn't affect the hash.
  # (method_name is encoded in the path, so no need to hash it.)
  params$method_name <- NULL
  # Sort for stability across R versions / platforms
  params <- params[order(names(params))]
  substr(digest::digest(params, algo = "xxhash64"), 1, 8)
}
