#------------------------------------------------------
# Block-Structured Gold Standard Network Generation
# Implements the simulation design from Section 5 of
# the DiNA Tutorial manuscript.
#------------------------------------------------------

library(igraph)
library(MASS)

#------------------------------------------------------
# Function 1: Generate a single block's graph structure
#------------------------------------------------------

generate_block_graph <- function(p_block,
                                  topology = "random",
                                  within_density_random = 0.15,
                                  m_pa = 3) {
  # `within_density_random` controls Erdős–Rényi edge probability and
  # is ignored for `topology = "scale_free"` (where density is set by
  # the BA preferential-attachment parameter `m_pa`). See
  # NOTES_network_density.md for why this parameter was renamed.
  if (topology == "random") {
    g <- sample_gnp(p_block, within_density_random)
  } else if (topology == "scale_free") {
    g <- sample_pa(n = p_block, power = 1, m = m_pa, directed = FALSE)
  } else {
    stop("topology must be 'random' or 'scale_free' (AR(1) banded ",
         "is handled by generate_gold_standard_ar1)")
  }
  return(g)
}

#------------------------------------------------------
# Function 2: Assemble B blocks into a full network
#------------------------------------------------------

assemble_block_network <- function(block_graphs,
                                    between_block_prob = 0.01,
                                    p = NULL) {

  B <- length(block_graphs)
  block_sizes <- sapply(block_graphs, vcount)
  if (is.null(p)) p <- sum(block_sizes)
  stopifnot(sum(block_sizes) == p)

  # Build block membership vector
  block_membership <- rep(seq_len(B), times = block_sizes)

  # Initialize full adjacency matrix
  adj_matrix <- matrix(0L, nrow = p, ncol = p)

  # Place each block's adjacency on the diagonal
  offset <- 0
  for (b in seq_len(B)) {
    idx <- (offset + 1):(offset + block_sizes[b])
    block_adj <- as.matrix(as_adjacency_matrix(block_graphs[[b]]))
    adj_matrix[idx, idx] <- block_adj
    offset <- offset + block_sizes[b]
  }

  # Add sparse between-block edges
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      if (block_membership[i] != block_membership[j]) {
        if (runif(1) < between_block_prob) {
          adj_matrix[i, j] <- 1L
          adj_matrix[j, i] <- 1L
        }
      }
    }
  }

  graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")

  return(list(
    adj_matrix = adj_matrix,
    graph = graph,
    block_membership = block_membership,
    block_sizes = block_sizes
  ))
}

#------------------------------------------------------
# Function 3: Assign edge weights
#------------------------------------------------------

assign_edge_weights <- function(adj_matrix,
                                 wt_method = "empirical",
                                 empirical_dist = NULL,
                                 unif_range = c(0.035, 0.75),
                                 abs_min_filter = 0.035) {

  p <- nrow(adj_matrix)
  # Count upper-triangle edges
  edge_idx <- which(adj_matrix[upper.tri(adj_matrix)] == 1)
  n_edges <- length(edge_idx)

  if (n_edges == 0) {
    return(matrix(0, nrow = p, ncol = p))
  }

  if (wt_method == "empirical") {
    if (is.null(empirical_dist)) {
      stop("empirical_dist must be provided when wt_method = 'empirical'")
    }
    # Filter to non-trivial weights
    emp <- empirical_dist[abs(empirical_dist$mids) > abs_min_filter, ]
    # Replace zero densities with minimum nonzero density
    min_dens <- min(emp$density[emp$density > 0])
    emp$density[emp$density == 0] <- min_dens

    if (nrow(emp) < 10) {
      warning("Fewer than 10 empirical bins after filtering; falling back to uniform weights")
      wt_method <- "uniform"
    } else {
      weights <- sample(emp$mids, size = n_edges, replace = TRUE, prob = emp$density)
    }
  }

  if (wt_method == "uniform") {
    a <- unif_range[1]
    b <- unif_range[2]
    samp_pos <- runif(ceiling(n_edges / 2), min = a, max = b)
    samp_neg <- runif(ceiling(n_edges / 2), min = -b, max = -a)
    weights <- sample(c(samp_pos, samp_neg), size = n_edges)
  }

  # Build symmetric weight matrix
  W <- matrix(0, nrow = p, ncol = p)
  ut <- upper.tri(W)
  ut_vals <- adj_matrix[ut]
  # Place weights at edge positions in upper triangle
  ut_weights <- rep(0, length(ut_vals))
  ut_weights[ut_vals == 1] <- weights
  W[ut] <- ut_weights
  # Mirror to lower triangle
  W <- W + t(W)

  return(W)
}

#------------------------------------------------------
# Function 4: Build precision matrix
#
# Two methods available:
#   "eigenvalue_boost" (default) — adds |min_eig|*1.01 + epsilon
#       to diagonal. Simple but inflates diagonals uniformly,
#       producing higher condition numbers (~12-27 for p=100).
#   "row_sum_normalize" — divides each off-diagonal entry by
#       1.5 * its row sum, then symmetrizes. Follows Danaher
#       et al. (2014). Produces well-conditioned matrices (~3-5)
#       with adaptive per-node diagonals.
#------------------------------------------------------

build_precision_matrix <- function(W, epsilon = 0.1,
                                    method = c("eigenvalue_boost", "row_sum_normalize"),
                                    d_scale = 0.6) {
  method <- match.arg(method)

  if (method == "eigenvalue_boost") {
    return(build_precision_eigenvalue_boost(W, epsilon))
  } else {
    return(build_precision_row_sum(W, d_scale = d_scale))
  }
}

#--- Eigenvalue boost (original method) ---

build_precision_eigenvalue_boost <- function(W, epsilon = 0.1) {
  p <- nrow(W)
  Theta <- W
  diag(Theta) <- 0

  # Find minimum eigenvalue of the off-diagonal matrix
  min_eig <- min(eigen(Theta, symmetric = TRUE, only.values = TRUE)$values)

  # Shift the diagonal by |min_eig| * 1.01 + epsilon
  diag(Theta) <- abs(min_eig) * 1.01 + epsilon

  # Final verification
  eig_vals <- eigen(Theta, symmetric = TRUE, only.values = TRUE)$values
  if (any(eig_vals <= 0)) {
    stop("Precision matrix is not positive definite after eigenvalue boosting.")
  }

  return(Theta)
}

#--- Row-sum normalization (Danaher et al. 2014) ---

build_precision_row_sum <- function(W, d_scale = NULL) {
  p <- nrow(W)
  A <- W
  diag(A) <- 0

  # Step 1: Divide each off-diagonal element by 1.5 * sum(|off-diag in row|)
  for (i in 1:p) {
    row_sum <- sum(abs(A[i, ]))
    if (row_sum > 0) {
      A[i, ] <- A[i, ] / (1.5 * row_sum)
    }
  }

  # Step 2: Symmetrize and set diagonal to 1
  A <- (A + t(A)) / 2
  diag(A) <- 1

  # Verify PD (row-sum normalization guarantees diagonal dominance,
  # but check anyway for safety)
  eigs <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigs <= 0)) {
    # Safety boost if needed (shouldn't happen with 1.5x divisor)
    diag(A) <- diag(A) + abs(min(eigs)) + 0.01
    warning("Row-sum normalized matrix required additional diagonal boost for PD")
  }

  # Optional: Danaher-style covariance scaling (d_scale)
  # Scales the covariance so that off-diagonal correlations are attenuated
  if (!is.null(d_scale)) {
    Sigma_raw <- solve(A)
    d_sqrt <- sqrt(diag(Sigma_raw))
    Sigma <- Sigma_raw / outer(d_sqrt, d_sqrt)
    Sigma <- Sigma * d_scale
    diag(Sigma) <- 1
    Sigma <- (Sigma + t(Sigma)) / 2
    A <- solve(Sigma)
    A <- (A + t(A)) / 2
  }

  return(A)
}

#------------------------------------------------------
# Function 5: Apply perturbation to create Group 2
#------------------------------------------------------

apply_perturbation <- function(adj_matrix_g1,
                                weight_matrix_g1,
                                block_membership,
                                setting = "A",
                                wt_method = "empirical",
                                empirical_dist = NULL,
                                unif_range = c(0.035, 0.75),
                                abs_min_filter = 0.035,
                                scale_factor = 1.5) {

  # Null setting: no perturbation (negative control)
  if (setting == "null") {
    return(list(
      adj_matrix_g2 = adj_matrix_g1,
      weight_matrix_g2 = weight_matrix_g1
    ))
  }

  adj_g2 <- adj_matrix_g1
  wt_g2  <- weight_matrix_g1

  block1_idx <- which(block_membership == 1)

  # ----- Setting C: purely activational (Block-1 magnitude shift) -----
  # The network STRUCTURE is identical in g1 and g2 (same edges,
  # same non-edges). The signal comes purely from a magnitude
  # change in Block-1 precision values, applied post-construction
  # in generate_gold_standard() (analogous to Setting B's Block-2
  # scaling). This contrasts with Setting A (delete) and
  # Setting B (delete + scale): Setting C is the "no structural
  # change at all, just shifted strength" perturbation.
  #
  # Here in apply_perturbation() we return g1 unchanged in
  # adjacency/weight space; the Block-1 Theta scaling happens
  # downstream.
  if (setting == "C") {
    return(list(
      adj_matrix_g2    = adj_matrix_g1,
      weight_matrix_g2 = weight_matrix_g1
    ))
  }

  # ----- Settings A and B (the original behavior) -----
  # Both delete all within-Block-1 edges at the weight / adjacency
  # level. Setting B's additional "scale Block 2 by scale_factor"
  # perturbation used to be applied at the weight level here
  # (resample block-2 weights, multiply by scale_factor, put back)
  # but that is completely cancelled by row-sum normalization in
  # build_precision_matrix(). The scaling is now applied
  # post-construction, directly on Theta, in
  # generate_gold_standard(); see the "Setting B: block-2 precision
  # scaling" block there. This change makes the perturbation
  # semantics invariant to the precision-construction choice.
  adj_g2[block1_idx, block1_idx] <- 0L
  wt_g2 [block1_idx, block1_idx] <- 0

  return(list(
    adj_matrix_g2    = adj_g2,
    weight_matrix_g2 = wt_g2
  ))
}

#------------------------------------------------------
# Function 6: Compute gold standard differential truths
#------------------------------------------------------

compute_gold_standard_truths <- function(Theta_g1, Theta_g2,
                                          adj_g1, adj_g2,
                                          diff_threshold = 0.01) {

  p <- nrow(Theta_g1)

  # Partial correlation: rho_ij = -Theta_ij / sqrt(Theta_ii * Theta_jj)
  pcor_from_precision <- function(Theta) {
    d <- sqrt(diag(Theta))
    Pcor <- -Theta / outer(d, d)
    diag(Pcor) <- 1
    return(Pcor)
  }

  Pcor_g1 <- pcor_from_precision(Theta_g1)
  Pcor_g2 <- pcor_from_precision(Theta_g2)

  # Delta_A (primary): magnitude-threshold differential. An edge is
  # differential iff the true precision difference exceeds
  # diff_threshold. This matches the primary detection rule in
  # calc_metrics (|theta_hat^(1) - theta_hat^(2)| > tau_delta),
  # ensuring the ground truth and the method's detection rule
  # operate on the same quantity. Captures both structural changes
  # (edges deleted, theta -> 0) and magnitude changes (edges with
  # different weights across groups).
  diff_theta <- abs(Theta_g1 - Theta_g2)
  Delta_A <- (diff_theta > diff_threshold) * 1L
  diag(Delta_A) <- 0L

  # Delta_A_XOR (secondary, legacy): pure structural XOR of
  # thresholded adjacencies. Flags only edges that exist in one
  # group but not the other — misses magnitude perturbations like
  # Setting B's block-2 scaling. Retained for the diff_XOR_*
  # diagnostic columns in calc_metrics.
  Delta_A_XOR <- (adj_g1 != adj_g2) * 1L
  diag(Delta_A_XOR) <- 0L

  # Union of edges across groups
  U <- ((adj_g1 == 1) | (adj_g2 == 1)) * 1L
  diag(U) <- 0L

  # Delta_P: partial correlation differences
  Delta_P <- Pcor_g1 - Pcor_g2
  diag(Delta_P) <- 0

  # Edge-restricted: only at union of edges
  Delta_P_edge <- Delta_P * U

  # Global: all entries (no masking)
  Delta_P_global <- Delta_P

  return(list(
    Delta_A        = Delta_A,        # PRIMARY: magnitude-threshold truth
    Delta_A_XOR    = Delta_A_XOR,    # secondary: structural XOR
    Delta_P_edge   = Delta_P_edge,
    Delta_P_global = Delta_P_global,
    Pcor_g1        = Pcor_g1,
    Pcor_g2        = Pcor_g2,
    diff_threshold = diff_threshold  # for reference / audit
  ))
}

#------------------------------------------------------
# Function 7: Generate complete gold standard
#------------------------------------------------------

generate_gold_standard <- function(p = 100,
                                    B = 5,
                                    p_block = 20,
                                    topology = "random",
                                    within_density_random = 0.15,
                                    between_prob = 0.01,
                                    setting = "A",
                                    wt_method = "empirical",
                                    empirical_dist = NULL,
                                    unif_range = c(0.035, 0.75),
                                    abs_min_filter = 0.035,
                                    epsilon = 0.1,
                                    precision_method = "row_sum_normalize",
                                    d_scale = NULL,
                                    scale_factor = 1.5,
                                    m_pa = 1,
                                    ar1_rho = 0.4,
                                    ar1_bandwidth = 1,
                                    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # AR(1) banded topology routes to its own generator below.
  if (identical(topology, "ar1_banded")) {
    return(generate_gold_standard_ar1(
      p = p, setting = setting,
      rho = ar1_rho, bandwidth = ar1_bandwidth,
      scale_factor = scale_factor,
      seed = NULL  # already consumed above
    ))
  }

  stopifnot(p == B * p_block)

  # 1. Generate B block graphs
  block_graphs <- lapply(seq_len(B), function(b) {
    generate_block_graph(p_block, topology, within_density_random, m_pa)
  })

  # Record realized within-block densities
  realized_densities <- sapply(block_graphs, function(g) {
    edge_density(g)
  })

  # 2. Assemble into full network
  net <- assemble_block_network(block_graphs, between_prob, p)

  # 3. Assign edge weights for G1
  W_g1 <- assign_edge_weights(
    adj_matrix = net$adj_matrix,
    wt_method = wt_method,
    empirical_dist = empirical_dist,
    unif_range = unif_range,
    abs_min_filter = abs_min_filter
  )

  # 4. Build precision matrix for G1
  Theta_g1 <- build_precision_matrix(W_g1, epsilon,
                                      method = precision_method,
                                      d_scale = d_scale)

  # 5. Apply perturbation to create G2
  g2 <- apply_perturbation(
    adj_matrix_g1 = net$adj_matrix,
    weight_matrix_g1 = W_g1,
    block_membership = net$block_membership,
    setting = setting,
    wt_method = wt_method,
    empirical_dist = empirical_dist,
    unif_range = unif_range,
    abs_min_filter = abs_min_filter,
    scale_factor = scale_factor
  )

  # 6. Build precision matrix for G2
  Theta_g2 <- build_precision_matrix(g2$weight_matrix_g2, epsilon,
                                      method = precision_method,
                                      d_scale = d_scale)

  # 6b. Setting B: block-2 precision scaling.
  #
  # Apply the scale_factor directly to the block-2 x block-2
  # OFF-DIAGONAL entries of Theta_g2 (diagonals untouched). Scaling
  # the off-diagonals by c scales the partial correlations at those
  # edges by ~c, giving a clean, method-agnostic magnitude
  # perturbation visible to every estimator (JGL / SSJGL / multiGGM
  # via precision; DiNGO via pcor).
  #
  # Doing this POST-construction bypasses row-sum normalization,
  # which would otherwise cancel the scaling (Theta = W / (1.5 *
  # rowsum(|W|)), so scaling W within a block scales both numerator
  # and denominator and leaves Theta ~unchanged).
  #
  # Positive-definiteness: off-diagonals within a small block stay
  # much smaller than the (approximately unit) diagonal, so a
  # 1.5x (or similar moderate) scale keeps the matrix PD. We
  # verify via chol() and, on failure, boost the diagonal to
  # recover PD and emit a warning.
  if (setting == "B") {
    block2_idx <- which(net$block_membership == 2)
    if (length(block2_idx) > 1) {
      sub <- Theta_g2[block2_idx, block2_idx]
      diag_sub <- diag(sub)
      sub <- sub * scale_factor           # scale whole block
      diag(sub) <- diag_sub               # restore diagonal (keep as-is)
      # Write scaled sub-block back and resymmetrize
      Theta_g2[block2_idx, block2_idx] <- sub
      Theta_g2 <- (Theta_g2 + t(Theta_g2)) / 2

      # Positive-definiteness check + safety boost
      pd_ok <- tryCatch({ chol(Theta_g2); TRUE },
                         error = function(e) FALSE)
      if (!pd_ok) {
        eig_min <- min(eigen(Theta_g2, symmetric = TRUE,
                              only.values = TRUE)$values)
        diag(Theta_g2) <- diag(Theta_g2) + abs(eig_min) + 0.01
        warning(sprintf(
          "Setting B: Theta_g2 lost PD after block-2 off-diagonal scaling (scale_factor = %g). Diagonal boosted by %.4f to recover PD.",
          scale_factor, abs(eig_min) + 0.01
        ))
      }
    }
  }

  # 6c. Setting C: Block-1 precision scaling.
  #
  # Applies scale_factor to the Block-1 x Block-1 OFF-DIAGONAL
  # entries of Theta_g2 (diagonals untouched). This produces a
  # purely activational signal: the network structure is
  # identical in g1 and g2 (same adjacency, same non-edges), but
  # the strength of Block-1 conditional dependences shifts by
  # ~scale_factor in g2. Methods that look at |theta_g1 -
  # theta_g2| see a graded signal localized to Block-1's existing
  # edges, with no help from edge-presence/absence cues.
  #
  # Like Setting B, this is applied post-construction (after
  # build_precision_matrix) so that row-sum normalization does
  # not cancel the scaling. Same PD safeguard.
  if (setting == "C") {
    block1_idx <- which(net$block_membership == 1)
    if (length(block1_idx) > 1) {
      sub <- Theta_g2[block1_idx, block1_idx]
      diag_sub <- diag(sub)
      sub <- sub * scale_factor           # scale whole block
      diag(sub) <- diag_sub               # restore diagonal (keep as-is)
      Theta_g2[block1_idx, block1_idx] <- sub
      Theta_g2 <- (Theta_g2 + t(Theta_g2)) / 2

      pd_ok <- tryCatch({ chol(Theta_g2); TRUE },
                         error = function(e) FALSE)
      if (!pd_ok) {
        eig_min <- min(eigen(Theta_g2, symmetric = TRUE,
                              only.values = TRUE)$values)
        diag(Theta_g2) <- diag(Theta_g2) + abs(eig_min) + 0.01
        warning(sprintf(
          "Setting C: Theta_g2 lost PD after block-1 off-diagonal scaling (scale_factor = %g). Diagonal boosted by %.4f to recover PD.",
          scale_factor, abs(eig_min) + 0.01
        ))
      }
    }
  }

  # 7. Compute covariance matrices
  Sigma_g1 <- solve(Theta_g1)
  Sigma_g2 <- solve(Theta_g2)
  # TODO do we need to check for symmetry/positive definiteness of these covariance matrices? Like in Danaher?

  # 8. Compute gold standard truths
  truths <- compute_gold_standard_truths(
    Theta_g1, Theta_g2,
    net$adj_matrix, g2$adj_matrix_g2
  )

  return(list(
    adj = list(g1 = net$adj_matrix, g2 = g2$adj_matrix_g2),
    weights = list(g1 = W_g1, g2 = g2$weight_matrix_g2),
    omega = list(g1 = Theta_g1, g2 = Theta_g2),
    sigma = list(g1 = Sigma_g1, g2 = Sigma_g2),
    truth = truths,
    block_membership = net$block_membership,
    params = list(
      p = p, B = B, p_block = p_block,
      topology = topology,
      within_density_random = within_density_random,
      between_prob = between_prob,
      setting = setting,
      wt_method = wt_method,
      epsilon = epsilon,
      precision_method = precision_method,
      d_scale = d_scale,
      scale_factor = scale_factor,
      m_pa = m_pa,
      seed = seed,
      realized_block_densities = realized_densities
    )
  ))
}

#------------------------------------------------------
# Data sampling function
#------------------------------------------------------

gen_dat <- function(gold, n, transform = NULL) {
  p <- nrow(gold$omega$g1)
  # Return raw MVN samples WITHOUT post-hoc mean-centering.
  #
  # NOTE: Previously this function subtracted each group's sample
  # column means (dat1 <- sweep(dat1, 2, colMeans(dat1))). That was
  # added out of caution but caused a major silent bug in DiNGO:
  #   (1) The data is already drawn from zero-mean MVN, so sample
  #       means are near zero by construction.
  #   (2) The sample covariance cov() is invariant to pre-centering
  #       (it subtracts means internally), so JGL/SSJGL/multiGGM
  #       results are unaffected.
  #   (3) DiNGO's Step 2 regresses the pooled-precision residuals
  #       on the group indicator x. Forcing each group's sample
  #       mean to exactly zero destroys the linear relationship
  #       between x and the data, so the regression coefficients
  #       on x are effectively zero and R1 = R2 at machine
  #       precision — the "R1 = R2 bug" we debugged for months.
  # Removing the centering step fixes DiNGO without affecting any
  # other method.
  #
  # OPTIONAL `transform` argument: apply a monotonic non-linear
  # coordinate-wise transformation to simulate non-Gaussian data.
  # Currently supports:
  #   "saturation" : f(x) = x^3 / (1 + x^2)  — sigmoid-like
  #                  compression mimicking RNA-seq / mass-spec
  #                  saturation and Hill-kinetic activation
  # TODO: add "zero_inflated" for scRNA-seq-like drop-out
  #       (Poisson/negative binomial with excess zeros)
  dat1 <- mvrnorm(n = n, mu = rep(0, p), Sigma = gold$sigma$g1)
  dat2 <- mvrnorm(n = n, mu = rep(0, p), Sigma = gold$sigma$g2)
  if (!is.null(transform)) {
    if (transform == "saturation") {
      # Monotonic sigmoid-like transform. Preserves rank dependence
      # but produces non-Gaussian marginals. See 09 report section
      # "Non-Gaussian setting" for biological interpretation.
      sat_fn <- function(x) x^3 / (1 + x^2)
      dat1 <- sat_fn(dat1)
      dat2 <- sat_fn(dat2)
    } else {
      stop(sprintf("Unknown transform: '%s'", transform))
    }
  }
  return(list(dat1 = dat1, dat2 = dat2))
}

#------------------------------------------------------
# AR(1) banded gold-standard generator (NON-modular)
#
# Builds Theta with banded off-diagonals decaying geometrically:
#   Theta[i,j] = rho^|i-j|  for 1 <= |i-j| <= bandwidth
#   Theta[i,j] = 0          for |i-j| > bandwidth
# Then boosts the diagonal so the matrix is PD.
#
# Perturbations operate on a contiguous "block" of nodes 1..p_block_ar1
# (default 25), mirroring the modular-perturbation semantics so the
# differential signal is comparable across topology families:
#   "A"  : zero out off-diagonals among nodes 1..p_block_ar1
#   "B"  : delete first block + scale next block off-diagonals
#   "C"  : scale first block off-diagonals by scale_factor
# Block size 25 matches p_block of the modular settings.
#------------------------------------------------------

generate_gold_standard_ar1 <- function(p = 250,
                                        setting = "A",
                                        rho = 0.4,
                                        bandwidth = 1,
                                        p_block_ar1 = 25,
                                        scale_factor = 1.5,
                                        diff_threshold = 0.01,
                                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(p_block_ar1 <= p, bandwidth >= 1, abs(rho) < 1)

  # Build Theta_g1
  Theta_g1 <- matrix(0, p, p)
  for (k in seq_len(bandwidth)) {
    idx <- seq_len(p - k)
    Theta_g1[cbind(idx, idx + k)] <- rho^k
    Theta_g1[cbind(idx + k, idx)] <- rho^k
  }
  diag(Theta_g1) <- 1
  # PD safety boost
  eig_min <- min(eigen(Theta_g1, symmetric = TRUE,
                        only.values = TRUE)$values)
  if (eig_min <= 0) {
    diag(Theta_g1) <- diag(Theta_g1) + abs(eig_min) + 0.01
  }

  # Adjacency from non-zero off-diagonals
  adj_g1 <- (Theta_g1 != 0) * 1L
  diag(adj_g1) <- 0L

  # Perturbation -> Theta_g2 + adj_g2
  Theta_g2 <- Theta_g1
  adj_g2   <- adj_g1
  pert_block <- seq_len(p_block_ar1)

  if (setting == "A") {
    # Zero out off-diagonals within the first block
    Theta_g2[pert_block, pert_block] <-
      Theta_g1[pert_block, pert_block] * diag(length(pert_block))
    diag(Theta_g2)[pert_block] <- diag(Theta_g1)[pert_block]
    adj_g2[pert_block, pert_block] <- 0L
  } else if (setting == "B") {
    # Delete first block; scale next block off-diagonals
    Theta_g2[pert_block, pert_block] <-
      Theta_g1[pert_block, pert_block] * diag(length(pert_block))
    diag(Theta_g2)[pert_block] <- diag(Theta_g1)[pert_block]
    adj_g2[pert_block, pert_block] <- 0L
    block2 <- seq_len(p_block_ar1) + p_block_ar1
    block2 <- block2[block2 <= p]
    if (length(block2) > 1) {
      sub <- Theta_g2[block2, block2]
      diag_sub <- diag(sub)
      sub <- sub * scale_factor
      diag(sub) <- diag_sub
      Theta_g2[block2, block2] <- sub
    }
  } else if (setting == "C") {
    # Scale first block off-diagonals only (purely magnitude)
    sub <- Theta_g2[pert_block, pert_block]
    diag_sub <- diag(sub)
    sub <- sub * scale_factor
    diag(sub) <- diag_sub
    Theta_g2[pert_block, pert_block] <- sub
  } else if (setting != "null") {
    stop(sprintf("Unknown AR(1) setting: %s", setting))
  }

  Theta_g2 <- (Theta_g2 + t(Theta_g2)) / 2
  # PD safety for g2
  eig_min_g2 <- min(eigen(Theta_g2, symmetric = TRUE,
                           only.values = TRUE)$values)
  if (eig_min_g2 <= 0) {
    diag(Theta_g2) <- diag(Theta_g2) + abs(eig_min_g2) + 0.01
  }

  truths <- compute_gold_standard_truths(
    Theta_g1, Theta_g2, adj_g1, adj_g2,
    diff_threshold = diff_threshold
  )

  list(
    adj     = list(g1 = adj_g1, g2 = adj_g2),
    weights = list(g1 = Theta_g1, g2 = Theta_g2),  # AR(1): W = Theta
    omega   = list(g1 = Theta_g1, g2 = Theta_g2),
    sigma   = list(g1 = solve(Theta_g1), g2 = solve(Theta_g2)),
    truth   = truths,
    block_membership = rep(1L, p),  # non-modular -> single block
    params  = list(
      p = p, setting = setting,
      topology = "ar1_banded",
      ar1_rho = rho, ar1_bandwidth = bandwidth,
      p_block_ar1 = p_block_ar1,
      scale_factor = scale_factor,
      seed = seed
    )
  )
}
