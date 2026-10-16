# Gold-standard network generation for the on-demand simulator.
#
# A lean port of the DiNA simulation generator
# (inspo_code/01_generate_gold_standard_networks.R), restricted to the
# block-structured Barabasi-Albert family the grid was built from. Semantics
# match the original; the only generalization is that `p` need not be a
# multiple of the block size.

# User-facing perturbation -> DiNA setting letter.
.PERT_SETTING <- c(structural = "A", both = "B", magnitude = "C")

# Block sizes for p features. Matches the DiNA registry at every gridded p
# (25-node blocks; p < 50 uses two equal blocks, as the p = 20 setting did);
# other p absorb the remainder into the last block.
.block_sizes <- function(p, p_block = 25L) {
  if (p < 2L * p_block) {
    h <- p %/% 2L
    return(c(h, p - h))
  }
  B <- max(2L, as.integer(round(p / p_block)))
  sizes <- rep(p_block, B)
  sizes[B] <- sizes[B] + (p - sum(sizes))
  sizes
}

# Draw edge weights for the upper-triangle edges of `adj` from the empirical
# distribution (histogram mids/density), filtering near-zero bins as DiNA did.
.draw_weights <- function(adj, weights, abs_min_filter = 0.035) {
  p <- nrow(adj)
  ut <- upper.tri(adj)
  is_edge <- adj[ut] == 1L
  W <- matrix(0, p, p)
  if (!any(is_edge)) return(W)
  emp <- weights[abs(weights$mids) > abs_min_filter, ]
  emp$density[emp$density == 0] <- min(emp$density[emp$density > 0])
  w_ut <- numeric(sum(ut))
  w_ut[is_edge] <- sample(emp$mids, sum(is_edge), replace = TRUE,
                          prob = emp$density)
  W[ut] <- w_ut
  W + t(W)
}

# Row-sum normalization (Danaher et al. 2014): off-diagonals divided by
# 1.5 x row |sum|, symmetrized, unit diagonal, PD safety boost.
.row_sum_precision <- function(W) {
  A <- W
  diag(A) <- 0
  rs <- rowSums(abs(A))
  nz <- rs > 0
  A[nz, ] <- A[nz, ] / (1.5 * rs[nz])
  A <- (A + t(A)) / 2
  diag(A) <- 1
  .make_pd(A)
}

.make_pd <- function(Theta) {
  ok <- tryCatch({ chol(Theta); TRUE }, error = function(e) FALSE)
  if (!ok) {
    eig_min <- min(eigen(Theta, symmetric = TRUE, only.values = TRUE)$values)
    diag(Theta) <- diag(Theta) + abs(eig_min) + 0.01
  }
  Theta
}

# Scale the off-diagonal entries of Theta within `idx` by `factor`.
.scale_block <- function(Theta, idx, factor) {
  if (length(idx) < 2L) return(Theta)
  sub <- Theta[idx, idx]
  d <- diag(sub)
  sub <- sub * factor
  diag(sub) <- d
  Theta[idx, idx] <- sub
  .make_pd((Theta + t(Theta)) / 2)
}

# Generate one gold standard: two precision/covariance matrices and the
# differential truth (upper-triangle logical vector).
#   p            : number of features
#   perturbation : "structural" | "magnitude" | "both"
#   family       : data-generating family (defaults to the shipped grid's)
.generate_gold <- function(p, perturbation, family = howbig_family,
                           weights = howbig_weights) {
  setting <- .PERT_SETTING[[perturbation]]
  sizes <- .block_sizes(p, family$p_block)
  B <- length(sizes)
  memb <- rep(seq_len(B), sizes)
  between_prob <- 0.01 * (5 / B)

  # 1-2. BA blocks on the diagonal + sparse between-block edges
  adj <- matrix(0L, p, p)
  off <- 0L
  for (b in seq_len(B)) {
    g <- igraph::sample_pa(sizes[b], power = 1, m = family$m_pa,
                           directed = FALSE)
    idx <- off + seq_len(sizes[b])
    adj[idx, idx] <- igraph::as_adjacency_matrix(g, sparse = FALSE)
    off <- off + sizes[b]
  }
  ut <- upper.tri(adj)
  cross <- ut & outer(memb, memb, "!=")
  add <- cross & (matrix(stats::runif(p * p), p, p) < between_prob)
  adj[add] <- 1L
  adj <- pmax(adj, t(adj))

  # 3-4. weights and group-1 precision
  W1 <- .draw_weights(adj, weights)
  Theta1 <- .row_sum_precision(W1)

  # 5-6. group-2: delete block-1 edges (A, B) before construction...
  W2 <- W1
  b1 <- which(memb == 1L)
  if (setting %in% c("A", "B")) W2[b1, b1] <- 0
  Theta2 <- .row_sum_precision(W2)
  # ...and apply magnitude scaling after it, so row-sum normalization
  # cannot cancel it (B: block 2; C: block 1).
  if (setting == "B") Theta2 <- .scale_block(Theta2, which(memb == 2L),
                                             family$scale_factor)
  if (setting == "C") Theta2 <- .scale_block(Theta2, b1, family$scale_factor)

  # 7-8. covariances and magnitude-threshold differential truth
  delta <- abs(Theta1 - Theta2) > family$diff_threshold
  list(p = p,
       sigma = list(solve(Theta1), solve(Theta2)),
       delta = delta[ut])
}

# Draw one dataset (list of two n_g x p matrices) from a gold standard.
.gen_data <- function(gold, n) {
  mu <- rep(0, gold$p)
  list(MASS::mvrnorm(n[1], mu, gold$sigma[[1]]),
       MASS::mvrnorm(n[2], mu, gold$sigma[[2]]))
}
