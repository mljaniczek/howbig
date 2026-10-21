# Simulate power for a design

Runs an on-demand simulation for a single design cell: generate networks
and data under the grid's data-generating family, fit the method, and
score its FDR-controlled differential calls. Use this directly to answer
designs the grid does not cover (other `p`, unequal group sizes, other
FDR targets), or via `power_for(..., simulate = TRUE)`.

## Usage

``` r
simulate_power(
  p,
  n,
  method,
  perturbation = "structural",
  fdr = 0.05,
  reps = 20L,
  n_perm = 30L,
  cores = 1L,
  seed = 1L,
  max_minutes = 60,
  verbose = TRUE
)
```

## Arguments

- p:

  Number of features.

- n:

  Per-group sample size; a single number, or `c(n1, n2)` for unequal
  groups.

- method:

  Estimation method label. See
  [`howbig_coverage()`](https://mljaniczek.github.io/howbig/reference/howbig_coverage.md)
  for the available methods.

- perturbation:

  Type of differential signal: `"structural"` (edges appear/disappear),
  `"magnitude"` (edge weights change), or `"both"`.

- fdr:

  Target false-discovery rate for the differential calls.

- reps:

  Number of simulated replicates (fresh network + data each).

- n_perm:

  Number of label permutations for the null (ignored for multiGGM).

- cores:

  Parallel worker processes. Workers are fresh R sessions (PSOCK), not
  forks, so they are safe with multithreaded BLAS such as Apple
  Accelerate; `howbig` must be installed for `cores > 1`.

- seed:

  Base seed; results are reproducible for a given seed regardless of
  `cores`.

- max_minutes:

  Refuse to start if the estimated runtime exceeds this many minutes.
  The estimate comes from per-fit timings of the production grid runs;
  raise this to run long designs deliberately.

- verbose:

  Print progress messages.

## Value

A `howbig_power` object with `provenance = "simulated"`.

## Details

Inference matches how the grid was scored. Frequentist methods are
called against a pooled label-permutation null (one representative
network, `n_perm` permutations with fresh data, Benjamini-Hochberg);
multiGGM uses its posterior probabilities directly (direct-posterior
FDR).

Simulation is slow: each replicate and each permutation is a full method
fit. Start with small `reps` / `n_perm` to gauge runtime.

## Examples

``` r
if (FALSE) { # \dontrun{
simulate_power(p = 100, n = 500, method = "JGL_BIC", reps = 10)
} # }
```
