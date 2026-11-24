# Power (sensitivity) for a given design

Estimates the power of a differential network analysis at a given sample
size and feature count, where **power is the sensitivity** – the
fraction of truly differential edges recovered – at a
false-discovery-rate target. The achieved FDR is reported alongside,
because several methods are not calibrated in every regime.

## Usage

``` r
power_for(
  p,
  n,
  method,
  perturbation = "structural",
  fdr = .HOWBIG_FDR_TARGET,
  scenario = "default",
  simulate = FALSE,
  ...
)
```

## Arguments

- p:

  Number of features.

- n:

  Per-group sample size; a single number, or `c(n1, n2)` for unequal
  groups. Group 1 is the reference group (its network is intact; group 2
  carries the change), so order matters.

- method:

  Estimation method label. See
  [`howbig_methods()`](https://mljaniczek.github.io/howbig/reference/howbig_methods.md).

- perturbation:

  Type of differential signal: `"structural"` (edges appear/disappear),
  `"magnitude"` (edge strengths change), or `"both"`.

- fdr:

  Target false-discovery rate. The grid is scored at `0.05`; other
  targets require simulation.

- scenario:

  Assumption about the networks and data, for checking how sensitive an
  answer is to them. `"default"` (block-structured Barabasi-Albert tree
  networks, Gaussian data, empirical edge weights) is the main grid. The
  others are in the grid only at p = 100, n = 1000 with a structural
  change: `"scale_free"` (denser Barabasi-Albert blocks),
  `"erdos_renyi"` (random blocks), `"ar1_banded"` (non-modular banded
  network), `"censored_5"`, `"censored_25"`, `"censored_50"` (that share
  of each feature's values set to zero, as at a detection limit), and
  `"uniform_weights"`.

- simulate:

  If `TRUE`, answer off-grid designs by on-demand simulation (slow;
  default scenario and simulatable methods only). If `FALSE` (default),
  they are flagged instead.

- ...:

  Passed to
  [`simulate_power()`](https://mljaniczek.github.io/howbig/reference/simulate_power.md)
  when simulating (e.g. `reps`, `n_perm`, `cores`, `seed`).

## Value

An object of class `howbig_power`: a list with the sensitivity estimate
and 80% band, achieved FDR, AUPR, the number of truly differential
edges, provenance, and any notes.

## Details

Answers come from a precomputed simulation grid where possible. Queries
that hit a grid cell are exact; queries on a swept arm are interpolated.
Anything else (other designs, other unequal group sizes, FDR targets
other than the grid's 0.05) is flagged `simulate_required`, and answered
by an on-demand simulation when `simulate = TRUE` (see
[`simulate_power()`](https://mljaniczek.github.io/howbig/reference/simulate_power.md)).

## Examples

``` r
power_for(p = 250, n = 1000, method = "JGL_BIC")
#> howbig power estimate  [in-grid]
#>   JGL_BIC | structural perturbation | p = 250, n = 1000
#>   Truly differential edges: ~50 of 31,125 pairs (0.16%)
#>   Sensitivity @ FDR 0.05:  0.60   (80% band 0.50–0.70)
#>   Achieved FDR:            0.059 ✓
#>   AUPR (threshold-free):   0.68
#>   based on 100 reps
power_for(p = 250, n = 1500, method = "JGL_BIC")          # interpolated
#> howbig power estimate  [interpolated]
#>   JGL_BIC | structural perturbation | p = 250, n = 1500
#>   Truly differential edges: ~50 of 31,125 pairs (0.16%)
#>   Sensitivity @ FDR 0.05:  0.66   (80% band 0.56–0.77)
#>   Achieved FDR:            0.063 ✓
#>   AUPR (threshold-free):   0.73
power_for(p = 100, n = c(500, 1500), method = "JGL_BIC")  # unequal groups
#> howbig power estimate  [in-grid]
#>   JGL_BIC | structural perturbation | p = 100, n = 500 / 1500
#>   Truly differential edges: ~45 of 4,950 pairs (0.92%)
#>   Sensitivity @ FDR 0.05:  0.62   (80% band 0.51–0.75)
#>   Achieved FDR:            0.047 ✓
#>   AUPR (threshold-free):   0.75
#>   based on 100 reps
power_for(p = 100, n = 1000, method = "JGL_BIC", scenario = "censored_25")
#> howbig power estimate  [in-grid]
#>   JGL_BIC | structural perturbation | p = 100, n = 1000 | scenario: censored_25
#>   Truly differential edges: ~45 of 4,950 pairs (0.92%)
#>   Sensitivity @ FDR 0.05:  0.65   (80% band 0.53–0.76)
#>   Achieved FDR:            0.092 ✓
#>   AUPR (threshold-free):   0.75
#>   based on 100 reps
power_for(p = 150, n = 1500, method = "JGL_BIC")          # off-grid: flagged
#> howbig power estimate  [OFF-GRID]
#>   JGL_BIC | structural perturbation | p = 150, n = 1500
#>   Not available from the grid.
#> note: not answerable from the grid: (p = 150, n = 1500) is off the grid
#>         for the structural perturbation.
#> note: rerun with simulate = TRUE to simulate this design (slow; see
#>         ?simulate_power).
```
