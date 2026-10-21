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
  simulate = FALSE,
  ...
)
```

## Arguments

- p:

  Number of features.

- n:

  Per-group sample size; a single number, or `c(n1, n2)` for unequal
  groups (always simulated).

- method:

  Estimation method label. See
  [`howbig_coverage()`](https://mljaniczek.github.io/howbig/reference/howbig_coverage.md)
  for the available methods.

- perturbation:

  Type of differential signal: `"structural"` (edges appear/disappear),
  `"magnitude"` (edge weights change), or `"both"`.

- fdr:

  Target false-discovery rate. The grid is scored at `0.05`; other
  targets require simulation.

- simulate:

  If `TRUE`, answer off-grid designs by on-demand simulation (slow). If
  `FALSE` (default), they are flagged instead.

- ...:

  Passed to
  [`simulate_power()`](https://mljaniczek.github.io/howbig/reference/simulate_power.md)
  when simulating (e.g. `reps`, `n_perm`, `cores`, `seed`).

## Value

An object of class `howbig_power`: a list with the sensitivity estimate
and 80% band, achieved FDR, AUPR, provenance, and any notes.

## Details

Answers come from a precomputed simulation grid where possible. Queries
that hit a grid cell are exact; queries on a swept arm are interpolated.
Anything else (other designs, unequal group sizes, FDR targets other
than the grid's 0.05) is flagged `simulate_required`, and answered by an
on-demand simulation when `simulate = TRUE` (see
[`simulate_power()`](https://mljaniczek.github.io/howbig/reference/simulate_power.md)).

## Examples

``` r
power_for(p = 250, n = 1000, method = "JGL_BIC")
#> howbig power estimate  [in-grid]
#>   JGL_BIC | structural perturbation | p = 250, n = 1000
#>   Sensitivity @ FDR 0.05:  0.60   (80% band 0.50–0.70)
#>   Achieved FDR:            0.059 ✓
#>   AUPR (threshold-free):   0.68
#>   based on 100 reps
power_for(p = 250, n = 1500, method = "JGL_BIC")   # interpolated
#> howbig power estimate  [interpolated]
#>   JGL_BIC | structural perturbation | p = 250, n = 1500
#>   Sensitivity @ FDR 0.05:  0.66   (80% band 0.56–0.77)
#>   Achieved FDR:            0.063 ✓
#>   AUPR (threshold-free):   0.73
power_for(p = 150, n = 1500, method = "JGL_BIC")   # off-grid: flagged
#> howbig power estimate  [OFF-GRID]
#>   JGL_BIC | structural perturbation | p = 150, n = 1500
#>   Not available from the grid.
#> note: not answerable from the grid: (p = 150, n = 1500) is off the grid
#>         for the structural perturbation.
#> note: rerun with simulate = TRUE to simulate this design (slow; see
#>         ?simulate_power).
```
