# Required sample size for a target power

Inverts
[`power_for()`](https://mljaniczek.github.io/howbig/reference/power_for.md):
finds the smallest per-group sample size `n` that reaches a target power
(sensitivity) at a feature count `p`.

## Usage

``` r
n_for(
  p,
  power = 0.8,
  method,
  perturbation = "structural",
  fdr = .HOWBIG_FDR_TARGET,
  simulate = FALSE,
  n_range = c(100, 5000),
  n_points = 5L,
  ...
)
```

## Arguments

- p:

  Number of features.

- power:

  Target sensitivity in (0, 1).

- method:

  Estimation method label.

- perturbation:

  Type of differential signal; see
  [`power_for()`](https://mljaniczek.github.io/howbig/reference/power_for.md).

- fdr:

  Target false-discovery rate. The grid is scored at `0.05`; other
  targets require simulation.

- simulate:

  If `TRUE`, answer off-grid designs by simulation (slow: `n_points`
  full
  [`simulate_power()`](https://mljaniczek.github.io/howbig/reference/simulate_power.md)
  runs).

- n_range, n_points:

  Per-group sample sizes searched when simulating.

- ...:

  Passed to
  [`simulate_power()`](https://mljaniczek.github.io/howbig/reference/simulate_power.md)
  (e.g. `reps`, `n_perm`, `cores`).

## Value

An object of class `howbig_n`: the required `n`, provenance, the `n`
range searched, the power curve used, and any notes.

## Details

From the grid, the inversion uses the precomputed sample-size arm, which
exists at `p = 250` for the `structural` perturbation at FDR 0.05. Other
designs are flagged `simulate_required`; with `simulate = TRUE` they are
answered by simulating power at `n_points` log-spaced sample sizes
across `n_range`, smoothing the curve to be monotone, and inverting it.

## Examples

``` r
n_for(p = 250, power = 0.8, method = "JGL_BIC")
#> howbig sample-size estimate  [interpolated]
#>   JGL_BIC | structural perturbation | p = 250 | target power = 0.80 @ FDR 0.05
#>   Required per-group n:  4512
```
