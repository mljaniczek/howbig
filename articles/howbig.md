# Get started with howbig

``` r

library(howbig)
```

## The question

Differential network analysis compares two or more groups (say, cases
and controls) and asks which edges in a graphical model differ between
them. `howbig` answers the planning question that comes first: **how
many samples do I need to find those edges?**

### What “power” means here

In a differential network, almost every possible edge is *not*
differential. With 500 features there are about 125,000 feature pairs,
and in the simulations behind `howbig` only around 50 of them differ. A
method that calls nothing differential still has a false positive rate
near zero, so FPR-based summaries say little.

`howbig` defines **power as sensitivity**: the fraction of truly
differential edges a method detects while controlling the false
discovery rate (FDR) at a target, 5% by default. It always reports the
FDR the method **actually achieved** next to it, because not every
method controls FDR in every setting.

## What you need to know

| Input | Argument | Typical knowledge |
|----|----|----|
| Number of features | `p` | Usually known approximately |
| Samples per group | `n` | A budget range; `c(n1, n2)` for unequal groups |
| Kind of change | `perturbation` | `"structural"` (edges appear/disappear), `"magnitude"` (strengths change), or `"both"` |
| FDR target | `fdr` | 0.05 by default |
| Method | `method` | The method you plan to use |
| Assumptions | `scenario` | Optional: network shape and data quality, for checking how much they matter |

``` r

howbig_methods()
#>        method                                               description
#> 1     JGL_BIC           Joint graphical lasso, fused penalty, BIC-tuned
#> 2 JGL_BIC_GGL           Joint graphical lasso, group penalty, BIC-tuned
#> 3       SSJGL       Spike-and-slab joint graphical lasso, fused penalty
#> 4   SSJGL_GGL       Spike-and-slab joint graphical lasso, group penalty
#> 5    multiGGM Bayesian joint GGM (stochastic search variable selection)
#> 6      InvCov               Unregularized inverse covariance (baseline)
#> 7        mGHS                     Bayesian multiple graphical horseshoe
#> 8  DiffNetFDR          Partial-correlation difference test (DiffNetFDR)
#> 9       ONDSA        Debiased (FastGGM-based) differential network test
#>                inference simulate
#> 1  permutation null + BH     TRUE
#> 2  permutation null + BH     TRUE
#> 3  permutation null + BH     TRUE
#> 4  permutation null + BH     TRUE
#> 5          posterior FDR     TRUE
#> 6  permutation null + BH     TRUE
#> 7          posterior FDR    FALSE
#> 8 built-in FDR procedure    FALSE
#> 9 built-in FDR procedure    FALSE
```

All methods were tuned the way they were in the simulation study; see
the [methodology
article](https://mljaniczek.github.io/howbig/articles/methodology.md).
Six of them can also be simulated for designs the grid does not cover;
`mGHS`, `DiffNetFDR` and `ONDSA` are available from the grid only.

## Power for a design

``` r

power_for(p = 250, n = 1000, method = "JGL_BIC")
#> howbig power estimate  [in-grid]
#>   JGL_BIC | structural perturbation | p = 250, n = 1000
#>   Truly differential edges: ~50 of 31,125 pairs (0.16%)
#>   Sensitivity @ FDR 0.05:  0.60   (80% band 0.50–0.70)
#>   Achieved FDR:            0.059 ✓
#>   AUPR (threshold-free):   0.68
#>   based on 100 reps
```

The output reads:

- **Truly differential edges**: how many edges differ between the groups
  in the simulated networks, out of all feature pairs. Power is measured
  against these.
- **Sensitivity @ FDR 0.05**: mean power across simulated replicates,
  with the 10th-90th percentile band. The band shows how much power
  varies from one network to another, not uncertainty in the mean.
- **Achieved FDR**: the false discovery rate the method actually
  reached; ✓ when within twice the target, ⚠ otherwise.
- **AUPR**: area under the precision-recall curve for the method’s
  ranking of edges, a threshold-free measure of how well it separates
  differential edges from the rest.

### Calibration is part of the answer

Change the kind of signal to a pure strength change and the same method
misses its FDR target. Its sensitivity looks respectable, but about one
in six of the edges it calls differential are not:

``` r

power_for(p = 100, n = 1000, method = "JGL_BIC", perturbation = "magnitude")
#> howbig power estimate  [in-grid]
#>   JGL_BIC | magnitude perturbation | p = 100, n = 1000
#>   Truly differential edges: ~24 of 4,950 pairs (0.48%)
#>   Sensitivity @ FDR 0.05:  0.64   (80% band 0.37–0.84)
#>   Achieved FDR:            0.168 ⚠
#>   AUPR (threshold-free):   0.74
#>   based on 100 reps
#> note: achieved FDR 0.17 exceeds target 0.05 — this method is not
#>         calibrated here.
```

## Sample size for a target power

[`n_for()`](https://mljaniczek.github.io/howbig/reference/n_for.md)
inverts the calculation:

``` r

n_for(p = 250, power = 0.8, method = "JGL_BIC")
#> howbig sample-size estimate  [interpolated]
#>   JGL_BIC | structural perturbation | p = 250 | target power = 0.80 @ FDR 0.05
#>   Required per-group n:  4512
```

The curves behind that answer show how slowly power grows with sample
size, and that no method is best everywhere. The Bayesian methods
(`multiGGM`, and `mGHS`, which tracks it closely) lead at small samples
but level off near 0.6; the joint graphical lasso overtakes them from
about n = 1000.

![Line chart of sensitivity at FDR 0.05 against per-group sample size
(log scale, 100 to 5000) at p = 250, for five methods. multiGGM starts
highest at 0.30 for n = 100 but levels off near 0.63. JGL_BIC rises from
0.08 to 0.81, crossing multiGGM near n = 1000. SSJGL rises from 0.08 to
0.73, DiffNetFDR from 0.05 to 0.71, and InvCov from 0 to 0.69. A dashed
line marks the 0.8 target, which only JGL_BIC crosses, near n =
4500.](howbig_files/figure-html/power-curves-1.png)

The same numbers as a table, for all nine methods:

|             | n = 100 | n = 250 | n = 500 | n = 1000 | n = 2000 | n = 5000 |
|:------------|--------:|--------:|--------:|---------:|---------:|---------:|
| DiffNetFDR  |    0.05 |    0.17 |    0.29 |     0.44 |     0.56 |     0.71 |
| InvCov      |    0.00 |    0.01 |    0.07 |     0.33 |     0.51 |     0.69 |
| JGL_BIC     |    0.08 |    0.28 |    0.42 |     0.60 |     0.70 |     0.81 |
| JGL_BIC_GGL |    0.04 |    0.25 |    0.37 |     0.55 |     0.65 |     0.79 |
| mGHS        |    0.25 |    0.41 |    0.48 |     0.57 |     0.60 |     0.62 |
| multiGGM    |    0.30 |    0.40 |    0.49 |     0.58 |     0.61 |     0.63 |
| ONDSA       |    0.02 |    0.15 |    0.26 |     0.43 |     0.56 |     0.71 |
| SSJGL       |    0.08 |    0.26 |    0.37 |     0.52 |     0.64 |     0.73 |
| SSJGL_GGL   |    0.03 |    0.23 |    0.35 |     0.52 |     0.62 |     0.74 |

Sensitivity at FDR 0.05, p = 250, structural change {.table}

## Where the grid has answers

The grid covers two **arms** that cross at p = 250, n = 1000: sample
size varies at p = 250, and feature count varies at n = 1000. At p = 100
there are also small-sample cells (n = 50 and 100) and an **anchor** at
n = 1000 that carries everything else: the three perturbation types,
unequal groups (500 / 1500), and the assumption scenarios.

![Scatter plot of grid coverage, feature count p against samples per
group n on a log scale. Points form a cross: a vertical arm at p = 250
for n from 100 to 5000, and a horizontal arm at n = 1000 for p from 50
to 500, meeting at p = 250, n = 1000. At p = 100 there are extra cells
at n = 50 and n = 100, and the anchor at p = 100, n = 1000 is ringed: it
holds the perturbation types, unequal groups and assumption
scenarios.](howbig_files/figure-html/coverage-1.png)

A query **on an arm** between grid points is interpolated:

``` r

power_for(p = 150, n = 1000, method = "JGL_BIC")
#> howbig power estimate  [interpolated]
#>   JGL_BIC | structural perturbation | p = 150, n = 1000
#>   Truly differential edges: ~47 of 11,175 pairs (0.42%)
#>   Sensitivity @ FDR 0.05:  0.65   (80% band 0.54–0.77)
#>   Achieved FDR:            0.041 ✓
#>   AUPR (threshold-free):   0.76
```

Interpolation stays on the arms because that is where it has been
validated: re-simulating this design gives the same 0.65. Anywhere else,
`howbig` declines to guess and says why:

``` r

power_for(p = 150, n = 1500, method = "JGL_BIC")
#> howbig power estimate  [OFF-GRID]
#>   JGL_BIC | structural perturbation | p = 150, n = 1500
#>   Not available from the grid.
#> note: not answerable from the grid: (p = 150, n = 1500) is off the grid
#>         for the structural perturbation.
#> note: rerun with simulate = TRUE to simulate this design (slow; see
#>         ?simulate_power).
```

The same happens for unequal group sizes other than the gridded 500 /
1500 (in that order: group 1 is the reference group) and for FDR targets
other than 0.05.

## How sensitive is the answer to your assumptions?

The default grid assumes block-structured, sparse networks with Gaussian
data. The `scenario` argument shows how the answer at the anchor changes
under other assumptions, for any method:

``` r

power_for(p = 100, n = 1000, method = "multiGGM", scenario = "censored_50")
#> howbig power estimate  [in-grid]
#>   multiGGM | structural perturbation | p = 100, n = 1000 | scenario: censored_50
#>   Truly differential edges: ~45 of 4,950 pairs (0.92%)
#>   Sensitivity @ FDR 0.05:  0.01   (80% band 0.00–0.03)
#>   Achieved FDR:            0.000 ✓
#>   AUPR (threshold-free):   0.61
#>   based on 100 reps
```

Across methods, some assumptions matter a lot:

|                 | JGL_BIC | SSJGL | multiGGM | mGHS | DiffNetFDR |
|:----------------|--------:|------:|---------:|-----:|-----------:|
| default         |    0.70 |  0.64 |     0.65 | 0.64 |       0.56 |
| uniform_weights |    0.69 |  0.64 |     0.66 | 0.65 |       0.57 |
| scale_free      |    0.55 |  0.51 |     0.50 | 0.48 |       0.41 |
| erdos_renyi     |    0.57 |  0.52 |     0.52 | 0.49 |       0.43 |
| ar1_banded      |    0.11 |  0.06 |     0.07 | 0.03 |       0.05 |
| censored_5      |    0.68 |  0.60 |     0.56 | 0.54 |       0.53 |
| censored_25     |    0.65 |  0.40 |     0.19 | 0.17 |       0.47 |
| censored_50     |    0.58 |  0.18 |     0.01 | 0.00 |       0.41 |

Sensitivity at FDR 0.05 by scenario (p = 100, n = 1000, structural
change) {.table}

- **Censoring** (a share of each feature’s values set to zero, as at a
  detection limit) is the biggest threat. At 50% the Bayesian methods
  find almost nothing, while the joint graphical lasso keeps most of its
  power.
- **Denser or random network blocks** (`scale_free`, `erdos_renyi`) cost
  every method roughly 0.1-0.2.
- **Banded, non-modular networks** (`ar1_banded`): no method recovers
  the differences, and the Bayesian methods also miss their FDR target
  there.
- **Edge-weight distribution** (`uniform_weights`) barely matters.

## Simulating your design

`simulate = TRUE`, or
[`simulate_power()`](https://mljaniczek.github.io/howbig/reference/simulate_power.md)
directly, simulates the exact design. It generates fresh networks and
data the way the grid was built, fits the method, and scores its
FDR-controlled calls against the truth. Here is a small run with unequal
groups and a 10% FDR target. It uses the fast `InvCov` baseline and few
replicates so this page builds quickly:

``` r

simulate_power(p = 60, n = c(150, 350), method = "InvCov", fdr = 0.10,
               reps = 10, n_perm = 10, seed = 1, verbose = FALSE)
#> howbig power estimate  [simulated]
#>   InvCov | structural perturbation | p = 60, n = 150 / 350
#>   Truly differential edges: ~44 of 1,770 pairs (2.51%)
#>   Sensitivity @ FDR 0.10:  0.04   (80% band 0.00–0.09)
#>   Monte Carlo SE of mean:  0.018
#>   Achieved FDR:            0.133 ✓
#>   AUPR (threshold-free):   0.21
#>   based on 10 reps, 10-permutation null; 0.0 min
```

For real use, keep the defaults (`reps = 20`, `n_perm = 30`) or raise
them, and watch the Monte Carlo SE: it is the uncertainty in the mean
sensitivity. Simulation uses the default scenario.

### Runtime

Every simulated replicate and every permutation for the null
distribution is a full model fit, and fit costs vary enormously by
method:

| Method     | p = 100 | p = 250  | p = 500    |
|------------|---------|----------|------------|
| `InvCov`   | seconds | seconds  | seconds    |
| `JGL_*`    | ~6 s    | ~1.5 min | ~15 min    |
| `SSJGL_*`  | ~30 s   | ~5 min   | ~30-50 min |
| `multiGGM` | ~4 min  | ~45 min  | ~7 hours   |

*Median time per fit in the production simulation runs.*

Before starting,
[`simulate_power()`](https://mljaniczek.github.io/howbig/reference/simulate_power.md)
prints an estimated runtime and refuses to start if it exceeds
`max_minutes` (60 by default). Use `cores` to run fits in parallel
worker processes:

``` r

simulate_power(p = 150, n = 1500, method = "JGL_BIC", cores = 4)
#> howbig: 50 JGL_BIC fits, estimated ~5 min on 4 core(s).
```

Simulated results are cached for the session, so asking again is
instant. `n_for(..., simulate = TRUE)` searches several sample sizes, so
it costs several simulations.

## Assumptions to keep in mind

- Default networks are built from 25-feature blocks, each a sparse
  Barabási-Albert tree, with Gaussian data; signal is confined to one or
  two blocks. Real biology may be harder or easier; the scenarios above
  give a sense of how much.
- Methods use fixed tuning choices (BIC for JGL, a fixed prior ladder
  for SSJGL, fixed MCMC settings for the Bayesian methods). Different
  choices could change power.
