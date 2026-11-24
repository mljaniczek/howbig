# Methods available in howbig

Lists the estimation methods in the grid, how each one decides which
edges are differential at the FDR target, and whether
[`simulate_power()`](https://mljaniczek.github.io/howbig/reference/simulate_power.md)
can run it for designs the grid does not cover.

## Usage

``` r
howbig_methods()
```

## Value

A data frame with columns `method`, `description`, `inference` and
`simulate`.

## Examples

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
