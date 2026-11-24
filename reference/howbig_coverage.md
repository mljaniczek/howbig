# Coverage of the built-in grid

Returns the precomputed cells `howbig` can answer directly, before any
interpolation or simulation: one row per method, scenario, perturbation,
number of features, and pair of group sizes.

## Usage

``` r
howbig_coverage(method = NULL, scenario = NULL)
```

## Arguments

- method:

  Optional method label to filter to.

- scenario:

  Optional scenario to filter to (see
  [`power_for()`](https://mljaniczek.github.io/howbig/reference/power_for.md)).

## Value

A data frame of grid cells with mean sensitivity, achieved FDR, AUPR,
the mean number of truly differential edges, the number of feature
pairs, and the replicate count behind each cell.

## Examples

``` r
head(howbig_coverage("JGL_BIC"))
#>    method    scenario perturbation   p   n1   n2 sensitivity_mean
#> 1 JGL_BIC  ar1_banded   structural 100 1000 1000        0.1141667
#> 2 JGL_BIC censored_25   structural 100 1000 1000        0.6472810
#> 3 JGL_BIC  censored_5   structural 100 1000 1000        0.6759562
#> 4 JGL_BIC censored_50   structural 100 1000 1000        0.5782838
#> 5 JGL_BIC     default         both 100 1000 1000        0.6773248
#> 6 JGL_BIC     default    magnitude 100 1000 1000        0.6395833
#>   achieved_fdr_mean aupr_mean n_diff_edges n_pairs n_reps
#> 1        0.02188347 0.5101047        24.00    4950    100
#> 2        0.09230322 0.7458427        45.49    4950    100
#> 3        0.03949096 0.7837786        45.49    4950    100
#> 4        0.08012585 0.6979188        45.49    4950    100
#> 5        0.10458554 0.7780510        70.23    4950    100
#> 6        0.16805053 0.7377327        24.00    4950    100
howbig_coverage("JGL_BIC", scenario = "erdos_renyi")
#>    method    scenario perturbation   p   n1   n2 sensitivity_mean
#> 1 JGL_BIC erdos_renyi   structural 100 1000 1000        0.5680751
#>   achieved_fdr_mean aupr_mean n_diff_edges n_pairs n_reps
#> 1        0.03751718 0.7277479        66.48    4950    100
```
