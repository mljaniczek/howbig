# Coverage of the built-in grid

Returns the precomputed cells `howbig` can answer directly (before any
interpolation or simulation), one row per (method, p, n, perturbation).

## Usage

``` r
howbig_coverage(method = NULL)
```

## Arguments

- method:

  Optional method label to filter to.

## Value

A data frame of grid cells with mean sensitivity, achieved FDR and AUPR,
and the replicate count behind each cell.
