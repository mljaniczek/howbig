# howbig

**How big a sample do you need to find what changed in a network?**

howbig gives power and sample-size calculations for differential network
analysis: comparing two groups’ Gaussian graphical models to find the
edges that differ. It is a network-native, method-aware analogue to
[pwr](https://cran.r-project.org/package=pwr).

Its answers are empirical. They come from a large simulation study that
ran seven estimation methods across sample sizes, feature counts, and
kinds of differential signal, plus an on-demand simulator for designs
the study didn’t cover. That evidence shows things a formula would hide:

- **Required sample sizes are large.** With 250 features, the joint
  graphical lasso needs roughly 4,500 samples per group to recover 80%
  of the differential edges.
- **Method choice matters, and no method wins everywhere.**
- **Some methods lose FDR control under some kinds of signal.** howbig
  reports the false discovery rate a method actually achieves, next to
  its power, so that assumption is visible.

## Installation

``` R
# install.packages("remotes")
remotes::install_github("mljaniczek/howbig")
```

Looking up and interpolating the grid needs nothing else. To *simulate*
a method you also need its package; `dependencies = TRUE` installs them
all (several come from GitHub):

``` R
remotes::install_github("mljaniczek/howbig", dependencies = TRUE)
```

## Example

You need four things: roughly how many features you’ll measure (`p`), a
sample size per group (`n`), the method you plan to use, and what kind
of change you expect between groups.

``` R
library(howbig)

power_for(p = 250, n = 1000, method = "JGL_BIC")
#> howbig power estimate  [in-grid]
#>   JGL_BIC | structural perturbation | p = 250, n = 1000
#>   Sensitivity @ FDR 0.05:  0.60   (80% band 0.50–0.70)
#>   Achieved FDR:            0.059 ✓
#>   AUPR (threshold-free):   0.68
#>   based on 100 reps
```

*Power* here is **sensitivity**: the fraction of truly differential
edges detected while controlling the false discovery rate at 5%. To go
the other way, ask for the sample size that reaches a target power:

``` R
n_for(p = 250, power = 0.8, method = "JGL_BIC")
#> howbig sample-size estimate  [interpolated]
#>   JGL_BIC | structural perturbation | p = 250 | target power = 0.80 @ FDR 0.05
#>   Required per-group n:  4512
```

If you expect edge **strengths** to change rather than edges to appear
or disappear, the picture changes. This method stops controlling its
false discovery rate:

``` R
power_for(p = 250, n = 1000, method = "JGL_BIC", perturbation = "magnitude")
#> howbig power estimate  [in-grid]
#>   JGL_BIC | magnitude perturbation | p = 250, n = 1000
#>   Sensitivity @ FDR 0.05:  0.52   (80% band 0.33–0.71)
#>   Achieved FDR:            0.237 ⚠
#>   AUPR (threshold-free):   0.60
#>   based on 100 reps
#> note: achieved FDR 0.24 exceeds target 0.05 — this method is not
#>         calibrated here.
```

## How answers are produced

Every answer states where it came from:

| Provenance | Meaning | Speed |
|----|----|----|
| `in-grid` | A precomputed simulation cell (typically 100 replicates) | instant |
| `interpolated` | Between grid cells along a swept axis | instant |
| `simulated` | A fresh simulation of your exact design | minutes to hours |
| `OFF-GRID` | Not answerable from the grid; howbig says why | \- |

Designs the grid doesn’t cover are flagged, not guessed:

``` R
power_for(p = 150, n = 1500, method = "JGL_BIC")
#> howbig power estimate  [OFF-GRID]
#>   JGL_BIC | structural perturbation | p = 150, n = 1500
#>   Not available from the grid.
#> note: not answerable from the grid: (p = 150, n = 1500) is off the grid
#>         for the structural perturbation.
#> note: rerun with simulate = TRUE to simulate this design (slow; see
#>         ?simulate_power).
```

Rerunning with `simulate = TRUE` simulates that design, as do unequal
group sizes (`n = c(600, 1400)`) and FDR targets other than 5%. Each
simulation estimates its runtime before it starts and refuses to start
runs longer than `max_minutes`.

## Learn more

- [Get
  started](https://mljaniczek.github.io/howbig/articles/howbig.html): a
  walkthrough of power, sample size, and simulation.
- [Methodology](https://mljaniczek.github.io/howbig/articles/methodology.html):
  how the networks are generated, how each method is scored, and how
  howbig was validated.

## Status

howbig is experimental. The grid it ships with is **interim**: it was
simulated on Barabási-Albert *tree* networks. It will be replaced by a
grid on denser scale-free networks (mean degree ~3.8, matching Danaher
et al. 2014), and numbers will shift when that lands.
