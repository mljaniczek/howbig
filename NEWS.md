# howbig 0.0.0.9000

* `power_for()` and `n_for()` answer power and sample-size questions from a
  precomputed simulation grid, interpolating along its swept axes.
* `simulate_power()` simulates designs the grid does not cover: any number of
  features, unequal group sizes, and any FDR target. It estimates runtime
  before starting and refuses runs longer than `max_minutes`.
* Every answer reports how it was produced (in-grid, interpolated, simulated,
  or off-grid), how many edges truly differ, and flags methods that miss
  their FDR target.
* The grid comes from the November 2026 production run on Barabási-Albert
  tree networks: nine methods (`howbig_methods()`), a sample-size arm down to
  n = 100, small-sample cells at p = 100, and an anchor (p = 100, n = 1000)
  with all perturbation types and unequal groups (500 / 1500).
* New `scenario` argument: check an answer at the anchor under other
  assumptions (denser, random, or banded networks; censoring at a detection
  limit; uniform edge weights).
* multiGGM simulations now call differential edges with the grid's
  magnitude rule, P(|Θ₁ − Θ₂| > 0.1).
