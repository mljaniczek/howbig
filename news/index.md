# Changelog

## howbig 0.0.0.9000

- [`power_for()`](https://mljaniczek.github.io/howbig/reference/power_for.md)
  and
  [`n_for()`](https://mljaniczek.github.io/howbig/reference/n_for.md)
  answer power and sample-size questions from a precomputed simulation
  grid, interpolating along its swept axes.
- [`simulate_power()`](https://mljaniczek.github.io/howbig/reference/simulate_power.md)
  simulates designs the grid does not cover: any number of features,
  unequal group sizes, and any FDR target. It estimates runtime before
  starting and refuses runs longer than `max_minutes`.
- Every answer reports how it was produced (in-grid, interpolated,
  simulated, or off-grid) and flags methods that miss their FDR target.
- The shipped grid is interim (Barabási-Albert tree networks) and will
  be replaced by the scale-free family.
