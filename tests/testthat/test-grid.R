test_that("in-grid cells are returned exactly", {
  r <- power_for(250, 1000, "JGL_BIC")
  expect_s3_class(r, "howbig_power")
  expect_equal(r$provenance, "in_grid")
  expect_equal(r$n_reps, 100)
  expect_true(r$sensitivity > 0 && r$sensitivity < 1)
})

test_that("on-arm queries interpolate between their neighbours", {
  lo <- power_for(250, 1000, "JGL_BIC")$sensitivity
  hi <- power_for(250, 2000, "JGL_BIC")$sensitivity
  mid <- power_for(250, 1500, "JGL_BIC")
  expect_equal(mid$provenance, "interpolated")
  expect_true(mid$sensitivity > lo && mid$sensitivity < hi)
  expect_equal(power_for(150, 1000, "JGL_BIC")$provenance, "interpolated")
})

test_that("off-grid designs are flagged with the reason, not guessed", {
  r <- power_for(150, 1500, "JGL_BIC")
  expect_equal(r$provenance, "simulate_required")
  expect_true(is.na(r$sensitivity))
  expect_match(paste(r$notes, collapse = " "), "off the grid")

  expect_match(paste(power_for(250, c(500, 1500), "JGL_BIC")$notes, collapse = " "),
               "unequal group sizes")
  expect_match(paste(power_for(250, 1000, "JGL_BIC", fdr = 0.1)$notes, collapse = " "),
               "FDR target")
  # magnitude perturbation exists only at the anchor
  expect_equal(power_for(250, 1000, "JGL_BIC", "magnitude")$provenance, "in_grid")
  expect_equal(power_for(250, 1500, "JGL_BIC", "magnitude")$provenance,
               "simulate_required")
})

test_that("miscalibration and metric holes are surfaced", {
  mag <- power_for(250, 1000, "JGL_BIC", "magnitude")
  expect_true(mag$achieved_fdr > 0.1)
  expect_match(paste(mag$notes, collapse = " "), "not calibrated")

  hole <- power_for(500, 1000, "JGL_BIC")
  expect_true(is.na(hole$sensitivity) && is.finite(hole$aupr))
})

test_that("n_for inverts the sample-size arm monotonically", {
  n70 <- n_for(250, 0.70, "JGL_BIC")
  n80 <- n_for(250, 0.80, "JGL_BIC")
  expect_true(n70$n < n80$n)
  expect_equal(n_for(100, 0.8, "JGL_BIC")$provenance, "simulate_required")
  expect_true(is.na(n_for(250, 0.99, "JGL_BIC")$n))
})

test_that("bad inputs error clearly", {
  expect_error(power_for(250, 1000, "NOPE"), "must be one of")
  expect_error(power_for(250, 1000, "JGL_BIC", "sideways"))
})
