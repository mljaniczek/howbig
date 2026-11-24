test_that("in-grid cells are returned exactly", {
  r <- power_for(250, 1000, "JGL_BIC")
  expect_s3_class(r, "howbig_power")
  expect_equal(r$provenance, "in_grid")
  expect_gte(r$n_reps, 50)
  expect_true(r$sensitivity > 0 && r$sensitivity < 1)
})

test_that("results report how many edges truly differ", {
  r <- power_for(500, 1000, "JGL_BIC")
  expect_equal(r$n_pairs, 500 * 499 / 2)
  expect_true(r$n_diff_edges > 10 && r$n_diff_edges < 200)
  expect_output(print(r), "Truly differential edges")
})

test_that("on-arm queries interpolate between their neighbours", {
  lo <- power_for(250, 1000, "JGL_BIC")$sensitivity
  hi <- power_for(250, 2000, "JGL_BIC")$sensitivity
  mid <- power_for(250, 1500, "JGL_BIC")
  expect_equal(mid$provenance, "interpolated")
  expect_true(mid$sensitivity > lo && mid$sensitivity < hi)
  expect_equal(power_for(150, 1000, "JGL_BIC")$provenance, "interpolated")
})

test_that("unequal groups are answered exactly, in order, or flagged", {
  r <- power_for(100, c(500, 1500), "JGL_BIC")
  expect_equal(r$provenance, "in_grid")
  expect_equal(r$n, c(500L, 1500L))
  expect_lt(r$sensitivity, power_for(100, 1000, "JGL_BIC")$sensitivity)
  # group 1 is the reference group: the reversed split is a different design
  rev <- power_for(100, c(1500, 500), "JGL_BIC")
  expect_equal(rev$provenance, "simulate_required")
  expect_match(paste(rev$notes, collapse = " "), "unequal group sizes")
})

test_that("assumption scenarios exist only at the anchor", {
  r <- power_for(100, 1000, "JGL_BIC", scenario = "censored_25")
  expect_equal(r$provenance, "in_grid")
  expect_equal(r$scenario, "censored_25")
  off <- power_for(250, 1000, "JGL_BIC", scenario = "erdos_renyi")
  expect_equal(off$provenance, "simulate_required")
  expect_match(paste(off$notes, collapse = " "), "only at p = 100")
  expect_error(power_for(250, 1000, "JGL_BIC", scenario = "erdos_renyi",
                         simulate = TRUE), "only scenario")
  expect_error(power_for(100, 1000, "JGL_BIC", scenario = "sideways"))
})

test_that("off-grid designs are flagged with the reason, not guessed", {
  r <- power_for(150, 1500, "JGL_BIC")
  expect_equal(r$provenance, "simulate_required")
  expect_true(is.na(r$sensitivity))
  expect_match(paste(r$notes, collapse = " "), "off the grid")
  expect_match(paste(power_for(250, 1000, "JGL_BIC", fdr = 0.1)$notes,
                     collapse = " "), "FDR target")
  # perturbation types are at the p = 100 anchor
  expect_equal(power_for(100, 1000, "JGL_BIC", "magnitude")$provenance, "in_grid")
  expect_equal(power_for(250, 1000, "JGL_BIC", "magnitude")$provenance,
               "simulate_required")
})

test_that("grid-only methods answer from the grid but refuse to simulate", {
  expect_equal(power_for(100, 1000, "mGHS")$provenance, "in_grid")
  off <- power_for(150, 1500, "mGHS")
  expect_match(paste(off$notes, collapse = " "), "grid only")
  expect_error(simulate_power(100, 1000, "mGHS", verbose = FALSE), "grid only")
  m <- howbig_methods()
  expect_setequal(m$method[m$simulate],
                  c("JGL_BIC", "JGL_BIC_GGL", "SSJGL", "SSJGL_GGL", "multiGGM",
                    "InvCov"))
})

test_that("thin cells are flagged", {
  thin <- power_for(500, 1000, "multiGGM")   # 2 replicates in the grid
  expect_match(paste(thin$notes, collapse = " "), "only 2 replicates")
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
