test_that("block sizes match the DiNA registry and cover any p", {
  expect_equal(.block_sizes(250), rep(25L, 10))
  expect_equal(.block_sizes(50), c(25L, 25L))
  expect_equal(.block_sizes(20), c(10L, 10L))
  for (p in c(37, 130, 333)) expect_equal(sum(.block_sizes(p)), p)
})

test_that("generator reproduces the grid's differential prevalence", {
  # DiNA: ~3.2% differential pairs at p = 50, ~0.04% at p = 500
  set.seed(1)
  prev50 <- mean(replicate(5, mean(.generate_gold(50, "structural")$delta)))
  prev500 <- mean(.generate_gold(500, "structural")$delta)
  expect_equal(prev50, 0.032, tolerance = 0.25)
  expect_equal(prev500, 0.0004, tolerance = 0.3)
  # magnitude-only changes no structure but still yields differential truth
  expect_gt(sum(.generate_gold(250, "magnitude")$delta), 0)
})

test_that("InvCov simulation matches its grid cell within Monte Carlo error", {
  grid <- power_for(100, 1000, "InvCov")
  sim <- simulate_power(100, 1000, "InvCov", reps = 12, n_perm = 12,
                        seed = 3, verbose = FALSE)
  expect_equal(sim$provenance, "simulated")
  expect_lt(abs(sim$sensitivity - grid$sensitivity), 4 * sim$sensitivity_mcse + 0.03)
  expect_lt(sim$achieved_fdr, 0.10)
})

test_that("simulation is reproducible and cached", {
  a <- simulate_power(60, 400, "InvCov", reps = 4, n_perm = 4, seed = 9,
                      verbose = FALSE)
  t <- system.time(b <- simulate_power(60, 400, "InvCov", reps = 4, n_perm = 4,
                                       seed = 9, verbose = FALSE))["elapsed"]
  expect_identical(a, b)
  expect_lt(t, 0.5)
})

test_that("power_for(simulate = TRUE) answers unequal groups and other FDRs", {
  r <- power_for(60, c(300, 700), "InvCov", fdr = 0.1, simulate = TRUE,
                 reps = 4, n_perm = 4, verbose = FALSE)
  expect_equal(r$provenance, "simulated")
  expect_equal(r$n, c(300L, 700L))
  expect_equal(r$fdr_target, 0.1)
})

test_that("parallel workers reproduce the serial result exactly", {
  # Workers are separate R processes that load the *installed* package, which
  # is the version under test only inside R CMD check (interactively it may
  # be a stale copy).
  skip_if(Sys.getenv("_R_CHECK_PACKAGE_NAME_") == "", "runs under R CMD check")
  skip_on_cran()
  a <- simulate_power(60, 300, "InvCov", reps = 4, n_perm = 4, seed = 21,
                      cores = 1, verbose = FALSE)
  rm(list = ls(.howbig_cache), envir = .howbig_cache)
  b <- simulate_power(60, 300, "InvCov", reps = 4, n_perm = 4, seed = 21,
                      cores = 2, verbose = FALSE)
  expect_equal(a$sensitivity, b$sensitivity)
  expect_equal(a$aupr, b$aupr)
})

test_that("runtime guard refuses designs that would take too long", {
  skip_if_not_installed("multiGGM")
  expect_error(simulate_power(250, 1000, "multiGGM", max_minutes = 60),
               "exceeds max_minutes")
})
