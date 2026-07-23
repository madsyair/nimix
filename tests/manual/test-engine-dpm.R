# These tests compile NIMBLE models (C++), so they are slow and skipped on CRAN
# and when nimble is unavailable. They form the Layer-3 smoke test of the
# engineering harness.

test_that("DPM smoke test: compiles and runs a short chain", {
  skip_on_cran()
  skip_if_not_installed("nimble")

  set.seed(1)
  y <- c(rnorm(60, -4, 1), rnorm(60, 4, 1))
  fit <- nimixClust(y, K_max = 6,
                     mcmcControl = list(niter = 600, nburnin = 100, thin = 1),
                     seed = 1, verbose = FALSE)
  expect_s4_class(fit, "FitResult")
  expect_equal(nrow(fit@clusterAllocation), nrow(fit@mcmcSamples))
  expect_equal(ncol(fit@clusterAllocation), length(y))
  expect_true(all(fit@Kposterior >= 1L))
  expect_true(all(fit@Kposterior <= 6L))   # cannot exceed truncation
})

test_that("empty components: K_max >> K_true runs without numeric error", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  # the test plan explicit scenario
  set.seed(2)
  y <- c(rnorm(50, -5, 1), rnorm(50, 5, 1))   # K_true = 2
  expect_error(
    nimixClust(y, K_max = 10,
                mcmcControl = list(niter = 400, nburnin = 100),
                seed = 2, verbose = FALSE),
    NA
  )
})
