# Parallel chains: each worker compiles its own model (fork-safe) and results
# match the sequential path in the label-invariant summary.

test_that("parallel chains recover the same parameters as sequential", {
  skip_on_cran()
  skip_on_os("windows")            # forking only
  set.seed(1)
  y <- c(rnorm(100, -3, 1), rnorm(100, 3, 1))

  nimixClearCache()
  fs <- nimixClust(y, K = 2, method = "fixedk",
                   mcmcControl = list(niter = 1200, nburnin = 500,
                                      nchains = 2, parallel = FALSE), seed = 5)
  nimixClearCache()
  fp <- nimixClust(y, K = 2, method = "fixedk",
                   mcmcControl = list(niter = 1200, nburnin = 500, nchains = 2,
                                      parallel = TRUE, ncores = 2), seed = 5)

  expect_equal(nrow(fs@mcmcSamples), nrow(fp@mcmcSamples))
  mus <- sort(relabel(fs)@relabeled$summary$mu_mean)
  mup <- sort(relabel(fp)@relabeled$summary$mu_mean)
  # both recover the two centres near -3 and 3
  expect_lt(abs(mup[1] - (-3)), 0.6)
  expect_lt(abs(mup[2] - 3), 0.6)
  # sequential and parallel agree to Monte Carlo error
  expect_lt(max(abs(mus - mup)), 0.75)
})
