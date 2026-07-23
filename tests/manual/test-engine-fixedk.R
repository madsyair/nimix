# Tests for the fixed-K finite-mixture engine. Argument guards are fast (pure
# R); the recovery test is skipped unless nimble is installed and not on CRAN.

test_that("K / K_max are not interchangeable", {
  y <- rnorm(20)
  # fixedk needs K, rejects K_max
  expect_error(nimixClust(y, method = "fixedk"), "needs the number of components")
  expect_error(nimixClust(y, K = 2, K_max = 5, method = "fixedk"),
               "Use K \\(not K_max\\)")
  # dpm rejects K
  expect_error(nimixClust(y, K = 2, method = "dpm"), "Use K_max \\(not K\\)")
  # same for the regression entry point
  df <- data.frame(y = rnorm(20), x = rnorm(20))
  expect_error(nimixReg(y ~ x, df, method = "fixedk"),
               "needs the number of components")
  expect_error(nimixReg(y ~ x, df, K = 2, method = "dpm"), "Use K_max")
})

test_that("FixedKEngine constructor validates dirichletConc", {
  expect_s4_class(FixedKEngine(), "FixedKEngine")
  expect_equal(FixedKEngine(2)@dirichletConc, 2)
  expect_error(FixedKEngine(0), "positive")
  expect_error(FixedKEngine(c(1, 2)), "positive scalar")
})

test_that("fixed-K finite mixture recovers two univariate components", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  set.seed(7)
  y <- c(rnorm(70, -4, 1), rnorm(70, 4, 1))
  fit <- nimixClust(y, K = 2, method = "fixedk",
                    mcmcControl = list(niter = 1500, nburnin = 600),
                    seed = 7, verbose = FALSE)
  expect_identical(fit@engineUsed, "fixedk")
  fit <- relabel(fit)
  means <- sort(fit@relabeled$summary$mu_mean)
  expect_true(means[1] < -2 && means[2] > 2)
})
