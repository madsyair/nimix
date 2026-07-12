# Tests for the univariate Student-t component. Prior/validation/density are
# pure R; recovery is nimble-gated.

test_that("Student-t defaultPrior is data-scaled and validates df", {
  set.seed(1); y <- rt(200, df = 4)
  spec <- StudentTUvSpec()
  pr <- defaultPrior(spec, y)
  expect_equal(pr$mu0, mean(y))
  expect_true(pr$muSd > 0 && pr$aTau > 0 && pr$bTau > 0)
  expect_equal(pr$df, 4)
  expect_silent(validateParams(spec, pr))
  expect_error(defaultPrior(spec, y, control = list(df = 2)), "df must exceed 2")
})

test_that("Student-t componentDensity is location-scale t", {
  d <- nimix:::componentDensity(StudentTUvSpec(), df = 5)
  expect_equal(d(2, list(mu = 1, tau = 1 / 4, df = 5)),
               stats::dt((2 - 1) / 2, df = 5) / 2, tolerance = 1e-10)
})

test_that("Student-t DPM recovers two heavy-tailed clusters", {
  skip_on_cran(); skip_if_not_installed("nimble")
  set.seed(5); y <- c(rt(70, df = 4) - 5, rt(70, df = 4) + 5)
  fit <- nimixClust(y, K_max = 8, distribution = "studentt", method = "dpm",
                    prior = list(df = 4),
                    mcmcControl = list(niter = 1500, nburnin = 600),
                    seed = 5, verbose = FALSE)
  fit <- relabel(fit)
  mu <- sort(fit@relabeled$summary$mu_mean)
  expect_true(mu[1] < -2 && mu[length(mu)] > 2)
})
