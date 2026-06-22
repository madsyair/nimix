# Discrete count components.

test_that("Poisson prior is data-scaled; Binomial needs size", {
  set.seed(1); y <- rpois(50, 4)
  pr <- defaultPrior(PoissonSpec(), y)
  expect_equal(pr$a0 / pr$b0, mean(y), tolerance = 1e-8)   # E[lambda] = mean(y)
  expect_silent(validateParams(PoissonSpec(), pr))
  expect_error(defaultPrior(BinomialSpec(), rbinom(20, 10, .3)), "size")
  prb <- defaultPrior(BinomialSpec(), rbinom(50, 10, .3), control = list(size = 10))
  expect_equal(prb$size, 10L)
})

test_that("Poisson and Binomial recover two components", {
  skip_on_cran(); skip_if_not_installed("nimble")
  set.seed(2); y <- c(rpois(80, 2), rpois(80, 12))
  fp <- nimixClust(y, K = 2, distribution = "poisson", method = "fixedk",
                   mcmcControl = list(niter = 1200, nburnin = 500),
                   seed = 2, verbose = FALSE)
  fp <- relabel(fp)
  lam <- sort(fp@relabeled$summary$lambda_mean)
  expect_true(lam[1] < 5 && lam[2] > 8)

  set.seed(3); yb <- c(rbinom(80, 20, 0.2), rbinom(80, 20, 0.75))
  fb <- nimixClust(yb, K = 2, distribution = "binomial", method = "fixedk",
                   prior = list(size = 20),
                   mcmcControl = list(niter = 1200, nburnin = 500),
                   seed = 3, verbose = FALSE)
  fb <- relabel(fb)
  pp <- sort(fb@relabeled$summary$prob_mean)
  expect_true(pp[1] < 0.4 && pp[2] > 0.6)
})
