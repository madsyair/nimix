# GLM mixture-of-regression: link functions, priors, recovery.

test_that("GLM reg specs are regression specs with correct inverse links", {
  expect_true(isRegressionSpec(PoissonRegSpec()))
  expect_true(isRegressionSpec(BinomialRegSpec()))
  expect_false(isRegressionSpec(NormalUvSpec()))
  expect_equal(linkInv(PoissonRegSpec(), log(3)), 3, tolerance = 1e-8)   # log link
  expect_equal(linkInv(BinomialRegSpec(), 0, prior = list(size = 10)), 5,
               tolerance = 1e-8)                                          # logit link
  expect_error(defaultPrior(BinomialRegSpec(), rbinom(20, 10, .3),
                            control = list(X = cbind(1, rnorm(20)))), "size")
})

test_that("nimixReg routes GLM distributions and rejects others", {
  d <- data.frame(x = rnorm(20), y = rpois(20, 3))
  # studentt / normalgamma residual families are available (v0.4.1); an unknown
  # family is the one that must error.
  expect_error(nimixReg(y ~ x, d, K_max = 4, distribution = "gamma"),
               "not available")
})

test_that("Poisson and Binomial regression mixtures recover slopes", {
  skip_on_cran(); skip_if_not_installed("nimble")
  set.seed(1); x1 <- rnorm(100); x2 <- rnorm(100)
  dP <- data.frame(x = c(x1, x2),
                   y = c(rpois(100, exp(0.6 + x1)), rpois(100, exp(0.6 - x2))))
  fp <- nimixReg(y ~ x, dP, K = 2, distribution = "poisson", method = "fixedk",
                 mcmcControl = list(niter = 3000, nburnin = 1200),
                 seed = 1, verbose = FALSE)
  fp <- relabel(fp); sl <- sort(fp@relabeled$summary$x)
  expect_true(sl[1] < -0.5 && sl[2] > 0.5)

  set.seed(2); xb <- rnorm(200)
  yb <- rbinom(200, 20, plogis(c(2 * xb[1:100], -2 * xb[101:200])))
  fb <- nimixReg(y ~ x, data.frame(x = xb, y = yb), K = 2,
                 distribution = "binomial", method = "fixedk",
                 prior = list(size = 20),
                 mcmcControl = list(niter = 3000, nburnin = 1200),
                 seed = 2, verbose = FALSE)
  fb <- relabel(fb); slb <- sort(fb@relabeled$summary$x)
  expect_true(slb[1] < -1 && slb[2] > 1)
})
