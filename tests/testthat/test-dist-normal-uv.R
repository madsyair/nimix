test_that("defaultPrior is data-scaled (Section 9.2)", {
  spec <- NormalUvSpec()
  set.seed(1)
  y <- rnorm(500, mean = 10, sd = 3)
  pr <- defaultPrior(spec, y)
  expect_equal(pr$mu0, mean(y), tolerance = 1e-8)
  # prior mean of s2 = s0 / (nu0 - 1) should be near var(y)
  expect_equal(pr$s0 / (pr$nu0 - 1), var(y), tolerance = 1e-6)
  # kappa0 = 1 / cLoc^2
  expect_equal(pr$kappa0, 1 / (pr$cLoc^2), tolerance = 1e-8)
})

test_that("defaultPrior honours control overrides", {
  spec <- NormalUvSpec()
  pr <- defaultPrior(spec, rnorm(100), control = list(cLoc = 4, nu0 = 5))
  expect_equal(pr$cLoc, 4)
  expect_equal(pr$nu0, 5)
  expect_error(defaultPrior(spec, rnorm(100), control = list(nu0 = 2)))
})

test_that("validateParams rejects degenerate priors", {
  spec <- NormalUvSpec()
  good <- list(mu0 = 0, kappa0 = 0.25, nu0 = 3, s0 = 2)
  expect_true(validateParams(spec, good))
  expect_error(validateParams(spec, list(mu0 = 0, kappa0 = 0.25, nu0 = 3)))
  expect_error(validateParams(spec, modifyList(good, list(nu0 = 2))))
  expect_error(validateParams(spec, modifyList(good, list(kappa0 = -1))))
})

test_that("simulateParams returns finite (mu, s2) of correct length", {
  spec <- NormalUvSpec()
  pr <- defaultPrior(spec, rnorm(200))
  set.seed(7)
  sp <- simulateParams(spec, pr, nClust = 5)
  expect_length(sp$mu, 5)
  expect_length(sp$s2, 5)
  expect_true(all(is.finite(sp$mu)))
  expect_true(all(sp$s2 > 0))
})

test_that("buildModelCode returns code + monitors for NormalUv x DPM", {
  mc <- buildModelCode(NormalUvSpec(), DPMEngine(), n = 50L, L = 8L)
  expect_true(inherits(mc$code, "nimbleCode") ||
              inherits(mc$code, "{") || is.call(mc$code))
  expect_setequal(mc$monitors, c("xi", "muTilde", "s2Tilde", "alpha"))
  expect_identical(unname(mc$paramNodes[["mu"]]), "muTilde")
})

test_that("componentDensity evaluates a Gaussian", {
  f <- componentDensity(NormalUvSpec())
  expect_equal(f(0, list(mu = 0, s2 = 1)), dnorm(0, 0, 1))
})
