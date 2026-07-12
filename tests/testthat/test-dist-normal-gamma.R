# Tests for the Normal-Gamma (scale-mixture Student-t) component. The key check
# is that the implemented marginal equals the analytic Student-t density.

test_that("Normal-Gamma marginal equals the analytic Student-t (numeric)", {
  mu <- 1.5; s2 <- 2.3; df <- 5
  mix <- function(y) stats::integrate(function(w)
    stats::dnorm(y, mu, sqrt(s2 / w)) * stats::dgamma(w, df / 2, df / 2),
    0, Inf)$value
  tdn <- function(y) stats::dt((y - mu) / sqrt(s2), df) / sqrt(s2)
  ys <- c(-3, -1, 0, 1.5, 3, 6)
  expect_equal(vapply(ys, mix, numeric(1)), vapply(ys, tdn, numeric(1)),
               tolerance = 1e-6)
  # and the spec's componentDensity matches the analytic t
  d <- nimix:::componentDensity(NormalGammaUvSpec(), df = df)
  expect_equal(vapply(ys, function(y) d(y, list(mu = mu, s2 = s2, df = df)),
                      numeric(1)),
               vapply(ys, tdn, numeric(1)), tolerance = 1e-10)
})

test_that("Normal-Gamma inherits the NIG prior and adds df", {
  set.seed(1); y <- rt(200, df = 5)
  pr <- defaultPrior(NormalGammaUvSpec(), y)
  expect_true(all(c("mu0", "kappa0", "nu0", "s0", "df") %in% names(pr)))
  expect_error(defaultPrior(NormalGammaUvSpec(), y, control = list(df = 2)),
               "df must exceed 2")
  expect_true(methods::is(NormalGammaUvSpec(), "NormalUvSpec"))   # inheritance
})

test_that("Normal-Gamma DPM recovers two heavy-tailed clusters", {
  skip_on_cran(); skip_if_not_installed("nimble")
  set.seed(6); y <- c(rt(70, df = 4) - 5, rt(70, df = 4) + 5)
  fit <- nimixClust(y, K_max = 8, distribution = "normalgamma", method = "dpm",
                    prior = list(df = 4),
                    mcmcControl = list(niter = 1500, nburnin = 600),
                    seed = 6, verbose = FALSE)
  fit <- relabel(fit)
  mu <- sort(fit@relabeled$summary$mu_mean)
  expect_true(mu[1] < -2 && mu[length(mu)] > 2)
})
