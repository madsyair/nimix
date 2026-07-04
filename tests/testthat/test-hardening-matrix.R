# Hardening matrix coverage (v0.9.0). Representative distribution x engine
# combinations that had no dedicated test anywhere before the 9.x-D hardening
# pass. The FULL matrix (every released family x both engines x 3 seeds) lives
# in inst/harness/run_recovery_suite.R (harness Layer 4); these tests pin the
# previously untested pairings so a regression cannot reopen the gap silently.

test_that("Student-t clustering works under the fixed-K engine", {
  skip_on_cran()
  set.seed(102)
  y <- c(-4 + rt(80, df = 5), 4 + rt(80, df = 5))
  fit <- relabel(nimixClust(y, K = 2, distribution = "studentt",
                            method = "fixedk",
                            mcmcControl = list(niter = 2500, nburnin = 1000),
                            seed = 7))
  mu <- sort(fit@relabeled$summary$mu_mean)
  expect_equal(fit@relabeled$modalK, 2L)
  expect_lt(abs(mu[1] - (-4)), 1.2)
  expect_lt(abs(mu[2] - 4), 1.2)
})

test_that("multivariate Normal-Gamma clustering works under the fixed-K engine", {
  skip_on_cran()
  set.seed(104)
  Y <- rbind(-4 + matrix(rt(100, df = 6), ncol = 2),
              4 + matrix(rt(100, df = 6), ncol = 2))
  fit <- relabel(nimixClust(Y, K = 2, distribution = "normalgamma",
                            method = "fixedk",
                            mcmcControl = list(niter = 2500, nburnin = 1000),
                            seed = 7))
  mu1 <- sort(fit@relabeled$summary$mu_1)
  expect_equal(fit@relabeled$modalK, 2L)
  expect_lt(abs(mu1[1] - (-4)), 1.5)
  expect_lt(abs(mu1[2] - 4), 1.5)
})

test_that("Poisson clustering works under the DPM engine (dominant components)", {
  skip_on_cran()
  set.seed(105)
  y <- c(rpois(100, 3), rpois(100, 15))
  fit <- relabel(nimixClust(y, K_max = 8, distribution = "poisson",
                            method = "dpm",
                            mcmcControl = list(niter = 3000, nburnin = 1000),
                            seed = 7))
  # On discrete counts the DPM posterior K is diffuse and small transient
  # components are expected even when recovery is clean (Miller & Harrison
  # 2013, NIPS): assess the DOMINANT components (weight >= 0.1).
  sm  <- fit@relabeled$summary
  dom <- sm[sm$weight >= 0.1, , drop = FALSE]
  lam <- sort(dom$lambda_mean)
  expect_equal(nrow(dom), 2L)
  expect_lt(abs(lam[1] - 3) / 3, 0.35)
  expect_lt(abs(lam[2] - 15) / 15, 0.35)
})

test_that("mixture regression works under the fixed-K engine", {
  skip_on_cran()
  set.seed(107)
  n <- 200; x <- runif(n, -3, 3)
  grp <- rep(1:2, c(160, 40))
  y <- ifelse(grp == 1, 2 * x, -2 * x) + rnorm(n, 0, 0.7)
  fit <- relabel(nimixReg(y ~ x, data.frame(y = y, x = x), K = 2,
                          method = "fixedk",
                          mcmcControl = list(niter = 3000, nburnin = 1000),
                          seed = 7))
  slopes <- sort(fit@relabeled$summary[["x"]])
  expect_lt(slopes[1], -1)
  expect_gt(slopes[2], 1)
})
