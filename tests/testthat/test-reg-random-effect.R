# Random intercept in the mixture-of-regressions (F5 increment 1).
# Design requirements from the F4 gate, asserted here:
#  * sum-to-zero parameterisation (free b had cor(beta0, mean(b)) = -0.979
#    and min ESS 25/2500; the constraint restored 205 with recovery intact)
#  * the P1 NIG sampler extended with the RE offset (production inprod form
#    gets no conjugacy detection) -- the scale-equivariance lock must still
#    hold with RE active.

.simRE <- function(seed = 31L) {
  set.seed(seed)
  G <- 12L; npg <- 20L; n <- G * npg
  g <- rep(seq_len(G), each = npg)
  b <- rnorm(G, 0, 0.8)
  x <- runif(n, -2, 2)
  zc <- rep(1:2, length.out = n)
  beta <- rbind(c(1, 2), c(-1, -2))
  y <- beta[zc, 1] + beta[zc, 2] * x + b[g] + rnorm(n, 0, 0.5)
  list(df = data.frame(y = y, x = x, g = factor(g)), b = b, G = G)
}

test_that("random intercept recovers b, tau, and component coefficients", {
  skip_on_cran()
  d <- .simRE()
  f <- nimixReg(y ~ x, d$df, random = ~ g, K = 2, method = "fixedk",
                mcmcControl = list(niter = 3500, nburnin = 1400), seed = 1)
  fr <- relabel(f)
  s <- fr@relabeled$summary
  expect_lt(max(abs(sort(s$x) - c(-2, 2))), 0.4)
  S <- f@mcmcSamples
  bh <- colMeans(S[, paste0("b[", seq_len(d$G), "]")])
  expect_gt(cor(bh, d$b - mean(d$b)), 0.9)
  # sum-to-zero enforced exactly, draw by draw
  expect_lt(max(abs(rowSums(S[, paste0("b[", seq_len(d$G), "]")]))), 1e-10)
  expect_lt(abs(mean(S[, "tauRE"]) - 0.8), 0.4)
})

test_that("scale-equivariance (the P1 lock) still holds with RE active", {
  skip_on_cran()
  d <- .simRE()
  f1 <- nimixReg(y ~ x, d$df, random = ~ g, K = 2, method = "fixedk",
                 mcmcControl = list(niter = 3000, nburnin = 1200), seed = 1)
  d2 <- d$df; d2$x <- d2$x * 1000
  f2 <- nimixReg(y ~ x, d2, random = ~ g, K = 2, method = "fixedk",
                 mcmcControl = list(niter = 3000, nburnin = 1200), seed = 1)
  s1 <- relabel(f1)@relabeled$summary
  s2 <- relabel(f2)@relabeled$summary
  expect_lt(max(abs(sort(s1$x) - sort(s2$x * 1000))), 0.15)
})

test_that("mixing meets the gate bar (sum-to-zero, not the free ridge)", {
  skip_on_cran()
  d <- .simRE()
  f <- nimixReg(y ~ x, d$df, random = ~ g, K = 2, method = "fixedk",
                mcmcControl = list(niter = 3500, nburnin = 1400), seed = 1)
  S <- f@mcmcSamples
  keep <- c("betaTilde[1, 1]", "betaTilde[1, 2]", "b[1]", "tauRE")
  ess <- coda::effectiveSize(coda::as.mcmc(S[, keep]))
  # free parameterisation measured 25; sum-to-zero 205-238 -- demand > 100
  expect_gt(min(ess), 100)
})

test_that("plain fits are untouched and guards fire", {
  skip_on_cran()
  d <- .simRE()
  # without random=: model has no b/tauRE nodes at all
  f0 <- nimixReg(y ~ x, d$df, K = 2, method = "fixedk",
                 mcmcControl = list(niter = 600, nburnin = 250), seed = 1)
  expect_false(any(c("b[1]", "tauRE") %in% colnames(f0@mcmcSamples)))

  df <- d$df
  expect_error(nimixReg(y ~ x, df, random = ~ g, K_max = 5, method = "dpm"),
               "fixedk")
  expect_error(nimixReg(y ~ x, df, random = ~ g + x, K = 2,
                        method = "fixedk"), "exactly one")
  expect_error(nimixReg(y ~ x, df, random = ~ zz, K = 2, method = "fixedk"),
               "not found")
  df$g2 <- factor(rep(1:2, 120))
  expect_error(nimixReg(y ~ x, df, random = ~ g2, K = 2, method = "fixedk"),
               "at least 3")
  expect_error(nimixReg(y ~ x, df, random = ~ g, K = 2, method = "fixedk",
                        distribution = "studentt"), "normal")
})
