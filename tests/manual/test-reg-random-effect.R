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
  # tau estimates the spread of the group offsets. With only G = 12 groups the
  # REALIZED sd of the drawn b's is itself a lottery (5th-95th percentile
  # 0.51-1.07 when the population sigma is 0.8), so the diagnostic check is
  # against the realized spread -- comparing to the nominal 0.8 tests the draw,
  # not the estimator. Both are asserted, the nominal one only loosely.
  expect_lt(abs(mean(S[, "tauRE"]) - sd(d$b - mean(d$b))), 0.25)
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

test_that("unbalanced groups across components do not break the random intercept", {
  skip_on_cran()
  # The F4 gate used groups balanced across components and flagged the
  # unbalanced case as an untested risk: groups that sit almost entirely in
  # one component could let b_g absorb that component's intercept. This makes
  # the risk an executable check -- groups 1-6 are ~90% component 1, groups
  # 7-12 ~90% component 2.
  set.seed(81)
  G <- 12L; npg <- 20L; n <- G * npg
  g <- rep(seq_len(G), each = npg)
  b <- rnorm(G, 0, 0.8)
  x <- runif(n, -2, 2)
  pc1 <- ifelse(g <= 6L, 0.9, 0.1)
  zc <- ifelse(runif(n) < pc1, 1L, 2L)
  beta <- rbind(c(1, 2), c(-1, -2))
  y <- beta[zc, 1] + beta[zc, 2] * x + b[g] + rnorm(n, 0, 0.5)
  df <- data.frame(y = y, x = x, gr = factor(g))

  f <- nimixReg(y ~ x, df, random = ~ gr, K = 2, method = "fixedk",
                mcmcControl = list(niter = 3000, nburnin = 1200), seed = 1)
  s <- relabel(f)@relabeled$summary
  # component coefficients survive: the offsets did NOT eat the intercepts
  expect_lt(max(abs(sort(s$x) - c(-2, 2))), 0.4)
  expect_lt(max(abs(sort(s[["(Intercept)"]]) - c(-1, 1))), 0.4)
  S <- f@mcmcSamples
  bh <- colMeans(S[, paste0("b[", seq_len(G), "]")])
  expect_gt(cor(bh, b - mean(b)), 0.9)
  # tau tracks the REALIZED spread (see the note in the recovery test above)
  expect_lt(abs(mean(S[, "tauRE"]) - sd(b - mean(b))), 0.25)
  # mixing does not collapse under imbalance
  ess <- coda::effectiveSize(coda::as.mcmc(S[, c("betaTilde[1, 1]", "b[1]",
                                                 "tauRE")]))
  expect_gt(min(ess), 100)
})

test_that("random-effect priors are data-scaled in the RESPONSE (regression test)", {
  skip_on_cran()
  # v1.3.0 shipped fixed bounds dunif(0.01, 5) on tauRE. With the response on
  # a large scale the needed tau sat far above the ceiling and the offsets
  # silently collapsed: at y x1000 the needed tauRE was 771 and
  # cor(b_hat, truth) fell from 0.992 to 0.091. The P1 scale-equivariance lock
  # never caught it because it only rescales PREDICTORS. This closes that gap.
  d <- .simRE()
  df <- d$df
  df$y <- df$y * 1000
  f <- nimixReg(y ~ x, df, random = ~ g, K = 2, method = "fixedk",
                mcmcControl = list(niter = 2000, nburnin = 800), seed = 1)
  bh <- colMeans(f@mcmcSamples[, paste0("b[", seq_len(d$G), "]")])
  bt <- (d$b - mean(d$b)) * 1000
  expect_gt(cor(bh, bt), 0.9)                  # was 0.091 before the fix
  expect_lt(abs(mean(f@mcmcSamples[, "tauRE"]) / sd(bt) - 1), 0.3)
})

test_that("random slope recovers group slopes and component coefficients", {
  skip_on_cran()
  set.seed(101)
  G <- 12L; npg <- 25L; n <- G * npg
  g <- rep(seq_len(G), each = npg)
  b <- rnorm(G, 0, 0.6); sg <- rnorm(G, 0, 0.3)
  x <- runif(n, -2, 2)
  zc <- rep(1:2, length.out = n)
  beta <- rbind(c(1, 2), c(-1, -2))
  y <- beta[zc, 1] + (beta[zc, 2] + sg[g]) * x + b[g] + rnorm(n, 0, 0.5)
  df <- data.frame(y = y, x = x, g = factor(g))

  f <- nimixReg(y ~ x, df, random = ~ x | g, K = 2, method = "fixedk",
                mcmcControl = list(niter = 3000, nburnin = 1200), seed = 1)
  S <- f@mcmcSamples
  expect_true(all(c("sRE[1]", "tauSlope", "b[1]", "tauRE") %in% colnames(S)))

  # Component slopes absorb mean(s_g) -- the sum-to-zero semantics, exactly as
  # the sum-to-zero constraint predicted (measured -2.14/1.88 against a prediction of
  # -2.13/1.87). Test the documented estimand, not the raw truth.
  s <- relabel(f)@relabeled$summary
  expect_lt(max(abs(sort(s$x) - (c(-2, 2) + mean(sg)))), 0.4)

  sh <- colMeans(S[, paste0("sRE[", seq_len(G), "]")])
  bh <- colMeans(S[, paste0("b[", seq_len(G), "]")])
  expect_gt(cor(sh, sg - mean(sg)), 0.9)
  expect_gt(cor(bh, b - mean(b)), 0.9)

  # both offsets sum to zero exactly, draw by draw
  expect_lt(max(abs(rowSums(S[, paste0("sRE[", seq_len(G), "]")]))), 1e-10)
  expect_lt(max(abs(rowSums(S[, paste0("b[", seq_len(G), "]")]))), 1e-10)
  # tau tracks the realized spread (not the population sigma)
  expect_lt(abs(mean(S[, "tauSlope"]) - sd(sg - mean(sg))), 0.15)
  ess <- coda::effectiveSize(coda::as.mcmc(S[, c("betaTilde[1, 2]", "sRE[1]",
                                                 "tauSlope")]))
  expect_gt(min(ess), 100)          # gate: free parameterisation gave 52
})

test_that("random slope is scale-equivariant in the predictor", {
  skip_on_cran()
  set.seed(101)
  G <- 12L; npg <- 25L; n <- G * npg
  g <- rep(seq_len(G), each = npg)
  b <- rnorm(G, 0, 0.6); sg <- rnorm(G, 0, 0.3)
  x <- runif(n, -2, 2)
  zc <- rep(1:2, length.out = n)
  beta <- rbind(c(1, 2), c(-1, -2))
  y <- beta[zc, 1] + (beta[zc, 2] + sg[g]) * x + b[g] + rnorm(n, 0, 0.5)
  df <- data.frame(y = y, x = x, g = factor(g))
  df2 <- df; df2$x <- df2$x * 1000

  f2 <- nimixReg(y ~ x, df2, random = ~ x | g, K = 2, method = "fixedk",
                 mcmcControl = list(niter = 2000, nburnin = 800), seed = 1)
  # tauSlope has units y/x, so it must shrink by the same factor -- with the
  # old fixed bounds it pinned at 0.01 instead of reaching 0.00029.
  expect_lt(abs(mean(f2@mcmcSamples[, "tauSlope"]) * 1000 /
                sd(sg - mean(sg)) - 1), 0.3)
  sh <- colMeans(f2@mcmcSamples[, paste0("sRE[", seq_len(G), "]")])
  expect_gt(cor(sh, (sg - mean(sg)) / 1000), 0.9)
})

test_that("random-slope guards fire and simpler paths stay untouched", {
  skip_on_cran()
  d <- .simRE()
  df <- d$df
  # slope variable must be a term of the fixed-effects formula
  expect_error(nimixReg(y ~ 1, df, random = ~ x | g, K = 2, method = "fixedk"),
               "term of the fixed")
  # only one slope variable
  df$w <- rnorm(nrow(df))
  expect_error(nimixReg(y ~ x + w, df, random = ~ x + w | g, K = 2,
                        method = "fixedk"), "exactly one random-slope")
  # constant slope variable carries no information
  df$const <- 1
  expect_error(nimixReg(y ~ const, df, random = ~ const | g, K = 2,
                        method = "fixedk"), "constant|term of the fixed")
  # intercept-only path gains no slope nodes (three-variant sampler safety)
  f1 <- nimixReg(y ~ x, d$df, random = ~ g, K = 2, method = "fixedk",
                 mcmcControl = list(niter = 600, nburnin = 250), seed = 1)
  expect_false(any(c("sRE[1]", "tauSlope") %in% colnames(f1@mcmcSamples)))
})

test_that("random effects work with heavy-tailed (student-t) residuals", {
  skip_on_cran()
  # Extending RE to this family needed no sampler work: unlike normal-reg it
  # keeps NIMBLE's default samplers (the inherited NIG Gibbs step is only
  # valid for Gaussian errors -- see test-reg-sampler-dispatch.R), so the
  # offsets just enter the linear predictor.
  set.seed(131)
  G <- 12L; npg <- 22L; n <- G * npg
  g <- rep(seq_len(G), each = npg)
  b <- rnorm(G, 0, 0.8)
  x <- runif(n, -2, 2)
  zc <- rep(1:2, length.out = n)
  beta <- rbind(c(1, 2), c(-1, -2))
  y <- beta[zc, 1] + beta[zc, 2] * x + b[g] + rt(n, 4) * 0.4
  df <- data.frame(y = y, x = x, gr = factor(g))

  f <- nimixReg(y ~ x, df, random = ~ gr, K = 2, method = "fixedk",
                distribution = "studentt",
                mcmcControl = list(niter = 3000, nburnin = 1200), seed = 1)
  S <- f@mcmcSamples
  expect_true(all(c("b[1]", "tauRE") %in% colnames(S)))
  s <- relabel(f)@relabeled$summary
  expect_lt(max(abs(sort(s$x) - c(-2, 2))), 0.4)
  bh <- colMeans(S[, paste0("b[", seq_len(G), "]")])
  expect_gt(cor(bh, b - mean(b)), 0.9)
  expect_lt(max(abs(rowSums(S[, paste0("b[", seq_len(G), "]")]))), 1e-10)
})

test_that("random slopes work with heavy-tailed residuals", {
  skip_on_cran()
  set.seed(141)
  G <- 12L; npg <- 25L; n <- G * npg
  g <- rep(seq_len(G), each = npg)
  b <- rnorm(G, 0, 0.6); sg <- rnorm(G, 0, 0.3)
  x <- runif(n, -2, 2)
  zc <- rep(1:2, length.out = n)
  beta <- rbind(c(1, 2), c(-1, -2))
  y <- beta[zc, 1] + (beta[zc, 2] + sg[g]) * x + b[g] + rt(n, 4) * 0.4
  df <- data.frame(y = y, x = x, gr = factor(g))

  f <- nimixReg(y ~ x, df, random = ~ x | gr, K = 2, method = "fixedk",
                distribution = "studentt",
                mcmcControl = list(niter = 3000, nburnin = 1200), seed = 1)
  S <- f@mcmcSamples
  expect_true(all(c("sRE[1]", "tauSlope") %in% colnames(S)))
  # same sum-to-zero semantics: component slopes absorb mean(s_g)
  s <- relabel(f)@relabeled$summary
  expect_lt(max(abs(sort(s$x) - (c(-2, 2) + mean(sg)))), 0.4)
  sh <- colMeans(S[, paste0("sRE[", seq_len(G), "]")])
  expect_gt(cor(sh, sg - mean(sg)), 0.9)
})
