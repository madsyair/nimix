# Markov-switching regression (Hamilton 1989): the coefficients and error
# variance switch with a latent first-order Markov regime.
#
# Two things about this path are worth stating up front, because both are
# deliberate and both look like omissions:
#
#  * There is NO conjugate NIG sampler here. customizeSamplers for
#    NormalRegSpec returns early unless a "z" node exists, and the HMM
#    marginalises the allocations away, so NIMBLE's defaults are used. That is
#    required, not lazy: the conjugate update conditions on the allocations,
#    and with them integrated out the coefficient full conditional is no
#    longer Normal-Inverse-Gamma. Installing it anyway is precisely the
#    silent-wrong-sampler bug fixed for Student-t regression.
#  * StudentTRegSpec and NormalGammaRegSpec `contains` NormalRegSpec, so they
#    INHERIT the HMM model code added here. They must be refused. The guard
#    matches on class(spec)[1] rather than is(), which is what makes that
#    work -- asserted below.

test_that("markov-switching regression recovers regime-specific coefficients", {
  skip_on_cran()
  set.seed(99)
  Tn <- 150L
  P <- rbind(c(.95, .05), c(.07, .93))
  z <- integer(Tn); z[1] <- 1L
  for (t in 2:Tn) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  x <- rnorm(Tn)
  B <- rbind(c(2, 1.5), c(-2, -1.5))     # regime 1 rises, regime 2 falls
  y <- B[z, 1] + B[z, 2] * x + rnorm(Tn, 0, 0.6)

  f <- nimixReg(y ~ x, data.frame(y = y, x = x), K = 2, method = "hmm",
                mcmcControl = list(niter = 400, nburnin = 200), seed = 1)
  r <- relabel(f)@relabeled$summary
  o <- order(r[["(Intercept)"]])
  expect_lt(max(abs(r[["(Intercept)"]][o] - c(-2, 2))), 0.4)
  expect_lt(max(abs(r[["x"]][o] - c(-1.5, 1.5))), 0.4)

  zv <- viterbiPath(f)
  expect_gt(max(mean(zv == z), mean((3L - zv) == z)), 0.9)
  expect_identical(length(zv), Tn)
  expect_true(all(zv %in% 1:2))          # K=2 regimes, not the DPM truncation
})

test_that("the hmm regression path refuses distributions without a kernel yet", {
  # The guard now admits normal, studentt, normalgamma, poisson, binomial --
  # each with its OWN kernel. It still refuses multivariate-response specs,
  # which have no HMM kernel; the exact-class check keeps them from silently
  # dispatching to an inherited univariate one.
  set.seed(1)
  df <- data.frame(y1 = rnorm(60), y2 = rnorm(60), x = rnorm(60),
                   g = rep(1:3, 20))
  expect_error(
    nimixReg(cbind(y1, y2) ~ x, df, K = 2, method = "hmm",
             distribution = "normal"),
    "currently supports|multivariate|hmm")
  # and the ordinary contract guards
  df <- data.frame(y = rnorm(60), x = rnorm(60), g = rep(1:3, 20))
  expect_error(nimixReg(y ~ x, df, K = 1, method = "hmm"), "K >= 2")
  expect_error(nimixReg(y ~ x, df, method = "hmm"), "regimes K")
  expect_error(nimixReg(y ~ x, df, K = 2, K_max = 5, method = "hmm"),
               "Use K \\(not K_max\\)")
  expect_error(nimixReg(y ~ x, df, K = 2, method = "hmm", random = ~ g),
               "fixedk")
})

test_that("the hmm regression path uses NIMBLE's default samplers, not the NIG conjugate one", {
  # Asserting the mechanism, not just the outcome: with the allocations
  # marginalised there is no "z" node, so customizeSamplers must decline to
  # install anything. If a future change installs the conjugate sampler here
  # the chain would target the wrong stationary distribution silently -- this
  # test is the tripwire.
  set.seed(1)
  Tn <- 60L
  x <- rnorm(Tn); y <- rnorm(Tn)
  spec <- getDistribution("normal-reg")
  mc <- nimix:::buildModelCode(spec, nimix:::HMMEngine(transConc = 1),
                               n = Tn, L = 2L)
  expect_false("z" %in% all.vars(mc$code))
  expect_identical(mc$allocNode, "zFFBS")     # decoded post-hoc, not sampled
  expect_true("P" %in% all.vars(mc$code))     # the transition matrix is
})

test_that("markov-switching POISSON regression recovers log-rate coefficients", {
  skip_on_cran()
  # The count counterpart of the Gaussian MS regression: a log-link Poisson
  # whose coefficients switch with the regime. Own kernel (dpois(exp(Xbeta)))
  # rather than the Gaussian one, and no error-variance parameter.
  set.seed(7)
  Tn <- 300L
  P <- rbind(c(.95, .05), c(.07, .93))
  z <- integer(Tn); z[1] <- 1L
  for (t in 2:Tn) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  x <- rnorm(Tn)
  B <- rbind(c(1.0, 0.8), c(2.2, -0.5))
  y <- rpois(Tn, exp(B[z, 1] + B[z, 2] * x))

  f <- nimixReg(y ~ x, data.frame(y = y, x = x), K = 2, method = "hmm",
                distribution = "poisson",
                mcmcControl = list(niter = 1500, nburnin = 600), seed = 1)
  r <- relabel(f)@relabeled$summary
  o <- order(r[["(Intercept)"]])
  expect_lt(max(abs(r[["(Intercept)"]][o] - c(1.0, 2.2))), 0.4)
  expect_lt(max(abs(r[["x"]][o] - c(0.8, -0.5))), 0.4)

  zv <- viterbiPath(f)
  expect_gt(max(mean(zv == z), mean((3L - zv) == z)), 0.9)
  expect_true(all(zv %in% 1:2))
})

test_that("forecasting a markov-switching Poisson regression yields counts", {
  skip_on_cran()
  set.seed(7)
  Tn <- 220L; H <- 6L
  P <- rbind(c(.95, .05), c(.07, .93))
  z <- integer(Tn + H); z[1] <- 1L
  for (t in 2:(Tn + H)) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  x <- rnorm(Tn + H); B <- rbind(c(1.0, 0.8), c(2.2, -0.5))
  y <- rpois(Tn + H, exp(B[z, 1] + B[z, 2] * x))
  df <- data.frame(y = y[1:Tn], x = x[1:Tn])
  nd <- data.frame(x = x[(Tn + 1):(Tn + H)])

  f <- nimixReg(y ~ x, df, K = 2, method = "hmm", distribution = "poisson",
                mcmcControl = list(niter = 800, nburnin = 300), seed = 1)
  fc <- nimixForecast(f, h = H, newdata = nd, draws = 150)
  expect_identical(nrow(fc$summary), H)
  expect_true(all(fc$draws == round(fc$draws) & fc$draws >= 0))  # counts
  expect_true(all(abs(rowSums(fc$regime) - 1) < 1e-8))
})

test_that("markov-switching BINOMIAL regression recovers logit coefficients", {
  skip_on_cran()
  # Proportion counterpart: a logit-link Binomial whose coefficients switch
  # with the regime, known number of trials via prior$size. Own kernel
  # (dbinom(size, plogis(Xbeta))).
  set.seed(11)
  Tn <- 300L; size <- 30
  P <- rbind(c(.95, .05), c(.07, .93))
  z <- integer(Tn); z[1] <- 1L
  for (t in 2:Tn) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  x <- rnorm(Tn)
  B <- rbind(c(-1.0, 1.2), c(1.5, -0.8))
  y <- rbinom(Tn, size, plogis(B[z, 1] + B[z, 2] * x))

  f <- nimixReg(y ~ x, data.frame(y = y, x = x), K = 2, method = "hmm",
                distribution = "binomial", prior = list(size = size),
                mcmcControl = list(niter = 1500, nburnin = 600), seed = 1)
  r <- relabel(f)@relabeled$summary
  o <- order(r[["(Intercept)"]])
  expect_lt(max(abs(r[["(Intercept)"]][o] - c(-1.0, 1.5))), 0.4)
  expect_lt(max(abs(r[["x"]][o] - c(1.2, -0.8))), 0.4)

  zv <- viterbiPath(f)
  expect_gt(max(mean(zv == z), mean((3L - zv) == z)), 0.9)
  expect_true(all(zv %in% 1:2))

  # forecast draws are valid proportions of `size`
  fc <- nimixForecast(f, h = 5L, newdata = data.frame(x = rnorm(5)),
                      draws = 150)
  expect_true(all(fc$draws >= 0 & fc$draws <= size &
                  fc$draws == round(fc$draws)))
})

test_that("markov-switching STUDENT-T regression recovers coefficients under heavy tails", {
  skip_on_cran()
  # Heavy-tail counterpart of the Gaussian MS regression: location and scale
  # switch with the regime, errors are Student-t (df fixed). Own kernel
  # (dt_nonstandard), NOT the inherited Gaussian one -- the guard refuses that
  # inheritance (9.29/9.41), and this OWN method runs instead.
  set.seed(5)
  Tn <- 300L
  P <- rbind(c(.95, .05), c(.07, .93))
  z <- integer(Tn); z[1] <- 1L
  for (t in 2:Tn) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  x <- rnorm(Tn)
  B <- rbind(c(1.5, 1.2), c(-1.5, -0.9))
  y <- B[z, 1] + B[z, 2] * x + rt(Tn, df = 4) * 0.6

  f <- nimixReg(y ~ x, data.frame(y = y, x = x), K = 2, method = "hmm",
                distribution = "studentt",
                mcmcControl = list(niter = 1500, nburnin = 600), seed = 1)
  r <- relabel(f)@relabeled$summary
  o <- order(r[["(Intercept)"]])
  expect_lt(max(abs(r[["(Intercept)"]][o] - c(-1.5, 1.5))), 0.4)
  expect_lt(max(abs(r[["x"]][o] - c(-0.9, 1.2))), 0.4)
  expect_gt(max(mean(viterbiPath(f) == z), mean((3L - viterbiPath(f)) == z)), 0.9)
})

test_that("normal-gamma HMM regression shares the Student-t kernel (sibling, explicit method)", {
  skip_on_cran()
  # NormalGammaRegSpec is a sibling of StudentTRegSpec, not a subclass, so it
  # needs its OWN buildModelCode method pointing at the same kernel -- and it
  # must recover the same regimes.
  set.seed(5)
  Tn <- 250L
  P <- rbind(c(.95, .05), c(.07, .93))
  z <- integer(Tn); z[1] <- 1L
  for (t in 2:Tn) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  x <- rnorm(Tn)
  B <- rbind(c(1.5, 1.2), c(-1.5, -0.9))
  y <- B[z, 1] + B[z, 2] * x + rt(Tn, 4) * 0.6

  f <- nimixReg(y ~ x, data.frame(y = y, x = x), K = 2, method = "hmm",
                distribution = "normalgamma",
                mcmcControl = list(niter = 1200, nburnin = 500), seed = 1)
  r <- relabel(f)@relabeled$summary
  o <- order(r[["(Intercept)"]])
  expect_lt(max(abs(r[["(Intercept)"]][o] - c(-1.5, 1.5))), 0.5)
  expect_gt(max(mean(viterbiPath(f) == z), mean((3L - viterbiPath(f)) == z)), 0.9)

  # forecast draws are finite (heavy-tailed but not degenerate)
  fc <- nimixForecast(f, h = 5L, newdata = data.frame(x = rnorm(5)), draws = 150)
  expect_true(all(is.finite(fc$draws)))
})

test_that("neo-normal families run under the HMM engine (generated kernels)", {
  skip_on_cran()
  # Markov-switching skewed/heavy-tail regression. The forward kernel, its
  # BUGSdist registration, and the three S4 methods are all generated from the
  # family's kernel base + shape list, so every family gains an HMM variant.
  mkY <- function(rng, seed = 11) {
    set.seed(seed); Tn <- 260L; x <- rnorm(Tn)
    z <- integer(Tn); z[1] <- 1L
    for (t in 2:Tn) z[t] <- if (runif(1) < 0.95) z[t - 1] else 3L - z[t - 1]
    b <- rbind(c(1.5, 1.2), c(-1.5, -0.8))
    list(y = vapply(1:Tn, function(t) b[z[t], 1] + b[z[t], 2] * x[t], 0) + rng(Tn),
         x = x)
  }
  cases <- list(
    msnburr  = list(rng = function(n) nimix:::rmsnburr(n, 0, 0.6, 2),
                    cols = c("sigma_mean", "alpha_mean")),
    sep      = list(rng = function(n) nimix:::rsep(n, 0, 0.6, 2),
                    cols = c("sigma_mean", "nu_mean")),
    gmsnburr = list(rng = function(n) nimix:::rgmsnburr(n, 0, 0.6, 2, 1.5),
                    cols = c("sigma_mean", "alpha_mean", "theta_mean"))
  )
  for (nm in names(cases)) {
    cs <- cases[[nm]]; d <- mkY(cs$rng)
    f <- nimixReg(d$y ~ d$x, data.frame(y = d$y, x = d$x), K = 2,
                  method = "hmm", distribution = nm,
                  mcmcControl = list(niter = 600, nburnin = 250), seed = 1)
    s <- relabel(f)@relabeled$summary
    expect_true(all(cs$cols %in% names(s)), info = nm)
    o <- order(s[["(Intercept)"]])
    expect_lt(max(abs(s[["(Intercept)"]][o] - c(-1.5, 1.5))), 0.6)
  }
})
