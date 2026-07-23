# The scale prior on the component error variance.
#
# defaultPrior centres the InvGamma prior on s2 at the residual variance of a
# GLOBAL OLS fit, which ignores the mixture. For separated components that
# measures the between-component spread as much as the within-component one,
# so s2 is biased upward. This file anchors the documented numbers (the
# caveat-as-executable-test pattern) and locks the override that lets a user
# who knows better say so.
#
# Measured in the study behind these tests: the bias runs ~9x at 25
# observations per component and ~1.0x by 5000; it grows with component
# separation (1.0x at +-0.5 intercepts, 4.7x at +-6); the coefficients and the
# DPM's K recovery are unaffected. The conservatism is deliberate -- it also
# regularises against splitting when K is over-specified -- so the default
# stays and the escape hatch is explicit.

.simSeparated <- function(seed = 5L, n = 300L) {
  set.seed(seed)
  x <- runif(n, -2, 2)
  zc <- rep(1:2, length.out = n)
  y <- c(3, -3)[zc] + c(2, -2)[zc] * x + rnorm(n, 0, 0.5)   # true s2 = 0.25
  data.frame(y = y, x = x)
}

.meanS2 <- function(fit)
  mean(c(mean(fit@mcmcSamples[, "s2Tilde[1]"]),
         mean(fit@mcmcSamples[, "s2Tilde[2]"])))

test_that("the prior scale can be set directly (s2Guess and s0)", {
  y <- rnorm(80)
  X <- stats::model.matrix(~ rnorm(80))
  spec <- getDistribution("normal-reg")

  d <- nimix:::defaultPrior(spec, y, control = list(X = X, s2Guess = 0.3))
  expect_equal(d$s0 / (d$nu0 - 1), 0.3)          # s2Guess IS the prior mean
  d2 <- nimix:::defaultPrior(spec, y, control = list(X = X, s0 = 0.6))
  expect_equal(d2$s0, 0.6)                        # s0 is the raw scale

  # the heavy-tailed family delegates via callNextMethod, so it inherits this
  dt <- nimix:::defaultPrior(getDistribution("student-t-reg"), y,
                             control = list(X = X, s2Guess = 0.3))
  expect_equal(dt$s0 / (dt$nu0 - 1), 0.3)
  expect_equal(dt$df, 4)

  expect_error(nimix:::defaultPrior(spec, y,
                                    control = list(X = X, s2Guess = -1)),
               "positive scalar")
  expect_error(nimix:::defaultPrior(spec, y, control = list(X = X, s0 = 0)),
               "positive scalar")
})

test_that("the default prior is conservative on separated components: the caveat", {
  skip_on_cran()
  # This test documents a deliberate design choice, not a defect. The default
  # over-estimates s2 for well-separated components; supplying the real scale
  # roughly halves the error; the slopes never move either way.
  df <- .simSeparated()
  mc <- list(niter = 3000, nburnin = 1200)

  fd <- nimixReg(y ~ x, df, K = 2, method = "fixedk", mcmcControl = mc,
                 seed = 1, verbose = FALSE)
  fs <- nimixReg(y ~ x, df, K = 2, method = "fixedk",
                 prior = list(s2Guess = 0.3), mcmcControl = mc, seed = 1,
                 verbose = FALSE)

  s2d <- .meanS2(fd)
  s2s <- .meanS2(fs)
  expect_gt(s2d / 0.25, 2)          # the default IS biased up (measured 3.1x)
  expect_lt(s2s, s2d)               # the informed prior pulls it back
  expect_lt(s2s / 0.25, 2)          # measured 1.4x

  # the coefficients are untouched by either choice -- that is why the bias
  # stayed invisible for so long
  for (f in list(fd, fs))
    expect_lt(max(abs(sort(relabel(f)@relabeled$summary$x) - c(-2, 2))), 0.2)
})

test_that("a tighter prior is safe with a correct K and risky when K is over-specified", {
  skip_on_cran()
  # This is why the default did not change. The conservatism regularises
  # against splitting, so tightening the scale is only safe when K is not
  # over-specified. Both halves are asserted; s2Guess is used to tighten,
  # since a relative multiplier on the (wrong-quantity) automatic scale was
  # rejected as a knob rather than a statement.
  df <- .simSeparated()
  mc <- list(niter = 1500, nburnin = 600)
  tight <- list(s2Guess = 1.5)          # ~1/10 of the automatic scale here
  nOcc <- function(f) as.integer(names(which.max(table(
    apply(f@clusterAllocation, 1L, function(v) length(unique(v)))))))

  k2a <- nimixReg(y ~ x, df, K = 2, method = "fixedk", mcmcControl = mc,
                  seed = 1, verbose = FALSE)
  k2b <- nimixReg(y ~ x, df, K = 2, method = "fixedk", prior = tight,
                  mcmcControl = mc, seed = 1, verbose = FALSE)
  expect_lt(.meanS2(k2b), .meanS2(k2a))       # tightening helps s2
  expect_identical(nOcc(k2a), 2L)
  expect_identical(nOcc(k2b), 2L)             # and costs nothing when K is right

  k4a <- nimixReg(y ~ x, df, K = 4, method = "fixedk", mcmcControl = mc,
                  seed = 1, verbose = FALSE)
  k4b <- nimixReg(y ~ x, df, K = 4, method = "fixedk", prior = tight,
                  mcmcControl = mc, seed = 1, verbose = FALSE)
  expect_identical(nOcc(k4a), 2L)             # default resists over-splitting
  expect_gt(nOcc(k4b), nOcc(k4a))             # tightened prior splits: the risk
})

test_that("the multivariate prior ellipse is anisotropic by default, and sigmaGuess fixes its shape", {
  # The multivariate analogue, and a different failure from the univariate
  # one. cov(data) is inflated only ALONG the direction separating the
  # components, so the prior mean of Sigma is the wrong SHAPE, not merely the
  # wrong size: measured 37.6x along the separation vector against 0.9x across
  # it, a condition number of 42.8 where the truth is a circle.
  #
  # The InverseWishart default survives this (df0 = d + 2 makes the prior worth
  # exactly one observation, for every d, so the posterior came back at 1.2x
  # with condition 1.4) -- which is why the default stays. What it cannot
  # express is shape, and sigmaGuess is that route.
  set.seed(9)
  n <- 400L; d <- 2L
  v <- c(1, 1) / sqrt(2)
  zc <- rep(1:2, length.out = n)
  mu <- rbind(3 * v, -3 * v)
  Y <- mu[zc, ] + matrix(rnorm(n * d, 0, 0.5), n, d)   # within-cov = 0.25 * I
  spec <- getDistribution("normal-mv")
  pmean <- function(p) p$S0 / (p$df0 - d - 1)
  cond <- function(M) {
    e <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
    max(e) / min(e)
  }

  a <- pmean(nimix:::defaultPrior(spec, Y))
  expect_gt(cond(a), 10)                       # the default ellipse IS skewed
  expect_gt(as.numeric(t(v) %*% a %*% v), 5)   # inflated along the separation

  # a scalar guess means "isotropic with this variance"
  b <- pmean(nimix:::defaultPrior(spec, Y, control = list(sigmaGuess = 0.25)))
  expect_equal(cond(b), 1)                     # a circle again
  expect_equal(as.numeric(t(v) %*% b %*% v), 0.25)

  # a matrix guess is taken as given
  cm <- pmean(nimix:::defaultPrior(spec, Y,
                                   control = list(sigmaGuess = diag(c(.25, .5)))))
  expect_equal(diag(cm), c(0.25, 0.5))

  expect_error(nimix:::defaultPrior(spec, Y, control = list(sigmaGuess = -1)),
               "positive scalar")
  expect_error(nimix:::defaultPrior(spec, Y,
                                    control = list(sigmaGuess = matrix(0, 2, 2))),
               "positive definite")
  expect_error(nimix:::defaultPrior(spec, Y,
                                    control = list(sigmaGuess = diag(3))),
               "2 x 2 matrix")
})

test_that("the multivariate REGRESSION prior is the univariate scale problem in matrix form", {
  # Two findings meet here. The mechanism is the univariate regression one
  # (9.30): .mvRegPrior centres Sigma on cov() of the GLOBAL OLS residuals,
  # which carry between-component variation as well as within. The severity
  # is the multivariate one (9.31): the InverseWishart default is worth one
  # observation, so the posterior mostly shrugs it off.
  #
  # Measured on two components differing ONLY in their coefficients, with
  # isotropic within-covariance 0.25 I: prior trace 22.4x truth and condition
  # 40 (a circle in truth), yet the fitted covariance came back at ~1.5x with
  # condition ~2. Compare the univariate regression case, ~2.4x at the same
  # per-component n against a much smaller prior error -- the IW weight is
  # what makes the difference.
  set.seed(7)
  nj <- 100L; n <- 2L * nj; d <- 2L
  x <- rnorm(n); X <- cbind(1, x)
  zc <- rep(1:2, each = nj)
  B1 <- rbind(c(-2, -2), c(1.5, 1.5))
  B2 <- rbind(c(2, 2), c(-1.5, -1.5))
  Y <- matrix(0, n, d)
  for (i in seq_len(n)) {
    B <- if (zc[i] == 1L) B1 else B2
    Y[i, ] <- X[i, ] %*% B + rnorm(d, 0, 0.5)
  }
  spec <- getDistribution("normal-mv-reg")
  pmean <- function(p) p$S0 / (p$df0 - d - 1)
  cond <- function(M) {
    e <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
    max(e) / min(e)
  }

  a <- pmean(nimix:::defaultPrior(spec, Y, control = list(X = X)))
  expect_gt(sum(diag(a)) / 0.5, 10)   # wrong size: the global residuals are inflated
  expect_gt(cond(a), 10)              # and wrong shape: skewed along the coefficient gap

  # a scalar guess means "isotropic residuals with this variance"
  b <- pmean(nimix:::defaultPrior(spec, Y,
                                  control = list(X = X, sigmaGuess = 0.25)))
  expect_equal(cond(b), 1)
  expect_equal(diag(b), rep(0.25, d))

  # the t and normal-gamma variants inherit it through callNextMethod()
  for (nm in c("student-t-mv-reg", "normal-gamma-mv-reg")) {
    p <- pmean(nimix:::defaultPrior(getDistribution(nm), Y,
                                    control = list(X = X, sigmaGuess = 0.25)))
    expect_equal(diag(p), rep(0.25, d))
  }

  expect_error(nimix:::defaultPrior(spec, Y,
                                    control = list(X = X, sigmaGuess = -1)),
               "positive scalar")
  expect_error(nimix:::defaultPrior(spec, Y,
                                    control = list(X = X,
                                                   sigmaGuess = matrix(0, 2, 2))),
               "positive definite")
  expect_error(nimix:::defaultPrior(spec, Y,
                                    control = list(X = X,
                                                   sigmaGuess = diag(3))),
               "2 x 2 matrix")
})
