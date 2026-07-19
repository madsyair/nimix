# brms-style prediction for mixtures of regressions. The three questions are
# genuinely different here, and the middle one is a trap.

.ppFit <- function(seed = 1L, n = 100L) {
  set.seed(seed)
  x <- rnorm(n); z <- rep(1:2, each = n / 2)
  y <- c(2, -2)[z] + c(1.5, -1.5)[z] * x + rnorm(n, 0, 0.5)
  list(fit = nimixReg(y ~ x, data.frame(y = y, x = x), K = 2,
                      method = "fixedk",
                      mcmcControl = list(niter = 600, nburnin = 250),
                      seed = 1),
       y = y, x = x)
}

test_that("posteriorLinpred returns each component's own line", {
  skip_on_cran()
  o <- .ppFit()
  nd <- data.frame(x = c(-1, 0, 1))
  lp <- posteriorLinpred(o$fit, nd, draws = 300)
  expect_identical(dim(lp), c(300L, 3L, 2L))     # draws x n x K
  m <- apply(lp, c(2, 3), mean)
  # the two components are the two crossing lines: +-(2 + 1.5x)
  hi <- m[, which.max(m[3, ])]; lo <- m[, which.min(m[3, ])]
  expect_lt(max(abs(hi - c(0.5, 2, 3.5))), 0.4)
  expect_lt(max(abs(lo - c(-0.5, -2, -3.5))), 0.4)
})

test_that("posteriorEpred averages the components away, which is the point of the warning", {
  skip_on_cran()
  o <- .ppFit()
  # With newdata the weights are the mixture weights, so E[Y|x] is the flat
  # line through the middle of two crossing ones -- measured 0.028/0.069/0.111
  # for data whose components have slopes +1.5 and -1.5. Correct, and useless.
  ep <- posteriorEpred(o$fit, data.frame(x = c(-1, 0, 1)), draws = 300)
  expect_identical(dim(ep), c(300L, 3L))
  expect_lt(max(abs(colMeans(ep))), 0.5)        # flat, near zero
  # In-sample the weights are the posterior ALLOCATIONS instead, and the same
  # function becomes useful: measured correlation 0.983 with y.
  ep2 <- posteriorEpred(o$fit, draws = 300)
  expect_gt(cor(colMeans(ep2), o$y), 0.9)
})

test_that("posteriorPredictive keeps the bimodality the expectation destroys", {
  skip_on_cran()
  o <- .ppFit()
  pp <- posteriorPredictive(o$fit, data.frame(x = c(-1, 0, 1)), draws = 300)
  expect_identical(dim(pp), c(300L, 3L))
  # at x = 1 the components sit at +3.5 and -3.5: the predictive must straddle
  # both, so its spread dwarfs the within-component sd of 0.5
  expect_gt(sd(pp[, 3]), 2)
  expect_gt(mean(pp[, 3] > 1), 0.2)
  expect_gt(mean(pp[, 3] < -1), 0.2)
})

test_that("posterior prediction guards its contract", {
  skip_on_cran()
  o <- .ppFit()
  expect_error(posteriorLinpred(o$fit, draws = 0), "integer >= 1")
  # an EXTRA column is fine -- model.matrix uses only what the formula names,
  # which is what brms does too. A MISSING one is the real error.
  expect_silent(posteriorLinpred(o$fit, newdata = data.frame(x = 1, w = 2),
                                 draws = 10L))
  # A MISSING predictor is the dangerous case: R resolves the name against the
  # formula's environment -- the fitted data -- so without a guard the caller
  # gets a silent full-length answer to a question they did not ask. Measured
  # before the fix: newdata of one row returned a 10 x 60 x 2 array.
  expect_error(posteriorLinpred(o$fit, newdata = data.frame(w = 2)),
               "missing predictor")
  expect_error(posteriorEpred(o$fit, newdata = data.frame(w = 2)),
               "missing predictor")
  expect_error(predict(o$fit, newdata = data.frame(w = 2)),
               "missing predictor")

  # clustering fits have no linear predictor to report
  set.seed(1)
  fk <- nimixClust(rnorm(60), K = 2, method = "fixedk",
                   mcmcControl = list(niter = 200, nburnin = 50), seed = 1)
  expect_error(posteriorLinpred(fk), "needs a regression fit")
})

test_that("epred with newdata is refused for hmm, where a weight is a function of time", {
  skip_on_cran()
  set.seed(99)
  Tn <- 120L
  P <- rbind(c(.95, .05), c(.07, .93))
  z <- integer(Tn); z[1] <- 1L
  for (t in 2:Tn) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  x <- rnorm(Tn); B <- rbind(c(2, 1.5), c(-2, -1.5))
  y <- B[z, 1] + B[z, 2] * x + rnorm(Tn, 0, 0.6)
  f <- nimixReg(y ~ x, data.frame(y = y, x = x), K = 2, method = "hmm",
                mcmcControl = list(niter = 300, nburnin = 150), seed = 1)

  # A future row has no decoded regime, and its regime probability depends on
  # how far ahead it is -- which is nimixForecast()'s job, not epred's.
  expect_error(posteriorEpred(f, data.frame(x = 1)), "nimixForecast")
  # but linpred needs no weights at all, so it works either way
  lp <- posteriorLinpred(f, data.frame(x = c(0, 1)), draws = 100)
  expect_identical(dim(lp), c(100L, 2L, 2L))
  # and in-sample epred works: the decoded regimes supply the weights
  ep <- posteriorEpred(f, draws = 100)
  expect_identical(dim(ep), c(100L, Tn))
  expect_gt(cor(colMeans(ep), y), 0.8)
})

test_that("linpred and epred respect the link, matching brms", {
  skip_on_cran()
  # The distinction brms draws: posterior_linpred is the linear predictor
  # (log-mu for a Poisson), posterior_epred is E[Y|X] on the response scale
  # (mu = exp(eta)). Getting these the same would be the classic GLM error.
  set.seed(1)
  n <- 150L; x <- rnorm(n)
  y <- rpois(n, exp(0.5 + 0.8 * x))
  f <- nimixReg(y ~ x, data.frame(y = y, x = x), K = 1,
                distribution = "poisson", method = "fixedk",
                mcmcControl = list(niter = 600, nburnin = 250), seed = 1)
  nd <- data.frame(x = c(-1, 0, 1))

  eta <- apply(posteriorLinpred(f, nd, draws = 200), c(2, 3), mean)[, 1]
  mu  <- colMeans(posteriorEpred(f, nd, draws = 200))

  # linpred is on the log scale: it should be roughly 0.5 + 0.8x, NOT exp of it
  expect_lt(max(abs(eta - (0.5 + 0.8 * c(-1, 0, 1)))), 0.4)
  # epred is E[Y|X] = exp(eta), the response mean
  expect_lt(max(abs(mu - exp(0.5 + 0.8 * c(-1, 0, 1)))), 0.6)
  # and they are genuinely different (the whole point)
  expect_gt(max(abs(mu - eta)), 1)

  # transform = TRUE turns linpred into the component response mean
  lpt <- apply(posteriorLinpred(f, nd, transform = TRUE, draws = 200),
               c(2, 3), mean)[, 1]
  expect_lt(max(abs(lpt - mu)), 0.3)
})

test_that("the predictive draws from the family, not a Gaussian around a transformed mean", {
  skip_on_cran()
  set.seed(1)
  n <- 150L; x <- rnorm(n)
  y <- rpois(n, exp(0.5 + 0.8 * x))
  f <- nimixReg(y ~ x, data.frame(y = y, x = x), K = 1,
                distribution = "poisson", method = "fixedk",
                mcmcControl = list(niter = 600, nburnin = 250), seed = 1)
  pp <- posteriorPredictive(f, data.frame(x = c(-1, 0, 1)), draws = 200)
  # a Poisson predictive draw is a count: non-negative integers
  expect_true(all(pp == round(pp) & pp >= 0))
  # and centred on E[Y|X] = exp(eta), not on eta
  expect_gt(mean(pp[, 3]), 2)
})

test_that("for a Normal fit link-awareness is a no-op (identity link)", {
  skip_on_cran()
  set.seed(2)
  n <- 120L; x <- rnorm(n); z <- rep(1:2, each = n / 2)
  y <- c(2, -2)[z] + c(1.5, -1.5)[z] * x + rnorm(n, 0, 0.5)
  g <- nimixReg(y ~ x, data.frame(y = y, x = x), K = 2, method = "fixedk",
                mcmcControl = list(niter = 500, nburnin = 200), seed = 1)
  nd <- data.frame(x = c(-1, 0, 1))
  lp  <- posteriorLinpred(g, nd, draws = 200)
  lpt <- posteriorLinpred(g, nd, transform = TRUE, draws = 200)
  expect_lt(max(abs(lp - lpt)), 1e-8)         # identity link
  expect_gt(cor(colMeans(posteriorEpred(g, draws = 200)), y), 0.9)
})

test_that("the predictive keeps each family's tails, not just the Gaussian one", {
  skip_on_cran()
  # epred is on the response scale (identity link for these, so = Xbeta), but
  # the PREDICTIVE must draw from the actual error law. For Student-t and
  # Normal-Gamma that is a scaled t, and a first cut drew Gaussian noise --
  # measured predictive kurtosis 3.14, indistinguishable from Normal, where a
  # t with df = 3 should be far heavier.
  set.seed(1)
  n <- 200L; x <- rnorm(n)
  y <- 1 + 2 * x + rt(n, df = 3) * 0.8
  kurt <- function(v) mean((v - mean(v))^4) / var(v)^2

  for (dist in c("studentt", "normalgamma")) {
    f <- nimixReg(y ~ x, data.frame(y = y, x = x), K = 1,
                  distribution = dist, method = "fixedk",
                  mcmcControl = list(niter = 500, nburnin = 200), seed = 1)
    # E[Y|X] = mu = Xbeta: identity link, so epred is unchanged
    expect_lt(abs(mean(posteriorEpred(f, data.frame(x = 0), draws = 300)) - 1),
              0.4)
    # but the predictive must be heavy-tailed: kurtosis well above 3
    pp <- posteriorPredictive(f, data.frame(x = 0), draws = 3000)
    expect_gt(kurt(pp[, 1]), 4)
  }
})

test_that("epred and the predictive agree in the mean for every univariate family", {
  skip_on_cran()
  # The predictive mean should track epred, whatever the family: a sanity
  # bridge between the two. Poisson (log link) is the interesting one -- if
  # epred forgot the link this would fail by a factor of exp().
  set.seed(3)
  n <- 200L; x <- rnorm(n)
  cases <- list(
    normal   = 1 + 0.5 * x + rnorm(n, 0, 0.5),
    studentt = 1 + 0.5 * x + rt(n, 5) * 0.5,
    poisson  = rpois(n, exp(0.3 + 0.4 * x))
  )
  for (dist in names(cases)) {
    f <- nimixReg(cases[[dist]] ~ x,
                  data.frame(y = cases[[dist]], x = x), K = 1,
                  distribution = dist, method = "fixedk",
                  mcmcControl = list(niter = 500, nburnin = 200), seed = 1)
    nd <- data.frame(x = c(-0.5, 0.5))
    ep <- colMeans(posteriorEpred(f, nd, draws = 400))
    pm <- colMeans(posteriorPredictive(f, nd, draws = 2000))
    expect_lt(max(abs(ep - pm) / pmax(abs(ep), 1)), 0.25)
  }
})
