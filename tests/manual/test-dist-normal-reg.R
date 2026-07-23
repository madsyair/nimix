# Tests for the Normal-linear regression component (v0.3.0). Prior/validation/
# density tests are pure R; the end-to-end recovery test (incl. the unbalanced
# 80/20 scenario required by the test plan) is skipped unless
# nimble is installed and not on CRAN.

test_that("defaultPrior builds a data-scaled NIG g-prior", {
  set.seed(1)
  x <- runif(100, -3, 3); X <- cbind(1, x); y <- 2 * x + rnorm(100)
  spec <- NormalRegSpec()
  pr <- defaultPrior(spec, y, control = list(X = X))
  expect_equal(pr$p, 2L)
  expect_length(pr$b0, 2L)
  expect_equal(dim(pr$B0), c(2L, 2L))
  expect_true(pr$nu0 > 2)
  expect_true(pr$s0 > 0)
  # default g-prior factor = n (unit information)
  expect_equal(pr$g, length(y))
  # B0 = g * (X'X)^{-1}
  expect_equal(pr$B0, length(y) * solve(crossprod(X)), tolerance = 1e-6)
  expect_silent(validateParams(spec, pr))
})

test_that("defaultPrior requires the design matrix", {
  expect_error(defaultPrior(NormalRegSpec(), rnorm(10)), "control\\$X")
})

test_that("validateParams enforces dimension and nu0 invariants", {
  spec <- NormalRegSpec()
  good <- list(b0 = c(0, 0), B0 = diag(2), nu0 = 3, s0 = 1, p = 2L)
  expect_silent(validateParams(spec, good))
  bad_nu <- good; bad_nu$nu0 <- 2
  expect_error(validateParams(spec, bad_nu), "nu0")
  bad_b <- good; bad_b$b0 <- c(0, 0, 0)
  expect_error(validateParams(spec, bad_b), "b0")
  bad_B <- good; bad_B$B0 <- matrix(c(1, 2, 2, 1), 2)   # not PD
  expect_error(validateParams(spec, bad_B), "positive definite")
})

test_that("simulateParams returns conformable beta and s2", {
  set.seed(2)
  pr <- list(b0 = c(0, 0), B0 = diag(2), nu0 = 5, s0 = 2, p = 2L)
  sp <- nimix:::simulateParams(NormalRegSpec(), pr, nClust = 4)
  expect_equal(dim(sp$beta), c(4L, 2L))
  expect_length(sp$s2, 4L)
  expect_true(all(sp$s2 > 0))
})

test_that("end-to-end DPM mixture of regressions recovers two slopes (80/20)", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  set.seed(11)
  # Unbalanced 80/20 mixture of two regression regimes (slope +2 vs -2),
  # the scenario required by the test plan
  n <- 200; n1 <- 160; n2 <- n - n1
  x <- runif(n, -3, 3)
  grp <- c(rep(1L, n1), rep(2L, n2))
  slope <- ifelse(grp == 1L, 2, -2)
  y <- slope * x + rnorm(n, 0, 0.6)
  fit <- nimixReg(y ~ x, data.frame(y = y, x = x), K_max = 8,
                  mcmcControl = list(niter = 3000, nburnin = 1000),
                  seed = 11, verbose = FALSE)
  modalK <- as.integer(names(sort(table(fit@Kposterior),
                                  decreasing = TRUE))[1])
  expect_true(modalK %in% c(2L, 3L))
  fit <- relabel(fit)
  # the two recovered slopes (column "x") should bracket 0 (one +, one -)
  slopes <- sort(fit@relabeled$summary[["x"]])
  expect_true(slopes[1] < 0 && slopes[length(slopes)] > 0)
})

test_that("FixedK regression is scale-equivariant in the predictors", {
  # Regression guarantee for the conjugate NIG Gibbs sampler. Before it,
  # betaTilde/s2Tilde fell back to RW_block in raw units (dynamic indexing
  # hides the linear-Gaussian structure from checkConjugacy), and rescaling a
  # predictor by 1000 shifted the recovered slope systematically (2.50 vs 2.01,
  # persisting at 6000 iterations). With the exact conditional the fit must be
  # invariant: slopes agree after undoing the rescale, up to MC noise.
  skip_on_cran()
  set.seed(3); n <- 200
  x <- runif(n, 0, 10); zc <- rbinom(n, 1, 0.5)
  y <- ifelse(zc == 1, 2 + 2 * x, -2 - 2 * x) + rnorm(n, 0, 1)
  run <- function(xx) {
    f <- relabel(nimixReg(y ~ x, data = data.frame(y = y, x = xx), K = 2,
                          method = "fixedk",
                          mcmcControl = list(niter = 1500, nburnin = 600),
                          seed = 1))
    s <- f@relabeled$summary; o <- order(s$x)
    list(sl = s$x[o], ic = s$`(Intercept)`[o])
  }
  a <- run(x); b <- run(x * 1000)
  expect_lt(max(abs(a$sl - b$sl * 1000)), 0.1)
  expect_lt(max(abs(a$ic - b$ic)), 0.3)
  # and the answer itself is right
  expect_lt(max(abs(sort(a$sl) - c(-2, 2))), 0.25)
})

test_that("conjugate NIG sampler replaces RW on the FixedK regression path", {
  skip_on_cran()
  set.seed(3); n <- 120
  x <- runif(n, 0, 10); zc <- rbinom(n, 1, 0.5)
  y <- ifelse(zc == 1, 2 + 2 * x, -2 - 2 * x) + rnorm(n, 0, 1)
  X <- cbind(1, x)
  sp <- getDistribution("normal-reg")
  pr <- defaultPrior(sp, y, control = list(X = X))
  mc <- buildModelCode(sp, new("FixedKEngine", dirichletConc = 1), n = n, L = 2)
  cn <- nimix:::buildConstants(sp, pr, n)
  cn$K <- 2; cn$alphaVec <- rep(1, 2); cn$p <- ncol(X); cn$X <- X
  ini <- nimix:::componentInits(sp, pr, y, 2)
  m <- suppressMessages(nimble::nimbleModel(
    mc$code, constants = cn, data = list(y = y),
    inits = c(list(z = ini$alloc, weights = c(.5, .5)), ini$params),
    calculate = FALSE))
  conf <- suppressMessages(nimble::configureMCMC(m, print = FALSE))
  customizeSamplers(sp, conf, m)
  nm <- vapply(conf$getSamplers(), function(s) s$name, character(1))
  tg <- vapply(conf$getSamplers(), function(s) s$target[1], character(1))
  onBeta <- nm[grepl("betaTilde", tg)]
  expect_true(length(onBeta) > 0)
  expect_false(any(grepl("^RW", onBeta)))          # no RW left on beta
  expect_false(any(grepl("s2Tilde", tg) & grepl("^RW", nm)))  # nor on s2
})
