# Tests for the multi-chain convergence diagnostics (v0.9.0). The split-Rhat /
# ESS helpers are exercised without NIMBLE; a real multi-chain fit is gated
# behind skip_on_cran.

test_that("split-Rhat is ~1 for agreeing chains and flags disagreeing chains", {
  set.seed(1)
  agree <- lapply(1:4, function(i) rnorm(500))
  expect_lt(abs(.splitRhat(agree) - 1), 0.05)

  # chains centred far apart do not mix -> large Rhat
  disagree <- list(rnorm(500, -5), rnorm(500, 5), rnorm(500, -5), rnorm(500, 5))
  expect_gt(.splitRhat(disagree), 1.1)
})

test_that("split-Rhat returns NA when undefined (one chain, or constant)", {
  expect_true(is.na(.splitRhat(list(rnorm(500)))))          # single chain
  expect_true(is.na(.splitRhat(list(rep(2, 500), rep(2, 500)))))  # constant
})

test_that(".sumESS adds per-chain effective sizes and ignores degenerate chains", {
  set.seed(2)
  chains <- lapply(1:3, function(i) rnorm(400))
  s <- .sumESS(chains)
  expect_true(is.finite(s) && s > 0)
  expect_equal(.sumESS(list(rep(1, 100))), 0)               # constant -> 0
})

test_that(".multiChainDiag assembles K (and optional alpha) diagnostics", {
  set.seed(3)
  K <- lapply(1:3, function(i) sample(2:3, 400, replace = TRUE))
  d <- .multiChainDiag(K)
  expect_equal(d$nchains, 3L)
  expect_true(all(c("RhatK", "essK") %in% names(d)))
  expect_false("RhatAlpha" %in% names(d))
  d2 <- .multiChainDiag(K, lapply(1:3, function(i) rgamma(400, 2, 2)))
  expect_true(all(c("RhatAlpha", "essAlpha") %in% names(d2)))
})

test_that("a multi-chain fit reuses the compiled model and reports Rhat/ESS", {
  skip_on_cran()
  set.seed(1); y <- c(rnorm(60, -3), rnorm(60, 3))
  ce <- .nimixModelCache
  nimixClearCache(); ce$builds <- 0L
  f <- nimixClust(y, K_max = 8, method = "dpm",
                  mcmcControl = list(niter = 1500, nburnin = 500, nchains = 3),
                  seed = 7)
  expect_equal(ce$builds, 1L)                    # only chain 1 compiles
  expect_equal(nrow(f@mcmcSamples), 3L * 1000L)  # pooled draws
  d <- f@diagnostics
  expect_equal(d$nchains, 3L)
  expect_true(is.finite(d$RhatK) && d$RhatK < 1.1)
  expect_true(is.finite(d$essK) && d$essK > 0)
  nimixClearCache()
})
