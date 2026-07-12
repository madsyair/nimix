# Tests for the configurable dispersed-initialisation ratio (mcmcControl$initRatio).
# The pure-R pieces (validation, seed-count) run everywhere; the end-to-end
# threading check only builds model code (no NIMBLE compile), so it stays fast.

test_that(".resolveInitRatio returns the default and validates input", {
  expect_equal(nimix:::.resolveInitRatio(list()), 0.8)               # not supplied
  expect_equal(nimix:::.resolveInitRatio(list(niter = 10)), 0.8)     # other controls only
  expect_equal(nimix:::.resolveInitRatio(list(initRatio = 0.6)), 0.6)
  expect_equal(nimix:::.resolveInitRatio(list(initRatio = 0.9)), 0.9)

  # out of the open interval (0, 1) / wrong type -> informative error
  expect_error(nimix:::.resolveInitRatio(list(initRatio = 0)),    "initRatio")
  expect_error(nimix:::.resolveInitRatio(list(initRatio = 1)),    "initRatio")
  expect_error(nimix:::.resolveInitRatio(list(initRatio = 1.5)),  "initRatio")
  expect_error(nimix:::.resolveInitRatio(list(initRatio = -0.1)), "initRatio")
  expect_error(nimix:::.resolveInitRatio(list(initRatio = "x")),  "initRatio")
  expect_error(nimix:::.resolveInitRatio(list(initRatio = NA_real_)),    "initRatio")
  expect_error(nimix:::.resolveInitRatio(list(initRatio = c(0.5, 0.7))), "initRatio")
})

test_that("a high but valid initRatio is accepted with a headroom warning", {
  expect_warning(r <- nimix:::.resolveInitRatio(list(initRatio = 0.96)), "headroom")
  expect_equal(r, 0.96)
  expect_warning(nimix:::.resolveInitRatio(list(initRatio = 0.95)), "headroom")
  # just below the warn threshold stays silent
  expect_silent(nimix:::.resolveInitRatio(list(initRatio = 0.9)))
})

test_that(".initRatioArg extracts the threaded value or falls back to default", {
  expect_equal(nimix:::.initRatioArg(initRatio = 0.5), 0.5)
  expect_equal(nimix:::.initRatioArg(), 0.8)             # direct call (e.g. unit test)
  expect_equal(nimix:::.initRatioArg(other = 1), 0.8)
})

test_that("initRatio scales the number of seeded clusters", {
  set.seed(1)
  Y <- rbind(matrix(rnorm(150 * 2, mean = -3), ncol = 2),
             matrix(rnorm(150 * 2, mean =  3), ncol = 2))
  sp <- nimix:::NormalMvSpec()
  pr <- defaultPrior(sp, Y)
  seed_k <- function(ir)
    length(unique(nimix:::componentInits(sp, pr, Y, 20, initRatio = ir)$alloc))

  # n = 300 -> ceil(sqrt(n)) = 18 caps the high end; floor(ir * 20) drives the rest
  expect_equal(seed_k(0.5), 10L)               # floor(0.5 * 20)
  expect_equal(seed_k(0.8), 16L)               # default ratio
  expect_lt(seed_k(0.5), seed_k(0.8))          # lower ratio -> fewer seeds
  expect_equal(nimix:::componentInits(sp, pr, Y, 20)$alloc |>
                 unique() |> length(), 16L)    # default when arg omitted
})

test_that("an invalid initRatio is rejected end-to-end before any MCMC run", {
  skip_if_not_installed("nimble")
  set.seed(1)
  Y <- rbind(matrix(rnorm(40 * 2, mean = -3), ncol = 2),
             matrix(rnorm(40 * 2, mean =  3), ncol = 2))
  expect_error(
    nimixClust(Y, distribution = "normal", method = "dpm",
               mcmcControl = list(niter = 200, nburnin = 50, initRatio = 1.5)),
    "initRatio")
})
