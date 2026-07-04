# Tests for the compiled-model cache (compile-once-reuse, v0.9.0). The cache
# logic is exercised without NIMBLE; the end-to-end reuse == fresh guarantee is
# a heavy NIMBLE test gated behind skip_on_cran.

test_that("nimixClearCache empties the cache and reports the count", {
  ce <- .nimixModelCache
  ce$entries <- list(list(key = list(a = 1), compiled = NULL),
                     list(key = list(a = 2), compiled = NULL))
  expect_equal(nimixClearCache(), 2L)
  expect_length(.nimixModelCache$entries, 0L)
  expect_equal(nimixClearCache(), 0L)
})

test_that("structural key ignores data/inits but tracks code, constants, monitors, spec", {
  spec <- .selectClusterSpec("normal", FALSE, 1L)
  mc1 <- list(code = quote({ y[1] ~ dnorm(0, 1) }), monitors = c("mu", "v"))
  k1  <- .cacheKey(mc1, constants = list(n = 10L, L = 6L), spec = spec)
  k1b <- .cacheKey(mc1, constants = list(n = 10L, L = 6L), spec = spec)
  expect_identical(k1, k1b)                        # deterministic
  # different truncation -> different key (different compiled structure)
  k2 <- .cacheKey(mc1, constants = list(n = 10L, L = 8L), spec = spec)
  expect_false(identical(k1, k2))
  # monitor order does not matter
  mc2 <- list(code = mc1$code, monitors = c("v", "mu"))
  expect_identical(.cacheKey(mc2, list(n = 10L, L = 6L), spec)$monitors,
                   k1$monitors)
})

test_that("cache get/put round-trips and LRU-evicts beyond the cap", {
  nimixClearCache()
  key <- function(i) list(id = i)
  for (i in seq_len(.NIMIX_CACHE_MAX)) .cachePut(key(i), list(tag = i))
  # newest is at the front, a hit moves it to the front
  hit <- .cacheGet(key(1L))
  expect_equal(hit$tag, 1L)
  expect_equal(.nimixModelCache$entries[[1]]$key, key(1L))
  # exceeding the cap evicts the least-recently-used entry
  .cachePut(key(99L), list(tag = 99L))
  expect_length(.nimixModelCache$entries, .NIMIX_CACHE_MAX)
  expect_null(.cacheGet(key(2L)))                  # the LRU one is gone
  nimixClearCache()
})

test_that("reuse is bit-for-bit identical to a fresh compile; structure misses", {
  skip_on_cran()
  set.seed(1); y <- c(rnorm(50, -3), rnorm(50, 3))
  mc <- list(niter = 1200, nburnin = 400)
  ce <- .nimixModelCache
  nimixClearCache(); ce$builds <- 0L

  f1 <- nimixClust(y, K_max = 6, method = "dpm", mcmcControl = mc, seed = 7)
  b1 <- ce$builds
  f2 <- nimixClust(y, K_max = 6, method = "dpm", mcmcControl = mc, seed = 7)
  expect_equal(ce$builds, b1)                       # reuse: no recompile
  f3 <- nimixClust(y, K_max = 6, method = "dpm",
                   mcmcControl = utils::modifyList(mc, list(reuse = FALSE)),
                   seed = 7)
  expect_equal(ce$builds, b1 + 1L)                  # fresh: recompiled

  # reuse, fresh, and self-reuse all agree exactly for identical (data, seed)
  expect_equal(f2@mcmcSamples, f3@mcmcSamples)
  expect_equal(f1@mcmcSamples, f2@mcmcSamples)

  # a different structure (K_max) is a cache miss and recompiles
  bPrev <- ce$builds
  nimixClust(y, K_max = 12, method = "dpm", mcmcControl = mc, seed = 7)
  expect_equal(ce$builds, bPrev + 1L)
  nimixClearCache()
})
