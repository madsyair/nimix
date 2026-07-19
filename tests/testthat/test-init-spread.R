test_that("spread init separates by scale where k-means separates by location", {
  skip_on_cran()
  # The one case measured to defeat k-means seeding: two components with very
  # different variances (sd 0.3 vs 3) but overlapping means. k-means forces a
  # spatial boundary; spread bands by |y - median|, catching the scale split.
  set.seed(1)
  z <- rep(1:2, c(150, 150))
  y <- c(rnorm(150, 0, 0.3), rnorm(150, 0, 3))
  acc <- function(a, b) { tab <- table(a, b); sum(apply(tab, 2, max)) / length(b) }

  cl_spread <- nimix:::.initClusters(y, 2L, "spread")
  cl_kmeans <- nimix:::.initClusters(y, 2L, "kmeans")
  # spread's initial allocation is markedly better on this case
  expect_gt(acc(cl_spread, z), acc(cl_kmeans, z))
  expect_gt(acc(cl_spread, z), 0.85)
})

test_that("spread is accepted end-to-end and falls back to k-means for mv", {
  skip_on_cran()
  set.seed(1)
  y <- c(rnorm(120, 0, 0.3), rnorm(120, 0, 3))
  f <- nimixClust(y, K = 2, method = "fixedk", initMethod = "spread",
                  mcmcControl = list(niter = 500, nburnin = 200), seed = 1)
  expect_s4_class(f, "FitResult")

  # multivariate response: spread has no univariate analogue, so it must fall
  # back to k-means rather than erroring
  Y <- cbind(c(rnorm(80, -2), rnorm(80, 2)), rnorm(160))
  cl <- nimix:::.initClusters(Y, 2L, "spread")
  expect_false(is.null(cl))            # fell back, did not fail
  expect_length(cl, nrow(Y))
})

test_that("spread rejects nothing k-means would accept (validation)", {
  skip_on_cran()
  set.seed(2); y <- c(rnorm(60, -3), rnorm(60, 3))
  expect_s4_class(
    nimixClust(y, K = 2, method = "fixedk", initMethod = "spread",
               mcmcControl = list(niter = 300, nburnin = 100), seed = 1,
               verbose = FALSE),
    "FitResult")
})
