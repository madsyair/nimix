# clusterValidity(): thin adaptor over cluster/fpc on the Binder partition.
# The separated-vs-overlapping contrast below is the documented caveat made
# executable: overlapping components are a legitimate mixture model yet score
# low on geometric indices -- that is expected behaviour, not a defect.

test_that("well-separated clusters score high on all indices", {
  skip_on_cran()
  skip_if_not_installed("cluster")
  skip_if_not_installed("fpc")
  set.seed(7)
  y <- c(rnorm(60, -4, 0.7), rnorm(60, 4, 0.7))
  f <- nimixClust(y, K = 2, method = "fixedk",
                  mcmcControl = list(niter = 800, nburnin = 300), seed = 1)
  v <- clusterValidity(f)
  expect_named(v, c("silhouette", "dunn", "ch"))
  expect_gt(v[["silhouette"]], 0.5)
  expect_gt(v[["dunn"]], 0)
  expect_gt(v[["ch"]], 100)
})

test_that("overlapping components (a legitimate model) score LOW: the caveat", {
  skip_on_cran()
  skip_if_not_installed("cluster")
  skip_if_not_installed("fpc")
  set.seed(7)
  y <- c(rnorm(60, -0.7), rnorm(60, 0.7))
  f <- nimixClust(y, K = 2, method = "fixedk",
                  mcmcControl = list(niter = 800, nburnin = 300), seed = 1)
  v <- clusterValidity(f)
  # measured: silhouette 0.52, dunn 0 -- clearly below the separated case
  expect_lt(v[["silhouette"]], 0.7)
  expect_lt(v[["dunn"]], 0.5)
})

test_that("explicit partition equals the default binderPartition route", {
  skip_on_cran()
  skip_if_not_installed("cluster")
  skip_if_not_installed("fpc")
  set.seed(7)
  y <- c(rnorm(60, -4, 0.7), rnorm(60, 4, 0.7))
  f <- nimixClust(y, K = 2, method = "fixedk",
                  mcmcControl = list(niter = 600, nburnin = 250), seed = 1)
  bp <- binderPartition(f)
  expect_identical(clusterValidity(f),
                   clusterValidity(f, partition = bp$partition))
  # metric subsetting honours order and names
  vs <- clusterValidity(f, metrics = "silhouette")
  expect_named(vs, "silhouette")
})

test_that("regression fits and single-cluster partitions are refused", {
  skip_on_cran()
  set.seed(7)
  x <- runif(60, -2, 2)
  df <- data.frame(y = 2 * x + rnorm(60, 0, 0.5), x = x)
  fr <- nimixReg(y ~ x, df, K = 2, method = "fixedk",
                 mcmcControl = list(niter = 400, nburnin = 150), seed = 1)
  expect_error(clusterValidity(fr), "regressions")

  y <- rnorm(40)
  f1 <- nimixClust(y, K = 1, method = "fixedk",
                   mcmcControl = list(niter = 300, nburnin = 100), seed = 1)
  expect_error(clusterValidity(f1), "at least two")
})

test_that("partition length mismatch and bad dist are caught", {
  skip_on_cran()
  skip_if_not_installed("cluster")
  set.seed(7)
  y <- c(rnorm(40, -4), rnorm(40, 4))
  f <- nimixClust(y, K = 2, method = "fixedk",
                  mcmcControl = list(niter = 400, nburnin = 150), seed = 1)
  expect_error(clusterValidity(f, partition = c(1L, 2L)), "length")
  expect_error(clusterValidity(f, dist = matrix(0, 2, 2)), "dist object")
})
