# Label-free partition summaries: psm() and binderPartition().
# These complement relabel(): relabel() conditions on the modal number of
# clusters (necessarily discarding other iterations) to align component
# parameters; psm/binder use EVERY draw because pairwise co-clustering is
# invariant to labels and to K.

test_that("psm is a valid similarity matrix and uses all draws", {
  skip_on_cran()
  set.seed(2); y <- c(rnorm(60, -3), rnorm(60, 3))
  f <- nimixClust(y, method = "dpm", K_max = 8,
                  mcmcControl = list(niter = 1200, nburnin = 500), seed = 1)
  S <- psm(f)
  n <- length(y)
  expect_equal(dim(S), c(n, n))
  expect_equal(S, t(S))
  expect_true(all(diag(S) == 1))
  expect_true(all(S >= 0 & S <= 1))
  # sharp block structure on well-separated clusters
  within <- mean(S[1:60, 1:60][upper.tri(S[1:60, 1:60])])
  between <- mean(S[1:60, 61:120])
  expect_gt(within, 0.9)
  expect_lt(between, 0.1)
})

test_that("binderPartition recovers the truth and needs no relabelling", {
  skip_on_cran()
  set.seed(2); y <- c(rnorm(60, -3), rnorm(60, 3)); zt <- rep(1:2, each = 60)
  f <- nimixClust(y, method = "dpm", K_max = 8,
                  mcmcControl = list(niter = 1200, nburnin = 500), seed = 1)
  # deliberately NO relabel() call anywhere
  bp <- binderPartition(f)
  expect_equal(bp$K, 2L)
  acc <- max(mean(bp$partition == zt), mean(bp$partition == (3L - zt)))
  expect_gt(acc, 0.95)
  # labels are recoded 1..K
  expect_setequal(unique(bp$partition), seq_len(bp$K))
  # the selected draw's co-clustering really attains the reported score
  z <- f@clusterAllocation[bp$draw, ]
  D <- outer(z, z, "==") - bp$psm
  expect_equal(sum(D * D), bp$score)
})

test_that("psm/binder work under FixedK too, and psm can be precomputed", {
  skip_on_cran()
  set.seed(3); y <- c(rnorm(50, -3), rnorm(50, 3))
  f <- nimixClust(y, K = 2, method = "fixedk",
                  mcmcControl = list(niter = 500, nburnin = 200), seed = 1)
  S <- psm(f)
  bp1 <- binderPartition(f)
  bp2 <- binderPartition(f, S = S)
  expect_equal(bp1$partition, bp2$partition)
  expect_equal(bp1$K, 2L)
})
