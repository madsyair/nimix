# drawsArray() / ppcData(): interop with bayesplot via plain base-R shapes.
# bayesplot itself stays in Suggests; these tests exercise the adaptors'
# structure and -- crucially -- the safety guard that refuses per-component
# draws before relabel(). The guard is the point: R-hat computed on a raw
# muTilde trace under label switching looks valid and means nothing.

test_that("drawsArray(invariant) preserves per-chain structure", {
  skip_on_cran()
  set.seed(2); y <- c(rnorm(50, -3), rnorm(50, 3))
  f <- nimixClust(y, method = "dpm", K_max = 8,
                  mcmcControl = list(niter = 600, nburnin = 250, nchains = 2),
                  seed = 1)
  da <- drawsArray(f)
  expect_equal(length(dim(da)), 3L)
  expect_equal(dim(da)[2L], 2L)
  expect_true(all(c("K", "entropy", "alpha") %in% dimnames(da)[[3L]]))
  # chain reconstruction is exact, element by element
  cid <- f@diagnostics$chainId
  expect_equal(as.numeric(da[, 2L, "K"]),
               as.numeric(f@Kposterior)[cid == 2L])
  expect_true(all(is.finite(da[, , "entropy"])))
})

test_that("drawsArray refuses component draws before relabel(), then serves them", {
  skip_on_cran()
  set.seed(2); y <- c(rnorm(50, -3), rnorm(50, 3))
  f <- nimixClust(y, K = 2, method = "fixedk",
                  mcmcControl = list(niter = 400, nburnin = 150), seed = 1)
  expect_error(drawsArray(f, "components"), "label permutations")
  fr <- relabel(f)
  dc <- drawsArray(fr, "components")
  expect_equal(length(dim(dc)), 3L)
  expect_equal(dim(dc)[2L], 1L)   # chain dim honestly collapsed post-relabel
  expect_true(all(c("mu[1]", "mu[2]", "s2[1]", "weight[1]") %in%
                    dimnames(dc)[[3L]]))
  # relabelled means must match the summary the user already sees
  s <- fr@relabeled$summary
  expect_equal(unname(sort(colMeans(dc[, 1L, c("mu[1]", "mu[2]")]))),
               sort(s$mu_mean), tolerance = 1e-8)
})

test_that("ppcData returns bayesplot-shaped y/yrep, with mv margins", {
  skip_on_cran()
  set.seed(2); y <- c(rnorm(50, -3), rnorm(50, 3))
  f <- nimixClust(y, K = 2, method = "fixedk",
                  mcmcControl = list(niter = 400, nburnin = 150), seed = 1)
  pd <- ppcData(f, ndraws = 20)
  expect_equal(length(pd$y), 100L)
  expect_equal(ncol(pd$yrep), 100L)
  expect_true(nrow(pd$yrep) <= 20L)
  # replicates mirror the data, so overlay plots are meaningful
  expect_lt(abs(mean(pd$yrep) - mean(pd$y)), 0.5)
  # multivariate: margin extraction and bounds check
  set.seed(4)
  Y <- rbind(matrix(rnorm(150, -3), 50, 3), matrix(rnorm(150, 3), 50, 3))
  fm <- nimixClust(Y, K = 2, method = "fixedk", distribution = "normal-mv",
                   mcmcControl = list(niter = 400, nburnin = 150), seed = 1)
  pd2 <- ppcData(fm, ndraws = 15, margin = 2L)
  expect_equal(pd2$y, as.numeric(Y[, 2L]))
  expect_equal(ncol(pd2$yrep), 100L)
  expect_error(ppcData(fm, margin = 7L), "margin")
})

test_that("outputs are accepted by bayesplot when it is installed", {
  skip_on_cran()
  skip_if_not_installed("bayesplot")
  set.seed(2); y <- c(rnorm(50, -3), rnorm(50, 3))
  f <- relabel(nimixClust(y, K = 2, method = "fixedk",
                          mcmcControl = list(niter = 400, nburnin = 150),
                          seed = 1))
  expect_s3_class(bayesplot::mcmc_trace(drawsArray(f)), "ggplot")
  expect_s3_class(bayesplot::mcmc_dens(drawsArray(f, "components")), "ggplot")
  pd <- ppcData(f, ndraws = 20)
  expect_s3_class(bayesplot::ppc_dens_overlay(pd$y, pd$yrep), "ggplot")
})
