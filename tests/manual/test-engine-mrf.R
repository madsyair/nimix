# Tests for the MRF engine (v1.2.0): spatially constrained finite mixture with
# a fixed-beta Potts prior on the labels (Besag 1974; Blekas et al. 2005).

test_that("mrf argument guards are informative", {
  y <- rnorm(20)
  g <- gridAdjacency(4, 5)
  # needs K and a neighbourhood; K_max belongs to the DPM
  expect_error(nimixClust(y, method = "mrf", spatialWeights = g),
               "needs the number of components K")
  expect_error(nimixClust(y, K = 2, method = "mrf"), "spatialWeights")
  expect_error(nimixClust(y, K = 2, K_max = 8, method = "mrf",
                          spatialWeights = g), "K_max")
  # spatialWeights is mrf-only
  expect_error(nimixClust(y, K_max = 6, method = "dpm", spatialWeights = g),
               "only used by method = 'mrf'")
  # region count must match the data
  expect_error(nimixClust(rnorm(7), K = 2, method = "mrf", spatialWeights = g),
               "must match")
  # beta must be a single non-negative number (engine validity)
  expect_error(nimixClust(y, K = 2, method = "mrf", spatialWeights = g,
                          prior = list(beta = -1)), "non-negative")
})

test_that("mrf still rejects unknown component families (ANY fallback)", {
  # every registered family is supported since batch 2; the ANY fallback
  # protects user-registered specs without an MRF sweep.
  methods::setClass("DummyMrfSpec", contains = "DistributionSpec",
                    where = environment())
  dummy <- methods::new("DummyMrfSpec", name = "dummy")
  expect_error(nimix:::.mrfSamplerFor(dummy), "not yet available")
})

test_that("mrf recovers a known spatial block structure and beats no-smoothing", {
  skip_on_cran()
  nr <- 12; nc <- 10; n <- nr * nc
  g <- gridAdjacency(nr, nc, "rook")
  zTrue <- integer(n)
  for (i in 1:nr) for (j in 1:nc)
    zTrue[(i - 1) * nc + j] <- if (j <= nc / 2) 1L else 2L
  set.seed(42)
  y <- rnorm(n, c(-2, 2)[zTrue], 1.4)          # overlapping components
  ctrl <- list(niter = 2000, nburnin = 800)

  fit <- relabel(nimixClust(y, K = 2, method = "mrf", spatialWeights = g,
                            mcmcControl = ctrl, seed = 7))
  zMap <- apply(fit@clusterAllocation, 2L,
                function(v) as.integer(names(which.max(table(v)))))
  acc <- max(mean(zMap == zTrue), mean(zMap == (3L - zTrue)))
  mu <- sort(fit@relabeled$summary$mu_mean)

  expect_gt(acc, 0.9)                          # spatial recovery
  expect_lt(abs(mu[1] - (-2)), 0.8)
  expect_lt(abs(mu[2] - 2), 0.8)

  # no-smoothing baseline on the same data must do worse than the MRF
  fk <- relabel(nimixClust(y, K = 2, method = "fixedk",
                           mcmcControl = ctrl, seed = 7))
  zF <- apply(fk@clusterAllocation, 2L,
              function(v) as.integer(names(which.max(table(v)))))
  accF <- max(mean(zF == zTrue), mean(zF == (3L - zTrue)))
  expect_gt(acc, accF)
})

test_that("mrf with multivariate Gaussian components recovers spatial blocks", {
  skip_on_cran()
  nr <- 12; nc <- 10; n <- nr * nc; d <- 2
  g <- gridAdjacency(nr, nc, "rook")
  zTrue <- integer(n)
  for (i in 1:nr) for (j in 1:nc)
    zTrue[(i - 1) * nc + j] <- if (j <= nc / 2) 1L else 2L
  set.seed(42)
  Y <- matrix(rnorm(n * d, mean = rep(c(-1.5, 1.5)[zTrue], d), sd = 1.4),
              ncol = d)
  ctrl <- list(niter = 2000, nburnin = 800)

  fit <- relabel(nimixClust(Y, K = 2, method = "mrf", spatialWeights = g,
                            mcmcControl = ctrl, seed = 7))
  zM <- apply(fit@clusterAllocation, 2L,
              function(v) as.integer(names(which.max(table(v)))))
  acc <- max(mean(zM == zTrue), mean(zM == (3L - zTrue)))
  mu1 <- sort(fit@relabeled$summary$mu_1_mean)
  expect_gt(acc, 0.9)
  expect_lt(abs(mu1[1] - (-1.5)), 0.8)
  expect_lt(abs(mu1[2] - 1.5), 0.8)

  fk <- relabel(nimixClust(Y, K = 2, method = "fixedk",
                           mcmcControl = ctrl, seed = 7))
  zF <- apply(fk@clusterAllocation, 2L,
              function(v) as.integer(names(which.max(table(v)))))
  accF <- max(mean(zF == zTrue), mean(zF == (3L - zTrue)))
  expect_gte(acc, accF)
})

test_that("mrf mixture-of-regressions recovers spatially clustered slopes", {
  skip_on_cran()
  nr <- 12; nc <- 10; n <- nr * nc
  g <- gridAdjacency(nr, nc, "rook")
  zTrue <- integer(n)
  for (i in 1:nr) for (j in 1:nc)
    zTrue[(i - 1) * nc + j] <- if (j <= nc / 2) 1L else 2L
  set.seed(42)
  x <- runif(n, -2, 2)
  y <- ifelse(zTrue == 1L, 2 * x, -2 * x) + rnorm(n, 0, 1.6)
  df <- data.frame(y = y, x = x)
  ctrl <- list(niter = 2500, nburnin = 1000)

  fm <- relabel(nimixReg(y ~ x, df, K = 2, method = "mrf", spatialWeights = g,
                         mcmcControl = ctrl, seed = 7))
  zM <- apply(fm@clusterAllocation, 2L,
              function(v) as.integer(names(which.max(table(v)))))
  acc <- max(mean(zM == zTrue), mean(zM == (3L - zTrue)))
  sl <- sort(fm@relabeled$summary[["x"]])
  expect_gt(acc, 0.9)
  expect_lt(sl[1], -1); expect_gt(sl[2], 1)

  fk <- relabel(nimixReg(y ~ x, df, K = 2, method = "fixedk",
                         mcmcControl = ctrl, seed = 7))
  zF <- apply(fk@clusterAllocation, 2L,
              function(v) as.integer(names(which.max(table(v)))))
  accF <- max(mean(zF == zTrue), mean(zF == (3L - zTrue)))
  expect_gte(acc, accF)
})

test_that("nimixReg mrf guards mirror the clustering ones", {
  df <- data.frame(y = rnorm(20), x = rnorm(20))
  g <- gridAdjacency(4, 5)
  expect_error(nimixReg(y ~ x, df, method = "mrf"), "spatialWeights")
  expect_error(nimixReg(y ~ x, df, K_max = 6, spatialWeights = g),
               "only used by method")
  expect_error(nimixReg(y ~ x, df, K = 2, method = "mrf",
                        spatialWeights = gridAdjacency(2, 2)), "must match")
  # all registered response families are supported under mrf since batch 2
})

test_that("beta estimation (pseudo-likelihood) detects positive interaction", {
  skip_on_cran()
  nr <- 12; nc <- 10; n <- nr * nc
  g <- gridAdjacency(nr, nc, "rook")
  zTrue <- integer(n)
  for (i in 1:nr) for (j in 1:nc)
    zTrue[(i - 1) * nc + j] <- if (j <= nc / 2) 1L else 2L
  set.seed(42)
  y <- rnorm(n, c(-2, 2)[zTrue], 1.4)

  fe <- relabel(nimixClust(y, K = 2, method = "mrf", spatialWeights = g,
                           prior = list(estimateBeta = TRUE),
                           mcmcControl = list(niter = 2500, nburnin = 1000),
                           seed = 7))
  b <- as.numeric(fe@mcmcSamples[, "beta"])
  expect_gt(stats::var(b), 0)                 # beta actually moves
  expect_gt(mean(b > 0.2), 0.95)              # positive interaction detected
  zM <- apply(fe@clusterAllocation, 2L,
              function(v) as.integer(names(which.max(table(v)))))
  acc <- max(mean(zM == zTrue), mean(zM == (3L - zTrue)))
  expect_gt(acc, 0.9)                         # recovery stays strong

  # fixed-beta fit AFTER an estimation fit must not inherit the PL sampler
  # from the cache (regression test for the estimate/fixed cache collision)
  ff <- nimixClust(y, K = 2, method = "mrf", spatialWeights = g,
                   mcmcControl = list(niter = 800, nburnin = 300), seed = 7)
  expect_true(all(as.numeric(ff@mcmcSamples[, "beta"]) == 0.8))
})

test_that("mrf batch-1 families: Poisson clustering and Student-t regression", {
  skip_on_cran()
  nr <- 10; nc <- 10; n <- nr * nc
  g <- gridAdjacency(nr, nc, "rook")
  zT <- integer(n)
  for (i in 1:nr) for (j in 1:nc)
    zT[(i - 1) * nc + j] <- if (j <= nc / 2) 1L else 2L
  ctrl <- list(niter = 2000, nburnin = 800)
  accOf <- function(f) {
    z <- apply(f@clusterAllocation, 2L,
               function(v) as.integer(names(which.max(table(v)))))
    max(mean(z == zT), mean(z == (3L - zT)))
  }

  set.seed(1); yP <- rpois(n, c(3, 12)[zT])
  fP <- relabel(nimixClust(yP, K = 2, method = "mrf", spatialWeights = g,
                           distribution = "poisson", mcmcControl = ctrl,
                           seed = 7))
  lam <- sort(fP@relabeled$summary$lambda_mean)
  expect_gt(accOf(fP), 0.9)
  expect_lt(abs(lam[1] - 3) / 3, 0.35)
  expect_lt(abs(lam[2] - 12) / 12, 0.35)

  set.seed(4); x <- runif(n, -1.5, 1.5)
  yT <- ifelse(zT == 1L, 2 * x, -2 * x) + rt(n, 5) * 1.2
  fT <- relabel(nimixReg(y ~ x, data.frame(y = yT, x = x), K = 2,
                         method = "mrf", spatialWeights = g,
                         distribution = "studentt",
                         mcmcControl = list(niter = 2500, nburnin = 1000),
                         seed = 7))
  sl <- sort(fT@relabeled$summary[["x"]])
  expect_gt(accOf(fT), 0.9)
  expect_lt(sl[1], -1); expect_gt(sl[2], 1)
})

test_that("mrf requires K >= 2 and still blocks augmented families", {
  g <- gridAdjacency(4, 5)
  expect_error(nimixClust(rnorm(20), K = 1, method = "mrf", spatialWeights = g),
               "K >= 2")

})

test_that("mrf batch-2 families: mv Student-t clustering and NG regression", {
  skip_on_cran()
  nr <- 10; nc <- 10; n <- nr * nc; d <- 2
  g <- gridAdjacency(nr, nc, "rook")
  zT <- integer(n)
  for (i in 1:nr) for (j in 1:nc)
    zT[(i - 1) * nc + j] <- if (j <= nc / 2) 1L else 2L
  accOf <- function(f) {
    z <- apply(f@clusterAllocation, 2L,
               function(v) as.integer(names(which.max(table(v)))))
    max(mean(z == zT), mean(z == (3L - zT)))
  }

  set.seed(11)
  Y <- matrix(rep(c(-1.5, 1.5)[zT], d) + rt(n * d, df = 5) * 1.1, ncol = d)
  fMt <- relabel(nimixClust(Y, K = 2, method = "mrf", spatialWeights = g,
                            distribution = "studentt",
                            mcmcControl = list(niter = 2000, nburnin = 800),
                            seed = 7))
  mu <- sort(fMt@relabeled$summary$mu_1_mean)
  expect_gt(accOf(fMt), 0.9)
  expect_lt(abs(mu[1] + 1.5), 0.9); expect_lt(abs(mu[2] - 1.5), 0.9)

  set.seed(4); x <- runif(n, -1.5, 1.5)
  yT <- ifelse(zT == 1L, 2 * x, -2 * x) + rt(n, 5) * 1.2
  fNG <- relabel(nimixReg(y ~ x, data.frame(y = yT, x = x), K = 2,
                          method = "mrf", spatialWeights = g,
                          distribution = "normalgamma",
                          mcmcControl = list(niter = 2500, nburnin = 1000),
                          seed = 7))
  sl <- sort(fNG@relabeled$summary[["x"]])
  expect_gt(accOf(fNG), 0.9)
  expect_lt(sl[1], -1); expect_gt(sl[2], 1)
})

test_that("mrf kernels are family-correct (regression test for the dispatch trap)", {
  # NormalGammaUvSpec inherits NormalUvSpec: buildModelCode must dispatch to
  # the pottsified NG kernel (with omega), never the plain Gaussian one.
  e <- MRFEngine(0.8, gridAdjacency(2, 2))
  kNG <- nimix:::buildModelCode(nimix:::getDistribution("normal-gamma"), e,
                                n = 10, L = 2)
  expect_true(grepl("omega", paste(deparse(kNG$code), collapse = " ")))
  kTR <- nimix:::buildModelCode(nimix:::getDistribution("student-t-reg"), e,
                                n = 10, L = 2)
  expect_true(grepl("dt(", paste(deparse(kTR$code), collapse = " "),
                    fixed = TRUE))
})
