# Ferreira-Steel skew multivariate Normal (FS 2007), A = chol(Sigma), O = I.

test_that("dskewmvn reduces to the multivariate Normal at gamma = 1", {
  Sg <- matrix(c(2, .8, .8, 1.2), 2, 2)
  dmvn <- function(x, mu, S) {
    k <- length(x)
    as.numeric(-0.5 * (k * log(2 * pi) + determinant(S)$modulus +
                         t(x - mu) %*% solve(S) %*% (x - mu)))
  }
  set.seed(1)
  for (i in 1:10) {
    x <- rnorm(2, 0, 3)
    expect_lt(abs(dskewmvn(x, c(0, 0), Sg, c(1, 1), log = TRUE) -
                    dmvn(x, c(0, 0), Sg)), 1e-10)
  }
})

test_that("dskewmvn integrates to one and rskewmvn matches the density", {
  Sg <- matrix(c(2, .8, .8, 1.2), 2, 2)
  gr <- seq(-10, 10, 0.05)
  z <- outer(gr, gr, function(a, b) dskewmvn(cbind(a, b), c(0, 0), Sg,
                                             c(0.6, 1.8)))
  expect_lt(abs(sum(z) * 0.05^2 - 1), 1e-3)
  # with Sigma = I the margins are independent FS-skew Normals:
  # P(X_j > mu_j) = gamma_j^2 / (1 + gamma_j^2)
  set.seed(6)
  Y <- rskewmvn(2e4, c(0, 0), diag(2), c(2, 0.5))
  expect_lt(abs(mean(Y[, 1] > 0) - 0.8), 0.02)
  expect_lt(abs(mean(Y[, 2] > 0) - 0.2), 0.02)
})

test_that("m = 1 skew-mv-Normal agrees with the univariate FSSN", {
  xg <- seq(-4, 4, 0.5)
  for (g in c(0.5, 1, 2)) {
    a <- dskewmvn(matrix(xg, ncol = 1), 0, matrix(1), g, log = TRUE)
    b <- dfssn(xg, 0, 1, g, log = TRUE)
    expect_lt(max(abs(a - b)), 1e-10)
  }
})

test_that("compiled dSkewMvN_k equals the R reference", {
  skip_on_cran()
  suppressMessages(nimix:::.nimixEnsureMSNBurr())
  cD <- nimble::compileNimble(get("dSkewMvN_k", envir = globalenv()))
  Sg <- matrix(c(2, .8, .8, 1.2), 2, 2)
  set.seed(7); err <- 0
  for (i in 1:20) {
    x <- rnorm(2, 0, 3); g <- exp(rnorm(2, 0, .6))
    err <- max(err, abs(cD(x, c(0, 0), Sg, g, log = 1) -
                          dskewmvn(x, c(0, 0), Sg, g, log = TRUE)))
  }
  expect_lt(err, 1e-8)
})

test_that("skew-mv-Normal mixture recovers structure (FixedK + DPM)", {
  skip_on_cran()
  set.seed(8)
  Sg <- matrix(c(1, .4, .4, .8), 2, 2)
  Y <- rbind(rskewmvn(120, c(-3, -3), Sg, c(0.6, 1.6)),
             rskewmvn(120, c( 3,  3), Sg, c(1.6, 0.6)))
  zt <- rep(1:2, each = 120)
  f <- relabel(nimixClust(Y, K = 2, method = "fixedk",
                          distribution = "skewnormal-mv",
                          mcmcControl = list(niter = 2200, nburnin = 900),
                          seed = 3))
  s <- f@relabeled$summary; o <- order(s$mu_1_mean)
  expect_lt(abs(s$mu_1_mean[o][1] - (-3)), 0.8)
  expect_lt(abs(s$mu_1_mean[o][2] - 3), 0.8)
  # per-dimension skewness recovered with the right direction
  expect_lt(s$gamma_1_mean[o][1], 1); expect_gt(s$gamma_1_mean[o][2], 1)
  expect_gt(s$gamma_2_mean[o][1], 1); expect_lt(s$gamma_2_mean[o][2], 1)
  z <- apply(f@clusterAllocation, 2L,
             function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(z == zt), mean(z == (3L - zt))), 0.9)

  fd <- relabel(nimixClust(Y, method = "dpm", distribution = "skewnormal-mv",
                           K_max = 6,
                           mcmcControl = list(niter = 1500, nburnin = 600),
                           seed = 3))
  expect_equal(fd@relabeled$modalK, 2L)
})
