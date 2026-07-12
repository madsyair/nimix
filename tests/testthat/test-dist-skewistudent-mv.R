# Ferreira-Steel skew multivariate independent-Student (FS 2007, Sec 5.2).
# Closed-form: no lambda augmentation, so no data-augmentation mixing penalty.

test_that("dskewmvit reduces correctly (nu -> Inf, gamma = 1, m = 1)", {
  Sg <- matrix(c(2, .8, .8, 1.2), 2, 2)
  xg <- c(0.7, -1.3)
  # nu -> Inf recovers the skew multivariate Normal
  expect_lt(abs(dskewmvit(xg, c(0, 0), Sg, c(0.6, 1.8), c(1e6, 1e6), log = TRUE) -
                  dskewmvn(xg, c(0, 0), Sg, c(0.6, 1.8), log = TRUE)), 1e-3)
  # gamma = 1 with Sigma = I is a product of independent Student margins
  expect_lt(abs(dskewmvit(xg, c(0, 0), diag(2), c(1, 1), c(5, 8), log = TRUE) -
                  (dt(xg[1], 5, log = TRUE) + dt(xg[2], 8, log = TRUE))), 1e-10)
  # m = 1 is the univariate fsst
  xv <- seq(-4, 4, 0.5)
  for (g in c(0.6, 1, 1.7))
    expect_lt(max(abs(dskewmvit(matrix(xv, ncol = 1), 0, matrix(1), g, 7,
                                log = TRUE) -
                        dfsst(xv, 0, 1, g, 7, log = TRUE))), 1e-10)
})

test_that("dskewmvit integrates to one and rskewmvit matches it", {
  Sg <- matrix(c(2, .8, .8, 1.2), 2, 2)
  gr <- seq(-25, 25, 0.05)
  z <- outer(gr, gr, function(a, b)
    dskewmvit(cbind(a, b), c(0, 0), Sg, c(0.6, 1.8), c(6, 9)))
  expect_lt(abs(sum(z) * 0.05^2 - 1), 5e-3)
  set.seed(9)
  Y <- rskewmvit(3e4, c(0, 0), diag(2), c(2, 0.5), c(8, 8))
  expect_lt(abs(mean(Y[, 1] > 0) - 0.8), 0.02)   # gamma^2/(1+gamma^2)
  expect_lt(abs(mean(Y[, 2] > 0) - 0.2), 0.02)
})

test_that("compiled dSkewMvIT_k equals the R reference", {
  skip_on_cran()
  suppressMessages(nimix:::.nimixEnsureMSNBurr())
  cD <- nimble::compileNimble(get("dSkewMvIT_k", envir = globalenv()))
  Sg <- matrix(c(2, .8, .8, 1.2), 2, 2)
  set.seed(10); err <- 0
  for (i in 1:20) {
    x <- rnorm(2, 0, 3); g <- exp(rnorm(2, 0, .6))
    nu <- c(3 + rexp(1, .2), 3 + rexp(1, .2))
    err <- max(err, abs(cD(x, c(0, 0), Sg, g, nu, log = 1) -
                          dskewmvit(x, c(0, 0), Sg, g, nu, log = TRUE)))
  }
  expect_lt(err, 1e-8)
})

test_that("skew-mv-IStudent mixture recovers structure (FixedK + DPM)", {
  skip_on_cran()
  set.seed(11)
  Sg <- matrix(c(1, .4, .4, .8), 2, 2)
  Y <- rbind(rskewmvit(120, c(-4, -4), Sg, c(0.6, 1.6), c(6, 6)),
             rskewmvit(120, c( 4,  4), Sg, c(1.6, 0.6), c(6, 6)))
  zt <- rep(1:2, each = 120)
  f <- relabel(nimixClust(Y, K = 2, method = "fixedk",
                          distribution = "skewistudent-mv",
                          mcmcControl = list(niter = 2200, nburnin = 900),
                          seed = 3))
  s <- f@relabeled$summary; o <- order(s$mu_1_mean)
  expect_lt(abs(s$mu_1_mean[o][1] - (-4)), 1.0)
  expect_lt(abs(s$mu_1_mean[o][2] - 4), 1.0)
  expect_lt(s$gamma_1_mean[o][1], 1); expect_gt(s$gamma_1_mean[o][2], 1)
  expect_gt(s$gamma_2_mean[o][1], 1); expect_lt(s$gamma_2_mean[o][2], 1)
  # nu is weakly identified; only the truncation bound is asserted
  expect_true(all(s$nu_1_mean > 2)); expect_true(all(s$nu_2_mean > 2))
  z <- apply(f@clusterAllocation, 2L,
             function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(z == zt), mean(z == (3L - zt))), 0.9)

  # dCRP must cluster mu (vector), Sigma (matrix), gamma and nu (element-wise,
  # truncated stochastic) simultaneously
  fd <- relabel(nimixClust(Y, method = "dpm", distribution = "skewistudent-mv",
                           K_max = 6,
                           mcmcControl = list(niter = 1500, nburnin = 600),
                           seed = 3))
  expect_equal(fd@relabeled$modalK, 2L)
})

test_that("MRF engine supports both skew multivariate families", {
  skip_on_cran()
  g <- gridAdjacency(8, 8, "rook"); nn <- 64
  zb <- integer(nn)
  for (i in 1:8) for (j in 1:8) zb[(i - 1) * 8 + j] <- if (j <= 4) 1L else 2L
  Sg <- matrix(c(0.6, 0.2, 0.2, 0.5), 2, 2)

  set.seed(12)
  Y <- matrix(NA, nn, 2)
  for (i in 1:nn) Y[i, ] <- if (zb[i] == 1L)
    rskewmvn(1, c(-3, -3), Sg, c(0.6, 1.6)) else
    rskewmvn(1, c(3, 3), Sg, c(1.6, 0.6))
  fm <- relabel(nimixClust(Y, K = 2, method = "mrf", spatialWeights = g,
                           distribution = "skewnormal-mv",
                           mcmcControl = list(niter = 1200, nburnin = 500),
                           seed = 3))
  z <- apply(fm@clusterAllocation, 2L,
             function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(z == zb), mean(z == (3L - zb))), 0.9)

  set.seed(13)
  Y2 <- matrix(NA, nn, 2)
  for (i in 1:nn) Y2[i, ] <- if (zb[i] == 1L)
    rskewmvit(1, c(-4, -4), Sg, c(0.6, 1.6), c(7, 7)) else
    rskewmvit(1, c(4, 4), Sg, c(1.6, 0.6), c(7, 7))
  fm2 <- relabel(nimixClust(Y2, K = 2, method = "mrf", spatialWeights = g,
                            distribution = "skewistudent-mv",
                            mcmcControl = list(niter = 1200, nburnin = 500),
                            seed = 3))
  z2 <- apply(fm2@clusterAllocation, 2L,
              function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(z2 == zb), mean(z2 == (3L - zb))), 0.9)
})
