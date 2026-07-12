# FS skew mv independent-Student with the orthogonal factor O estimated via the
# Householder angle theta (m = 2).

test_that("dskewmvito reduces correctly (nu -> Inf, theta = 0)", {
  Sg <- matrix(c(2, .8, .8, 1.2), 2, 2); x <- c(0.7, -1.3)
  # nu -> Inf recovers the skew-Normal-O family at the same theta
  for (th in c(-0.3, 0, 0.3))
    expect_lt(abs(dskewmvito(x, c(0, 0), Sg, c(0.6, 1.8), c(1e6, 1e6), th,
                             log = TRUE) -
                    dskewmvno(x, c(0, 0), Sg, c(0.6, 1.8), th, log = TRUE)),
              1e-3)
  # theta = 0 gives O = diag(1, -1): equals skewistudent-mv with 1/gamma_2
  set.seed(5)
  for (i in 1:10) {
    xx <- rnorm(2, 0, 3); g <- exp(rnorm(2, 0, .6))
    expect_lt(abs(dskewmvito(xx, c(0, 0), Sg, g, c(5, 9), 0, log = TRUE) -
                    dskewmvit(xx, c(0, 0), Sg, c(g[1], 1 / g[2]), c(5, 9),
                              log = TRUE)), 1e-10)
  }
})

test_that("dskewmvito integrates to one", {
  Sg <- matrix(c(2, .8, .8, 1.2), 2, 2)
  for (th in c(-0.3, 0, 0.3)) {
    gr <- seq(-30, 30, 0.05)
    z <- outer(gr, gr, function(a, b)
      dskewmvito(cbind(a, b), c(0, 0), Sg, c(0.6, 1.8), c(6, 9), th))
    expect_lt(abs(sum(z) * 0.05^2 - 1), 5e-3)
  }
})

test_that("theta stays identified under symmetry (independent-t is not spherical)", {
  # Contrast with skewnormal-mv-o, where gamma = 1 makes the density
  # theta-invariant because spherical Normal errors carry no direction.
  Sg <- matrix(c(2, .8, .8, 1.2), 2, 2); x <- c(0.7, -1.3)
  base <- dskewmvito(x, c(0, 0), Sg, c(1, 1), c(6, 6), 0, log = TRUE)
  dev <- max(vapply(c(-0.35, 0.35), function(th)
    abs(dskewmvito(x, c(0, 0), Sg, c(1, 1), c(6, 6), th, log = TRUE) - base),
    numeric(1)))
  expect_gt(dev, 1e-3)                        # theta still matters at gamma = 1
  # letting nu -> Inf restores sphericity and hence theta-invariance
  base2 <- dskewmvito(x, c(0, 0), Sg, c(1, 1), c(1e7, 1e7), 0, log = TRUE)
  dev2 <- max(vapply(c(-0.35, 0.35), function(th)
    abs(dskewmvito(x, c(0, 0), Sg, c(1, 1), c(1e7, 1e7), th, log = TRUE) - base2),
    numeric(1)))
  expect_lt(dev2, 1e-5)
})

test_that("compiled dSkewMvITO_k equals the R reference", {
  skip_on_cran()
  suppressMessages(nimix:::.nimixEnsureMSNBurr())
  cK <- nimble::compileNimble(get("dSkewMvITO_k", envir = globalenv()))
  Sg <- matrix(c(2, .8, .8, 1.2), 2, 2)
  set.seed(3); err <- 0
  for (i in 1:20) {
    x <- rnorm(2, 0, 3); g <- exp(rnorm(2, 0, .6))
    nu <- c(3 + rexp(1, .2), 3 + rexp(1, .2)); th <- runif(1, -pi/8, pi/8)
    err <- max(err, abs(cK(x, c(0, 0), Sg, g, nu, th, log = 1) -
                          dskewmvito(x, c(0, 0), Sg, g, nu, th, log = TRUE)))
  }
  expect_lt(err, 1e-8)
})

test_that("skewistudent-mv-o recovers theta and partition", {
  skip_on_cran()
  set.seed(24)
  Sg <- matrix(c(1, .3, .3, .8), 2, 2)
  Y <- rbind(rskewmvito(400, c(-5, -5), Sg, c(0.4, 2.5), c(7, 7),  0.30),
             rskewmvito(400, c( 5,  5), Sg, c(2.5, 0.4), c(7, 7), -0.25))
  zt <- rep(1:2, each = 400)
  f <- relabel(nimixClust(Y, K = 2, method = "fixedk",
                          distribution = "skewistudent-mv-o",
                          mcmcControl = list(niter = 3000, nburnin = 1300),
                          seed = 3))
  s <- f@relabeled$summary; o <- order(s$mu_1_mean)
  # sign and magnitude of the angles; the mirror mode makes tight coverage
  # unrealistic at this n, so we assert direction plus a tolerance
  expect_gt(s$theta_mean[o][1], 0)
  expect_lt(s$theta_mean[o][2], 0)
  expect_lt(abs(s$theta_mean[o][1] - 0.30), 0.12)
  expect_lt(abs(s$theta_mean[o][2] - (-0.25)), 0.12)
  expect_true(all(s$nu_1_mean > 2), all(s$nu_2_mean > 2))
  z <- apply(f@clusterAllocation, 2L,
             function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(z == zt), mean(z == (3L - zt))), 0.9)
})

test_that("both -o families run under DPM; univariate input rejected", {
  skip_on_cran()
  set.seed(24); Sg <- matrix(c(1, .3, .3, .8), 2, 2)
  Y <- rbind(rskewmvito(200, c(-6, -6), Sg, c(0.4, 2.5), c(7, 7),  0.30),
             rskewmvito(200, c( 6,  6), Sg, c(2.5, 0.4), c(7, 7), -0.25))
  fd <- relabel(nimixClust(Y, method = "dpm",
                           distribution = "skewistudent-mv-o", K_max = 6,
                           mcmcControl = list(niter = 1500, nburnin = 600),
                           seed = 3))
  expect_equal(fd@relabeled$modalK, 2L)

  # Pre-C.4 d = 3 was refused with "m = 2 only"; it now ROUTES to the
  # general-m spec (full routing + recovery asserted in
  # test-skew-mv-o-general.R). Univariate input must still be rejected.
  expect_error(nimixClust(rnorm(30), K = 2, method = "fixedk",
                          distribution = "skewistudent-mv-o"),
               "multivariate")
})

test_that("MRF engine supports both estimated-O families", {
  skip_on_cran()
  g <- gridAdjacency(8, 8, "rook"); nn <- 64
  zb <- integer(nn)
  for (i in 1:8) for (j in 1:8) zb[(i - 1) * 8 + j] <- if (j <= 4) 1L else 2L
  Sg <- matrix(c(0.6, 0.2, 0.2, 0.5), 2, 2)

  set.seed(30)
  Y <- matrix(NA, nn, 2)
  for (i in 1:nn) Y[i, ] <- if (zb[i] == 1L)
    rskewmvno(1, c(-3, -3), Sg, c(0.5, 2.0),  0.30) else
    rskewmvno(1, c( 3,  3), Sg, c(2.0, 0.5), -0.25)
  fm <- relabel(nimixClust(Y, K = 2, method = "mrf", spatialWeights = g,
                           distribution = "skewnormal-mv-o",
                           mcmcControl = list(niter = 1200, nburnin = 500),
                           seed = 3))
  z <- apply(fm@clusterAllocation, 2L,
             function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(z == zb), mean(z == (3L - zb))), 0.9)

  set.seed(31)
  Y2 <- matrix(NA, nn, 2)
  for (i in 1:nn) Y2[i, ] <- if (zb[i] == 1L)
    rskewmvito(1, c(-4, -4), Sg, c(0.5, 2.0), c(7, 7),  0.30) else
    rskewmvito(1, c( 4,  4), Sg, c(2.0, 0.5), c(7, 7), -0.25)
  fm2 <- relabel(nimixClust(Y2, K = 2, method = "mrf", spatialWeights = g,
                            distribution = "skewistudent-mv-o",
                            mcmcControl = list(niter = 1200, nburnin = 500),
                            seed = 3))
  z2 <- apply(fm2@clusterAllocation, 2L,
              function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(z2 == zb), mean(z2 == (3L - zb))), 0.9)
})
