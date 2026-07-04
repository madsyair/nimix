# MSNBurr / MSNBurr-IIa neo-normal families (Iriawan 2000; Choir 2020):
# numerical-stability guarantees + mixture integration.

test_that("R-level densities are numerically stable and correct", {
  # integrate to one across a wide alpha grid (incl. tiny alpha)
  for (a in c(0.05, 0.5, 1, 5, 100)) {
    IM <- integrate(function(y) dmsnburr(y, 1, 2, a), -Inf, Inf,
                    rel.tol = 1e-9)$value
    I2 <- integrate(function(y) dmsnburr2a(y, 1, 2, a), -Inf, Inf,
                    rel.tol = 1e-9)$value
    expect_lt(abs(IM - 1), 1e-5)
    expect_lt(abs(I2 - 1), 1e-5)
  }
  # finite log-densities far into the tails (stability guarantee)
  ext <- c(-1000, -100, 100, 1000)
  expect_true(all(is.finite(dmsnburr(ext, 1, 2, 0.5, log = TRUE))))
  expect_true(all(is.finite(dmsnburr2a(ext, 1, 2, 0.5, log = TRUE))))
  # alpha = 1 is exactly logistic(mu, sigma/omega)
  om1 <- exp(nimix:::.log_omega_msnburr(1))
  xg <- seq(-8, 8, 0.5)
  expect_lt(max(abs(dmsnburr(xg, 1, 2, 1, log = TRUE) -
                    dlogis(xg, 1, 2 / om1, log = TRUE))), 1e-8)
  # MSNBurr and MSNBurr-IIa are mirror images: f_IIa(x) = f(2mu - x)
  expect_lt(max(abs(dmsnburr2a(xg, 1, 2, 0.6) -
                    dmsnburr(2 * 1 - xg, 1, 2, 0.6))), 1e-8)
})

test_that("NIMBLE densities equal the R references (stability preserved)", {
  skip_on_cran()
  suppressMessages(nimix:::.nimixEnsureMSNBurr())
  cdM <- nimble::compileNimble(get("dMSNBurr_k", envir = globalenv()))
  cd2 <- nimble::compileNimble(get("dMSNBurr2a_k", envir = globalenv()))
  grid <- expand.grid(x = c(-500, -8, 0, 8, 500), a = c(0.05, 0.6, 1, 80))
  vM <- mapply(function(x, a) cdM(x, 1, 2, a, log = 1), grid$x, grid$a)
  v2 <- mapply(function(x, a) cd2(x, 1, 2, a, log = 1), grid$x, grid$a)
  expect_lt(max(abs(vM - dmsnburr(grid$x, 1, 2, grid$a, log = TRUE))), 1e-8)
  expect_lt(max(abs(v2 - dmsnburr2a(grid$x, 1, 2, grid$a, log = TRUE))), 1e-8)
  expect_true(all(is.finite(vM)) && all(is.finite(v2)))
})

test_that("MSNBurr mixtures recover a known skew-component structure", {
  skip_on_cran()
  set.seed(7)
  y <- c(rmsnburr(160, -2.5, 0.9, 0.5), rmsnburr(140, 2.5, 0.9, 3))
  zTrue <- rep(1:2, c(160, 140))
  f <- relabel(nimixClust(y, K = 2, method = "fixedk", distribution = "msnburr",
                          mcmcControl = list(niter = 4000, nburnin = 1500),
                          seed = 3))
  s <- f@relabeled$summary; o <- order(s$mu_mean)
  expect_lt(abs(s$mu_mean[o][1] - (-2.5)), 0.6)
  expect_lt(abs(s$mu_mean[o][2] - 2.5), 0.6)
  zM <- apply(f@clusterAllocation, 2L,
              function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(zM == zTrue), mean(zM == (3L - zTrue))), 0.9)

  # posterior predictive check on the well-specified model should not flag
  pc <- ppCheck(f, nrep = 120, statistics = c("mean", "sd"))
  expect_true(all(pc$ppp > 0.05 & pc$ppp < 0.95))
})

test_that("MSNBurr-IIa works across all three engines", {
  skip_on_cran()
  set.seed(8)
  y <- c(rmsnburr2a(150, -2.5, 0.9, 3), rmsnburr2a(150, 2.5, 0.9, 0.5))
  zT <- rep(1:2, each = 150)
  ff <- relabel(nimixClust(y, K = 2, method = "fixedk",
                           distribution = "msnburr2a",
                           mcmcControl = list(niter = 3000, nburnin = 1200),
                           seed = 3))
  zf <- apply(ff@clusterAllocation, 2L,
              function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(zf == zT), mean(zf == (3L - zT))), 0.9)

  g <- gridAdjacency(10, 10, "rook"); nn <- 100
  zb <- integer(nn)
  for (i in 1:10) for (j in 1:10) zb[(i - 1) * 10 + j] <- if (j <= 5) 1L else 2L
  set.seed(9)
  ym <- ifelse(zb == 1L, rmsnburr(nn, -2, 0.8, 0.5), rmsnburr(nn, 2, 0.8, 0.5))
  fm <- relabel(nimixClust(ym, K = 2, method = "mrf", spatialWeights = g,
                           distribution = "msnburr",
                           mcmcControl = list(niter = 2000, nburnin = 800),
                           seed = 3))
  zm <- apply(fm@clusterAllocation, 2L,
              function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(zm == zb), mean(zm == (3L - zb))), 0.9)
})

test_that("predict() returns a valid predictive density for MSNBurr", {
  skip_on_cran()
  set.seed(7)
  y <- c(rmsnburr(120, -2.5, 0.9, 0.5), rmsnburr(120, 2.5, 0.9, 3))
  f <- relabel(nimixClust(y, K = 2, method = "fixedk",
                          distribution = "msnburr",
                          mcmcControl = list(niter = 1500, nburnin = 600),
                          seed = 3))
  pr <- predict(f)
  expect_s3_class(pr, "data.frame")
  expect_true(all(c("x", "density") %in% names(pr)))
  expect_true(all(pr$density >= 0) && all(is.finite(pr$density)))
})
