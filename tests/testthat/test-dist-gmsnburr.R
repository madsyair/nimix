# GMSNBurr generalized neo-normal family (Iriawan 2000; Choir 2020).

test_that("GMSNBurr density is stable, integrates to one, and reduces correctly", {
  for (a in c(0.3, 1, 3)) for (th in c(0.5, 1, 4)) {
    I <- integrate(function(y) dgmsnburr(y, 1, 2, a, th), -Inf, Inf,
                   rel.tol = 1e-9)$value
    expect_lt(abs(I - 1), 1e-5)
  }
  ext <- c(-1000, -100, 100, 1000)
  expect_true(all(is.finite(dgmsnburr(ext, 1, 2, 2, 3, log = TRUE))))
  xg <- seq(-6, 6, 0.5)
  # theta = 1 -> MSNBurr; alpha = 1 -> MSNBurr-IIa
  expect_lt(max(abs(dgmsnburr(xg, 0, 1, 2, 1) - dmsnburr(xg, 0, 1, 2))), 1e-10)
  expect_lt(max(abs(dgmsnburr(xg, 0, 1, 1, 2) - dmsnburr2a(xg, 0, 1, 2))), 1e-10)
})

test_that("GMSNBurr NIMBLE density equals the R reference", {
  skip_on_cran()
  suppressMessages(nimix:::.nimixEnsureMSNBurr())
  cdG <- nimble::compileNimble(get("dGMSNBurr_k", envir = globalenv()))
  grid <- expand.grid(x = c(-8, -1, 0, 0.7, 8), a = c(0.5, 1, 3),
                      th = c(0.5, 1, 3))
  v <- mapply(function(x, a, th) cdG(x, 1, 2, a, th, log = 1),
              grid$x, grid$a, grid$th)
  expect_lt(max(abs(v - dgmsnburr(grid$x, 1, 2, grid$a, grid$th, log = TRUE))),
            1e-8)
  expect_true(all(is.finite(v)))
})

test_that("GMSNBurr mixtures recover skew structure across engines", {
  skip_on_cran()
  set.seed(21)
  y <- c(rgmsnburr(150, -3, 0.9, 1, 3), rgmsnburr(150, 3, 0.9, 3, 1))
  zt <- rep(1:2, each = 150)
  f <- relabel(nimixClust(y, K = 2, method = "fixedk",
                          distribution = "gmsnburr",
                          mcmcControl = list(niter = 3000, nburnin = 1200),
                          seed = 3))
  s <- f@relabeled$summary; o <- order(s$mu_mean)
  expect_lt(abs(s$mu_mean[o][1] - (-3)), 0.7)
  expect_lt(abs(s$mu_mean[o][2] - 3), 0.7)
  z <- apply(f@clusterAllocation, 2L,
             function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(z == zt), mean(z == (3L - zt))), 0.9)

  fd <- relabel(nimixClust(y, method = "dpm", distribution = "gmsnburr",
                           K_max = 6,
                           mcmcControl = list(niter = 2000, nburnin = 800),
                           seed = 3))
  expect_equal(fd@relabeled$modalK, 2L)
})
