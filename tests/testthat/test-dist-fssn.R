# FSSN Fernandez-Steel skew Normal (Batch B; FS 1998, Choir 2020).

test_that("FSSN density stable, integrates to one, reduces to Normal", {
  for (a in c(0.3, 0.5, 1, 2, 3)) {
    expect_lt(abs(integrate(function(y) dfssn(y, 1, 2, a), -Inf, Inf,
                            rel.tol = 1e-9)$value - 1), 1e-5)
  }
  xg <- seq(-6, 6, 0.5)
  expect_lt(max(abs(dfssn(xg, 0, 1, 1) - dnorm(xg, 0, 1))), 1e-8)  # alpha=1
  expect_true(all(is.finite(dfssn(c(-1e4, 1e4), 1, 2, 2, log = TRUE))))
})

test_that("FSSN quantile vectorises without recycling warnings", {
  expect_silent(qfssn(c(.2, .7, .3, .9), c(-3, -3, 3, 3), c(1, 1, .8, .8),
                      c(.5, .5, 2, 2)))
})

test_that("FSSN NIMBLE density equals the R reference (scalar params)", {
  skip_on_cran()
  suppressMessages(nimix:::.nimixEnsureMSNBurr())
  cdF <- nimble::compileNimble(get("dFSSN_k", envir = globalenv()))
  err <- 0
  for (a in c(0.5, 1, 2, 3)) for (x in c(-8, -1, 0, 2, 8))
    err <- max(err, abs(cdF(x, 1, 2, a, log = 1) - dfssn(x, 1, 2, a, log = TRUE)))
  expect_lt(err, 1e-8)
})

test_that("FSSN mixture recovers skew components (FixedK + DPM)", {
  skip_on_cran()
  set.seed(41)
  y <- c(rfssn(150, -4, 1, 0.5), rfssn(150, 4, 1, 2))
  zt <- rep(1:2, each = 150)
  f <- relabel(nimixClust(y, K = 2, method = "fixedk", distribution = "fssn",
                          mcmcControl = list(niter = 2500, nburnin = 1000),
                          seed = 3))
  s <- f@relabeled$summary; o <- order(s$mu_mean)
  expect_lt(abs(s$mu_mean[o][1] - (-4)), 0.8)
  expect_lt(abs(s$mu_mean[o][2] - 4), 0.8)
  # left component alpha < 1, right component alpha > 1
  expect_lt(s$alpha_mean[o][1], 1)
  expect_gt(s$alpha_mean[o][2], 1)
  z <- apply(f@clusterAllocation, 2L,
             function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(z == zt), mean(z == (3L - zt))), 0.9)
})
