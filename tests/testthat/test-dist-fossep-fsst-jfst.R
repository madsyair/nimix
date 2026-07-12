# FOSSEP, FSST and JFST -- the four-parameter Batch B neo-normal families.
# (Choir 2020; Fernandez & Steel 1998, 1999; Jones & Faddy 2003.)

test_that("FOSSEP/FSST/JFST densities integrate to one and stay finite", {
  for (a in c(0.5, 1, 2)) for (th in c(1.5, 2, 3)) {
    expect_lt(abs(integrate(function(y) dfossep(y, 1, 2, a, th), -Inf, Inf,
                            rel.tol = 1e-9)$value - 1), 1e-5)
  }
  for (a in c(0.5, 1, 2)) for (nu in c(3, 6, 15)) {
    expect_lt(abs(integrate(function(y) dfsst(y, 1, 2, a, nu), -Inf, Inf,
                            rel.tol = 1e-9)$value - 1), 1e-5)
  }
  for (a in c(1, 3, 5)) for (th in c(1, 3, 5)) {
    expect_lt(abs(integrate(function(y) djfst(y, 1, 2, a, th), -Inf, Inf,
                            rel.tol = 1e-9)$value - 1), 1e-5)
  }
  expect_true(all(is.finite(dfossep(c(-1e4, 1e4), 1, 2, 2, 3, log = TRUE))))
  expect_true(all(is.finite(dfsst(c(-1e4, 1e4), 1, 2, 2, 6, log = TRUE))))
  expect_true(all(is.finite(djfst(c(-1e4, 0, 1e4), 1, 2, 3, 5, log = TRUE))))
})

test_that("known reductions hold exactly", {
  xg <- seq(-6, 6, 0.5)
  # FSST with alpha = 1 is the (standardised) symmetric Student-t
  expect_lt(max(abs(dfsst(xg, 0, 1, 1, 6) - dt(xg, 6))), 1e-8)
  # JFST with alpha = theta is symmetric about mu
  expect_lt(max(abs(djfst(xg, 0, 1, 4, 4) - djfst(-xg, 0, 1, 4, 4))), 1e-10)
})

test_that("quantiles vectorise without recycling warnings", {
  u <- c(.2, .7, .3, .9)
  mu <- c(-3, -3, 3, 3); sg <- c(1, 1, .8, .8)
  expect_silent(qfossep(u, mu, sg, c(.5, .5, 2, 2), c(2, 2, 3, 3)))
  expect_silent(qfsst(u, mu, sg, c(.5, .5, 2, 2), c(5, 5, 8, 8)))
  expect_silent(qjfst(u, mu, sg, c(3, 3, 5, 5), c(5, 5, 3, 3)))
  # vector call equals element-wise scalar calls
  vs <- vapply(seq_along(u),
               function(i) qfsst(u[i], mu[i], sg[i], c(.5, .5, 2, 2)[i],
                                 c(5, 5, 8, 8)[i]), numeric(1))
  expect_lt(max(abs(qfsst(u, mu, sg, c(.5, .5, 2, 2), c(5, 5, 8, 8)) - vs)),
            1e-10)
  # r-functions with per-observation parameters (the ppCheck pattern)
  expect_silent(rfsst(4, mu, sg, c(.5, .5, 2, 2), c(5, 5, 8, 8)))
  expect_true(all(is.finite(rfsst(4, mu, sg, c(.5, .5, 2, 2), c(5, 5, 8, 8)))))
})

test_that("NIMBLE densities equal the R reference (scalar params)", {
  skip_on_cran()
  suppressMessages(nimix:::.nimixEnsureMSNBurr())
  cdF <- nimble::compileNimble(get("dFOSSEP_k", envir = globalenv()))
  cdS <- nimble::compileNimble(get("dFSST_k", envir = globalenv()))
  cdJ <- nimble::compileNimble(get("dJFST_k", envir = globalenv()))
  eF <- 0; eS <- 0; eJ <- 0
  for (x in c(-8, -1, 0, 2, 8)) {
    for (a in c(0.5, 1, 2)) for (th in c(1.5, 2, 3))
      eF <- max(eF, abs(cdF(x, 1, 2, a, th, log = 1) -
                          dfossep(x, 1, 2, a, th, log = TRUE)))
    for (a in c(0.5, 1, 2)) for (nu in c(3, 6, 15))
      eS <- max(eS, abs(cdS(x, 1, 2, a, nu, log = 1) -
                          dfsst(x, 1, 2, a, nu, log = TRUE)))
    for (a in c(1, 3, 5)) for (th in c(1, 3, 5))
      eJ <- max(eJ, abs(cdJ(x, 1, 2, a, th, log = 1) -
                          djfst(x, 1, 2, a, th, log = TRUE)))
  }
  expect_lt(eF, 1e-8); expect_lt(eS, 1e-8); expect_lt(eJ, 1e-8)
})

test_that("FOSSEP mixture recovers skew components (FixedK + DPM)", {
  skip_on_cran()
  set.seed(51)
  y <- c(rfossep(150, -4, 1, 0.6, 2), rfossep(150, 4, 1, 1.8, 2.5))
  zt <- rep(1:2, each = 150)
  f <- relabel(nimixClust(y, K = 2, method = "fixedk", distribution = "fossep",
                          mcmcControl = list(niter = 2500, nburnin = 1000),
                          seed = 3))
  s <- f@relabeled$summary; o <- order(s$mu_mean)
  expect_lt(abs(s$mu_mean[o][1] - (-4)), 0.9)
  expect_lt(abs(s$mu_mean[o][2] - 4), 0.9)
  expect_lt(s$alpha_mean[o][1], 1)          # left component skewed left
  expect_gt(s$alpha_mean[o][2], 1)          # right component skewed right
  z <- apply(f@clusterAllocation, 2L,
             function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(z == zt), mean(z == (3L - zt))), 0.9)

  fd <- relabel(nimixClust(y, method = "dpm", distribution = "fossep",
                           K_max = 8,
                           mcmcControl = list(niter = 1800, nburnin = 700),
                           seed = 3))
  expect_equal(fd@relabeled$modalK, 2L)
})

test_that("FSST mixture recovers heavy-tailed skew components", {
  skip_on_cran()
  set.seed(61)
  y <- c(rfsst(150, -5, 1, 0.6, 8), rfsst(150, 5, 1, 1.7, 8))
  zt <- rep(1:2, each = 150)
  f <- relabel(nimixClust(y, K = 2, method = "fixedk", distribution = "fsst",
                          mcmcControl = list(niter = 2000, nburnin = 800),
                          seed = 3))
  s <- f@relabeled$summary; o <- order(s$mu_mean)
  expect_lt(abs(s$mu_mean[o][1] - (-5)), 1.0)
  expect_lt(abs(s$mu_mean[o][2] - 5), 1.0)
  # nu is only weakly identified, so we only require it to stay above the
  # truncation bound rather than to recover the simulating value.
  expect_true(all(s$nu_mean > 2))
  z <- apply(f@clusterAllocation, 2L,
             function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(z == zt), mean(z == (3L - zt))), 0.9)

  # nuTilde must be stochastic, else dCRP refuses to cluster it
  fd <- relabel(nimixClust(y, method = "dpm", distribution = "fsst", K_max = 8,
                           mcmcControl = list(niter = 1500, nburnin = 600),
                           seed = 3))
  expect_equal(fd@relabeled$modalK, 2L)
})

test_that("JFST mixture recovers components (FixedK + DPM)", {
  skip_on_cran()
  set.seed(71)
  y <- c(rjfst(150, -5, 1, 2, 6), rjfst(150, 5, 1, 6, 2))
  zt <- rep(1:2, each = 150)
  f <- relabel(nimixClust(y, K = 2, method = "fixedk", distribution = "jfst",
                          mcmcControl = list(niter = 2500, nburnin = 1000),
                          seed = 3))
  s <- f@relabeled$summary; o <- order(s$mu_mean)
  expect_lt(abs(s$mu_mean[o][1] - (-5)), 1.2)
  expect_lt(abs(s$mu_mean[o][2] - 5), 1.2)
  z <- apply(f@clusterAllocation, 2L,
             function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(z == zt), mean(z == (3L - zt))), 0.85)

  fd <- relabel(nimixClust(y, method = "dpm", distribution = "jfst", K_max = 8,
                           mcmcControl = list(niter = 1500, nburnin = 600),
                           seed = 3))
  expect_equal(fd@relabeled$modalK, 2L)
})
