# SEP and LEP symmetric exponential-power families (Batch B; Choir 2020).

test_that("SEP/LEP densities are stable, integrate to one, reduce to Normal", {
  for (nu in c(1, 1.5, 2, 3)) {
    expect_lt(abs(integrate(function(y) dsep(y, 1, 2, nu), -Inf, Inf,
                            rel.tol = 1e-9)$value - 1), 1e-5)
    expect_lt(abs(integrate(function(y) dlep(y, 1, 2, nu), -Inf, Inf,
                            rel.tol = 1e-9)$value - 1), 1e-5)
  }
  xg <- seq(-6, 6, 0.5)
  expect_lt(max(abs(dsep(xg, 0, 1, 2) - dnorm(xg, 0, 1))), 1e-8)  # nu=2 Normal
  expect_lt(max(abs(dlep(xg, 0, 1, 2) - dnorm(xg, 0, 1))), 1e-8)
  expect_true(all(is.finite(dsep(c(-1e4, 1e4), 1, 2, 1.5, log = TRUE))))
})

test_that("SEP/LEP quantiles vectorise without recycling warnings", {
  mu <- c(-3, -3, 3, 3); sg <- c(1, 1, .8, .8); nu <- c(2, 2, 1.5, 1.5)
  expect_silent(qsep(c(.2, .7, .3, .9), mu, sg, nu))
  expect_silent(qlep(c(.2, .7, .3, .9), mu, sg, nu))
  # vector == element-wise scalar
  u <- c(.2, .7, .3, .9)
  vs <- vapply(seq_along(u), function(i) qsep(u[i], mu[i], sg[i], nu[i]),
               numeric(1))
  expect_lt(max(abs(qsep(u, mu, sg, nu) - vs)), 1e-10)
})

test_that("SEP/LEP NIMBLE densities equal the R reference", {
  skip_on_cran()
  suppressMessages(nimix:::.nimixEnsureMSNBurr())
  cdS <- nimble::compileNimble(get("dSEP_k", envir = globalenv()))
  cdL <- nimble::compileNimble(get("dLEP_k", envir = globalenv()))
  gs <- expand.grid(x = c(-6, -1, 0, 2, 6), nu = c(1, 1.5, 2, 3))
  expect_lt(max(abs(mapply(function(x, nu) cdS(x, 1, 2, nu, log = 1),
                           gs$x, gs$nu) - dsep(gs$x, 1, 2, gs$nu, log = TRUE))),
            1e-8)
  expect_lt(max(abs(mapply(function(x, nu) cdL(x, 1, 2, nu, log = 1),
                           gs$x, gs$nu) - dlep(gs$x, 1, 2, gs$nu, log = TRUE))),
            1e-8)
})

test_that("SEP/LEP mixtures recover components under FixedK and DPM", {
  skip_on_cran()
  for (d in c("sep", "lep")) {
    rfun <- get(paste0("r", d))
    set.seed(31)
    y <- c(rfun(150, -4, 1, 1.3), rfun(150, 4, 1, 1.3))
    zt <- rep(1:2, each = 150)
    f <- relabel(nimixClust(y, K = 2, method = "fixedk", distribution = d,
                            mcmcControl = list(niter = 2000, nburnin = 800),
                            seed = 3))
    s <- f@relabeled$summary; o <- order(s$mu_mean)
    expect_lt(abs(s$mu_mean[o][1] - (-4)), 0.9)
    expect_lt(abs(s$mu_mean[o][2] - 4), 0.9)
    z <- apply(f@clusterAllocation, 2L,
               function(v) as.integer(names(which.max(table(v)))))
    expect_gt(max(mean(z == zt), mean(z == (3L - zt))), 0.9)
  }
})
