# The Fernandez-Steel families (fssn, fsst, fossep) must share ONE skewness
# convention: the exported `alpha` IS the FS skewness gamma, so
#   P(X > mu) = alpha^2 / (1 + alpha^2),
# i.e. alpha > 1 skews right and alpha = 1 is symmetric. This test locks that
# guarantee across d/p/q/r and the compiled NIMBLE kernels.

test_that("FS families share the alpha == gamma convention (density)", {
  xg <- seq(-5, 5, 0.25)
  fsN <- function(x, g) { s <- ifelse(x < 0, x * g, x / g)
                          (2 / (g + 1 / g)) * dnorm(s) }
  fsT <- function(x, g, nu) { s <- ifelse(x < 0, x * g, x / g)
                              (2 / (g + 1 / g)) * dt(s, nu) }
  for (a in c(0.4, 1, 2.5)) {
    expect_lt(max(abs(dfssn(xg, 0, 1, a) - fsN(xg, a))), 1e-10)
    expect_lt(max(abs(dfsst(xg, 0, 1, a, 7) - fsT(xg, a, 7))), 1e-10)
  }
})

test_that("right-tail mass equals alpha^2/(1+alpha^2) in every FS family", {
  for (a in c(0.4, 0.5, 1, 2, 2.5)) {
    target <- a^2 / (1 + a^2)
    expect_lt(abs((1 - pfssn(0, 0, 1, a)) - target), 1e-8)
    expect_lt(abs((1 - pfsst(0, 0, 1, a, 7)) - target), 1e-8)
    expect_lt(abs((1 - pfossep(0, 0, 1, a, 2)) - target), 1e-8)
  }
})

test_that("quantiles invert the CDF and RNG matches the density", {
  u <- c(.1, .4, .75, .95)
  for (a in c(0.5, 2)) {
    expect_lt(max(abs(pfssn(qfssn(u, 1, 2, a), 1, 2, a) - u)), 1e-8)
    expect_lt(max(abs(pfsst(qfsst(u, 1, 2, a, 7), 1, 2, a, 7) - u)), 1e-8)
    expect_lt(max(abs(pfossep(qfossep(u, 1, 2, a, 2), 1, 2, a, 2) - u)), 1e-8)
  }
  set.seed(2)
  for (a in c(0.5, 2)) {
    target <- a^2 / (1 + a^2)
    expect_lt(abs(mean(rfssn(5e4, 0, 1, a) > 0) - target), 0.01)
    expect_lt(abs(mean(rfsst(5e4, 0, 1, a, 20) > 0) - target), 0.01)
    expect_lt(abs(mean(rfossep(5e4, 0, 1, a, 2) > 0) - target), 0.01)
  }
})

test_that("compiled NIMBLE kernels follow the harmonised convention", {
  skip_on_cran()
  suppressMessages(nimix:::.nimixEnsureMSNBurr())
  cdN <- nimble::compileNimble(get("dFSSN_k", envir = globalenv()))
  cdT <- nimble::compileNimble(get("dFSST_k", envir = globalenv()))
  eN <- 0; eT <- 0
  for (a in c(0.4, 1, 2.5)) for (x in c(-8, -1, 0, 2, 8)) {
    eN <- max(eN, abs(cdN(x, 1, 2, a, log = 1) - dfssn(x, 1, 2, a, log = TRUE)))
    eT <- max(eT, abs(cdT(x, 1, 2, a, 6, log = 1) - dfsst(x, 1, 2, a, 6, log = TRUE)))
  }
  expect_lt(eN, 1e-8); expect_lt(eT, 1e-8)
})
