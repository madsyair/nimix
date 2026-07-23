test_that("MSNBurr regression recovers a random intercept", {
  skip_on_cran()
  # Neo-normal random effects (F10 prototype). The random intercept enters the
  # location just as it does for a Gaussian -- neo-normal is location-scale --
  # with the same sum-to-zero parameterisation. Only the emission line differs.
  set.seed(5)
  n <- 240L; x <- rnorm(n); grp <- rep(1:6, length.out = n)
  b_grp <- rnorm(6, 0, 0.9); b_grp <- b_grp - mean(b_grp)
  y <- 1 + 0.5 * x + b_grp[grp] + nimix:::rmsnburr(n, 0, 0.5, 2)
  f <- nimixReg(y ~ x, data.frame(y = y, x = x, grp = factor(grp)), K = 1,
                method = "fixedk", distribution = "msnburr", random = ~ grp,
                mcmcControl = list(niter = 1000, nburnin = 400), seed = 1)
  sm <- f@mcmcSamples
  bhat <- colMeans(sm[, grep("^b\\[", colnames(sm)), drop = FALSE])
  expect_length(bhat, 6L)
  expect_gt(cor(bhat, b_grp), 0.9)
})

test_that("MSNBurr regression recovers a random intercept AND slope", {
  skip_on_cran()
  set.seed(7)
  n <- 300L; x <- rnorm(n); grp <- rep(1:6, length.out = n)
  b_grp <- rnorm(6, 0, 0.8); b_grp <- b_grp - mean(b_grp)
  s_grp <- rnorm(6, 0, 0.5); s_grp <- s_grp - mean(s_grp)
  y <- 1 + 0.5 * x + b_grp[grp] + s_grp[grp] * x + nimix:::rmsnburr(n, 0, 0.5, 2)
  f <- nimixReg(y ~ x, data.frame(y = y, x = x, grp = factor(grp)), K = 1,
                method = "fixedk", distribution = "msnburr", random = ~ x | grp,
                mcmcControl = list(niter = 1000, nburnin = 450), seed = 1)
  sm <- f@mcmcSamples
  bhat <- colMeans(sm[, grep("^b\\[", colnames(sm)), drop = FALSE])
  shat <- colMeans(sm[, grep("^sRE\\[", colnames(sm)), drop = FALSE])
  expect_gt(cor(bhat, b_grp), 0.9)
  expect_gt(cor(shat, s_grp), 0.9)
})

test_that("random intercept generalises across all neo-normal families", {
  skip_on_cran()
  # Every family gains RE through the same generic path (.neoRegFixedKCode gets
  # re = TRUE, .neoRegREConstants appends the group constants). Covers 2- and
  # 3-shape families and the log-normal / truncated-Gamma shape priors.
  set.seed(9)
  n <- 200L; x <- rnorm(n); grp <- rep(1:6, length.out = n)
  bg <- rnorm(6, 0, 0.8); bg <- bg - mean(bg)
  gens <- list(
    msnburr2a = function(m) nimix:::rmsnburr2a(m, 0, 0.5, 2),
    sep       = function(m) nimix:::rsep(m, 0, 0.5, 2),
    gmsnburr  = function(m) nimix:::rgmsnburr(m, 0, 0.5, 2, 1.5),
    fsst      = function(m) nimix:::rfsst(m, 0, 0.5, 1.5, 5)
  )
  for (nm in names(gens)) {
    y <- 1 + 0.5 * x + bg[grp] + gens[[nm]](n)
    f <- nimixReg(y ~ x, data.frame(y = y, x = x, grp = factor(grp)), K = 1,
                  method = "fixedk", distribution = nm, random = ~ grp,
                  mcmcControl = list(niter = 700, nburnin = 300), seed = 1)
    sm <- f@mcmcSamples
    bhat <- colMeans(sm[, grep("^b\\[", colnames(sm)), drop = FALSE])
    expect_gt(cor(bhat, bg), 0.9)
  }
})

test_that("random effects work under the DPM engine for neo-normal families", {
  skip_on_cran()
  # DPM allocates observations to a random number of clusters, but the external
  # grouping factor is fixed; both act on iid exchangeable observations, so the
  # random effect sits cleanly on the location just as under fixed-K.
  set.seed(17)
  n <- 220L; x <- rnorm(n); grp <- rep(1:6, length.out = n)
  bg <- rnorm(6, 0, 0.7); bg <- bg - mean(bg)
  for (nm in c("msnburr", "sep", "gmsnburr")) {
    gen <- switch(nm,
      msnburr  = nimix:::rmsnburr(n, 0, 0.5, 2),
      sep      = nimix:::rsep(n, 0, 0.5, 2),
      gmsnburr = nimix:::rgmsnburr(n, 0, 0.5, 2, 1.5))
    y <- 1 + 0.5 * x + bg[grp] + gen
    f <- nimixReg(y ~ x, data.frame(y = y, x = x, grp = factor(grp)),
                  method = "dpm", distribution = nm, random = ~ grp,
                  mcmcControl = list(niter = 800, nburnin = 350), seed = 1)
    sm <- f@mcmcSamples
    bhat <- colMeans(sm[, grep("^b\\[", colnames(sm)), drop = FALSE])
    expect_gt(cor(bhat, bg), 0.9)
  }
})

test_that("random effects are refused under the HMM engine (no exchangeable grouping)", {
  skip_on_cran()
  # A single HMM time series has no exchangeable grouping, so a random effect is
  # unidentified and confounds with the regime transitions. The guard refuses it
  # with a message pointing at panel-HMM as the prerequisite.
  set.seed(1)
  n <- 60L; x <- rnorm(n); grp <- rep(1:3, length.out = n); y <- x + rnorm(n)
  expect_error(
    nimixReg(y ~ x, data.frame(y = y, x = x, grp = factor(grp)), K = 2,
             method = "hmm", distribution = "msnburr", random = ~ grp,
             mcmcControl = list(niter = 50, nburnin = 20), seed = 1),
    regexp = "HMM is excluded|panel")
})
