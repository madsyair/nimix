test_that("Poisson regression supports a random intercept (GLMM, RE in the link)", {
  skip_on_cran()
  # The genuine GLMM case: the random effect enters INSIDE the log link, so the
  # group offset is multiplicative on the count scale. tauRE is on the link
  # scale, so its bounds are fixed (not sd(y)) -- a large count would otherwise
  # give an absurdly wide ceiling.
  set.seed(11)
  n <- 240L; x <- rnorm(n); grp <- rep(1:6, length.out = n)
  b_grp <- rnorm(6, 0, 0.5); b_grp <- b_grp - mean(b_grp)
  y <- rpois(n, exp(0.5 + 0.4 * x + b_grp[grp]))
  f <- nimixReg(y ~ x, data.frame(y = y, x = x, grp = factor(grp)), K = 1,
                method = "fixedk", distribution = "poisson", random = ~ grp,
                mcmcControl = list(niter = 900, nburnin = 400), seed = 1)
  sm <- f@mcmcSamples
  bhat <- colMeans(sm[, grep("^b\\[", colnames(sm)), drop = FALSE])
  expect_length(bhat, 6L)
  expect_gt(cor(bhat, b_grp), 0.85)
})

test_that("Binomial regression supports a random intercept (GLMM, logit link)", {
  skip_on_cran()
  set.seed(13)
  n <- 260L; x <- rnorm(n); grp <- rep(1:6, length.out = n); sz <- 10L
  b_grp <- rnorm(6, 0, 0.6); b_grp <- b_grp - mean(b_grp)
  y <- rbinom(n, sz, plogis(0.2 + 0.5 * x + b_grp[grp]))
  f <- nimixReg(y ~ x, data.frame(y = y, x = x, grp = factor(grp)), K = 1,
                method = "fixedk", distribution = "binomial", random = ~ grp,
                prior = list(size = sz),
                mcmcControl = list(niter = 900, nburnin = 400), seed = 1)
  sm <- f@mcmcSamples
  bhat <- colMeans(sm[, grep("^b\\[", colnames(sm)), drop = FALSE])
  expect_gt(cor(bhat, b_grp), 0.9)
})

test_that("Poisson GLMM supports a random slope too", {
  skip_on_cran()
  set.seed(15)
  n <- 300L; x <- rnorm(n); grp <- rep(1:6, length.out = n)
  bg <- rnorm(6, 0, 0.4); bg <- bg - mean(bg)
  sg <- rnorm(6, 0, 0.3); sg <- sg - mean(sg)
  y <- rpois(n, exp(0.4 + 0.3 * x + bg[grp] + sg[grp] * x))
  f <- nimixReg(y ~ x, data.frame(y = y, x = x, grp = factor(grp)), K = 1,
                method = "fixedk", distribution = "poisson", random = ~ x | grp,
                mcmcControl = list(niter = 1000, nburnin = 450), seed = 1)
  sm <- f@mcmcSamples
  bh <- colMeans(sm[, grep("^b\\[", colnames(sm)), drop = FALSE])
  sh <- colMeans(sm[, grep("^sRE\\[", colnames(sm)), drop = FALSE])
  expect_gt(cor(bh, bg), 0.85)
  expect_gt(cor(sh, sg), 0.85)
})
