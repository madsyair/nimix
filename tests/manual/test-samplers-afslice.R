# AF_slice default for correlated-parameter univariate families.
# Measured motivation (this container, n = 300, K = 2, 2500 iters):
# fssn min ESS 12 -> 622 (x52), gmsnburr 46 -> 417; RW_block made fssn WORSE
# (12 -> 7), so the block-sampler escalation ladder mattered.

test_that("AF_slice blocks replace per-parameter samplers on fssn (FixedK)", {
  skip_on_cran()
  suppressMessages(nimix:::.nimixEnsureMSNBurr())   # kernels are lazily defined
  set.seed(9)
  y <- c(rfssn(80, -4, 1, 0.5), rfssn(80, 4, 1, 2))
  sp <- getDistribution("fssn")
  pr <- defaultPrior(sp, y)
  mc <- buildModelCode(sp, new("FixedKEngine", dirichletConc = 1), n = 160, L = 2)
  cn <- nimix:::buildConstants(sp, pr, 160); cn$K <- 2; cn$alphaVec <- rep(1, 2)
  ini <- nimix:::componentInits(sp, pr, y, 2)
  m <- suppressMessages(nimble::nimbleModel(
    mc$code, constants = cn, data = list(y = y),
    inits = c(list(z = ini$alloc, weights = c(.5, .5)), ini$params),
    calculate = FALSE))
  conf <- suppressMessages(nimble::configureMCMC(m, print = FALSE))
  customizeSamplers(sp, conf, m)
  nm <- vapply(conf$getSamplers(), function(s) s$name, character(1))
  tg <- lapply(conf$getSamplers(), function(s) s$target)
  onMu <- nm[vapply(tg, function(t) any(grepl("muTilde", t)), logical(1))]
  expect_true(all(grepl("AF_slice", onMu)))
  # each block groups the component's mu, sigma, alpha together
  blk <- tg[[which(vapply(tg, function(t) any(grepl("muTilde\\[1\\]", t)),
                          logical(1)))[1]]]
  expect_setequal(blk, c("muTilde[1]", "sigmaTilde[1]", "alphaTilde[1]"))
})

test_that("fssn mixing improves by an order of magnitude and recovery holds", {
  skip_on_cran()
  set.seed(9)
  y <- c(rfssn(150, -4, 1, 0.5), rfssn(150, 4, 1, 2))
  f <- relabel(nimixClust(y, K = 2, method = "fixedk", distribution = "fssn",
                          mcmcControl = list(niter = 2000, nburnin = 800),
                          seed = 1))
  S <- f@mcmcSamples
  ess <- coda::effectiveSize(coda::as.mcmc(
    S[, grep("Tilde", colnames(S)), drop = FALSE]))
  # baseline min ESS was 12/1500 draws; require a clear improvement, with
  # slack for MC noise (measured 622)
  expect_gt(min(ess), 150)
  s <- f@relabeled$summary; o <- order(s$mu_mean)
  expect_lt(max(abs(s$mu_mean[o] - c(-4, 4))), 0.8)
  expect_lt(max(abs(sort(s$alpha_mean) - c(0.5, 2))), 0.5)
})

test_that("AF_slice defaults do not break the DPM path", {
  skip_on_cran()
  set.seed(9)
  y <- c(rfssn(100, -4, 1, 0.5), rfssn(100, 4, 1, 2))
  f <- nimixClust(y, method = "dpm", K_max = 8,
                  mcmcControl = list(niter = 1000, nburnin = 400), seed = 1)
  Kp <- apply(f@clusterAllocation, 1L, function(v) length(unique(v)))
  expect_equal(as.integer(names(sort(table(Kp), decreasing = TRUE))[1]), 2L)
})
