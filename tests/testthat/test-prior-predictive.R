test_that("priorPredictive passes a sane prior and flags a bad one", {
  # No MCMC here: parameters come from simulateParams, data from the prior --
  # cheap enough for the automatic tier.
  set.seed(5)
  y <- c(rnorm(80, -2, .6), rnorm(120, 3, .8))
  pp <- priorPredictive(y, K = 2, distribution = "normal", nsim = 150, seed = 1)
  expect_s3_class(pp, "nimixPriorPred")
  expect_equal(sum(pp$pTail < 0.05), 0)          # data-scaled prior covers data
  # A prior centred at 100 with tiny spread cannot reach data centred near 1
  pb <- priorPredictive(y, K = 2, distribution = "normal", nsim = 150, seed = 1,
                        prior = list(mu0 = 100, kappa0 = 100))
  expect_lt(pb$pTail["mean"], 0.05)              # mean flagged as unreachable
})

test_that("priorPredictive works across component families (nClust contract)", {
  # Regression guard for the simulateParams signature unification: methods
  # must take nClust (the generic's name), or positional dispatch loses it.
  set.seed(5)
  y <- c(rnorm(60, -2, .6), rnorm(60, 3, .8))
  for (d in c("msnburr", "sep", "fsst", "studentt", "gmsnburr")) {
    pp <- priorPredictive(y, K = 2, distribution = d, nsim = 40, seed = 1)
    expect_equal(dim(pp$sim), c(40L, 5L))
    expect_true(all(is.finite(pp$pTail)))
  }
})

test_that("priorPredictive print and plot methods work", {
  set.seed(5)
  y <- rnorm(100)
  pp <- priorPredictive(y, K = 2, distribution = "normal", nsim = 40, seed = 1)
  expect_output(print(pp), "Prior predictive check")
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, 400, 300)
  plot(pp)
  grDevices::dev.off()
  expect_gt(file.info(tmp)$size, 1000)
})
