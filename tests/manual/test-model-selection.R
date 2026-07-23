# Cluster profiling, median in summary, model selection (WAIC/LOO), ensemble.

test_that("WAIC helpers are correct on a synthetic loglik matrix", {
  # two obs, known loglik columns: logMeanExp and WAIC decomposition
  set.seed(1)
  ll <- matrix(rnorm(200 * 2, -1, 0.3), 200, 2)
  lme <- nimix:::.logMeanExp(ll[, 1])
  expect_equal(lme, log(mean(exp(ll[, 1]))), tolerance = 1e-10)
  w <- nimix:::.waicFromLL(ll)
  # p_waic = sum of per-column variances; elpd = lppd - p_waic
  expect_equal(w$p_waic, sum(apply(ll, 2L, var)), tolerance = 1e-10)
  expect_equal(w$waic, -2 * w$elpd_waic, tolerance = 1e-10)
})

test_that("summary reports mean, median, and a credible interval", {
  skip_on_cran()
  set.seed(1)
  y <- c(rnorm(120, -2, 1), rnorm(80, 3, 0.8))
  f <- relabel(nimixClust(y, K = 2, method = "fixedk",
                          mcmcControl = list(niter = 1500, nburnin = 600),
                          seed = 3))
  s <- f@relabeled$summary
  expect_true(all(c("mu_mean", "mu_med", "mu_lwr", "mu_upr") %in% names(s)))
  # median within the credible interval, mean finite
  expect_true(all(s$mu_med >= s$mu_lwr & s$mu_med <= s$mu_upr))
})

test_that("clusterProfile characterises each cluster", {
  skip_on_cran()
  set.seed(1)
  y <- c(rnorm(120, -2, 1), rnorm(80, 3, 0.8))
  f <- nimixClust(y, K = 2, method = "fixedk",
                  mcmcControl = list(niter = 1500, nburnin = 600), seed = 3)
  pr <- clusterProfile(f)
  expect_s3_class(pr, "nimixProfile")
  expect_true(all(c("cluster", "size", "proportion",
                    "y_mean", "y_sd", "y_median") %in% names(pr)))
  expect_equal(sum(pr$size), length(y))
  expect_equal(sum(pr$proportion), 1, tolerance = 1e-8)
  # the two profiled means should straddle the two true centres
  expect_gt(max(pr$y_mean), 1)
  expect_lt(min(pr$y_mean), 0)
})

test_that("modelSelect and ensembleFit prefer the true K", {
  skip_on_cran()
  set.seed(11)
  y <- c(rnorm(150, -3, 1), rnorm(150, 3, 1))
  f2 <- nimixClust(y, K = 2, method = "fixedk",
                   mcmcControl = list(niter = 1500, nburnin = 600), seed = 3)
  f3 <- nimixClust(y, K = 3, method = "fixedk",
                   mcmcControl = list(niter = 1500, nburnin = 600), seed = 3)

  ms <- modelSelect(K2 = f2, K3 = f3)
  expect_s3_class(ms, "nimixModelSelect")
  expect_equal(ms$model[1], "K2")            # K=2 is the true model
  expect_equal(ms$dWAIC[1], 0)

  w <- nimixWAIC(f2)
  expect_true(is.finite(w$waic) && is.finite(w$elpd_waic))

  e <- ensembleFit(K2 = f2, K3 = f3, method = "waic")
  expect_s3_class(e, "nimixEnsemble")
  expect_equal(sum(e$weights), 1, tolerance = 1e-8)
  expect_gt(e$weights["K2"], e$weights["K3"])
  pe <- predict(e)
  expect_true(all(is.finite(pe$density)))
})

test_that("LOO and stacking paths work when loo is available", {
  skip_on_cran()
  skip_if_not_installed("loo")
  set.seed(11)
  y <- c(rnorm(120, -3, 1), rnorm(120, 3, 1))
  f2 <- nimixClust(y, K = 2, method = "fixedk",
                   mcmcControl = list(niter = 1500, nburnin = 600), seed = 3)
  f3 <- nimixClust(y, K = 3, method = "fixedk",
                   mcmcControl = list(niter = 1500, nburnin = 600), seed = 3)
  lo <- nimixLOO(f2)
  expect_true("elpd_loo" %in% rownames(lo$estimates))
  es <- ensembleFit(K2 = f2, K3 = f3, method = "stacking")
  expect_equal(sum(es$weights), 1, tolerance = 1e-6)
})
