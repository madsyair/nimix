# Convergence suite (Vehtari et al. 2021) + posterior predictive workflow.

test_that("bulk/tail ESS and folded Rhat match theory on known series", {
  set.seed(1)
  ch <- list(rnorm(1000), rnorm(1000))
  expect_gt(nimix:::.bulkESS(ch), 1500)          # iid: ~ total draws
  expect_gt(nimix:::.tailESS(ch), 1000)
  expect_lt(abs(nimix:::.foldedRhat(ch) - 1), 0.02)

  # scale disagreement: folded Rhat must flag what location Rhat misses
  chS <- list(rnorm(1000, 0, 1), rnorm(1000, 0, 3))
  expect_lt(nimix:::.splitRhat(chS), 1.02)
  expect_gt(nimix:::.foldedRhat(chS), 1.1)

  # strong autocorrelation: bulk ESS far below the draw count
  ar <- function(n, rho) { x <- numeric(n)
    for (i in 2:n) x[i] <- rho * x[i - 1] + rnorm(1); x }
  set.seed(2)
  chA <- list(ar(2000, 0.9), ar(2000, 0.9))
  expect_lt(nimix:::.bulkESS(chA), 600)
})

test_that("allocation entropy is exact on known partitions", {
  zM <- rbind(c(1L, 1L, 2L, 2L), c(1L, 1L, 1L, 1L))
  h <- nimix:::.allocEntropy(zM, 2L)
  expect_equal(h[1], log(2))
  expect_equal(h[2], 0)
})

test_that("ppCheck flags misfit and passes a well-specified model", {
  skip_on_cran()
  set.seed(5)
  y <- c(rnorm(120, -2, 1), rnorm(80, 2, 1))
  f <- nimixClust(y, K = 2, method = "fixedk",
                  mcmcControl = list(niter = 2000, nburnin = 800, nchains = 2),
                  seed = 3)
  # diagnostics carry the Vehtari functional table incl. allocation entropy
  ft <- f@diagnostics$functionals
  expect_true(is.data.frame(ft) && "entropy" %in% ft$functional)
  pc <- ppCheck(f, nrep = 120, statistics = c("mean", "sd", "skew"))
  expect_true(all(pc$ppp > 0.05 & pc$ppp < 0.95))

  # Poisson mixture on overdispersed counts: the tail must be flagged
  set.seed(9); yod <- rnbinom(200, mu = 8, size = 2)
  fp <- nimixClust(yod, K = 2, method = "fixedk", distribution = "poisson",
                   mcmcControl = list(niter = 1500, nburnin = 600), seed = 3)
  pcp <- ppCheck(fp, nrep = 120, statistics = c("max"))
  expect_lt(pcp$ppp[pcp$statistic == "max"], 0.05)
})
