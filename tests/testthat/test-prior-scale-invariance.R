# Prior scale invariance across engines.
#
# nimix assigns priors from the data scale, so multiplying the response by a
# constant must multiply the location estimates by the same constant and leave
# everything else alone. This file locks that invariant per engine.
#
# It exists because v1.3.0 shipped a violation that no gate caught: the
# random-effect prior used fixed bounds dunif(0.01, 5), so at y x1000 the
# needed tauRE (771) sat far above the ceiling and cor(b_hat, truth) collapsed
# from 0.992 to 0.091 -- silently. The pre-existing scale-equivariance lock
# missed it because it rescales PREDICTORS only, never the response. The
# lesson generalised: an invariance lock must exercise every direction of the
# invariance it claims to protect.

.simTwoClusters <- function(seed = 7L) {
  set.seed(seed)
  c(rnorm(120, -3, 0.7), rnorm(120, 3, 0.9))
}

.locs <- function(fit, sc) sort(relabel(fit)@relabeled$summary$mu_mean) / sc

test_that("fixedk clustering is exactly scale-equivariant in the response", {
  skip_on_cran()
  y <- .simTwoClusters()
  mc <- list(niter = 1500, nburnin = 600)
  f1 <- nimixClust(y, K = 2, method = "fixedk", mcmcControl = mc, seed = 1,
                   verbose = FALSE)
  f2 <- nimixClust(y * 1000, K = 2, method = "fixedk", mcmcControl = mc,
                   seed = 1, verbose = FALSE)
  # Conjugate samplers on data-scaled priors reproduce the draws exactly, so
  # this is an equality, not an approximation.
  expect_equal(.locs(f1, 1), .locs(f2, 1000), tolerance = 1e-6)
})

test_that("dpm clustering is exactly scale-equivariant in the response", {
  skip_on_cran()
  y <- .simTwoClusters()
  mc <- list(niter = 1500, nburnin = 600)
  f1 <- nimixClust(y, K_max = 6, method = "dpm", mcmcControl = mc, seed = 1,
                   verbose = FALSE)
  f2 <- nimixClust(y * 1000, K_max = 6, method = "dpm", mcmcControl = mc,
                   seed = 1, verbose = FALSE)
  expect_equal(.locs(f1, 1), .locs(f2, 1000), tolerance = 1e-6)
})

test_that("hmm is scale-equivariant up to Monte Carlo error", {
  skip_on_cran()
  # The HMM path uses adaptive RW samplers, whose proposal scale starts at
  # NIMBLE's default of 1 regardless of the data scale. Short chains therefore
  # show a transient adaptation difference (measured 0.08 at 2500 iterations);
  # with adequate burnin the two agree to within MCSE (measured 0.001 against
  # an MCSE of 0.0013 at 8000). This asserts the converged behaviour.
  set.seed(11)
  P <- rbind(c(.95, .05), c(.10, .90)); mu <- c(-2, 2)
  z <- integer(250); z[1] <- 1L
  for (t in 2:250) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  y <- rnorm(250, mu[z], c(.6, .8)[z])
  mc <- list(niter = 6000, nburnin = 2500)
  f1 <- nimixClust(y, K = 2, method = "hmm", mcmcControl = mc, seed = 1,
                   verbose = FALSE)
  f2 <- nimixClust(y * 1000, K = 2, method = "hmm", mcmcControl = mc, seed = 1,
                   verbose = FALSE)
  expect_lt(max(abs(.locs(f1, 1) - .locs(f2, 1000))), 0.05)
})
