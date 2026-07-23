# Recovery test: simulate data from a
# mixture with KNOWN K and parameters, then check that the posterior of K
# concentrates near the truth and that RELABELLED component estimates recover
# the true parameters. Label switching makes raw summaries meaningless, so we
# must relabel first.

test_that("recovery: two well-separated Gaussians (post-relabel)", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  skip_if_not_installed("label.switching")

  trueMu <- c(-5, 5); trueSd <- c(1, 1); nEach <- 120
  recovered_ok <- logical(3)
  modalK <- integer(3)

  for (s in 1:3) {                       # >= 3 seeds
    set.seed(100 + s)
    y <- c(rnorm(nEach, trueMu[1], trueSd[1]),
           rnorm(nEach, trueMu[2], trueSd[2]))
    fit <- nimixClust(y, K_max = 8,
                       mcmcControl = list(niter = 3000, nburnin = 1000),
                       seed = s, verbose = FALSE)
    fit <- relabel(fit)
    modalK[s] <- fit@relabeled$modalK
    est <- sort(fit@relabeled$summary$mu_mean)
    # relabelled means should be close to the true (sorted) means
    recovered_ok[s] <- isTRUE(all(abs(est - sort(trueMu)) < 1.0)) &&
                       fit@relabeled$modalK == 2L
  }

  # modal K should be 2 in the clear majority of seeds
  expect_gte(sum(modalK == 2L), 2L)
  expect_gte(sum(recovered_ok), 2L)
})

test_that("relabel conditions on modal K and returns aligned summaries", {
  skip_on_cran()
  skip_if_not_installed("nimble")

  set.seed(11)
  y <- c(rnorm(80, -4, 1), rnorm(80, 4, 1))
  fit <- nimixClust(y, K_max = 6,
                     mcmcControl = list(niter = 1500, nburnin = 500),
                     seed = 11, verbose = FALSE)
  fit <- relabel(fit)
  expect_equal(nrow(fit@relabeled$summary), fit@relabeled$modalK)
  expect_true(all(fit@relabeled$summary$s2_mean > 0))
  # mixing weights sum to ~1
  expect_equal(sum(fit@relabeled$summary$weight), 1, tolerance = 0.05)
})
