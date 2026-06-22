# Unit tests for the multivariate Gaussian spec (v0.2.0). The prior / validation
# / density tests are pure R and run everywhere; the end-to-end recovery test is
# skipped unless nimble is installed and not on CRAN.

test_that("defaultPrior is data-scaled and well-formed", {
  set.seed(1)
  Y <- rbind(matrix(rnorm(100, -3), ncol = 2),
             matrix(rnorm(100,  3), ncol = 2))
  spec <- NormalMvSpec()
  pr <- defaultPrior(spec, Y)
  expect_equal(pr$d, 2L)
  expect_length(pr$mu0, 2L)
  expect_equal(dim(pr$S0), c(2L, 2L))
  expect_true(pr$df0 > pr$d + 1)              # Section 9.3
  # E[Sigma] under InvWishart(S0, df0) = S0 / (df0 - d - 1) should match cov(Y).
  expect_equal(pr$S0 / (pr$df0 - pr$d - 1), stats::cov(Y),
               tolerance = 1e-8)
  expect_silent(validateParams(spec, pr))
})

test_that("validateParams enforces the dimension and df0 invariants", {
  spec <- NormalMvSpec()
  good <- list(mu0 = c(0, 0), kappa0 = 0.25, df0 = 4, S0 = diag(2), d = 2L)
  expect_silent(validateParams(spec, good))
  # df0 must exceed d + 1 (Section 9.3)
  bad_df <- good; bad_df$df0 <- 3
  expect_error(validateParams(spec, bad_df), "df0")
  # mu0 length must equal d
  bad_mu <- good; bad_mu$mu0 <- c(0, 0, 0)
  expect_error(validateParams(spec, bad_mu), "mu0")
  # S0 must be d x d and positive definite
  bad_R <- good; bad_R$S0 <- matrix(c(1, 2, 2, 1), 2)   # not PD
  expect_error(validateParams(spec, bad_R), "positive definite")
})

test_that("componentDensity matches a direct multivariate-normal computation", {
  spec <- NormalMvSpec()
  f <- componentDensity(spec)
  mu <- c(1, -1); Sig <- matrix(c(2, 0.5, 0.5, 1), 2)
  x <- c(0.3, 0.2)
  # direct formula
  d <- 2
  quad <- t(x - mu) %*% solve(Sig) %*% (x - mu)
  ref <- as.numeric((2 * pi)^(-d / 2) * det(Sig)^(-0.5) * exp(-0.5 * quad))
  expect_equal(f(x, list(mu = mu, Sigma = Sig)), ref, tolerance = 1e-8)
})

test_that("simulateParams returns conformable mu and Sigma", {
  set.seed(2)
  spec <- NormalMvSpec()
  pr <- list(mu0 = c(0, 0), kappa0 = 0.25, df0 = 6, S0 = diag(2), d = 2L)
  sp <- simulateParams(spec, pr, nClust = 3)
  expect_equal(dim(sp$mu), c(3L, 2L))
  expect_equal(dim(sp$Sigma), c(2L, 2L, 3L))
  # each covariance is symmetric positive definite
  for (j in 1:3)
    expect_true(min(eigen(sp$Sigma[, , j], symmetric = TRUE,
                          only.values = TRUE)$values) > 0)
})

test_that("end-to-end multivariate DPM recovers two well-separated clusters", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  set.seed(7)
  n <- 120
  Y <- rbind(matrix(rnorm(n, mean = -3), ncol = 2),
             matrix(rnorm(n, mean =  3), ncol = 2))
  fit <- nimixClust(Y, K_max = 8,
                    mcmcControl = list(niter = 3000, nburnin = 1000),
                    seed = 7, verbose = FALSE)
  modalK <- as.integer(names(sort(table(fit@Kposterior),
                                  decreasing = TRUE))[1])
  expect_true(modalK %in% c(2L, 3L))   # 2 true; allow mild over-clustering
  fit <- relabel(fit)
  expect_equal(nrow(fit@relabeled$summary), modalK)
})
