# Univariate heavy-tail mixture regression (Student-t / Normal-Gamma errors).

test_that("heavy-tail reg specs inherit NormalRegSpec and validate df", {
  expect_true(methods::is(StudentTRegSpec(), "NormalRegSpec"))
  expect_true(methods::is(NormalGammaRegSpec(), "NormalRegSpec"))
  expect_true(isRegressionSpec(StudentTRegSpec()))
  expect_true(isRegressionSpec(NormalRegSpec()))           # v0.4.1 bug fix
  X <- cbind(1, rnorm(20))
  expect_error(defaultPrior(StudentTRegSpec(), rnorm(20),
                            control = list(X = X, df = 2)), "df")
})

test_that("nimixReg routes heavy-tail distributions", {
  d <- data.frame(x = rnorm(20), y = rnorm(20))
  expect_silent(spec <- getDistribution("normal-gamma-reg"))
})

test_that("Student-t and Normal-Gamma regression recover slopes", {
  skip_on_cran(); skip_if_not_installed("nimble")
  set.seed(1); x <- rnorm(200)
  y <- c(2 * x[1:100], -2 * x[101:200]) + c(rt(100, 4), rt(100, 4)) * 0.8
  d <- data.frame(x = x, y = y)
  for (dist in c("studentt", "normalgamma")) {
    fit <- nimixReg(y ~ x, d, K = 2, distribution = dist, method = "fixedk",
                    prior = list(df = 4),
                    mcmcControl = list(niter = 3000, nburnin = 1200),
                    seed = 1, verbose = FALSE)
    fit <- relabel(fit); sl <- sort(fit@relabeled$summary$x)
    expect_true(sl[1] < -1 && sl[2] > 1, info = dist)
  }
})
