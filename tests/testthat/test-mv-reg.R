# Multivariate-response mixture regression.

test_that("mv-reg specs are regression specs; heavy-tail inherit Normal mv-reg", {
  expect_true(isRegressionSpec(NormalMvRegSpec()))
  expect_true(methods::is(StudentTMvRegSpec(), "NormalMvRegSpec"))
  expect_true(methods::is(NormalGammaMvRegSpec(), "NormalMvRegSpec"))
  X <- cbind(1, rnorm(20)); Y <- matrix(rnorm(40), 20, 2)
  pr <- defaultPrior(NormalMvRegSpec(), Y, control = list(X = X))
  expect_equal(pr$d, 2L); expect_equal(pr$p, 2L)
  expect_true(pr$df0 > pr$d + 1)
})

test_that("multivariate-response regression recovers coefficient matrices", {
  skip_on_cran(); skip_if_not_installed("nimble")
  set.seed(1); n <- 200; x <- rnorm(n)
  B1 <- rbind(c(1, 0.5), c(3, 1)); B2 <- rbind(c(1, 0.5), c(-3, -1))
  X <- cbind(1, x); E <- matrix(rnorm(n * 2), n, 2) * 0.7
  Y <- rbind(X[1:100, ] %*% B1, X[101:200, ] %*% B2) + E
  dat <- data.frame(y1 = Y[, 1], y2 = Y[, 2], x = x)
  for (dist in c("normal", "studentt", "normalgamma")) {
    fit <- nimixReg(cbind(y1, y2) ~ x, dat, K = 2, distribution = dist,
                    method = "fixedk",
                    prior = if (dist == "normal") list() else list(df = 5),
                    mcmcControl = list(niter = 2500, nburnin = 1000),
                    seed = 1, verbose = FALSE)
    fit <- relabel(fit); sl1 <- sort(fit@relabeled$summary[["x:y1"]])
    expect_true(sl1[1] < -1.5 && sl1[2] > 1.5, info = dist)
  }
})

test_that("predict returns one fitted column per response", {
  skip_on_cran(); skip_if_not_installed("nimble")
  set.seed(2); n <- 120; x <- rnorm(n)
  Y <- cbind(2 * x, -x) + matrix(rnorm(n * 2), n, 2) * 0.5
  dat <- data.frame(y1 = Y[, 1], y2 = Y[, 2], x = x)
  fit <- nimixReg(cbind(y1, y2) ~ x, dat, K = 1, distribution = "normal",
                  method = "fixedk",
                  mcmcControl = list(niter = 800, nburnin = 300),
                  seed = 2, verbose = FALSE)
  pp <- predict(fit, newdata = data.frame(x = c(-1, 0, 1)))
  expect_true(all(c(".fitted.y1", ".fitted.y2") %in% colnames(pp)))
})
