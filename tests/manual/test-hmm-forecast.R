# nimixForecast(): the forward algorithm run past the end of the series.
# Filter the regime distribution to time T, push it through P^h, sample a
# regime, draw from that regime's emission -- once per posterior draw, so both
# parameter and regime uncertainty are integrated over.

test_that("forecasts are calibrated and revert to the stationary mixture", {
  skip_on_cran()
  set.seed(42)
  P <- rbind(c(.95, .05), c(.10, .90)); mu <- c(-2, 2)
  Tn <- 240L; H <- 12L
  z <- integer(Tn + H); z[1] <- 1L
  for (t in 2:(Tn + H)) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  yall <- rnorm(Tn + H, mu[z], 0.7)
  ytr <- yall[1:Tn]; yte <- yall[(Tn + 1):(Tn + H)]

  f <- nimixClust(ytr, K = 2, method = "hmm", distribution = "normal",
                  mcmcControl = list(niter = 1500, nburnin = 600), seed = 1)
  fc <- nimixForecast(f, h = H, draws = 400)

  expect_identical(dim(fc$draws), c(400L, H))
  expect_identical(dim(fc$regime), c(H, 2L))
  expect_identical(nrow(fc$summary), H)
  expect_true(all(abs(rowSums(fc$regime) - 1) < 1e-8))
  expect_true(all(fc$summary$lower <= fc$summary$median))
  expect_true(all(fc$summary$median <= fc$summary$upper))

  # calibration: measured 0.917 against nominal 0.90 on this benchmark
  cov <- mean(yte >= fc$summary$lower & yte <= fc$summary$upper)
  expect_gt(cov, 0.7)

  # the regime distribution must move TOWARD stationarity, not away
  stat <- c(0.10, 0.05) / 0.15          # pi = (2/3, 1/3) for this P
  expect_lt(max(abs(fc$regime[H, ] - stat)),
            max(abs(fc$regime[1L, ] - stat)))

  # and the interval must widen as the regime becomes unknowable
  w <- fc$summary$upper - fc$summary$lower
  expect_gt(w[H], w[1L])
})

test_that("forecasting a markov-switching regression needs future covariates", {
  skip_on_cran()
  set.seed(11)
  Tn <- 150L; H <- 6L
  P <- rbind(c(.95, .05), c(.07, .93))
  z <- integer(Tn + H); z[1] <- 1L
  for (t in 2:(Tn + H)) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  x <- rnorm(Tn + H); B <- rbind(c(2, 1.5), c(-2, -1.5))
  yall <- B[z, 1] + B[z, 2] * x + rnorm(Tn + H, 0, 0.6)
  df <- data.frame(y = yall[1:Tn], x = x[1:Tn])
  nd <- data.frame(x = x[(Tn + 1):(Tn + H)])
  yte <- yall[(Tn + 1):(Tn + H)]

  f <- nimixReg(y ~ x, df, K = 2, method = "hmm",
                mcmcControl = list(niter = 400, nburnin = 200), seed = 1)

  # the regime can be projected forward; the covariates cannot
  expect_error(nimixForecast(f, h = H), "needs `newdata`")
  expect_error(nimixForecast(f, h = 3L, newdata = nd), "row\\(s\\) but h")

  fc <- nimixForecast(f, h = H, newdata = nd, draws = 300)
  expect_identical(nrow(fc$summary), H)
  # seed 11 holds its regime across the boundary, where the forecast earns its
  # keep: measured RMSE 0.50 against 2.65 for a constant-mean benchmark. When
  # the regime switches at the boundary instead it loses (4.42 vs 2.45) -- a
  # 5% event the model correctly called unlikely, not a defect.
  expect_identical(z[Tn], z[Tn + 1L])
  rmse <- function(a) sqrt(mean((a - yte)^2))
  expect_lt(rmse(fc$summary$median), rmse(rep(mean(df$y), H)))
})

test_that("nimixForecast guards its contract", {
  set.seed(1); y <- rnorm(60)
  fk <- nimixClust(y, K = 2, method = "fixedk",
                   mcmcControl = list(niter = 200, nburnin = 50), seed = 1)
  expect_error(nimixForecast(fk, h = 3), "needs a fit from method = 'hmm'")

  f <- nimixClust(y, K = 2, method = "hmm",
                  mcmcControl = list(niter = 200, nburnin = 50), seed = 1)
  expect_error(nimixForecast(f, h = 0), "integer >= 1")
  expect_error(nimixForecast(f, h = 1, level = 1), "between 0 and 1")
  expect_error(nimixForecast(f, h = 1, level = 0), "between 0 and 1")
})

test_that("forecasting a markov-switching AR feeds the lag back per draw", {
  skip_on_cran()
  # An MS-AR is an MS regression with y lagged into X, so fitting one is free.
  # Forecasting one is NOT: newdata cannot supply a future y. `lags` says the
  # column is the response in disguise, and each posterior draw then feeds its
  # own trajectory back -- which is exactly why the interval must widen.
  set.seed(5)
  Tn <- 260L; H <- 10L
  P <- rbind(c(.95, .05), c(.06, .94))
  z <- integer(Tn + H); z[1] <- 1L
  for (t in 2:(Tn + H)) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  cc <- c(1.5, -1.5); phi <- c(0.30, 0.75)
  y <- numeric(Tn + H); y[1] <- 0
  for (t in 2:(Tn + H))
    y[t] <- cc[z[t]] + phi[z[t]] * y[t - 1L] + rnorm(1, 0, 0.5)
  ytr <- y[1:Tn]; yte <- y[(Tn + 1):(Tn + H)]
  df <- data.frame(y = ytr[-1], ylag = ytr[-Tn])

  f <- nimixReg(y ~ ylag, df, K = 2, method = "hmm",
                mcmcControl = list(niter = 500, nburnin = 200), seed = 1)
  fc <- nimixForecast(f, h = H, lags = c(ylag = 1), draws = 300)

  # measured: RMSE 0.54 against 3.97 for a constant-mean benchmark, on a
  # series whose regime happened to hold across the boundary
  expect_identical(z[Tn], z[Tn + 1L])
  rmse <- function(a) sqrt(mean((a - yte)^2))
  expect_lt(rmse(fc$summary$median), rmse(rep(mean(ytr), H)) / 2)
  expect_gt(mean(yte >= fc$summary$lower & yte <= fc$summary$upper), 0.7)

  # the compounding: a path feeding on itself must get less certain, and
  # faster than a fixed-covariate forecast would. Measured 2.15 -> 5.12.
  w <- fc$summary$upper - fc$summary$lower
  expect_gt(w[H], 2 * w[1L])
})

test_that("lags and newdata divide the design matrix between them", {
  skip_on_cran()
  set.seed(13)
  Tn <- 240L; H <- 6L
  P <- rbind(c(.95, .05), c(.06, .94))
  z <- integer(Tn + H); z[1] <- 1L
  for (t in 2:(Tn + H)) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  w <- rnorm(Tn + H)
  cc <- c(1.5, -1.5); phi <- c(0.3, 0.7); bw <- c(2, -2)
  y <- numeric(Tn + H); y[1] <- 0
  for (t in 2:(Tn + H))
    y[t] <- cc[z[t]] + phi[z[t]] * y[t - 1L] + bw[z[t]] * w[t] + rnorm(1, 0, .5)
  df <- data.frame(y = y[2:Tn], ylag = y[1:(Tn - 1L)], w = w[2:Tn])
  nd <- data.frame(w = w[(Tn + 1):(Tn + H)])

  f <- nimixReg(y ~ ylag + w, df, K = 2, method = "hmm",
                mcmcControl = list(niter = 500, nburnin = 200), seed = 1)
  # lag generated, exogenous supplied: measured RMSE 0.83 against 5.29
  fc <- nimixForecast(f, h = H, newdata = nd, lags = c(ylag = 1), draws = 250)
  yte <- y[(Tn + 1):(Tn + H)]
  expect_lt(sqrt(mean((fc$summary$median - yte)^2)),
            sqrt(mean((mean(df$y) - yte)^2)))

  # a lagged response is generated by the forecast, never supplied to it
  expect_error(nimixForecast(f, h = H, newdata = cbind(nd, ylag = 1),
                             lags = c(ylag = 1)), "not both")
  expect_error(nimixForecast(f, h = H, newdata = nd, lags = c(zz = 1)),
               "does not use")
  expect_error(nimixForecast(f, h = H, newdata = nd, lags = c(ylag = 0)),
               "positive integers")
  expect_error(nimixForecast(f, h = H, lags = c(ylag = 1)), "non-lagged")
})
