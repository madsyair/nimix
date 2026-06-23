## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", eval = FALSE)

## -----------------------------------------------------------------------------
# library(nimix)
# 
# set.seed(1)
# n <- 250
# x <- runif(n, -3, 3)
# grp <- c(rep(1L, 200), rep(2L, 50)) # 80% / 20%
# y <- ifelse(grp == 1L, 2 * x, -2 * x) + rnorm(n, 0, 0.7)
# df <- data.frame(y = y, x = x)

## -----------------------------------------------------------------------------
# fit <- nimixReg(
#  y ~ x, data = df,
#  K_max = 8, method = "dpm",
#  mcmcControl = list(niter = 4000, nburnin = 1000),
#  verbose = FALSE
# )
# summary(fit)

## -----------------------------------------------------------------------------
# predict(fit, newdata = data.frame(x = c(-2, 0, 2)))
# 
# plot(fit, type = "fitted") # observed vs fitted

## -----------------------------------------------------------------------------
# fit2 <- nimixReg(
#  y ~ x, data = df,
#  K = 2, method = "fixedk",
#  mcmcControl = list(niter = 4000, nburnin = 1000),
#  verbose = FALSE
# )
# summary(fit2)

