## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", eval = FALSE)

## -----------------------------------------------------------------------------
# library(nimix)
# 
# set.seed(1)
# y <- c(rnorm(150, -4, 1), rnorm(150, 4, 1)) # two clusters, K_true = 2
# 
# fit <- nimixClust(
#  y, K_max = 10, distribution = "normal", method = "dpm",
#  mcmcControl = list(niter = 6000, nburnin = 2000)
# )
# fit

## -----------------------------------------------------------------------------
# summary(fit)

## -----------------------------------------------------------------------------
# plot(fit, type = "trace_raw") # zig-zags between levels = switching
# plot(fit, type = "trace_relabeled") # stable bands after relabelling

## -----------------------------------------------------------------------------
# plot(fit, type = "K")

## -----------------------------------------------------------------------------
# plot(fit, type = "density")
# pp <- predict(fit, newdata = seq(-8, 8, length.out = 200))
# head(pp)

