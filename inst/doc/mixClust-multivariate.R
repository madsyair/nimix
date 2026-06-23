## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", eval = FALSE)

## -----------------------------------------------------------------------------
# library(nimix)
# 
# set.seed(1)
# # two well-separated 2-D clusters, K_true = 2
# Y <- rbind(
#  matrix(rnorm(150 * 2, mean = -3), ncol = 2),
#  matrix(rnorm(150 * 2, mean = 3), ncol = 2)
# )
# 
# fit <- nimixClust(
#  Y, K_max = 10, distribution = "normal", method = "dpm",
#  mcmcControl = list(niter = 6000, nburnin = 2000)
# )
# fit

## -----------------------------------------------------------------------------
# summary(fit)

## -----------------------------------------------------------------------------
# plot(fit, type = "K") # posterior of the number of occupied clusters
# plot(fit, type = "cluster") # data coloured by MAP cluster (first two dims)

## -----------------------------------------------------------------------------
# plot(fit, type = "trace_raw")
# plot(fit, type = "trace_relabeled")

## -----------------------------------------------------------------------------
# newpts <- rbind(c(-3, -3), c(0, 0), c(3, 3))
# predict(fit, newdata = newpts)

