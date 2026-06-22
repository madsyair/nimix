## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", eval = FALSE)

## -----------------------------------------------------------------------------
# mu <- 1.5; s2 <- 2.3; df <- 5
# mix <- function(y) integrate(function(w)
#   dnorm(y, mu, sqrt(s2 / w)) * dgamma(w, df / 2, df / 2), 0, Inf)$value
# tdn <- function(y) dt((y - mu) / sqrt(s2), df) / sqrt(s2)
# ys <- c(-2, 0, 1.5, 4)
# data.frame(y = ys, scale_mixture = sapply(ys, mix), student_t = sapply(ys, tdn))
# # the two columns agree to numerical-integration error

## -----------------------------------------------------------------------------
# library(nimix)
# set.seed(1)
# y <- c(rt(120, df = 4) - 5, rt(120, df = 4) + 5)   # two heavy-tailed clusters
# 
# fit_t <- nimixClust(y, K_max = 8, distribution = "studentt", prior = list(df = 4),
#                     mcmcControl = list(niter = 4000, nburnin = 1000),
#                     verbose = FALSE)
# fit_ng <- nimixClust(y, K_max = 8, distribution = "normalgamma", prior = list(df = 4),
#                      mcmcControl = list(niter = 4000, nburnin = 1000),
#                      verbose = FALSE)
# summary(fit_t)
# summary(fit_ng)

