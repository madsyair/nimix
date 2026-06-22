# Mixtures of linear regressions

This vignette demonstrates `nimixReg`: a Bayesian mixture of Gaussian
linear regressions. Each component is a linear model
$`y \sim N(x^\top \beta_k, \sigma^2_k)`$ with a conjugate
Normal-Inverse-Gamma cluster prior. The number of components can be
inferred with a Dirichlet Process Mixture (`method = "dpm"`), or fixed
(`method = "fixedk"`).

> Chunks use `eval = FALSE` so the vignette builds quickly; run them
> interactively. Pass `verbose = FALSE` to silence NIMBLE’s compilation
> notes.

## A two-regime example

We simulate two regression regimes that share an intercept but have
opposite slopes, with an unbalanced 80/20 split between them.

``` r

library(nimix)

set.seed(1)
n <- 250
x <- runif(n, -3, 3)
grp <- c(rep(1L, 200), rep(2L, 50)) # 80% / 20%
y <- ifelse(grp == 1L, 2 * x, -2 * x) + rnorm(n, 0, 0.7)
df <- data.frame(y = y, x = x)
```

## Estimating the number of regimes (DPM)

``` r

fit <- nimixReg(
 y ~ x, data = df,
 K_max = 8, method = "dpm",
 mcmcControl = list(niter = 4000, nburnin = 1000),
 verbose = FALSE
)
summary(fit)
```

`summary` relabels the components first (the mixture likelihood is
invariant to permuting the component labels, so raw per-component
averages are not meaningful) and then reports, per component, the mixing
weight, the regression coefficients named after the design columns, and
the residual variance.

## Predictions

`predict` returns the posterior predictive mean $`E[y \mid x]`$,
averaging the component linear predictors with the posterior mixing
weights.

``` r

predict(fit, newdata = data.frame(x = c(-2, 0, 2)))

plot(fit, type = "fitted") # observed vs fitted
```

## Fixing the number of regimes

When the number of regimes is known or assumed, the finite-mixture
engine is a fast baseline. It uses a Dirichlet prior on the weights and
a categorical allocation, so there is no truncation to tune.

``` r

fit2 <- nimixReg(
 y ~ x, data = df,
 K = 2, method = "fixedk",
 mcmcControl = list(niter = 4000, nburnin = 1000),
 verbose = FALSE
)
summary(fit2)
```

## Notes

- Gating is constant: the mixing weights do not depend on the
  covariates. Covariate-dependent (concomitant) gating is a planned
  option.
- The residual-variance prior is data-scaled with a finite prior
  variance, which keeps the variance of a small component from
  collapsing toward zero.
- As with clustering, always read component estimates from `summary`
  (which relabels) rather than from the raw MCMC draws.
