# Bayesian mixture of linear regressions

Fit a mixture of Gaussian linear regressions. Each component is \\y \sim
N(x^\top \beta_k, \sigma^2_k)\\ with a conjugate Normal-Inverse-Gamma
cluster prior. The number of components can be inferred with a Dirichlet
Process Mixture (`method = "dpm"`, using `K_max`) or fixed
(`method = "fixedk"`, using `K`). Gating is constant: the mixing weights
do not depend on the covariates.

## Usage

``` r
nimixReg(
  formula,
  data,
  K = NULL,
  K_max = NULL,
  distribution = "normal",
  method = c("dpm", "fixedk", "rjmcmc"),
  gating = c("constant", "covariate"),
  prior = list(),
  mcmcControl = list(),
  initMethod = c("kmeans", "single"),
  seed = 1L,
  verbose = TRUE
)
```

## Arguments

- formula:

  A model formula, e.g. `y ~ x1 + x2`.

- data:

  A data frame containing the formula variables.

- K:

  Integer number of components for the finite mixture
  (`method = "fixedk"`). Required for that method.

- K_max:

  Integer truncation level for the DPM (`method = "dpm"`); a data-aware
  default of `min(10, floor(n / 5))` is used when missing.

- distribution:

  Component distribution. Currently `"normal"` (Gaussian linear
  component). Other GLM families are planned.

- method:

  Engine: `"dpm"` (default; estimate the number of components),
  `"fixedk"` (fixed `K`), or `"rjmcmc"` (planned).

- gating:

  Mixing-weight model: `"constant"` (default; weights do not depend on
  covariates). `"covariate"` (concomitant gating) is a planned opt-in
  and currently errors.

- prior:

  A named list of prior overrides passed to
  [`defaultPrior`](https://madsyair.github.io/nimix/reference/defaultPrior.md)
  (e.g. `g` for the g-prior factor, `nu0` for the InvGamma shape) plus,
  for the DPM, optional `concPrior = c(shape, rate)`, or, for the finite
  mixture, `dirichletConc`.

- mcmcControl:

  A named list with `niter`, `nburnin`, `thin`.

- initMethod:

  Initialisation: `"kmeans"` (default) or `"single"`.

- seed:

  Integer RNG seed.

- verbose:

  Logical; print NIMBLE configuration and progress. When `FALSE`,
  NIMBLE's compilation notes are silenced.

## Value

A
[`FitResult`](https://madsyair.github.io/nimix/reference/FitResult-class.md).
[`summary()`](https://rdrr.io/r/base/summary.html) reports relabelled
per-component regression coefficients and residual variances;
`predict(fit, newdata)` returns the posterior predictive mean;
`plot(fit, type = "fitted")` shows observed vs fitted.

## References

Hurn, M., Justel, A., & Robert, C.P. (2003). Estimating mixtures of
regressions. *JCGS*, 12(1), 55–79.
[doi:10.1198/1061860031329](https://doi.org/10.1198/1061860031329)

Grün, B., & Leisch, F. (2008). FlexMix version 2. *JSS*, 28(4), 1–35.
[doi:10.18637/jss.v028.i04](https://doi.org/10.18637/jss.v028.i04)

## Examples

``` r
# \donttest{
set.seed(1)
x <- runif(200, -3, 3)
grp <- rep(1:2, each = 100)
y <- ifelse(grp == 1, 2 * x, -2 * x) + rnorm(200, 0, 0.7)
df <- data.frame(y = y, x = x)

## number of regimes estimated (DPM)
fit <- nimixReg(y ~ x, df, K_max = 8,
                mcmcControl = list(niter = 2000, nburnin = 1000),
                verbose = FALSE)
summary(fit)
#> Relabelling MCMC output before summarising (label switching)...
#> nimix mixture summary (engine: dpm, distribution: normal-reg)
#> Observations: 200 (dimension d = 1)
#> Relabelling: ECR-ITERATIVE-1 conditioned on modal K = 2 (846 draws)
#> 
#> Posterior of number of occupied clusters:
#> 
#>     2     3     4     5 
#> 0.846 0.137 0.015 0.002 
#> 
#> Relabelled component estimates (posterior mean; CIs for univariate):
#>  component weight (Intercept)     x s2_mean s2_lwr s2_upr
#>          1  0.499    -0.07315  1.97    1.01  0.773   1.31
#>          2  0.501     0.00391 -1.97    1.08  0.825   1.39
#> 
#> Mixing diagnostic (single chain): ESS(alpha) = 736, ESS(#clusters) = 398
#> Note: cross-chain Rhat requires multiple chains (planned v0.9.0).
predict(fit, newdata = data.frame(x = c(-2, 0, 2)))
#>    x     .fitted
#> 1 -2 -0.01956009
#> 2  0 -0.03470138
#> 3  2 -0.04984267

## fixed number of regimes (finite mixture)
fit2 <- nimixReg(y ~ x, df, K = 2, method = "fixedk",
                 mcmcControl = list(niter = 2000, nburnin = 1000),
                 verbose = FALSE)
summary(fit2)
#> Relabelling MCMC output before summarising (label switching)...
#> Warning: Low effective sample size for the cluster count (ESS = 0 of 1000 draws): the chain may be mixing poorly across partitions. Consider a longer run or k-means initialisation.
#> nimix mixture summary (engine: fixedk, distribution: normal-reg)
#> Observations: 200 (dimension d = 1)
#> Relabelling: ECR-ITERATIVE-1 conditioned on modal K = 2 (1000 draws)
#> 
#> Posterior of number of occupied clusters:
#> 
#> 2 
#> 1 
#> 
#> Relabelled component estimates (posterior mean; CIs for univariate):
#>  component weight (Intercept)     x s2_mean s2_lwr s2_upr
#>          1  0.497    -0.06728  1.97   0.972  0.745   1.26
#>          2  0.503     0.00823 -1.97   1.069  0.786   1.45
#> 
#> Mixing diagnostic (single chain): ESS(alpha) = NA, ESS(#clusters) = 0
#> Note: cross-chain Rhat requires multiple chains (planned v0.9.0).
# }
```
