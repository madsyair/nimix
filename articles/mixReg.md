# Mixtures of linear regressions

[`nimixReg()`](https://madsyair.github.io/nimix/reference/nimixReg.md)
fits a Bayesian mixture of Gaussian linear regressions: each component
is $`y \sim N(x^\top \beta_k, \sigma^2_k)`$ with a conjugate
Normal-Inverse-Gamma cluster prior. The number of regimes can be
inferred (`method = "dpm"`) or fixed (`method = "fixedk"`).

This vignette uses the **Preston curve** – the classic
official-statistics relationship between national income and life
expectancy – on the `wdi2022` dataset
([`?wdi2022`](https://madsyair.github.io/nimix/reference/wdi2022.md)).

> MCMC chunks use `eval = FALSE` (CRAN time limits); printed results are
> from an actual run.

## The data

``` r

library(nimix)
#> Loading required package: nimble
#> nimble version 1.4.2 is loaded.
#> For more information on NIMBLE and a User Manual,
#> please visit https://R-nimble.org.
#> 
#> Attaching package: 'nimble'
#> The following object is masked from 'package:stats':
#> 
#>     simulate
#> The following object is masked from 'package:base':
#> 
#>     declare
data(wdi2022)
df <- data.frame(life = wdi2022$life_exp,
                 lgdp = log(wdi2022$gdp_pc))
plot(df$lgdp, df$life, xlab = "log GDP per capita",
     ylab = "life expectancy", main = "Preston curve, 2022")
```

![](mixReg_files/figure-html/data-1.png)

## Does the Preston curve have latent regimes?

``` r

fit <- nimixReg(life ~ lgdp, df, K_max = 6, method = "dpm",
                mcmcControl = list(niter = 5000, nburnin = 2000),
                seed = 1)
fit <- relabel(fit)
summary(fit)
```

An instructive result: on log income the DPM concentrates on a **single
regime** –

    #> modalK = 1
    #>   (Intercept)  lgdp
    #> 1      36.782 4.067      # life ~ 36.8 + 4.07 x log GDP pc

roughly four extra years of life expectancy per log-dollar of income.
The lesson is that a nonparametric mixture does **not** invent regimes:
when one regression line (with its residual spread) explains the data,
the posterior collapses to it. Estimated K always depends on the assumed
component form – heterogeneity here lives in the residuals, not in
distinct slopes.

Even when two components are *offered* (fixed K = 2), one simply
empties:

``` r

fit2 <- nimixReg(life ~ lgdp, df, K = 2, method = "fixedk",
                 mcmcControl = list(niter = 5000, nburnin = 2000),
                 seed = 1)
summary(relabel(fit2))
#> one occupied component; weight ~ 1 on the Preston line
```

## When regimes do exist

With genuinely regime-structured data the same call recovers them. A
simulated benchmark (two latent groups with slopes +2 / -2, mixed
80/20):

``` r

set.seed(107)
n <- 250; x <- runif(n, -3, 3)
grp <- rep(1:2, c(200, 50))
y <- ifelse(grp == 1, 2 * x, -2 * x) + rnorm(n, 0, 0.7)
fs <- nimixReg(y ~ x, data.frame(y = y, x = x), K_max = 8,
               mcmcControl = list(niter = 5000, nburnin = 2000), seed = 1)
summary(relabel(fs))    # recovers slopes near +2 and -2, weights near 0.8/0.2
```

## Prediction

``` r

newd <- data.frame(lgdp = log(c(1000, 10000, 50000)))
predict(fit, newdata = newd)
```

Heavy-tailed residuals (`distribution = "studentt"` / `"normalgamma"`)
and multivariate responses (`cbind(y1, y2) ~ x`) use the same interface.
