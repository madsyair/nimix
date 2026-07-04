# Multivariate mixture clustering with the DPM engine

This vignette demonstrates **multivariate** Gaussian mixture clustering
on the DPM engine. From the user’s point of view the only change from
the univariate case is that `data` is a numeric **matrix** (one row per
observation, one column per dimension). The example clusters countries
in the `wdi2022` official-statistics dataset
([`?wdi2022`](https://madsyair.github.io/nimix/reference/wdi2022.md)) on
two development indicators jointly.

> MCMC chunks use `eval = FALSE` (CRAN time limits); printed results are
> from an actual run.

## The data: income and longevity jointly

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
Y <- cbind(log_gdp  = log(wdi2022$gdp_pc),
           life_exp = wdi2022$life_exp)
plot(Y, main = "207 countries, 2022", xlab = "log GDP per capita",
     ylab = "life expectancy (years)")
```

![](mixClust-multivariate_files/figure-html/data-1.png)

## The model

Each component is a multivariate normal with a conjugate
Normal-Inverse-Wishart cluster base measure (`dmnorm(cov=)` +
`dinvwish`), data-scaled by default. NIMBLE’s native conjugate CRP
samplers are used – nimix deliberately does not reimplement them.

## Fit

``` r

fit <- nimixClust(Y, K_max = 8,
                  mcmcControl = list(niter = 4000, nburnin = 1500),
                  seed = 1)
fit <- relabel(fit)
summary(fit)
```

On the 2022 data the dominant structure is a clearly separated
**high-income / high-longevity** cluster plus overlapping
developing-country clusters:

    #> modalK = 4
    #>   component weight   mu_1   mu_2      # (log GDP pc, life expectancy)
    #> 3         3  0.309 10.363 80.302     # rich, long-lived cluster
    #> 4         4  0.444  8.087 69.402     # large developing cluster
    #> 1         1  0.204  8.570 71.788     # overlapping middle cluster
    #> 2         2  0.043  8.556 69.914     # small transient component

Overlapping components with small weights are expected in a DPM; report
the dominant components and, when a hard partition is needed, use the
MAP allocation:

``` r

alloc <- apply(fit@clusterAllocation, 2L, function(z)
  as.integer(names(which.max(table(z)))))
plot(Y, col = alloc, pch = 19,
     main = "MAP cluster allocation")
```

## Fixed-K comparison

``` r

fit3 <- nimixClust(Y, K = 3, method = "fixedk",
                   mcmcControl = list(niter = 4000, nburnin = 1500),
                   seed = 1)
summary(relabel(fit3))
```
