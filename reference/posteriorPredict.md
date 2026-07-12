# Posterior predictive replicates

Simulates replicated data sets from the fitted mixture, conditionally on
each retained draw's parameters and allocations, and returns them – the
input that graphical posterior predictive checking (for example
`bayesplot::ppc_dens_overlay(y, yrep)`) consumes.
[`ppCheck`](https://madsyair.github.io/nimix/reference/ppCheck.md)
computes summary-statistic tail probabilities from the same replicates;
this function exposes the replicates themselves.

## Usage

``` r
posteriorPredict(fit, ndraws = 100, seed = 1L)
```

## Arguments

- fit:

  A `FitResult` from a clustering fit.

- ndraws:

  Number of posterior draws to use (thinned evenly).

- seed:

  RNG seed for the replicate noise.

## Value

For univariate data, an `ndraws x n` matrix; for multivariate data, an
`ndraws x n x d` array. The attribute `"draws"` records which posterior
iterations were used.

## See also

[`ppCheck`](https://madsyair.github.io/nimix/reference/ppCheck.md) for
tail-probability summaries.
