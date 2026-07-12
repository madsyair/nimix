# Observed data and posterior predictive replicates for graphical PPC

Packages `y` and `yrep` in the shapes `bayesplot`'s `ppc_*` functions
consume, e.g. `bayesplot::ppc_dens_overlay(d$y, d$yrep)`.

## Usage

``` r
ppcData(fit, ndraws = 100, margin = 1L, seed = 1L)
```

## Arguments

- fit:

  A `FitResult` from a clustering fit.

- ndraws:

  Number of replicate draws (rows of `yrep`).

- margin:

  For multivariate fits, which data dimension to extract (bayesplot's
  PPC graphics are univariate). Ignored for univariate fits.

- seed:

  RNG seed passed to
  [`posteriorPredict`](https://madsyair.github.io/nimix/reference/posteriorPredict.md).

## Value

A list with `y` (numeric vector) and `yrep` (`ndraws x n` matrix).

## See also

[`posteriorPredict`](https://madsyair.github.io/nimix/reference/posteriorPredict.md),
[`ppCheck`](https://madsyair.github.io/nimix/reference/ppCheck.md).
