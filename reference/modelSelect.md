# Compare mixture models by predictive fit

Ranks several fitted mixtures by WAIC (native) and, when loo is
available, by PSIS-LOO. Models must be fitted to the same data. Use to
choose K, or to compare component families (e.g. Normal vs Student-t vs
MSNBurr) on the same data.

## Usage

``` r
modelSelect(..., maxDraws = 1000L)
```

## Arguments

- ...:

  Two or more clustering
  [`FitResult`](https://madsyair.github.io/nimix/reference/FitResult-class.md)
  objects, or a single named list of them.

- maxDraws:

  Cap on posterior draws used per fit. Default 1000.

## Value

A data.frame ordered best-first, with WAIC, elpd, and (if available) LOO
columns, plus `dWAIC` relative to the best model.
