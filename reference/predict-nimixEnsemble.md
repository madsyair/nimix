# Weighted predictive density from a mixture ensemble

Weighted predictive density from a mixture ensemble

## Usage

``` r
# S4 method for class 'nimixEnsemble'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  A `nimixEnsemble` from
  [`ensembleFit`](https://madsyair.github.io/nimix/reference/ensembleFit.md).

- newdata:

  Points at which to evaluate the density (univariate) or a matrix of
  rows (multivariate). Defaults to each model's own data grid.

- ...:

  Unused.

## Value

A data.frame of evaluation points and the ensemble-weighted density.
