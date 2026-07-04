# Profile the clusters of a fitted mixture

Assigns each observation to its maximum a posteriori (MAP) cluster – the
component it occupies most often across the retained MCMC draws – and
summarises the observed data within each cluster: cluster size,
proportion, and per-variable within-cluster mean, standard deviation,
and median. This characterises what each cluster *is* in terms of the
data, complementing [`summary`](https://rdrr.io/r/base/summary.html),
which reports the fitted component parameters, and
`plot(fit, type = "cluster")`, which shows the partition visually.

## Usage

``` r
clusterProfile(fit, variables = NULL)
```

## Arguments

- fit:

  A
  [`FitResult`](https://madsyair.github.io/nimix/reference/FitResult-class.md)
  from
  [`nimixClust`](https://madsyair.github.io/nimix/reference/nimixClust.md).
  Regression fits are also accepted: the response and each covariate are
  profiled per regime.

- variables:

  Optional character vector selecting which data columns to profile
  (multivariate / regression); default profiles all.

## Value

A data.frame with one row per occupied cluster and columns `cluster`,
`size`, `proportion`, followed by `<var>_mean`, `<var>_sd`,
`<var>_median` for each variable.

## Details

Cluster ids are recoded 1, 2, ... by descending size (largest cluster
first). Because mixture labels are not identified, this data-side
profile does not require
[`relabel`](https://madsyair.github.io/nimix/reference/relabel.md): it
summarises a partition, and a partition is label-invariant.

## See also

[`summary`](https://rdrr.io/r/base/summary.html),
`plot(fit, type = "cluster")`

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- nimixClust(y, K = 3, method = "fixedk")
clusterProfile(fit)
} # }
```
