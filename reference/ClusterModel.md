# Construct a ClusterModel

Construct a ClusterModel

## Usage

``` r
ClusterModel(data, distSpec, engine, Kmax, prior)
```

## Arguments

- data:

  Numeric vector (univariate) or numeric matrix (multivariate, one row
  per observation).

- distSpec:

  A
  [`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md).

- engine:

  An
  [`EngineConfig`](https://madsyair.github.io/nimix/reference/EngineConfig-class.md).

- Kmax:

  Integer truncation level.

- prior:

  Named list of prior hyperparameters.

## Value

A `ClusterModel`.
