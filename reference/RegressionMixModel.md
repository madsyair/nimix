# Construct a RegressionMixModel

Construct a RegressionMixModel

## Usage

``` r
RegressionMixModel(data, X, formula, distSpec, engine, Kmax, prior)
```

## Arguments

- data:

  Numeric response vector.

- X:

  Numeric design matrix.

- formula:

  The model formula.

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

A `RegressionMixModel`.
