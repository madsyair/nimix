# Is this a regression component spec?

Predicate the predict path uses to route to the regression branch.
Default `FALSE`; regression specs override it to `TRUE`.

## Usage

``` r
isRegressionSpec(spec, ...)

# S4 method for class 'DistributionSpec'
isRegressionSpec(spec, ...)

# S4 method for class 'PoissonRegSpec'
isRegressionSpec(spec, ...)

# S4 method for class 'BinomialRegSpec'
isRegressionSpec(spec, ...)

# S4 method for class 'NormalRegSpec'
isRegressionSpec(spec, ...)

# S4 method for class 'NormalMvRegSpec'
isRegressionSpec(spec, ...)
```

## Arguments

- spec:

  A
  [`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md).

- ...:

  Unused.

## Value

Logical scalar.

## Functions

- `isRegressionSpec(DistributionSpec)`: Default: not a regression spec.

- `isRegressionSpec(PoissonRegSpec)`: Poisson regression is a regression
  spec.

- `isRegressionSpec(BinomialRegSpec)`: Binomial regression is a
  regression spec.

- `isRegressionSpec(NormalRegSpec)`: Normal-linear regression is a
  regression spec.

- `isRegressionSpec(NormalMvRegSpec)`: Multivariate-response regression
  spec.
