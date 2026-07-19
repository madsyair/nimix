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

# S4 method for class 'MSNBurrRegSpec'
isRegressionSpec(spec)

# S4 method for class 'SEPRegSpec'
isRegressionSpec(spec)

# S4 method for class 'MSNBurr2aRegSpec'
isRegressionSpec(spec)

# S4 method for class 'FSSNRegSpec'
isRegressionSpec(spec)

# S4 method for class 'GMSNBurrRegSpec'
isRegressionSpec(spec)

# S4 method for class 'LEPRegSpec'
isRegressionSpec(spec)

# S4 method for class 'FSSTRegSpec'
isRegressionSpec(spec)

# S4 method for class 'FOSSEPRegSpec'
isRegressionSpec(spec)

# S4 method for class 'JFSTRegSpec'
isRegressionSpec(spec)
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

- `isRegressionSpec(MSNBurrRegSpec)`: MSNBurr regression is a regression
  spec.

- `isRegressionSpec(SEPRegSpec)`: SEP regression is a regression spec.

- `isRegressionSpec(MSNBurr2aRegSpec)`: MSNBurr-IIa regression is a
  regression spec.

- `isRegressionSpec(FSSNRegSpec)`: FSSN regression is a regression spec.

- `isRegressionSpec(GMSNBurrRegSpec)`: GMSNBurr regression is a
  regression spec.

- `isRegressionSpec(LEPRegSpec)`: LEP regression is a regression spec.

- `isRegressionSpec(FSSTRegSpec)`: FSST regression is a regression spec.

- `isRegressionSpec(FOSSEPRegSpec)`: FOSSEP regression is a regression
  spec.

- `isRegressionSpec(JFSTRegSpec)`: JFST regression is a regression spec.
