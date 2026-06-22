# Assemble distribution-specific NIMBLE constants for the DPM engine

Assemble distribution-specific NIMBLE constants for the DPM engine

## Usage

``` r
buildConstants(spec, prior, n, ...)

# S4 method for class 'NormalUvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'NormalMvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'NormalRegSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'PoissonRegSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'BinomialRegSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'StudentTRegSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'NormalGammaRegSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'NormalMvRegSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'StudentTMvRegSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'NormalGammaMvRegSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'StudentTUvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'NormalGammaUvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'StudentTMvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'NormalGammaMvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'PoissonSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'BinomialSpec'
buildConstants(spec, prior, n, ...)
```

## Arguments

- spec:

  A
  [`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md).

- prior:

  The prior list (from
  [`defaultPrior`](https://madsyair.github.io/nimix/reference/defaultPrior.md)).

- n:

  Integer number of observations.

- ...:

  Reserved for methods.

## Value

A named list of constants (excluding the concentration hyperprior, which
the engine appends).

## Functions

- `buildConstants(NormalUvSpec)`: Univariate Normal-Inverse-Gamma
  constants.

- `buildConstants(NormalMvSpec)`: Multivariate Normal-Inverse-Wishart
  constants (includes the dimension `d`, mean vector `mu0` and scale
  matrix `S0`).

- `buildConstants(NormalRegSpec)`: Regression constants (design matrix
  `X`, number of predictors `p`, and the NIG hyperparameters).

- `buildConstants(PoissonRegSpec)`: Poisson regression constants.

- `buildConstants(BinomialRegSpec)`: Binomial regression constants plus
  `size`.

- `buildConstants(StudentTRegSpec)`: NIG g-prior constants plus `df`.

- `buildConstants(NormalGammaRegSpec)`: NIG g-prior constants plus `df`.

- `buildConstants(NormalMvRegSpec)`: Multivariate-regression constants.

- `buildConstants(StudentTMvRegSpec)`: Multivariate-regression constants
  plus `df`.

- `buildConstants(NormalGammaMvRegSpec)`: Multivariate-regression
  constants plus `df`.

- `buildConstants(StudentTUvSpec)`: Student-t constants
  (location/precision prior + df).

- `buildConstants(NormalGammaUvSpec)`: Normal-Inverse-Gamma constants
  plus `df`.

- `buildConstants(StudentTMvSpec)`: Normal-Inverse-Wishart constants
  plus `df`.

- `buildConstants(NormalGammaMvSpec)`: Normal-Inverse-Wishart constants
  plus `df`.

- `buildConstants(PoissonSpec)`: Poisson Gamma-prior constants.

- `buildConstants(BinomialSpec)`: Binomial Beta-prior constants plus
  `size`.
