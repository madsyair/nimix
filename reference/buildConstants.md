# Assemble distribution-specific NIMBLE constants for the DPM engine

Assemble distribution-specific NIMBLE constants for the DPM engine

## Usage

``` r
buildConstants(spec, prior, n, ...)

# S4 method for class 'FSSNUvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'FOSSEPUvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'FSSTUvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'PoissonRegSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'BinomialRegSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'GMSNBurrUvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'NormalRegSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'StudentTRegSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'NormalGammaRegSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'JFSTUvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'SEPUvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'LEPUvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'MSNBurrUvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'MSNBurr2aUvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'NormalMvRegSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'StudentTMvRegSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'NormalGammaMvRegSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'NormalMvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'NormalGammaMvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'NormalUvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'NormalGammaUvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'PoissonSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'BinomialSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'SkewNormalMvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'SkewNormalMvOSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'SkewIStudentMvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'SkewIStudentMvOSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'SkewNormalMvOGenSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'SkewIStudentMvOGenSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'StudentTMvSpec'
buildConstants(spec, prior, n, ...)

# S4 method for class 'StudentTUvSpec'
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

- `buildConstants(FSSNUvSpec)`: FSSN constants.

- `buildConstants(FOSSEPUvSpec)`: FOSSEP constants.

- `buildConstants(FSSTUvSpec)`: FSST constants.

- `buildConstants(PoissonRegSpec)`: Poisson regression constants.

- `buildConstants(BinomialRegSpec)`: Binomial regression constants plus
  `size`.

- `buildConstants(GMSNBurrUvSpec)`: GMSNBurr constants.

- `buildConstants(NormalRegSpec)`: Regression constants (design matrix
  `X`, number of predictors `p`, and the NIG hyperparameters).

- `buildConstants(StudentTRegSpec)`: NIG g-prior constants plus `df`.

- `buildConstants(NormalGammaRegSpec)`: NIG g-prior constants plus `df`.

- `buildConstants(JFSTUvSpec)`: JFST constants.

- `buildConstants(SEPUvSpec)`: SEP constants.

- `buildConstants(LEPUvSpec)`: LEP constants.

- `buildConstants(MSNBurrUvSpec)`: MSNBurr constants.

- `buildConstants(MSNBurr2aUvSpec)`: MSNBurr-IIa constants.

- `buildConstants(NormalMvRegSpec)`: Multivariate-regression constants.

- `buildConstants(StudentTMvRegSpec)`: Multivariate-regression constants
  plus `df`.

- `buildConstants(NormalGammaMvRegSpec)`: Multivariate-regression
  constants plus `df`.

- `buildConstants(NormalMvSpec)`: Multivariate Normal-Inverse-Wishart
  constants (includes the dimension `d`, mean vector `mu0` and scale
  matrix `S0`).

- `buildConstants(NormalGammaMvSpec)`: Normal-Inverse-Wishart constants
  plus `df`.

- `buildConstants(NormalUvSpec)`: Univariate Normal-Inverse-Gamma
  constants.

- `buildConstants(NormalGammaUvSpec)`: Normal-Inverse-Gamma constants
  plus `df`.

- `buildConstants(PoissonSpec)`: Poisson Gamma-prior constants.

- `buildConstants(BinomialSpec)`: Binomial Beta-prior constants plus
  `size`.

- `buildConstants(SkewNormalMvSpec)`: Skew-mv-Normal constants.

- `buildConstants(SkewNormalMvOSpec)`: Skew-mv-Normal-O constants.

- `buildConstants(SkewIStudentMvSpec)`: Skew-mv-IStudent constants.

- `buildConstants(SkewIStudentMvOSpec)`: Skew-mv-IStudent-O constants.

- `buildConstants(SkewNormalMvOGenSpec)`: General-m constants (angle box
  included).

- `buildConstants(SkewIStudentMvOGenSpec)`: General-m constants (angle
  box included).

- `buildConstants(StudentTMvSpec)`: Normal-Inverse-Wishart constants
  plus `df`.

- `buildConstants(StudentTUvSpec)`: Student-t constants
  (location/precision prior + df).
