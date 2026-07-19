# Validate component parameters or a prior specification

Checks that a candidate parameter / prior list is internally consistent
for a given
[`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md).
Each concrete distribution must implement this.

## Usage

``` r
validateParams(spec, params, ...)

# S4 method for class 'FSSNUvSpec'
validateParams(spec, params, ...)

# S4 method for class 'FOSSEPUvSpec'
validateParams(spec, params, ...)

# S4 method for class 'FSSTUvSpec'
validateParams(spec, params, ...)

# S4 method for class 'PoissonRegSpec'
validateParams(spec, params, ...)

# S4 method for class 'BinomialRegSpec'
validateParams(spec, params, ...)

# S4 method for class 'GMSNBurrUvSpec'
validateParams(spec, params, ...)

# S4 method for class 'NormalRegSpec'
validateParams(spec, params, ...)

# S4 method for class 'JFSTUvSpec'
validateParams(spec, params, ...)

# S4 method for class 'SEPUvSpec'
validateParams(spec, params, ...)

# S4 method for class 'LEPUvSpec'
validateParams(spec, params, ...)

# S4 method for class 'MSNBurrUvSpec'
validateParams(spec, params, ...)

# S4 method for class 'MSNBurr2aUvSpec'
validateParams(spec, params, ...)

# S4 method for class 'NormalMvRegSpec'
validateParams(spec, params, ...)

# S4 method for class 'MSNBurrRegSpec'
validateParams(spec, params, ...)

# S4 method for class 'SEPRegSpec'
validateParams(spec, params, ...)

# S4 method for class 'MSNBurr2aRegSpec'
validateParams(spec, params, ...)

# S4 method for class 'FSSNRegSpec'
validateParams(spec, params, ...)

# S4 method for class 'GMSNBurrRegSpec'
validateParams(spec, params, ...)

# S4 method for class 'LEPRegSpec'
validateParams(spec, params, ...)

# S4 method for class 'FSSTRegSpec'
validateParams(spec, params, ...)

# S4 method for class 'FOSSEPRegSpec'
validateParams(spec, params, ...)

# S4 method for class 'JFSTRegSpec'
validateParams(spec, params, ...)

# S4 method for class 'NormalMvSpec'
validateParams(spec, params, ...)

# S4 method for class 'NormalGammaMvSpec'
validateParams(spec, params, ...)

# S4 method for class 'NormalUvSpec'
validateParams(spec, params, ...)

# S4 method for class 'NormalGammaUvSpec'
validateParams(spec, params, ...)

# S4 method for class 'PoissonSpec'
validateParams(spec, params, ...)

# S4 method for class 'BinomialSpec'
validateParams(spec, params, ...)

# S4 method for class 'SkewNormalMvSpec'
validateParams(spec, params, ...)

# S4 method for class 'SkewNormalMvOSpec'
validateParams(spec, params, ...)

# S4 method for class 'SkewIStudentMvSpec'
validateParams(spec, params, ...)

# S4 method for class 'SkewIStudentMvOSpec'
validateParams(spec, params, ...)

# S4 method for class 'SkewNormalMvOGenSpec'
validateParams(spec, params, ...)

# S4 method for class 'SkewIStudentMvOGenSpec'
validateParams(spec, params, ...)

# S4 method for class 'StudentTMvSpec'
validateParams(spec, params, ...)

# S4 method for class 'StudentTUvSpec'
validateParams(spec, params, ...)
```

## Arguments

- spec:

  A
  [`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md).

- params:

  A named list of parameters or hyperparameters to validate.

- ...:

  Reserved for methods.

## Value

Invisibly `TRUE` if valid; otherwise an error is raised.

## Functions

- `validateParams(FSSNUvSpec)`: FSSN hyperparameter checks.

- `validateParams(FOSSEPUvSpec)`: FOSSEP hyperparameter checks.

- `validateParams(FSSTUvSpec)`: FSST hyperparameter checks.

- `validateParams(PoissonRegSpec)`: Validate the coefficient prior.

- `validateParams(BinomialRegSpec)`: Validate the coefficient prior and
  `size`.

- `validateParams(GMSNBurrUvSpec)`: GMSNBurr hyperparameter checks.

- `validateParams(NormalRegSpec)`: Validate a Normal-Inverse-Gamma
  regression prior and enforce \\\dim(b_0) = \dim(B_0) = p\\ and \\nu_0
  \> 2\\.

- `validateParams(JFSTUvSpec)`: JFST hyperparameter checks.

- `validateParams(SEPUvSpec)`: SEP hyperparameter checks.

- `validateParams(LEPUvSpec)`: LEP hyperparameter checks.

- `validateParams(MSNBurrUvSpec)`: MSNBurr hyperparameter checks.

- `validateParams(MSNBurr2aUvSpec)`: MSNBurr-IIa hyperparameter checks.

- `validateParams(NormalMvRegSpec)`: Validate the
  multivariate-regression prior.

- `validateParams(MSNBurrRegSpec)`: MSNBurr regression prior validation.

- `validateParams(SEPRegSpec)`: SEP regression validation.

- `validateParams(MSNBurr2aRegSpec)`: MSNBurr-IIa regression validation.

- `validateParams(FSSNRegSpec)`: FSSN regression validation.

- `validateParams(GMSNBurrRegSpec)`: GMSNBurr regression validation.

- `validateParams(LEPRegSpec)`: LEP regression validation.

- `validateParams(FSSTRegSpec)`: FSST regression validation.

- `validateParams(FOSSEPRegSpec)`: FOSSEP regression validation.

- `validateParams(JFSTRegSpec)`: JFST regression validation.

- `validateParams(NormalMvSpec)`: Validate a Normal-Inverse-Wishart
  prior list and enforce the dimension invariant \\\dim(\mu_0) =
  \dim(S_0) = d\\ and \\df_0 \> d + 1\\ (.b encapsulation, ).

- `validateParams(NormalGammaMvSpec)`: Validate the NIW prior and the
  fixed `df`.

- `validateParams(NormalUvSpec)`: Validate a Normal-Inverse-Gamma prior
  list.

- `validateParams(NormalGammaUvSpec)`: Validate the NIG prior and the
  fixed `df`.

- `validateParams(PoissonSpec)`: Validate the Gamma rate prior.

- `validateParams(BinomialSpec)`: Validate the Beta prior and `size`.

- `validateParams(SkewNormalMvSpec)`: Skew-mv-Normal hyperparameter
  checks.

- `validateParams(SkewNormalMvOSpec)`: Checks the Householder bound.

- `validateParams(SkewIStudentMvSpec)`: Skew-mv-IStudent hyperparameter
  checks.

- `validateParams(SkewIStudentMvOSpec)`: Checks the Householder bound.

- `validateParams(SkewNormalMvOGenSpec)`: Checks the angle box.

- `validateParams(SkewIStudentMvOGenSpec)`: Checks the angle box.

- `validateParams(StudentTMvSpec)`: Validate the NIW prior and the fixed
  `df`.

- `validateParams(StudentTUvSpec)`: Validate the Student-t prior list.
