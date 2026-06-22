# Validate component parameters or a prior specification

Checks that a candidate parameter / prior list is internally consistent
for a given
[`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md).
Each concrete distribution must implement this.

## Usage

``` r
validateParams(spec, params, ...)

# S4 method for class 'NormalUvSpec'
validateParams(spec, params, ...)

# S4 method for class 'NormalMvSpec'
validateParams(spec, params, ...)

# S4 method for class 'NormalRegSpec'
validateParams(spec, params, ...)

# S4 method for class 'PoissonRegSpec'
validateParams(spec, params, ...)

# S4 method for class 'BinomialRegSpec'
validateParams(spec, params, ...)

# S4 method for class 'NormalMvRegSpec'
validateParams(spec, params, ...)

# S4 method for class 'StudentTUvSpec'
validateParams(spec, params, ...)

# S4 method for class 'NormalGammaUvSpec'
validateParams(spec, params, ...)

# S4 method for class 'StudentTMvSpec'
validateParams(spec, params, ...)

# S4 method for class 'NormalGammaMvSpec'
validateParams(spec, params, ...)

# S4 method for class 'PoissonSpec'
validateParams(spec, params, ...)

# S4 method for class 'BinomialSpec'
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

- `validateParams(NormalUvSpec)`: Validate a Normal-Inverse-Gamma prior
  list.

- `validateParams(NormalMvSpec)`: Validate a Normal-Inverse-Wishart
  prior list and enforce the dimension invariant \\\dim(\mu_0) =
  \dim(S_0) = d\\ and \\df_0 \> d + 1\\ (.b encapsulation, ).

- `validateParams(NormalRegSpec)`: Validate a Normal-Inverse-Gamma
  regression prior and enforce \\\dim(b_0) = \dim(B_0) = p\\ and \\nu_0
  \> 2\\.

- `validateParams(PoissonRegSpec)`: Validate the coefficient prior.

- `validateParams(BinomialRegSpec)`: Validate the coefficient prior and
  `size`.

- `validateParams(NormalMvRegSpec)`: Validate the
  multivariate-regression prior.

- `validateParams(StudentTUvSpec)`: Validate the Student-t prior list.

- `validateParams(NormalGammaUvSpec)`: Validate the NIG prior and the
  fixed `df`.

- `validateParams(StudentTMvSpec)`: Validate the NIW prior and the fixed
  `df`.

- `validateParams(NormalGammaMvSpec)`: Validate the NIW prior and the
  fixed `df`.

- `validateParams(PoissonSpec)`: Validate the Gamma rate prior.

- `validateParams(BinomialSpec)`: Validate the Beta prior and `size`.
