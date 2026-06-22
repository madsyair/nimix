# Build a data-scaled default prior for a distribution

Returns a named list of prior hyperparameters scaled to the observed
data, following the weakly-informative, data-scaled philosophy in
project knowledge (priors for location parameters must not be made
arbitrarily vague when `K_max` is large).

## Usage

``` r
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'NormalUvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'NormalMvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'NormalRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'PoissonRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'BinomialRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'StudentTRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'NormalGammaRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'NormalMvRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'StudentTMvRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'NormalGammaMvRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'StudentTUvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'NormalGammaUvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'StudentTMvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'NormalGammaMvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'PoissonSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'BinomialSpec'
defaultPrior(spec, data, control = list(), ...)
```

## Arguments

- spec:

  A
  [`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md).

- data:

  Numeric data used to scale the prior.

- control:

  A named list of user overrides (merged over the defaults).

- ...:

  Reserved for methods.

## Value

A named list of prior hyperparameters.

## Functions

- `defaultPrior(NormalUvSpec)`: Data-scaled Normal-Inverse-Gamma prior.

  Control overrides: `cLoc` (location spread multiplier; prior sd of
  `mu` ~ `cLoc * sd(data)`, default 2), `nu0` (InvGamma shape, default
  3, must exceed 2 for a finite prior variance), and `concPrior`.

- `defaultPrior(NormalMvSpec)`: Data-scaled Normal-Inverse-Wishart prior
  (multivariate).

  Control overrides: `cLoc` (mean-dispersion multiplier; prior
  covariance of `mu` is `Sigma / kappa0` with `kappa0 = 1 / cLoc^2`,
  default `cLoc = 2`) and `df0` (inverse-Wishart degrees of freedom,
  default `d + 2`, must exceed `d + 1` for a finite, non-singular prior
  covariance on empty components).

- `defaultPrior(NormalRegSpec)`: Data-scaled Normal-Inverse-Gamma
  g-prior for the regression component. Requires the design matrix in
  `control$X`.

  Control overrides: `g` (g-prior factor; prior covariance of `beta` is
  `s2 * g * solve(crossprod(X))`, default `g = n`, the unit-information
  prior) and `nu0` (InvGamma shape, default 3, must exceed 2 for a
  finite prior variance).

- `defaultPrior(PoissonRegSpec)`: g-prior on the coefficients.

- `defaultPrior(BinomialRegSpec)`: g-prior on the coefficients; needs
  `size`.

- `defaultPrior(StudentTRegSpec)`: NIG g-prior plus a fixed `df`
  (default 4, \> 2).

- `defaultPrior(NormalGammaRegSpec)`: NIG g-prior plus a fixed `df`
  (default 4, \> 2).

- `defaultPrior(NormalMvRegSpec)`: Inverse-Wishart on Sigma +
  matrix-normal coefficients.

- `defaultPrior(StudentTMvRegSpec)`: Inverse-Wishart + matrix-normal
  plus a fixed `df`.

- `defaultPrior(NormalGammaMvRegSpec)`: Inverse-Wishart + matrix-normal
  plus a fixed `df`.

- `defaultPrior(StudentTUvSpec)`: Data-scaled Normal location / Gamma
  precision prior for the Student-t component.

  Control overrides: `cLoc` (location spread multiplier, default 2),
  `df` (degrees of freedom, a fixed hyperparameter, default 4, must
  exceed 2 for a finite component variance).

- `defaultPrior(NormalGammaUvSpec)`: Normal-Inverse-Gamma prior plus a
  fixed `df` (degrees of freedom, default 4, must exceed 2 for a finite
  component variance).

- `defaultPrior(StudentTMvSpec)`: Normal-Inverse-Wishart prior plus a
  fixed `df` (default 5, must exceed 2 for a finite component
  covariance).

- `defaultPrior(NormalGammaMvSpec)`: Normal-Inverse-Wishart prior plus a
  fixed `df` (default 5, must exceed 2 for a finite component
  covariance).

- `defaultPrior(PoissonSpec)`: Data-scaled Gamma prior on the Poisson
  rate (`E[lambda] ~= mean(y)`).

- `defaultPrior(BinomialSpec)`: Data-scaled Beta prior on the success
  probability. Requires the number of trials in `control$size`.
