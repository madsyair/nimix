# Build a data-scaled default prior for a distribution

Returns a named list of prior hyperparameters scaled to the observed
data, following the weakly-informative, data-scaled philosophy in
project knowledge (priors for location parameters must not be made
arbitrarily vague when `K_max` is large).

## Usage

``` r
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'FSSNUvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'FOSSEPUvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'FSSTUvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'PoissonRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'BinomialRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'GMSNBurrUvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'NormalRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'StudentTRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'NormalGammaRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'JFSTUvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'SEPUvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'LEPUvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'MSNBurrUvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'MSNBurr2aUvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'NormalMvRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'StudentTMvRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'NormalGammaMvRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'MSNBurrRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'SEPRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'MSNBurr2aRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'FSSNRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'GMSNBurrRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'LEPRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'FSSTRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'FOSSEPRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'JFSTRegSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'NormalMvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'NormalGammaMvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'NormalUvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'NormalGammaUvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'PoissonSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'BinomialSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'SkewNormalMvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'SkewNormalMvOSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'SkewIStudentMvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'SkewIStudentMvOSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'SkewNormalMvOGenSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'SkewIStudentMvOGenSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'StudentTMvSpec'
defaultPrior(spec, data, control = list(), ...)

# S4 method for class 'StudentTUvSpec'
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

- `defaultPrior(FSSNUvSpec)`: Data-scaled FSSN prior (log-normal
  skewness).

- `defaultPrior(FOSSEPUvSpec)`: Data-scaled FOSSEP prior.

- `defaultPrior(FSSTUvSpec)`: Data-scaled FSST prior.

- `defaultPrior(PoissonRegSpec)`: g-prior on the coefficients.

- `defaultPrior(BinomialRegSpec)`: g-prior on the coefficients; needs
  `size`.

- `defaultPrior(GMSNBurrUvSpec)`: Data-scaled GMSNBurr prior.

- `defaultPrior(NormalRegSpec)`: Data-scaled Normal-Inverse-Gamma
  g-prior for the regression component. Requires the design matrix in
  `control$X`.

  Control overrides: `g` (g-prior factor; prior covariance of `beta` is
  `s2 * g * solve(crossprod(X))`, default `g = n`, the unit-information
  prior) and `nu0` (InvGamma shape, default 3, must exceed 2 for a
  finite prior variance).

- `defaultPrior(StudentTRegSpec)`: NIG g-prior plus a fixed `df`
  (default 4, \> 2).

- `defaultPrior(NormalGammaRegSpec)`: NIG g-prior plus a fixed `df`
  (default 4, \> 2).

- `defaultPrior(JFSTUvSpec)`: Data-scaled JFST prior (symmetric in
  alpha/theta).

- `defaultPrior(SEPUvSpec)`: Data-scaled SEP prior.

- `defaultPrior(LEPUvSpec)`: Data-scaled LEP prior.

- `defaultPrior(MSNBurrUvSpec)`: Data-scaled MSNBurr prior
  (location/scale/shape).

- `defaultPrior(MSNBurr2aUvSpec)`: Data-scaled MSNBurr-IIa prior.

- `defaultPrior(NormalMvRegSpec)`: Inverse-Wishart on Sigma +
  matrix-normal coefficients.

  Control overrides: `g` (the coefficient g-prior scale, default `n`)
  and `sigmaGuess`, the prior mean of a component's residual covariance
  – a positive scalar (read as isotropic) or a \\d \times d\\
  positive-definite matrix.

  `sigmaGuess` exists because the automatic reference is
  [`cov()`](https://rdrr.io/r/stats/cor.html) of the *global* OLS
  residuals, which carry the between-component variation as well as the
  within-component one. This is the multivariate face of the univariate
  problem behind
  [`nimixReg`](https://madsyair.github.io/nimix/reference/nimixReg.md)'s
  `s2Guess`. Measured on two components differing only in their
  coefficients, with isotropic within-covariance: the prior mean of
  Sigma was 22.4x the truth in trace, with condition number 40 where the
  truth is a circle – wrong size and wrong shape.

  It stays the default anyway, for the reason the clustering case does:
  `df0 = d + 2` makes the InverseWishart prior worth exactly one
  observation for every \\d\\ (the univariate InvGamma is worth four),
  so the fitted covariance came back at about 1.5x with condition ~2
  against a prior off by 22x. Reach for `sigmaGuess` when you know the
  residual scale – most often `sigmaGuess = s` for isotropic residuals
  of variance `s`. `"studentt"` and `"normalgamma"` inherit it.

- `defaultPrior(StudentTMvRegSpec)`: Inverse-Wishart + matrix-normal
  plus a fixed `df`.

- `defaultPrior(NormalGammaMvRegSpec)`: Inverse-Wishart + matrix-normal
  plus a fixed `df`.

- `defaultPrior(MSNBurrRegSpec)`: MSNBurr regression prior.

- `defaultPrior(SEPRegSpec)`: SEP regression prior.

- `defaultPrior(MSNBurr2aRegSpec)`: MSNBurr-IIa regression prior (shares
  the MSNBurr scale/shape prior).

- `defaultPrior(FSSNRegSpec)`: FSSN regression prior (log-normal
  skewness).

- `defaultPrior(GMSNBurrRegSpec)`: GMSNBurr regression prior.

- `defaultPrior(LEPRegSpec)`: LEP regression prior.

- `defaultPrior(FSSTRegSpec)`: FSST regression prior.

- `defaultPrior(FOSSEPRegSpec)`: FOSSEP regression prior.

- `defaultPrior(JFSTRegSpec)`: JFST regression prior.

- `defaultPrior(NormalMvSpec)`: Data-scaled Normal-Inverse-Wishart prior
  (multivariate).

  Control overrides: `cLoc` (mean-dispersion multiplier; prior
  covariance of `mu` is `Sigma / kappa0` with `kappa0 = 1 / cLoc^2`,
  default `cLoc = 2`) and `df0` (inverse-Wishart degrees of freedom,
  default `d + 2`, must exceed `d + 1` for a finite, non-singular prior
  covariance on empty components).

  `sigmaGuess` sets the prior mean of a component's covariance directly:
  a positive scalar (read as isotropic) or a \\d \times d\\
  positive-definite matrix. The default reference is `cov(data)`, the
  *global* covariance, which for separated components is inflated only
  along the direction that separates them – so the prior ellipse has the
  wrong shape, not merely the wrong size. Measured on isotropic
  components separated along a vector v: 37.6x along v against 0.9x
  across it, a condition number of 42.8 where the truth is a circle. The
  InverseWishart default absorbs this well (`df0 = d + 2` makes the
  prior worth exactly one observation for every \\d\\; the fitted
  covariance came back at 1.2x with condition 1.4), so it stays.
  `sigmaGuess` is for callers who know the within-component shape – most
  often `sigmaGuess = s` for an isotropic component of variance `s`.

- `defaultPrior(NormalGammaMvSpec)`: Normal-Inverse-Wishart prior plus a
  fixed `df` (default 5, must exceed 2 for a finite component
  covariance).

- `defaultPrior(NormalUvSpec)`: Data-scaled Normal-Inverse-Gamma prior.

  Control overrides: `cLoc` (location spread multiplier; prior sd of
  `mu` ~ `cLoc * sd(data)`, default 2), `nu0` (InvGamma shape, default
  3, must exceed 2 for a finite prior variance), and `concPrior`.

- `defaultPrior(NormalGammaUvSpec)`: Normal-Inverse-Gamma prior plus a
  fixed `df` (degrees of freedom, default 4, must exceed 2 for a finite
  component variance).

- `defaultPrior(PoissonSpec)`: Data-scaled Gamma prior on the Poisson
  rate (`E[lambda] ~= mean(y)`).

- `defaultPrior(BinomialSpec)`: Data-scaled Beta prior on the success
  probability. Requires the number of trials in `control$size`.

- `defaultPrior(SkewNormalMvSpec)`: Data-scaled skew-mv-Normal prior.

- `defaultPrior(SkewNormalMvOSpec)`: Adds the Householder angle bound to
  the mv prior.

- `defaultPrior(SkewIStudentMvSpec)`: Data-scaled skew-mv-IStudent
  prior.

- `defaultPrior(SkewIStudentMvOSpec)`: Adds the Householder angle bound.

- `defaultPrior(SkewNormalMvOGenSpec)`: Adds the Householder angle box.

- `defaultPrior(SkewIStudentMvOGenSpec)`: Adds the Householder angle
  box.

- `defaultPrior(StudentTMvSpec)`: Normal-Inverse-Wishart prior plus a
  fixed `df` (default 5, must exceed 2 for a finite component
  covariance).

- `defaultPrior(StudentTUvSpec)`: Data-scaled Normal location / Gamma
  precision prior for the Student-t component.

  Control overrides: `cLoc` (location spread multiplier, default 2),
  `df` (degrees of freedom, a fixed hyperparameter, default 4, must
  exceed 2 for a finite component variance).
