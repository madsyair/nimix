# Component density evaluator (R-level, for posterior predictive checks)

Returns a function `f(x, params)` giving the component density at `x`.
Used by [`predict()`](https://rdrr.io/r/stats/predict.html) to build
posterior predictive densities. For univariate specs `x` is a scalar and
`params` a list with scalar entries; for multivariate specs `x` is a
length-`d` vector and `params$mu` a vector, `params$Sigma` a `d x d`
matrix.

## Usage

``` r
componentDensity(spec, ...)

# S4 method for class 'FSSNUvSpec'
componentDensity(spec, ...)

# S4 method for class 'FOSSEPUvSpec'
componentDensity(spec, ...)

# S4 method for class 'FSSTUvSpec'
componentDensity(spec, ...)

# S4 method for class 'PoissonRegSpec'
componentDensity(spec, ...)

# S4 method for class 'BinomialRegSpec'
componentDensity(spec, size = NULL, ...)

# S4 method for class 'GMSNBurrUvSpec'
componentDensity(spec, ...)

# S4 method for class 'NormalRegSpec'
componentDensity(spec, ...)

# S4 method for class 'JFSTUvSpec'
componentDensity(spec, ...)

# S4 method for class 'SEPUvSpec'
componentDensity(spec, ...)

# S4 method for class 'LEPUvSpec'
componentDensity(spec, ...)

# S4 method for class 'MSNBurrUvSpec'
componentDensity(spec, ...)

# S4 method for class 'MSNBurr2aUvSpec'
componentDensity(spec, ...)

# S4 method for class 'NormalMvRegSpec'
componentDensity(spec, ...)

# S4 method for class 'NormalMvSpec'
componentDensity(spec, ...)

# S4 method for class 'NormalGammaMvSpec'
componentDensity(spec, df = 5, ...)

# S4 method for class 'NormalUvSpec'
componentDensity(spec, ...)

# S4 method for class 'NormalGammaUvSpec'
componentDensity(spec, df = 4, ...)

# S4 method for class 'PoissonSpec'
componentDensity(spec, ...)

# S4 method for class 'BinomialSpec'
componentDensity(spec, size = NULL, ...)

# S4 method for class 'SkewNormalMvSpec'
componentDensity(spec, ...)

# S4 method for class 'SkewNormalMvOSpec'
componentDensity(spec, ...)

# S4 method for class 'SkewIStudentMvSpec'
componentDensity(spec, ...)

# S4 method for class 'SkewIStudentMvOSpec'
componentDensity(spec, ...)

# S4 method for class 'SkewNormalMvOGenSpec'
componentDensity(spec, ...)

# S4 method for class 'SkewIStudentMvOGenSpec'
componentDensity(spec, ...)

# S4 method for class 'StudentTMvSpec'
componentDensity(spec, df = 5, ...)

# S4 method for class 'StudentTUvSpec'
componentDensity(spec, df = 4, ...)
```

## Arguments

- spec:

  A
  [`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md).

- ...:

  Reserved for methods.

## Value

A function of `(x, params)`.

## Functions

- `componentDensity(FSSNUvSpec)`: FSSN density closure (stable reference
  form).

- `componentDensity(FOSSEPUvSpec)`: FOSSEP density closure (stable
  reference form).

- `componentDensity(FSSTUvSpec)`: FSST density closure (stable reference
  form).

- `componentDensity(PoissonRegSpec)`: Poisson pmf at the fitted rate.

- `componentDensity(BinomialRegSpec)`: Binomial pmf at the fitted
  probability.

- `componentDensity(GMSNBurrUvSpec)`: GMSNBurr density closure (stable
  reference form).

- `componentDensity(NormalRegSpec)`: Gaussian density of a response
  given its linear predictor. `params` must carry `mu` (the fitted mean)
  and `s2`.

- `componentDensity(JFSTUvSpec)`: JFST density closure (stable reference
  form).

- `componentDensity(SEPUvSpec)`: SEP density closure (stable reference
  form).

- `componentDensity(LEPUvSpec)`: LEP density closure (stable reference
  form).

- `componentDensity(MSNBurrUvSpec)`: MSNBurr density closure (stable
  reference form).

- `componentDensity(MSNBurr2aUvSpec)`: MSNBurr-IIa density closure
  (stable form).

- `componentDensity(NormalMvRegSpec)`: Multivariate normal density of
  the residual.

- `componentDensity(NormalMvSpec)`: Multivariate normal density for
  predictive use.

- `componentDensity(NormalGammaMvSpec)`: Multivariate Student-t marginal
  density.

- `componentDensity(NormalUvSpec)`: Gaussian density for posterior
  predictive use.

- `componentDensity(NormalGammaUvSpec)`: Student-t marginal density
  (location `mu`, scale \\\sqrt{s2}\\, `df`).

- `componentDensity(PoissonSpec)`: Poisson pmf.

- `componentDensity(BinomialSpec)`: Binomial pmf.

- `componentDensity(SkewNormalMvSpec)`: Skew-mv-Normal density closure.

- `componentDensity(SkewNormalMvOSpec)`: Density closure with the
  Householder angle.

- `componentDensity(SkewIStudentMvSpec)`: Skew-mv-IStudent density
  closure.

- `componentDensity(SkewIStudentMvOSpec)`: Density closure with nu and
  the angle.

- `componentDensity(SkewNormalMvOGenSpec)`: Density closure with a
  general orthogonal factor.

- `componentDensity(SkewIStudentMvOGenSpec)`: Density closure, general O
  and per-margin nu.

- `componentDensity(StudentTMvSpec)`: Multivariate Student-t density.

- `componentDensity(StudentTUvSpec)`: Location-scale Student-t density.
