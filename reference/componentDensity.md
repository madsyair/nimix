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

# S4 method for class 'NormalUvSpec'
componentDensity(spec, ...)

# S4 method for class 'NormalMvSpec'
componentDensity(spec, ...)

# S4 method for class 'NormalRegSpec'
componentDensity(spec, ...)

# S4 method for class 'PoissonRegSpec'
componentDensity(spec, ...)

# S4 method for class 'BinomialRegSpec'
componentDensity(spec, size = NULL, ...)

# S4 method for class 'NormalMvRegSpec'
componentDensity(spec, ...)

# S4 method for class 'StudentTUvSpec'
componentDensity(spec, df = 4, ...)

# S4 method for class 'NormalGammaUvSpec'
componentDensity(spec, df = 4, ...)

# S4 method for class 'StudentTMvSpec'
componentDensity(spec, df = 5, ...)

# S4 method for class 'NormalGammaMvSpec'
componentDensity(spec, df = 5, ...)

# S4 method for class 'PoissonSpec'
componentDensity(spec, ...)

# S4 method for class 'BinomialSpec'
componentDensity(spec, size = NULL, ...)
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

- `componentDensity(NormalUvSpec)`: Gaussian density for posterior
  predictive use.

- `componentDensity(NormalMvSpec)`: Multivariate normal density for
  predictive use.

- `componentDensity(NormalRegSpec)`: Gaussian density of a response
  given its linear predictor. `params` must carry `mu` (the fitted mean)
  and `s2`.

- `componentDensity(PoissonRegSpec)`: Poisson pmf at the fitted rate.

- `componentDensity(BinomialRegSpec)`: Binomial pmf at the fitted
  probability.

- `componentDensity(NormalMvRegSpec)`: Multivariate normal density of
  the residual.

- `componentDensity(StudentTUvSpec)`: Location-scale Student-t density.

- `componentDensity(NormalGammaUvSpec)`: Student-t marginal density
  (location `mu`, scale \\\sqrt{s2}\\, `df`).

- `componentDensity(StudentTMvSpec)`: Multivariate Student-t density.

- `componentDensity(NormalGammaMvSpec)`: Multivariate Student-t marginal
  density.

- `componentDensity(PoissonSpec)`: Poisson pmf.

- `componentDensity(BinomialSpec)`: Binomial pmf.
