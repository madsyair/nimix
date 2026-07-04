# Build dispersed initial values (engine-agnostic)

Build dispersed initial values (engine-agnostic)

## Usage

``` r
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'PoissonRegSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'BinomialRegSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'GMSNBurrUvSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'NormalRegSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'NormalGammaRegSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'MSNBurrUvSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'MSNBurr2aUvSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'NormalMvRegSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'NormalGammaMvRegSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'NormalMvSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'NormalGammaMvSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'NormalUvSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'NormalGammaUvSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'PoissonSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'BinomialSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'StudentTUvSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)
```

## Arguments

- spec:

  A
  [`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md).

- prior:

  The prior list.

- data:

  The observed data.

- count:

  Integer number of component slots (K_max for the DPM, K for the finite
  mixture).

- initMethod:

  `"kmeans"` (default) or `"single"`.

- ...:

  Reserved for methods.

## Value

A list with `alloc` (an integer allocation vector) and `params` (a named
list of component-parameter initial values).

## Functions

- `componentInits(PoissonRegSpec)`: Global-GLM start with k-means
  allocation.

- `componentInits(BinomialRegSpec)`: Global-GLM start with k-means
  allocation.

- `componentInits(GMSNBurrUvSpec)`: Dispersed k-means start for
  GMSNBurr.

- `componentInits(NormalRegSpec)`: k-means-on-(predictors, response)
  start with local OLS coefficients per cluster (dispersed start).

- `componentInits(NormalGammaRegSpec)`: k-means start (inherited) plus
  unit `omega`.

- `componentInits(MSNBurrUvSpec)`: Dispersed k-means start for MSNBurr.

- `componentInits(MSNBurr2aUvSpec)`: Dispersed k-means start for
  MSNBurr-IIa.

- `componentInits(NormalMvRegSpec)`: Global multivariate-OLS start,
  k-means allocation.

- `componentInits(NormalGammaMvRegSpec)`: Global multivariate-OLS start
  plus unit `omega`.

- `componentInits(NormalMvSpec)`: k-means dispersed start for the
  multivariate DPM.

- `componentInits(NormalGammaMvSpec)`: k-means start (inherited) plus
  unit `omega`.

- `componentInits(NormalUvSpec)`: k-means dispersed start for univariate
  DPM.

  A k-means allocation gives a dispersed initial partition that shortens
  burn-in.

- `componentInits(NormalGammaUvSpec)`: k-means start (inherited) plus
  unit `omega`.

- `componentInits(PoissonSpec)`: k-means start on counts.

- `componentInits(BinomialSpec)`: k-means start on proportions.

- `componentInits(StudentTUvSpec)`: k-means dispersed start (location +
  precision).
