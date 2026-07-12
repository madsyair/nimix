# Build dispersed initial values (engine-agnostic)

Build dispersed initial values (engine-agnostic)

## Usage

``` r
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'FSSNUvSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'FOSSEPUvSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'FSSTUvSpec'
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

# S4 method for class 'JFSTUvSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'SEPUvSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'LEPUvSpec'
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

# S4 method for class 'SkewNormalMvSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'SkewNormalMvOSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'SkewIStudentMvSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'SkewIStudentMvOSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'SkewNormalMvOGenSpec'
componentInits(spec, prior, data, count, initMethod = "kmeans", ...)

# S4 method for class 'SkewIStudentMvOGenSpec'
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

- `componentInits(FSSNUvSpec)`: Dispersed k-means start for FSSN (alpha
  at 1).

- `componentInits(FOSSEPUvSpec)`: Dispersed k-means start for FOSSEP.

- `componentInits(FSSTUvSpec)`: Dispersed k-means start for FSST.

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

- `componentInits(JFSTUvSpec)`: Dispersed k-means start for JFST
  (symmetric start).

- `componentInits(SEPUvSpec)`: Dispersed k-means start for SEP.

- `componentInits(LEPUvSpec)`: Dispersed k-means start for LEP.

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

- `componentInits(SkewNormalMvSpec)`: Dispersed k-means start (gamma at
  1: symmetric).

- `componentInits(SkewNormalMvOSpec)`: k-means start; theta
  grid-initialised per cluster.

- `componentInits(SkewIStudentMvSpec)`: k-means start (gamma at 1, nu at
  8).

- `componentInits(SkewIStudentMvOSpec)`: k-means start; theta
  grid-initialised per cluster.

- `componentInits(SkewNormalMvOGenSpec)`: k-means start; angles at the
  box midpoint.

- `componentInits(SkewIStudentMvOGenSpec)`: k-means start; angles at the
  box midpoint.

- `componentInits(StudentTUvSpec)`: k-means dispersed start (location +
  precision).
