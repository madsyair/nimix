# Parse raw cluster-parameter traces from the MCMC sample matrix

Parse raw cluster-parameter traces from the MCMC sample matrix

## Usage

``` r
extractParamTraces(spec, samples, L, ...)

# S4 method for class 'PoissonRegSpec'
extractParamTraces(spec, samples, L, d = NULL, prior = NULL, ...)

# S4 method for class 'BinomialRegSpec'
extractParamTraces(spec, samples, L, d = NULL, prior = NULL, ...)

# S4 method for class 'GMSNBurrUvSpec'
extractParamTraces(spec, samples, L, ...)

# S4 method for class 'NormalRegSpec'
extractParamTraces(spec, samples, L, d = NULL, prior = NULL, ...)

# S4 method for class 'MSNBurrUvSpec'
extractParamTraces(spec, samples, L, ...)

# S4 method for class 'MSNBurr2aUvSpec'
extractParamTraces(spec, samples, L, ...)

# S4 method for class 'NormalMvRegSpec'
extractParamTraces(spec, samples, L, d = NULL, prior = NULL, ...)

# S4 method for class 'NormalMvSpec'
extractParamTraces(spec, samples, L, d = NULL, ...)

# S4 method for class 'NormalUvSpec'
extractParamTraces(spec, samples, L, ...)

# S4 method for class 'PoissonSpec'
extractParamTraces(spec, samples, L, ...)

# S4 method for class 'BinomialSpec'
extractParamTraces(spec, samples, L, ...)

# S4 method for class 'StudentTUvSpec'
extractParamTraces(spec, samples, L, ...)
```

## Arguments

- spec:

  A
  [`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md).

- samples:

  The MCMC sample matrix (iterations x monitored nodes).

- L:

  Integer truncation length.

- ...:

  Reserved for methods.

## Value

A named list of raw parameter traces (matrices/arrays) keyed by the
logical parameter names of the distribution.

## Functions

- `extractParamTraces(PoissonRegSpec)`: Parse coefficient traces.

- `extractParamTraces(BinomialRegSpec)`: Parse coefficient traces.

- `extractParamTraces(GMSNBurrUvSpec)`: Parse mu / sigma / alpha / theta
  traces.

- `extractParamTraces(NormalRegSpec)`: Parse betaTilde (L x p) and
  s2Tilde (L) traces; `prior$coefNames` (if present) labels the
  coefficients.

- `extractParamTraces(MSNBurrUvSpec)`: Parse mu / sigma / alpha traces.

- `extractParamTraces(MSNBurr2aUvSpec)`: Parse mu / sigma / alpha
  traces.

- `extractParamTraces(NormalMvRegSpec)`: Parse coefficient and
  covariance traces.

- `extractParamTraces(NormalMvSpec)`: Parse muTilde (L x d) and covTilde
  (L x d x d) traces into arrays.

- `extractParamTraces(NormalUvSpec)`: Parse muTilde / s2Tilde traces.

- `extractParamTraces(PoissonSpec)`: Parse lambda traces.

- `extractParamTraces(BinomialSpec)`: Parse prob traces.

- `extractParamTraces(StudentTUvSpec)`: Parse muTilde / tauTilde traces.
