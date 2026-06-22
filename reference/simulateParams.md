# Simulate component parameters from a prior (for inits / recovery tests)

Simulate component parameters from a prior (for inits / recovery tests)

## Usage

``` r
simulateParams(spec, prior, nClust, ...)

# S4 method for class 'NormalUvSpec'
simulateParams(spec, prior, nClust, ...)

# S4 method for class 'NormalMvSpec'
simulateParams(spec, prior, nClust, ...)

# S4 method for class 'NormalRegSpec'
simulateParams(spec, prior, nClust, ...)

# S4 method for class 'StudentTUvSpec'
simulateParams(spec, prior, nClust, ...)

# S4 method for class 'PoissonSpec'
simulateParams(spec, prior, nClust, ...)

# S4 method for class 'BinomialSpec'
simulateParams(spec, prior, nClust, ...)
```

## Arguments

- spec:

  A
  [`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md).

- prior:

  A prior list (typically from
  [`defaultPrior`](https://madsyair.github.io/nimix/reference/defaultPrior.md)).

- nClust:

  Integer number of components to simulate.

- ...:

  Reserved for methods.

## Value

A named list of simulated parameter vectors/matrices.

## Functions

- `simulateParams(NormalUvSpec)`: Draw (mu, s2) from the NIG prior.

- `simulateParams(NormalMvSpec)`: Draw (mu, Sigma) per cluster from the
  Normal-Inverse-Wishart prior. Returns `mu` (nClust x d) and `Sigma` (d
  x d x nClust).

- `simulateParams(NormalRegSpec)`: Draw (beta, s2) per cluster from the
  NIG prior. Returns `beta` (nClust x p) and `s2` (length nClust).

- `simulateParams(StudentTUvSpec)`: Draw (mu, tau) from the
  location/precision prior.

- `simulateParams(PoissonSpec)`: Draw lambda from the Gamma prior.

- `simulateParams(BinomialSpec)`: Draw prob from the Beta prior.
