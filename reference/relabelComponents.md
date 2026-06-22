# Permute cluster parameters and build the relabelled component summary

Called by
[`relabel`](https://madsyair.github.io/nimix/reference/relabel.md) after
the label-permutation matrix has been derived from the allocation
vectors (which is distribution-independent). The spec is responsible
only for permuting its own parameters and producing a tidy per-component
summary, so multivariate covariance handling stays inside
[`NormalMvSpec`](https://madsyair.github.io/nimix/reference/NormalMvSpec-class.md).

## Usage

``` r
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'NormalUvSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'NormalMvSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'NormalRegSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'PoissonRegSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'BinomialRegSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'NormalMvRegSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'StudentTUvSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'PoissonSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'BinomialSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)
```

## Arguments

- spec:

  A
  [`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md).

- paramTrace:

  The raw parameter trace list (from `extractParamTraces`).

- idx:

  Integer indices of retained (modal-K) iterations.

- occList:

  A list (length = number of retained iterations) of the sorted occupied
  cluster labels at each retained iteration.

- perms:

  The permutation matrix (retained iterations x modalK).

- modalK:

  Integer modal number of occupied clusters.

- weights:

  A retained-iterations x modalK matrix of already-permuted mixing
  weights.

- ...:

  Reserved for methods.

## Value

A named list with at least `summary` (a data.frame) plus permuted
parameter arrays cached on the `FitResult`.

## Functions

- `relabelComponents(NormalUvSpec)`: Permute (mu, s2) and summarise
  components.

- `relabelComponents(NormalMvSpec)`: Permute multivariate (mu, Sigma)
  and summarise.

- `relabelComponents(NormalRegSpec)`: Permute (beta, s2) and summarise
  the regression components (one coefficient column per predictor).

- `relabelComponents(PoissonRegSpec)`: Permute coefficients and
  summarise.

- `relabelComponents(BinomialRegSpec)`: Permute coefficients and
  summarise.

- `relabelComponents(NormalMvRegSpec)`: Permute coefficient matrices and
  summarise.

- `relabelComponents(StudentTUvSpec)`: Permute (mu, tau) and summarise;
  the reported scale is \\\sigma = \tau^{-1/2}\\.

- `relabelComponents(PoissonSpec)`: Permute lambda and summarise.

- `relabelComponents(BinomialSpec)`: Permute prob and summarise.
