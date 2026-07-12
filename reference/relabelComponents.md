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

# S4 method for class 'FSSNUvSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'FOSSEPUvSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'FSSTUvSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'PoissonRegSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'BinomialRegSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'GMSNBurrUvSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'NormalRegSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'JFSTUvSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'SEPUvSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'LEPUvSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'MSNBurrUvSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'MSNBurr2aUvSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'NormalMvRegSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'NormalMvSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'NormalUvSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'PoissonSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'BinomialSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'SkewNormalMvSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'SkewNormalMvOSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'SkewIStudentMvSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'SkewIStudentMvOSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'SkewNormalMvOGenSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'SkewIStudentMvOGenSpec'
relabelComponents(spec, paramTrace, idx, occList, perms, modalK, weights, ...)

# S4 method for class 'StudentTUvSpec'
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

- `relabelComponents(FSSNUvSpec)`: Permute mu / sigma / alpha,
  summarise.

- `relabelComponents(FOSSEPUvSpec)`: Permute mu / sigma / alpha / theta,
  summarise.

- `relabelComponents(FSSTUvSpec)`: Permute mu / sigma / alpha / nu,
  summarise.

- `relabelComponents(PoissonRegSpec)`: Permute coefficients and
  summarise.

- `relabelComponents(BinomialRegSpec)`: Permute coefficients and
  summarise.

- `relabelComponents(GMSNBurrUvSpec)`: Permute mu / sigma / alpha /
  theta, summarise.

- `relabelComponents(NormalRegSpec)`: Permute (beta, s2) and summarise
  the regression components (one coefficient column per predictor).

- `relabelComponents(JFSTUvSpec)`: Permute mu / sigma / alpha / theta,
  summarise.

- `relabelComponents(SEPUvSpec)`: Permute mu / sigma / nu, summarise.

- `relabelComponents(LEPUvSpec)`: Permute mu / sigma / nu, summarise.

- `relabelComponents(MSNBurrUvSpec)`: Permute mu / sigma / alpha and
  summarise.

- `relabelComponents(MSNBurr2aUvSpec)`: Permute mu / sigma / alpha and
  summarise.

- `relabelComponents(NormalMvRegSpec)`: Permute coefficient matrices and
  summarise.

- `relabelComponents(NormalMvSpec)`: Permute multivariate (mu, Sigma)
  and summarise.

- `relabelComponents(NormalUvSpec)`: Permute (mu, s2) and summarise
  components.

- `relabelComponents(PoissonSpec)`: Permute lambda and summarise.

- `relabelComponents(BinomialSpec)`: Permute prob and summarise.

- `relabelComponents(SkewNormalMvSpec)`: Permute mv (mu, Sigma, gamma)
  and summarise.

- `relabelComponents(SkewNormalMvOSpec)`: Permute mv params plus theta,
  summarise.

- `relabelComponents(SkewIStudentMvSpec)`: Permute mv (mu, Sigma, gamma,
  nu), summarise.

- `relabelComponents(SkewIStudentMvOSpec)`: Permute mv params plus nu
  and theta.

- `relabelComponents(SkewNormalMvOGenSpec)`: Permute components, then
  canonicalise each draw's orthogonal factor via FS restriction (8),
  adjusting gamma with it.

- `relabelComponents(SkewIStudentMvOGenSpec)`: Permute components, then
  canonicalise each draw's orthogonal factor, carrying gamma and nu with
  the permutation.

- `relabelComponents(StudentTUvSpec)`: Permute (mu, tau) and summarise;
  the reported scale is \\\sigma = \tau^{-1/2}\\.
