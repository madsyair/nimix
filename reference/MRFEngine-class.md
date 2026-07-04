# Markov random field engine (spatially constrained finite mixture)

Latent component labels follow a Potts model on the neighbourhood graph
of a
[`SpatialWeightSpec`](https://madsyair.github.io/nimix/reference/SpatialWeightSpec-class.md)
instead of being independent across observations: neighbouring regions
favour the same component, with fixed interaction strength `beta` (Potts
1952; Besag 1974; spatially variant finite mixtures, Blekas et al.
2005). `beta = 0` removes the spatial smoothing. Bayesian estimation of
`beta` (a hyperprior instead of a fixed value) is planned for a later
1.x release.

## Usage

``` r
MRFEngine(beta = 0.8, spatial, estimateBeta = FALSE, betaMax = 2)
```

## Arguments

- beta:

  Non-negative interaction strength (fixed value or start).

- spatial:

  A
  [`SpatialWeightSpec`](https://madsyair.github.io/nimix/reference/SpatialWeightSpec-class.md).

- estimateBeta:

  Logical; estimate beta by pseudo-likelihood Metropolis.

- betaMax:

  Upper bound of the uniform prior on beta.

## Slots

- `beta`:

  Non-negative spatial interaction strength: the fixed value when
  `estimateBeta = FALSE`, otherwise the chain's starting value.

- `spatial`:

  The
  [`SpatialWeightSpec`](https://madsyair.github.io/nimix/reference/SpatialWeightSpec-class.md)
  neighbourhood.

- `estimateBeta`:

  Logical; update `beta` by pseudo-likelihood Metropolis (Besag 1975)
  instead of holding it fixed.

- `betaMax`:

  Upper bound of the uniform prior on `beta`.

## References

Besag, J. (1974). Spatial interaction and the statistical analysis of
lattice systems. *JRSS B*, 36(2), 192–236.
