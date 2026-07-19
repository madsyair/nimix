# Spatial neighbourhood structure for spatially constrained mixtures

Represents the neighbourhood graph of `n` regions (or observations) as a
symmetric, zero-diagonal, non-negative weight matrix. It is the
structural ingredient of the spatially constrained mixture models
planned for the 1.x series, in which the latent component labels follow
a Markov random field on this graph rather than being independent across
observations (Besag 1974; spatially variant finite mixtures, Blekas et
al. 2005).

## Usage

``` r
# S4 method for class 'SpatialWeightSpec'
show(object)
```

## Arguments

- object:

  A `SpatialWeightSpec`.

## Details

A `SpatialWeightSpec` is intentionally independent of
[`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md):
the same neighbourhood structure can be paired with any registered
component distribution.

## Functions

- `show(SpatialWeightSpec)`: Compact printout: number of regions, edges,
  degree range, and weight type.

## Slots

- `edges`:

  Integer matrix, one row per undirected edge with endpoints `i < j`, in
  canonical column-major order (see
  [`getEdges`](https://madsyair.github.io/nimix/reference/getEdges.md)).
  This is the canonical representation since v1.5.0; the dense matrix is
  derived on demand and refused for large graphs.

- `edgeWeights`:

  Positive numeric weights, one per edge (1 for binary contiguity).

- `nNodes`:

  Integer number of regions.

- `regionIds`:

  Character vector of length `nNodes` naming the regions (defaults to
  row names of the matrix, or `"region1"`, ... when absent).

## References

Besag, J. (1974). Spatial interaction and the statistical analysis of
lattice systems. *Journal of the Royal Statistical Society B*, 36(2),
192–236.

Blekas, K., Likas, A., Galatsanos, N.P., & Lagaris, I.E. (2005). A
spatially constrained mixture model for image segmentation. *IEEE
Transactions on Neural Networks*, 16(2), 494–498.
[doi:10.1109/TNN.2004.841773](https://doi.org/10.1109/TNN.2004.841773)

Anselin, L. (1988). *Spatial Econometrics: Methods and Models*. Kluwer,
Dordrecht. (Queen/rook contiguity conventions, Ch. 3.)

## See also

[`spatialWeights`](https://madsyair.github.io/nimix/reference/spatialWeights.md),
[`gridAdjacency`](https://madsyair.github.io/nimix/reference/gridAdjacency.md),
[`neighborsOf`](https://madsyair.github.io/nimix/reference/neighborsOf.md)
