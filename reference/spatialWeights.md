# Construct a SpatialWeightSpec from an adjacency/weight matrix

Construct a SpatialWeightSpec from an adjacency/weight matrix

## Usage

``` r
spatialWeights(
  adjacency = NULL,
  regionIds = NULL,
  edges = NULL,
  nNodes = NULL,
  edgeWeights = NULL
)
```

## Arguments

- adjacency:

  Numeric n x n matrix: symmetric, zero diagonal, non-negative. Binary
  0/1 contiguity is the typical case. Give either this *or* `edges`, not
  both.

- regionIds:

  Optional character vector of n unique region names. Defaults to the
  matrix row names, or `"region1"`, ....

- edges:

  Alternative sparse input: a two-column matrix of node indices, one row
  per undirected edge (order and duplicates are normalised). Never
  allocates an \\n \times n\\ matrix, so it is the route for large
  graphs. Requires `nNodes`.

- nNodes:

  Number of nodes when constructing from `edges`.

- edgeWeights:

  Optional positive weights, one per row of `edges`; defaults to 1
  (binary contiguity).

## Value

A validated
[`SpatialWeightSpec`](https://madsyair.github.io/nimix/reference/SpatialWeightSpec-class.md).

## Examples

``` r
A <- matrix(0, 3, 3); A[1, 2] <- A[2, 1] <- 1; A[2, 3] <- A[3, 2] <- 1
sw <- spatialWeights(A, regionIds = c("west", "centre", "east"))
nRegions(sw)
#> [1] 3
neighborsOf(sw, "centre")
#> [1] "west" "east"
```
