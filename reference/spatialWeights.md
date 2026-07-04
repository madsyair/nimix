# Construct a SpatialWeightSpec from an adjacency/weight matrix

Construct a SpatialWeightSpec from an adjacency/weight matrix

## Usage

``` r
spatialWeights(adjacency, regionIds = NULL)
```

## Arguments

- adjacency:

  Numeric n x n matrix: symmetric, zero diagonal, non-negative. Binary
  0/1 contiguity is the typical case.

- regionIds:

  Optional character vector of n unique region names. Defaults to the
  matrix row names, or `"region1"`, ....

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
