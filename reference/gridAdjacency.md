# Rook/queen contiguity on a regular grid

Builds the
[`SpatialWeightSpec`](https://madsyair.github.io/nimix/reference/SpatialWeightSpec-class.md)
of a `nrow` x `ncol` regular lattice under rook (shared edge) or queen
(shared edge or corner) contiguity – the standard conventions of spatial
econometrics (Anselin 1988, Ch. 3). Regular grids with known block
structure are also the prescribed synthetic-graph setting for the
spatial recovery tests of the MRF engine planned for v0.6.0.

## Usage

``` r
gridAdjacency(nrow, ncol, contiguity = c("rook", "queen"))
```

## Arguments

- nrow, ncol:

  Grid dimensions (each \>= 1).

- contiguity:

  `"rook"` (default) or `"queen"`.

## Value

A
[`SpatialWeightSpec`](https://madsyair.github.io/nimix/reference/SpatialWeightSpec-class.md)
with regions named `"r<i>c<j>"` in row-major order.

## Examples

``` r
g <- gridAdjacency(3, 3, "queen")
neighborsOf(g, "r2c2")   # interior cell: 8 queen neighbours
#> [1] "r1c1" "r1c2" "r1c3" "r2c1" "r2c3" "r3c1" "r3c2" "r3c3"
```
