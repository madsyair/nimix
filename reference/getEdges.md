# Edge list of a spatial weight structure

The canonical sparse form: one row per undirected edge, endpoints as
integer node indices with `i < j`, ordered by second endpoint then first
(the column-major order of `which(upper.tri(A) & A > 0)` on the dense
form, so downstream constants are identical either way). Unlike
[`getAdjacency`](https://madsyair.github.io/nimix/reference/getAdjacency.md)
this never allocates an \\n \times n\\ matrix, so it works at any graph
size.

## Usage

``` r
getEdges(spec)

# S4 method for class 'SpatialWeightSpec'
getEdges(spec)
```

## Arguments

- spec:

  A
  [`SpatialWeightSpec`](https://madsyair.github.io/nimix/reference/SpatialWeightSpec-class.md).

## Value

Integer matrix with two columns.
