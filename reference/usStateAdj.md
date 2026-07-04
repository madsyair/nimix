# Contiguity of the 48 contiguous US states + DC (official derivation)

Symmetric binary adjacency matrix (49 x 49, postal-code dimnames)
between the contiguous United States and the District of Columbia,
derived from the U.S. Census Bureau's 2023 county adjacency file: two
states are adjacent when any pair of their counties is listed as
adjacent across the state line. The Census county file's conventions
(including some water adjacencies) are inherited as-is.

## Usage

``` r
usStateAdj
```

## Format

A 49 x 49 numeric 0/1 matrix; rows/columns ordered by state FIPS and
named by postal code. Use `spatialWeights(usStateAdj)` to obtain a
[`SpatialWeightSpec`](https://madsyair.github.io/nimix/reference/SpatialWeightSpec-class.md).

## Source

U.S. Census Bureau, county adjacency file,
`https://www2.census.gov/geo/docs/reference/county_adjacency/`

## Details

Derived on 2026-07-03 from `county_adjacency2023.txt`; 112 undirected
edges; e.g. Tennessee has the well-known 8 neighbours and Maine exactly
1.
