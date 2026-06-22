# Construct a fixed-K finite-mixture engine configuration

Construct a fixed-K finite-mixture engine configuration

## Usage

``` r
FixedKEngine(dirichletConc = 1)
```

## Arguments

- dirichletConc:

  Positive scalar concentration of the symmetric Dirichlet prior on the
  mixing weights. Defaults to `1` (uniform on the simplex).

## Value

A
[`FixedKEngine`](https://madsyair.github.io/nimix/reference/FixedKEngine-class.md)
object.

## Examples

``` r
eng <- FixedKEngine(dirichletConc = 1)
eng
#> An object of class "FixedKEngine"
#> Slot "dirichletConc":
#> [1] 1
#> 
#> Slot "name":
#> [1] "fixedk"
#> 
```
