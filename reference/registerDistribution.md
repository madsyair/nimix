# Register a component distribution

Adds a
[`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md)
to the registry under its `name` slot so it can be selected by name. New
built-in distributions (Student-t, Poisson/Binomial) are planned for
v0.4.0.

## Usage

``` r
registerDistribution(spec, overwrite = FALSE)
```

## Arguments

- spec:

  A
  [`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md)
  instance.

- overwrite:

  Logical; overwrite an existing entry of the same name?

## Value

Invisibly, the registered name.

## Examples

``` r
registerDistribution(NormalUvSpec(), overwrite = TRUE)
listDistributions()
#>  [1] "binomial"            "binomial-reg"        "normal"             
#>  [4] "normal-gamma"        "normal-gamma-mv"     "normal-gamma-mv-reg"
#>  [7] "normal-gamma-reg"    "normal-mv"           "normal-mv-reg"      
#> [10] "normal-reg"          "normal-uv"           "poisson"            
#> [13] "poisson-reg"         "student-t"           "student-t-mv"       
#> [16] "student-t-mv-reg"    "student-t-reg"      
```
