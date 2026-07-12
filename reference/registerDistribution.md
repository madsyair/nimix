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
#>  [1] "binomial"            "binomial-reg"        "fossep"             
#>  [4] "fssn"                "fsst"                "gmsnburr"           
#>  [7] "jfst"                "lep"                 "msnburr"            
#> [10] "msnburr2a"           "normal"              "normal-gamma"       
#> [13] "normal-gamma-mv"     "normal-gamma-mv-reg" "normal-gamma-reg"   
#> [16] "normal-mv"           "normal-mv-reg"       "normal-reg"         
#> [19] "normal-uv"           "poisson"             "poisson-reg"        
#> [22] "sep"                 "skewistudent-mv"     "skewistudent-mv-o"  
#> [25] "skewistudent-mv-og"  "skewnormal-mv"       "skewnormal-mv-o"    
#> [28] "skewnormal-mv-og"    "student-t"           "student-t-mv"       
#> [31] "student-t-mv-reg"    "student-t-reg"      
```
