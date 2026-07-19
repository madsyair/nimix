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
#>  [4] "fossep-reg"          "fssn"                "fssn-reg"           
#>  [7] "fsst"                "fsst-reg"            "gmsnburr"           
#> [10] "gmsnburr-reg"        "jfst"                "jfst-reg"           
#> [13] "lep"                 "lep-reg"             "msnburr"            
#> [16] "msnburr-reg"         "msnburr2a"           "msnburr2a-reg"      
#> [19] "normal"              "normal-gamma"        "normal-gamma-mv"    
#> [22] "normal-gamma-mv-reg" "normal-gamma-reg"    "normal-mv"          
#> [25] "normal-mv-reg"       "normal-reg"          "normal-uv"          
#> [28] "poisson"             "poisson-reg"         "sep"                
#> [31] "sep-reg"             "skewistudent-mv"     "skewistudent-mv-o"  
#> [34] "skewistudent-mv-og"  "skewnormal-mv"       "skewnormal-mv-o"    
#> [37] "skewnormal-mv-og"    "student-t"           "student-t-mv"       
#> [40] "student-t-mv-reg"    "student-t-reg"      
```
