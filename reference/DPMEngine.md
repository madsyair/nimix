# Construct a DPM engine configuration

Construct a DPM engine configuration

## Usage

``` r
DPMEngine(concPrior = c(2, 4))
```

## Arguments

- concPrior:

  Length-2 numeric `c(shape, rate)` for the Gamma prior on the DP
  concentration \\\alpha\\. Defaults to `c(2, 4)` (weakly informative,
  prior mean 0.5).

## Value

A
[`DPMEngine`](https://madsyair.github.io/nimix/reference/DPMEngine-class.md)
object.

## Examples

``` r
eng <- DPMEngine(concPrior = c(2, 4))
eng
#> An object of class "DPMEngine"
#> Slot "concPrior":
#> [1] 2 4
#> 
#> Slot "name":
#> [1] "dpm"
#> 
```
