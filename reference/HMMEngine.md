# Construct a hidden-Markov mixture engine configuration

Construct a hidden-Markov mixture engine configuration

## Usage

``` r
HMMEngine(transConc = 1)
```

## Arguments

- transConc:

  Positive scalar concentration for the symmetric Dirichlet prior on
  each transition-matrix row. Default `1`.

## Value

An
[`HMMEngine`](https://madsyair.github.io/nimix/reference/HMMEngine-class.md)
object.

## Examples

``` r
eng <- HMMEngine()
eng
#> An object of class "HMMEngine"
#> Slot "transConc":
#> [1] 1
#> 
#> Slot "name":
#> [1] "hmm"
#> 
```
