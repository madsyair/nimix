# Numerically Stable log(1 - exp(x))

Computes log(1 - exp(x)) for x \< 0. Used in upper-tail probability
computations.

## Usage

``` r
log1mexp(x)
```

## Arguments

- x:

  Numeric vector, must be \< 0.

## Value

Numeric vector of the same length as \`x\`.
