# Numerically stable log(exp(a) + exp(b))

Computes log(exp(a) + exp(b)) avoiding overflow. Equivalent to
\\\max(a,b) + \log(1 + \exp(-\|a - b\|))\\.

## Usage

``` r
log_sum_exp(a, b)
```

## Arguments

- a:

  Numeric vector.

- b:

  Numeric vector.

## Value

Numeric vector of the same length.
