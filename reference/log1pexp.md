# Numerically Stable log(1 + exp(x))

Computes log(1 + exp(x)) avoiding overflow for large positive x and
underflow for large negative x.

## Usage

``` r
log1pexp(x)
```

## Arguments

- x:

  Numeric vector.

## Value

Numeric vector of the same length as \`x\`.

## Details

Piecewise approximation for numerical stability:

- \\x \ge 33\\: returns \\x\\ (above this, \\e^x\\ overflows double
  precision; \\1 + e^x \approx e^x\\)

- \\x \le -37\\: returns \\e^x\\ (below this, \\\log(1 + \epsilon)
  \approx \epsilon\\ and \\1 + e^x\\ evaluates to exactly 1 in double
  precision)

- otherwise: returns
  [`log1p`](https://rdrr.io/r/base/Log.html)`(`[`exp`](https://rdrr.io/r/base/Log.html)`(x))`
