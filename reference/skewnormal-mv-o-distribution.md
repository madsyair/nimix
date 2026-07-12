# Ferreira-Steel skew multivariate Normal with estimated orthogonal factor

Density and random generation for the FS skew multivariate Normal in
which the orthogonal factor of \\A = OU\\ is estimated via the
Householder angle `theta`. Only `m = 2` is supported: `theta` must lie
in \\(-\pi/8, \pi/8)\\, which is exactly FS's identifiability
restriction (8).

## Usage

``` r
dskewmvno(x, mu, Sigma, gamma, theta, log = FALSE)

rskewmvno(n, mu, Sigma, gamma, theta)
```

## Arguments

- x:

  Numeric vector of length `d`, or an `n x d` matrix.

- mu:

  Numeric location vector (length `d`).

- Sigma:

  Positive-definite `d x d` scatter matrix.

- gamma:

  Positive numeric vector of per-dimension FS skewness parameters
  (length `d`).

- theta:

  Householder angle in \\(-\pi/8, \pi/8)\\.

- log:

  Logical; return the log-density?

- n:

  Integer number of draws.

## Value

`dskewmvno` a numeric vector; `rskewmvno` an `n x 2` matrix.

## References

Ferreira & Steel (2007), Statistica Sinica 17, 505–529.
