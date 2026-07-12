# Ferreira-Steel skew multivariate independent-Student with estimated O

Density and random generation for the FS skew multivariate
independent-Student in which the orthogonal factor of \\A = OU\\ is
estimated through the Householder angle `theta`. Bivariate only: `theta`
lies in \\(-\pi/8, \pi/8)\\, exactly FS's restriction (8).

## Usage

``` r
dskewmvito(x, mu, Sigma, gamma, nu, theta, log = FALSE)

rskewmvito(n, mu, Sigma, gamma, nu, theta)
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

- nu:

  Positive numeric vector of per-dimension degrees of freedom.

- theta:

  Householder angle in \\(-\pi/8, \pi/8)\\.

- log:

  Logical; return the log-density?

- n:

  Integer number of draws.

## Value

`dskewmvito` a numeric vector; `rskewmvito` an `n x 2` matrix.

## References

Ferreira & Steel (2007), Statistica Sinica 17, 505–529.
