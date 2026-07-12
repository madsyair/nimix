# Ferreira-Steel skew multivariate independent-Student

Density and RNG for the FS skew multivariate independent-Student:
FS-skew Student-t margins with per-dimension `nu`, transformed by the
upper-triangular Cholesky factor of `Sigma`. `gamma = 1` gives the
symmetric independent-Student; `nu -> Inf` recovers
[`dskewmvn`](https://madsyair.github.io/nimix/reference/skewnormal-mv-distribution.md).

## Usage

``` r
dskewmvit(x, mu, Sigma, gamma, nu, log = FALSE)

rskewmvit(n, mu, Sigma, gamma, nu)
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

- log:

  Logical; return the log-density?

- n:

  Integer number of draws.

## Value

`dskewmvit` numeric vector; `rskewmvit` an `n x d` matrix.

## References

Ferreira & Steel (2007), Statistica Sinica 17, 505–529.
