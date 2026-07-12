# Ferreira-Steel skew multivariate Normal

Density and random generation for the Ferreira-Steel skew multivariate
Normal with location `mu`, scatter `Sigma` (its upper-triangular
Cholesky factor is the FS transformation `A`), and per-dimension
skewness `gamma` in the shared Fernandez-Steel convention (`gamma = 1`
symmetric; `gamma > 1` skews right).

## Usage

``` r
dskewmvn(x, mu, Sigma, gamma, log = FALSE)

rskewmvn(n, mu, Sigma, gamma)
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

- log:

  Logical; return the log-density?

- n:

  Integer number of draws.

## Value

`dskewmvn` a numeric vector of (log-)densities; `rskewmvn` an `n x d`
matrix of draws.

## References

Ferreira, J. T. A. S. & Steel, M. F. J. (2007). A new class of skewed
multivariate distributions with applications to regression analysis.
Statistica Sinica 17, 505–529.
