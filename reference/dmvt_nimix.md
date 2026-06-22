# Multivariate Student-t log density (nimbleFunction)

User-defined NIMBLE distribution for a location-scale multivariate-t
with location `mu`, scale matrix `cov`, and degrees of freedom `df`.
Registered with NIMBLE when the package loads.

## Usage

``` r
dmvt_nimix(x, mu, cov, df, log = 0)

rmvt_nimix(n, mu, cov, df)
```

## Arguments

- x, mu:

  Numeric vectors (observation, location).

- cov:

  Scale matrix.

- df:

  Degrees of freedom.

- log:

  Return the log density?

- n:

  Number of draws (always 1 for `rmvt_nimix`).

## Value

A density (or log density), or a single draw.
