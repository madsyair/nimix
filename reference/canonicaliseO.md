# Canonical representative of an FS orthogonal factor

Maps `(O, gamma, nu)` to the unique signed row permutation satisfying
Ferreira & Steel's identifiability restriction (8). Exactly one such
representative exists. `Sigma` is invariant under the map.

## Usage

``` r
canonicaliseO(O, gamma, nu = NULL)
```

## Arguments

- O:

  Orthogonal matrix with determinant \\(-1)^{m+1}\\.

- gamma:

  Positive skewness vector, length `m`.

- nu:

  Optional positive degrees-of-freedom vector, length `m`.

## Value

A list with the canonical `O`, `gamma`, and (if supplied) `nu`. If no
representative is found (a measure-zero event, e.g. a zero in the first
column of `O`), the inputs are returned unchanged with
`canonical = FALSE`.

## References

Ferreira & Steel (2007), Statistica Sinica 17, 505–529.
