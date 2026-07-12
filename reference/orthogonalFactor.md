# Orthogonal factor from Householder angles

Builds \\O = O\_{\theta^m} \cdots O\_{\theta^2}\\ from a flat vector of
\\m(m-1)/2\\ Householder angles (Ferreira & Steel 2007, Appendix A).

## Usage

``` r
orthogonalFactor(theta, m)
```

## Arguments

- theta:

  Numeric vector of `m * (m - 1) / 2` angles.

- m:

  Dimension.

## Value

An `m x m` orthogonal matrix with determinant \\(-1)^{m+1}\\.

## References

Ferreira & Steel (2007), Statistica Sinica 17, 505–529.
