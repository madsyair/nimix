# Skew multivariate independent-Student mixture components (Ferreira-Steel)

Heavy-tailed multivariate component: FS-skew Student-t margins with
per-dimension degrees of freedom \\\nu_j\\ (truncated below at 2),
per-dimension skewness `gamma`, and \\A = \mathrm{chol}(\Sigma)\\ upper
triangular (orthogonal factor fixed at the identity; see
[`SkewNormalMvSpec-class`](https://madsyair.github.io/nimix/reference/SkewNormalMvSpec-class.md)).
`gamma = 1` is the symmetric independent-Student; \\\nu \to \infty\\
recovers the skew multivariate Normal. Non-conjugate.

## Usage

``` r
SkewIStudentMvSpec()
```

## References

Ferreira & Steel (2007), Statistica Sinica 17, 505–529.
