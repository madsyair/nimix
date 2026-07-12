# Skew multivariate Normal mixture components (Ferreira-Steel)

Multivariate component family \\\eta = A^T \epsilon + \mu\\ with
independent Fernandez-Steel skew-Normal margins for `eps` and `A` the
upper-triangular Cholesky factor of `Sigma` (the orthogonal factor of FS
Lemma 1 is fixed at the identity; see the file header for what that
implies). `gamma = 1` in every dimension recovers the multivariate
Normal exactly. Non-conjugate.

## Usage

``` r
SkewNormalMvSpec()
```

## References

Ferreira, J. T. A. S. & Steel, M. F. J. (2007). Statistica Sinica 17,
505–529.
