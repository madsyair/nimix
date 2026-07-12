# Skew multivariate Normal components with estimated O, general dimension

As
[`SkewNormalMvOSpec-class`](https://madsyair.github.io/nimix/reference/SkewNormalMvOSpec-class.md)
but for any \\m \ge 3\\: the orthogonal factor of \\A = OU\\ is
parameterised by \\m(m-1)/2\\ Householder angles, sampled on the FS
angle box, and FS's identifiability restriction (8) is applied as a
post-hoc canonicalisation of the posterior draws (see
[`canonicaliseO`](https://madsyair.github.io/nimix/reference/canonicaliseO.md)).
Reached via `distribution = "skewnormal-mv-o"` when the data have more
than two columns.

## Usage

``` r
SkewNormalMvOGenSpec()
```

## References

Ferreira & Steel (2007), Statistica Sinica 17, 505–529.
