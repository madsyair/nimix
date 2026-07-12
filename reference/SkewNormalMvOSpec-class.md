# Skew multivariate Normal components with estimated orthogonal factor

As
[`SkewNormalMvSpec-class`](https://madsyair.github.io/nimix/reference/SkewNormalMvSpec-class.md),
but the orthogonal factor of \\A = OU\\ is estimated through the
Householder angle `theta` with a uniform prior on \\(-\pi/8, \pi/8)\\
(FS restriction (8)). Bivariate data only. `theta` is identified only
through the skewness: at `gamma = 1` the density is invariant in
`theta`.

## Usage

``` r
SkewNormalMvOSpec()
```

## References

Ferreira & Steel (2007), Statistica Sinica 17, 505–529.
