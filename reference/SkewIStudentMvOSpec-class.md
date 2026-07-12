# Skew mv independent-Student components with estimated orthogonal factor

As
[`SkewIStudentMvSpec-class`](https://madsyair.github.io/nimix/reference/SkewIStudentMvSpec-class.md),
but the orthogonal factor of \\A = OU\\ is estimated through the
Householder angle `theta` with a uniform prior on \\(-\pi/8, \pi/8)\\
(FS restriction (8)). Bivariate data only. `theta` is identified only
through the skewness.

## Usage

``` r
SkewIStudentMvOSpec()
```

## References

Ferreira & Steel (2007), Statistica Sinica 17, 505–529.
