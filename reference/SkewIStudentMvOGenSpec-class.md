# Skew mv independent-Student components with estimated O, general dimension

As
[`SkewIStudentMvOSpec-class`](https://madsyair.github.io/nimix/reference/SkewIStudentMvOSpec-class.md)
but for any \\m \ge 3\\, using the general Householder parameterisation
and the canonicalisation of
[`canonicaliseO`](https://madsyair.github.io/nimix/reference/canonicaliseO.md).
Reached via `distribution = "skewistudent-mv-o"` for data with more than
two columns. Experimental; see
[`SkewNormalMvOGenSpec-class`](https://madsyair.github.io/nimix/reference/SkewNormalMvOGenSpec-class.md)
for how to read `gamma` and `O` after canonicalisation.

## Usage

``` r
SkewIStudentMvOGenSpec()
```

## References

Ferreira & Steel (2007), Statistica Sinica 17, 505–529.
