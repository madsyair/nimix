# FSST mixture components (Fernandez-Steel skew Student-t)

Four-parameter component with location `mu`, scale `sigma`, skewness
`alpha` (`alpha = 1` symmetric) and degrees of freedom `nu`, truncated
below at 2 so the component variance exists. Heavy tails for small `nu`.
Non-conjugate. Note the skew-\\t\\ pitfalls of Fernandez & Steel (1999):
very small `nu` weakens identifiability, and `nu` is typically only
weakly identified by the data.

## Usage

``` r
FSSTUvSpec()
```

## References

Fernandez, C. & Steel, M. F. J. (1999). Multivariate Student-t
regression models: Pitfalls and inference. Biometrika 86, 153–167.
Choir, A. S. (2020). The New Neo-Normal Distributions and their
Properties.
