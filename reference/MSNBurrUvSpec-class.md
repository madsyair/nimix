# MSNBurr mixture components (neo-normal, left-skew capable)

Univariate neo-normal component family with location `mu`, scale `sigma`
and shape `alpha` controlling skewness (Iriawan 2000; Choir 2020).
`alpha = 1` is the logistic distribution. Non-conjugate; NIMBLE assigns
adaptive samplers to the component parameters.

## Usage

``` r
MSNBurrUvSpec()
```

## References

Iriawan, N. (2000). Computationally Intensive Approaches to Inference in
Neo-Normal Linear Models. PhD Thesis, Curtin University of Technology.

Choir, A. S. (2020). The New Neo-Normal Distributions and their
Properties. Dissertation. Institut Teknologi Sepuluh Nopember.
