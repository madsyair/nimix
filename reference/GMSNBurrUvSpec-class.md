# GMSNBurr mixture components (generalized neo-normal)

Univariate neo-normal component family with location `mu`, scale
`sigma`, and two shape parameters `alpha`, `theta` governing skewness
(Iriawan 2000; Choir 2020). `theta = 1` is MSNBurr, `alpha = 1` is
MSNBurr-IIa, and `alpha = theta` is symmetric. Non-conjugate; NIMBLE
assigns adaptive samplers to the component parameters.

## Usage

``` r
GMSNBurrUvSpec()
```

## References

Choir, A. S. (2020). The New Neo-Normal Distributions and their
Properties. Dissertation. Institut Teknologi Sepuluh Nopember.
