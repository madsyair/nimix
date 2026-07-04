# MSNBurr-IIa Distribution

Density, distribution function, quantile function, and random generation
for the MSNBurr-IIa distribution.

## Usage

``` r
dmsnburr2a(x, mu = 0, sigma = 1, alpha = 1, log = FALSE)

pmsnburr2a(q, mu = 0, sigma = 1, alpha = 1, lower.tail = TRUE, log.p = FALSE)

qmsnburr2a(p, mu = 0, sigma = 1, alpha = 1, lower.tail = TRUE, log.p = FALSE)

rmsnburr2a(n, mu = 0, sigma = 1, alpha = 1)
```

## Arguments

- x, q:

  Numeric vector of quantiles.

- mu:

  Numeric. Location parameter (default = 0).

- sigma:

  Numeric. Scale parameter, must be positive (default = 1).

- alpha:

  Numeric. Shape parameter, must be positive (default = 1).

- log, log.p:

  Logical. If TRUE, probabilities/densities are given as log. Default is
  FALSE.

- lower.tail:

  Logical. If TRUE (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- p:

  Numeric vector of probabilities.

- n:

  Integer. Number of observations to generate.

## Value

`dmsnburr2a` gives the density, `pmsnburr2a` gives the distribution
function, `qmsnburr2a` gives the quantile function, and `rmsnburr2a`
generates random deviates.

## Details

The MSNBurr-IIa distribution with parameters \\\mu\\, \\\sigma\\, and
\\\alpha\\ has probability density function:

\$\$f(x\|\mu,\sigma,\alpha) = \frac{\omega}{\sigma}
\exp\left(\omega\frac{x-\mu}{\sigma}\right) \left(1 +
\frac{1}{\alpha}\exp\left(\omega\frac{x-\mu}{\sigma}\right)\right)^{-(\alpha+1)}\$\$

where \\\omega\\ is the same normalizing constant as MSNBurr. The
MSNBurr-IIa is a special case of the GMSNBurr distribution with
\\\alpha\_{\text{gmsnburr}} = 1\\ and \\\beta\_{\text{gmsnburr}} =
\alpha\\. It is right-skewed (positive skewness).

## References

Choir, A. S. (2020). The New Neo-Normal Distributions and their
Properties. Dissertation. Institut Teknologi Sepuluh Nopember.

## See also

`dmsnburr`

## Examples

``` r
dmsnburr2a(0, mu = 0, sigma = 1, alpha = 2)
#> [1] 0.3989423

x <- seq(-3, 10, by = 0.1)
plot(x, dmsnburr2a(x, alpha = 0.1), type = "l",
     ylab = "Density", main = "MSNBurr-IIa Densities")
lines(x, dmsnburr2a(x, alpha = 0.5), col = "red")
lines(x, dmsnburr2a(x, alpha = 1), col = "blue")
lines(x, dmsnburr2a(x, alpha = 5), col = "forestgreen")


pmsnburr2a(c(-2, 0, 2), mu = 0, sigma = 1, alpha = 1)
#> [1] 0.0394854 0.5000000 0.9605146
qmsnburr2a(c(0.025, 0.5, 0.975), alpha = 1)
#> [1] -2.295797  0.000000  2.295797

set.seed(123)
r <- rmsnburr2a(1000)
hist(r, breaks = 30, freq = FALSE, main = "MSNBurr-IIa Random Samples")
curve(dmsnburr2a(x), add = TRUE, col = "red", lwd = 2)

```
