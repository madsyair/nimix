# MSNBurr Distribution

Density, distribution function, quantile function, and random generation
for the MSNBurr (Modified to be Stable as Normal from Burr)
distribution.

## Usage

``` r
dmsnburr(x, mu = 0, sigma = 1, alpha = 1, log = FALSE)

pmsnburr(q, mu = 0, sigma = 1, alpha = 1, lower.tail = TRUE, log.p = FALSE)

qmsnburr(p, mu = 0, sigma = 1, alpha = 1, lower.tail = TRUE, log.p = FALSE)

rmsnburr(n, mu = 0, sigma = 1, alpha = 1)
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

`dmsnburr` gives the density, `pmsnburr` gives the distribution
function, `qmsnburr` gives the quantile function, and `rmsnburr`
generates random deviates.

## Details

The MSNBurr distribution with parameters \\\mu\\, \\\sigma\\, and
\\\alpha\\ has probability density function:

\$\$f(x\|\mu,\sigma,\alpha) = \frac{\omega}{\sigma}
\exp\left(-\omega\frac{x-\mu}{\sigma}\right) \left(1 +
\frac{1}{\alpha}\exp\left(-\omega\frac{x-\mu}{\sigma}\right)\right)^{-(\alpha+1)}\$\$

where \\-\infty \< x \< \infty\\, \\-\infty \< \mu \< \infty\\, \\\sigma
\> 0\\, \\\alpha \> 0\\, and \\\omega\\ is the normalizing constant:

\$\$\omega = \frac{1}{\sqrt{2\pi}} \left(1 +
\frac{1}{\alpha}\right)^{\alpha+1}\$\$

The MSNBurr is a special case of the GMSNBurr distribution with \\\theta
= 1\\. It is left-skewed (negative skewness).

## References

Iriawan, N. (2000). Computationally Intensive Approaches to Inference in
Neo-Normal Linear Models. PhD Thesis, Curtin University of Technology.

Choir, A. S. (2020). The New Neo-Normal Distributions and their
Properties. Dissertation. Institut Teknologi Sepuluh Nopember.

## See also

`dmsnburr2a`

## Examples

``` r
dmsnburr(0, mu = 0, sigma = 1, alpha = 1)
#> [1] 0.3989423

x <- seq(-10, 3, by = 0.1)
plot(x, dmsnburr(x, alpha = 0.1), type = "l",
     ylab = "Density", main = "MSNBurr Densities")
lines(x, dmsnburr(x, alpha = 0.5), col = "red")
lines(x, dmsnburr(x, alpha = 1), col = "blue")
lines(x, dmsnburr(x, alpha = 5), col = "forestgreen")


pmsnburr(c(-2, 0, 2), mu = 0, sigma = 1, alpha = 1)
#> [1] 0.0394854 0.5000000 0.9605146
qmsnburr(c(0.025, 0.5, 0.975), alpha = 1)
#> [1] -2.295797  0.000000  2.295797

set.seed(123)
r <- rmsnburr(1000)
hist(r, breaks = 30, freq = FALSE, main = "MSNBurr Random Samples")
curve(dmsnburr(x), add = TRUE, col = "red", lwd = 2)

```
