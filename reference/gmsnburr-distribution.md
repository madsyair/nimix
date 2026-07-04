# GMSNBurr Distribution

Density, distribution function, quantile function, and random generation
for the Generalized MSNBurr (GMSNBurr) distribution.

## Usage

``` r
dgmsnburr(x, mu = 0, sigma = 1, alpha = 1, theta = 1, log = FALSE)

pgmsnburr(
  q,
  mu = 0,
  sigma = 1,
  alpha = 1,
  theta = 1,
  lower.tail = TRUE,
  log.p = FALSE
)

qgmsnburr(
  p,
  mu = 0,
  sigma = 1,
  alpha = 1,
  theta = 1,
  lower.tail = TRUE,
  log.p = FALSE
)

rgmsnburr(n, mu = 0, sigma = 1, alpha = 1, theta = 1)
```

## Arguments

- x, q:

  Numeric vector of quantiles.

- mu:

  Numeric. Location parameter (default = 0).

- sigma:

  Numeric. Scale parameter, must be positive (default = 1).

- alpha:

  Numeric. First shape parameter, must be positive (default = 1).

- theta:

  Numeric. Second shape parameter, must be positive (default = 1).

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

`dgmsnburr` gives the density, `pgmsnburr` gives the distribution
function, `qgmsnburr` gives the quantile function, and `rgmsnburr`
generates random deviates.

## Details

The GMSNBurr distribution with parameters \\\mu\\, \\\sigma\\,
\\\alpha\\, and \\\theta\\ has probability density function:

\$\$f(x\|\mu,\sigma,\alpha,\theta) =
\frac{\omega}{B(\alpha,\beta)\sigma}
\left(\frac{\theta}{\alpha}\right)^\theta
\exp\left(-\theta\omega\frac{x-\mu}{\sigma}\right) \left(1 +
\frac{\theta}{\alpha}
\exp\left(-\omega\frac{x-\mu}{\sigma}\right)\right)^{-(\alpha+\theta)}\$\$

Special cases:

- \\\alpha = \theta\\: Symmetric distribution

- \\\alpha \< \theta\\: Left-skewed

- \\\alpha \> \theta\\: Right-skewed

- \\\theta = 1\\: Reduces to MSNBurr

- \\\alpha = 1\\: Reduces to MSNBurr-IIa

- \\\alpha = \theta \to \infty\\: Converges to \\\mathcal{N}(\mu,
  \sigma^2)\\

## References

Choir, A. S. (2020). The New Neo-Normal Distributions and their
Properties. Dissertation. Institut Teknologi Sepuluh Nopember.

Iriawan, N. (2000). Computationally Intensive Approaches to Inference in
Neo-Normal Linear Models. PhD Thesis, Curtin University of Technology.

## See also

\[MSNBurr()\], \[MSNBurr2a()\]

## Examples

``` r
dgmsnburr(0, mu = 0, sigma = 1, alpha = 1, theta = 1)
#> [1] 0.3989423

x <- seq(-4, 4, by = 0.1)
plot(x, dgmsnburr(x, alpha = 1, theta = 1), type = "l",
     ylab = "Density", main = "GMSNBurr Densities")
lines(x, dgmsnburr(x, alpha = 2, theta = 1), col = "red")
lines(x, dgmsnburr(x, alpha = 1, theta = 2), col = "blue")


pgmsnburr(c(-2, 0, 2))
#> [1] 0.0394854 0.5000000 0.9605146
qgmsnburr(c(0.025, 0.5, 0.975))
#> [1] -2.295797  0.000000  2.295797

set.seed(123)
r <- rgmsnburr(1000)
hist(r, breaks = 30, freq = FALSE, main = "GMSNBurr Random Samples")
curve(dgmsnburr(x), add = TRUE, col = "red", lwd = 2)

```
