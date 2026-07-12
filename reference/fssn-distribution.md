# Fernandez-Steel Skew Normal Distribution

Density, distribution function, quantile function, and random generation
for the Fernandez-Steel Skew Normal (FSSN) distribution.

## Usage

``` r
dfssn(x, mu = 0, sigma = 1, alpha = 1, log = FALSE)

pfssn(q, mu = 0, sigma = 1, alpha = 1, lower.tail = TRUE, log.p = FALSE)

qfssn(p, mu = 0, sigma = 1, alpha = 1, lower.tail = TRUE, log.p = FALSE)

rfssn(n, mu = 0, sigma = 1, alpha = 1)
```

## Arguments

- x, q:

  Numeric vector of quantiles.

- mu:

  Numeric. Location parameter (default = 0).

- sigma:

  Numeric. Scale parameter, must be positive (default = 1).

- alpha:

  Numeric. Skewness parameter, must be positive (default = 1,
  symmetric). Follows the Fernandez-Steel convention shared by all skew
  families in nimix: `alpha` is the FS skewness \\\gamma\\, so \\P(X \>
  \mu) = \alpha^2/(1 + \alpha^2)\\ and `alpha > 1` skews right. \\\alpha
  = 1\\ gives the standard normal distribution.

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

`dfssn` gives the density, `pfssn` gives the distribution function,
`qfssn` gives the quantile function, and `rfssn` generates random
deviates.

## Details

The Fernandez-Steel Skew Normal distribution is a special case of the
FOSSEP distribution with \\\theta = 2\\. It applies the Fernandez-Steel
skewing mechanism to a standard normal kernel.

\$\$f(y\|\mu,\sigma,\alpha) = \frac{2}{\sigma(\alpha + 1/\alpha)}
\phi(z/\alpha) \quad\text{if } y \< \mu\$\$

\$\$f(y\|\mu,\sigma,\alpha) = \frac{2}{\sigma(\alpha + 1/\alpha)}
\phi(\alpha z) \quad\text{if } y \ge \mu\$\$

where \\z = (y - \mu)/\sigma\\ and \\\phi\\ is the standard normal
density.

## References

Fernandez, C., Osiewalski, J., & Steel, M. F. (1995). Modeling and
inference with v-spherical distributions. Journal of the American
Statistical Association, 90(432), pp 1331-1340.

## Examples

``` r
dfssn(0, mu = 0, sigma = 1, alpha = 1)
#> [1] 0.3989423

x <- seq(-5, 5, by = 0.1)
plot(x, dfssn(x, alpha = 0.5), type = "l",
     ylab = "Density", main = "FSSN Densities")
lines(x, dfssn(x, alpha = 1), col = "red")
lines(x, dfssn(x, alpha = 2), col = "blue")


pfssn(c(-1, 0, 1), mu = 0, sigma = 1, alpha = 2)
#> [1] 0.009100053 0.200000000 0.506339938
qfssn(c(0.025, 0.5, 0.975), mu = 0, sigma = 1, alpha = 2)
#> [1] -0.7670603  0.9775528  4.3077494

set.seed(123)
r <- rfssn(1000, mu = 0, sigma = 1, alpha = 2)
hist(r, breaks = 30, freq = FALSE, main = "FSSN Random Samples")
curve(dfssn(x, alpha = 2), add = TRUE, col = "red", lwd = 2)

```
