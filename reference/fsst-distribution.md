# Fernandez-Steel Skew t Distribution

Density, distribution function, quantile function, and random generation
for the Fernandez-Steel Skew t (FSST) distribution.

## Usage

``` r
dfsst(x, mu = 0, sigma = 1, alpha = 1, nu = 5, log = FALSE)

pfsst(
  q,
  mu = 0,
  sigma = 1,
  alpha = 1,
  nu = 5,
  lower.tail = TRUE,
  log.p = FALSE
)

qfsst(
  p,
  mu = 0,
  sigma = 1,
  alpha = 1,
  nu = 5,
  lower.tail = TRUE,
  log.p = FALSE
)

rfsst(n, mu = 0, sigma = 1, alpha = 1, nu = 5)
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
  = 1\\ gives the standard Student-t distribution with \\\nu\\ degrees
  of freedom.

- nu:

  Numeric. Degrees of freedom, must be positive (default = 5). As \\\nu
  \to \infty\\, the FSST approaches the FSSN distribution.

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

`dfsst` gives the density, `pfsst` gives the distribution function,
`qfsst` gives the quantile function, and `rfsst` generates random
deviates.

## Details

The Fernandez-Steel Skew t distribution applies the Fernandez-Steel
skewing mechanism to a Student-t kernel with \\\nu\\ degrees of freedom.

\$\$f(y\|\mu,\sigma,\alpha,\nu) = \frac{2}{\sigma(\alpha + 1/\alpha)}
t\_\nu(z/\alpha) \quad\text{if } y \< \mu\$\$

\$\$f(y\|\mu,\sigma,\alpha,\nu) = \frac{2}{\sigma(\alpha + 1/\alpha)}
t\_\nu(\alpha z) \quad\text{if } y \ge \mu\$\$

where \\z = (y - \mu)/\sigma\\ and \\t\_\nu\\ is the standard Student-t
density with \\\nu\\ degrees of freedom.

The mean exists for \\\nu \> 1\\ and the variance exists for \\\nu \>
2\\.

## References

Fernandez, C., Osiewalski, J., & Steel, M. F. (1995). Modeling and
inference with v-spherical distributions. Journal of the American
Statistical Association, 90(432), pp 1331-1340.

## Examples

``` r
dfsst(0, mu = 0, sigma = 1, alpha = 1, nu = 5)
#> [1] 0.3796067

x <- seq(-5, 5, by = 0.1)
plot(x, dfsst(x, alpha = 1, nu = 5), type = "l",
     ylab = "Density", main = "FSST Densities")
lines(x, dfsst(x, alpha = 2, nu = 3), col = "red")
lines(x, dfsst(x, alpha = 0.5, nu = 10), col = "blue")


pfsst(c(-1, 0, 1), mu = 0, sigma = 1, alpha = 2, nu = 5)
#> [1] 0.0203879 0.2000000 0.4893609
qfsst(c(0.025, 0.5, 0.975), mu = 0, sigma = 1, alpha = 2, nu = 5)
#> [1] -0.9204588  1.0407951  5.9348394

set.seed(123)
r <- rfsst(1000, mu = 0, sigma = 1, alpha = 2, nu = 5)
hist(r, breaks = 30, freq = FALSE, main = "FSST Random Samples")
curve(dfsst(x, alpha = 2, nu = 5), add = TRUE, col = "red", lwd = 2)

```
