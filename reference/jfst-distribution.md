# Jones-Faddy Skew-t Distribution

Density, distribution function, quantile function, and random generation
for the Jones-Faddy Skew-t distribution.

## Usage

``` r
djfst(x, mu = 0, sigma = 1, alpha = 3, theta = 3, log = FALSE)

pjfst(
  q,
  mu = 0,
  sigma = 1,
  alpha = 3,
  theta = 3,
  lower.tail = TRUE,
  log.p = FALSE
)

qjfst(
  p,
  mu = 0,
  sigma = 1,
  alpha = 3,
  theta = 3,
  lower.tail = TRUE,
  log.p = FALSE
)

rjfst(n, mu = 0, sigma = 1, alpha = 3, theta = 3)
```

## Arguments

- x, q:

  Numeric vector of quantiles.

- mu:

  Numeric. Location parameter (default = 0).

- sigma:

  Numeric. Scale parameter, must be positive (default = 1).

- alpha:

  Numeric. Left tail parameter (\> 0), the Jones-Faddy \\a\\ parameter
  (default = 3). \\\alpha = \theta\\ gives symmetry.

- theta:

  Numeric. Right tail parameter (\> 0), the Jones-Faddy \\b\\ parameter
  (default = 3). \\\alpha = \theta\\ gives symmetry.

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

`djfst` gives the density, `pjfst` gives the distribution function,
`qjfst` gives the quantile function, and `rjfst` generates random
deviates.

## Details

The Jones-Faddy Skew-t distribution with parameters \\\mu\\, \\\sigma\\,
\\\alpha\\ (left tail), and \\\theta\\ (right tail) has density:

\$\$f(y\|\mu,\sigma,\alpha,\theta) = \frac{c}{\sigma} \left\[1 +
\frac{z}{\sqrt{\alpha+\theta+z^2}}\right\]^{\alpha+0.5} \left\[1 -
\frac{z}{\sqrt{\alpha+\theta+z^2}}\right\]^{\theta+0.5}\$\$

where \\z = (y-\mu)/\sigma\\, \\c = \left\[2^{\alpha+\theta-1}
\sqrt{\alpha+\theta} \\ B(\alpha,\theta)\right\]^{-1}\\.

\\\alpha \< \theta\\: left-skewed (heavier right tail). \\\alpha \>
\theta\\: right-skewed (heavier left tail). \\\alpha = \theta\\:
Student-\$t\$ distribution (symmetric, \\\nu = 2\alpha\\).

When \\\alpha = \theta\\, JFST reduces to the Student-\$t\$ distribution
with \\\nu = 2\alpha\\ degrees of freedom. As \\\alpha = \theta \to
\infty\\, JFST approaches the Normal distribution.

Moments exist only when certain conditions on \\\alpha\\ and \\\theta\\
are met: mean requires \\\alpha,\theta \> 0.5\\, variance requires
\\\alpha,\theta \> 1\\, skewness requires \\\alpha,\theta \> 1.5\\,
kurtosis requires \\\alpha,\theta \> 2\\.

## References

Jones, M.C. and Faddy, M.J. (2003). A skew extension of the t
distribution, with applications. Journal of the Royal Statistical
Society, Series B, 65, pp 159-174.

## Examples

``` r
djfst(0, mu = 0, sigma = 1, alpha = 3, theta = 3)
#> [1] 0.3827328

x <- seq(-5, 5, by = 0.1)
plot(x, djfst(x, alpha = 3, theta = 3), type = "l",
     ylab = "Density", main = "JFST Densities")
lines(x, djfst(x, alpha = 5, theta = 1), col = "red")
lines(x, djfst(x, alpha = 1, theta = 5), col = "blue")


pjfst(c(-2, 0, 2), mu = 0, sigma = 1, alpha = 3, theta = 3)
#> [1] 0.04621316 0.50000000 0.95378684
qjfst(c(0.025, 0.5, 0.975), mu = 0, sigma = 1, alpha = 3, theta = 3)
#> [1] -2.446912  0.000000  2.446912

set.seed(123)
r <- rjfst(1000, mu = 0, sigma = 1, alpha = 3, theta = 3)
hist(r, breaks = 30, freq = FALSE, main = "JFST Random Samples")
curve(djfst(x, alpha = 3, theta = 3), add = TRUE, col = "red", lwd = 2)

```
