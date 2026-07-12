# Fernandez-Osiewalski-Steel Skew Exponential Power Distribution

Density, distribution function, quantile function, and random generation
for the Fernandez-Osiewalski-Steel Skew Exponential Power (FOSSEP)
distribution.

## Usage

``` r
dfossep(x, mu = 0, sigma = 1, alpha = 2, theta = 2, log = FALSE)

pfossep(
  q,
  mu = 0,
  sigma = 1,
  alpha = 2,
  theta = 2,
  lower.tail = TRUE,
  log.p = FALSE
)

qfossep(
  p,
  mu = 0,
  sigma = 1,
  alpha = 2,
  theta = 2,
  lower.tail = TRUE,
  log.p = FALSE
)

rfossep(n, mu = 0, sigma = 1, alpha = 2, theta = 2)
```

## Arguments

- x, q:

  Numeric vector of quantiles.

- mu:

  Numeric. Location parameter (default = 0).

- sigma:

  Numeric. Scale parameter, must be positive (default = 1).

- alpha:

  Numeric. Skewness parameter, must be positive (default = 2). Follows
  the Fernandez-Steel convention shared by all skew families in nimix:
  `alpha` is the FS skewness \\\gamma\\, so \\P(X \> \mu) =
  \alpha^2/(1 + \alpha^2)\\ and `alpha > 1` skews right.

- theta:

  Numeric. Kurtosis/exponential-power parameter, must be positive
  (default = 2).

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

`dfossep` gives the density, `pfossep` gives the distribution function,
`qfossep` gives the quantile function, and `rfossep` generates random
deviates.

## Details

The Fernandez-Osiewalski-Steel Skew Exponential Power distribution with
parameters \\\mu\\, \\\sigma\\, \\\alpha\\, and \\\theta\\ has density:

\$\$f(y\|\mu,\sigma,\alpha,\theta) = \frac{c}{\sigma}
\exp\left(-\frac{1}{2} \|\alpha z\|^\theta\right) \quad\text{if } y \<
\mu\$\$

\$\$f(y\|\mu,\sigma,\alpha,\theta) = \frac{c}{\sigma}
\exp\left(-\frac{1}{2} \|z/\alpha\|^\theta\right) \quad\text{if } y \ge
\mu\$\$

where \\z = (y - \mu)/\sigma\\, \\c =
\frac{\alpha\theta}{(1+\alpha^2)\\2^{1/\theta}\\\Gamma(1/\theta)}\\.

When \\\theta = 2\\, it reduces to the Fernandez-Steel Skew Normal
(FSSN).

## References

Fernandez, C., Osiewalski, J., & Steel, M. F. (1995). Modeling and
inference with v-spherical distributions. Journal of the American
Statistical Association, 90(432), pp 1331-1340.

Rigby, R. A., Stasinopoulos, M. D., Heller, G. Z., & De Bastiani, F.
(2019). Distributions for Modeling Location, Scale, and Shape: Using
GAMLSS in R. CRC Press.

## Examples

``` r
dfossep(0, mu = 0, sigma = 1, alpha = 2, theta = 2)
#> [1] 0.3191538

x <- seq(-5, 5, by = 0.1)
plot(x, dfossep(x, alpha = 0.5, theta = 1.5), type = "l",
     ylab = "Density", main = "FOSSEP Densities")
lines(x, dfossep(x, alpha = 2, theta = 3), col = "red")


pfossep(c(-1, 0, 1), mu = 0, sigma = 1, alpha = 2, theta = 2)
#> [1] 0.009100053 0.200000000 0.506339938
qfossep(c(0.025, 0.5, 0.975), mu = 0, sigma = 1, alpha = 2, theta = 2)
#> [1] -0.7670603  0.9775528  4.3077494

set.seed(123)
r <- rfossep(1000, mu = 0, sigma = 1, alpha = 2, theta = 2)
hist(r, breaks = 30, freq = FALSE, main = "FOSSEP Random Samples")
curve(dfossep(x, alpha = 2, theta = 2), add = TRUE, col = "red", lwd = 2)

```
