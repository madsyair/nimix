# Lunetta Exponential Power Distribution

Density, distribution function, quantile function, and random generation
for the Lunetta Exponential Power (LEP) distribution.

## Usage

``` r
dlep(x, mu = 0, sigma = 1, nu = 2, log = FALSE)

plep(q, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE)

qlep(p, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE)

rlep(n, mu = 0, sigma = 1, nu = 2)
```

## Arguments

- x, q:

  Numeric vector of quantiles.

- mu:

  Numeric. Location parameter (default = 0).

- sigma:

  Numeric. Scale parameter, must be positive (default = 1).

- nu:

  Numeric. Shape parameter, must be positive (default = 2).

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

`dlep` gives the density, `plep` gives the distribution function, `qlep`
gives the quantile function, and `rlep` generates random deviates.

## Details

The LEP distribution with parameters \\\mu\\, \\\sigma\\, and \\\nu\\
has density:

\$\$f(x\|\mu,\sigma,\nu) = \frac{1}{2\nu^{1/\nu} \Gamma(1+1/\nu) \sigma}
\exp\left(-\frac{\|x-\mu\|^\nu}{\nu \sigma^\nu}\right)\$\$

where \\-\infty \< x \< \infty\\, \\-\infty \< \mu \< \infty\\, \\\sigma
\> 0\\, \\\nu \> 0\\.

The LEP is symmetric around \\\mu\\. When \\\nu = 2\\, it reduces to the
normal distribution. When \\\nu = 1\\, it becomes the Laplace (double
exponential) distribution.

## References

Lunetta, G. (1963). Di una generalizzazione dello schema della curva
normale. *Annali della Facolta di Economia e Commercio di Palermo*, 17,
237-244.

## Examples

``` r
dlep(0, mu = 0, sigma = 1, nu = 2)
#> [1] 0.3989423

x <- seq(-5, 5, by = 0.1)
plot(x, dlep(x, nu = 2), type = "l",
     ylab = "Density", main = "LEP Densities")
lines(x, dlep(x, nu = 1), col = "red")
lines(x, dlep(x, nu = 4), col = "blue")
legend("topright", c("nu=2 (Normal)", "nu=1 (Laplace)", "nu=4"),
       col = c("black", "red", "blue"), lty = 1)


plep(c(-2, 0, 2), mu = 0, sigma = 1, nu = 2)
#> [1] 0.02275013 0.50000000 0.97724987
qlep(c(0.025, 0.5, 0.975), mu = 0, sigma = 1, nu = 2)
#> [1] -1.959964  0.000000  1.959964

set.seed(123)
r <- rlep(1000, mu = 0, sigma = 1, nu = 2)
hist(r, breaks = 30, freq = FALSE, main = "LEP Random Samples")
curve(dlep(x, nu = 2), add = TRUE, col = "red", lwd = 2)

```
