# nimix

**nimix** is a Bayesian mixture-modelling package built on top of
[NIMBLE](https://r-nimble.org/). It provides mixture *clustering*
(univariate and multivariate) and mixtures *of regressions* (univariate
and multivariate response) through an extensible S4 `DistributionSpec`
contract, so new component families plug into every inference engine
without touching the engine code.

## Features

- **Clustering** — univariate and multivariate Gaussian, plus Student-t,
  Normal-Gamma (heavy-tailed), Poisson and Binomial components.
- **Mixture of regressions** (`nimixReg`) — Gaussian linear, Student-t
  and Normal-Gamma residuals, Poisson/Binomial GLM components, and a
  multivariate (matrix) response via `cbind(y1, y2) ~ x`.
- **Two inference engines**, selected with `method`:
  - `"dpm"` — a Dirichlet Process Mixture that *infers* the number of
    components (NIMBLE’s Chinese Restaurant Process), using `K_max`.
  - `"fixedk"` — a finite mixture with a fixed, known `K` (including the
    single-component `K = 1` baseline for model comparison).
- Post-hoc label-switching correction via the
  [`label.switching`](https://CRAN.R-project.org/package=label.switching)
  package, with [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`predict()`](https://rdrr.io/r/stats/predict.html) methods.

## Installation

`nimix` depends on `nimble`, which compiles model code, so you need a
working C/C++ toolchain (Rtools on Windows, Xcode CLT on macOS,
build-essential on Linux).

``` r

# install.packages("pak")
pak::pak("madsyair/nimix")

# or
# install.packages("remotes")
remotes::install_github("madsyair/nimix")
```

## Quick start

``` r

library(nimix)

## Univariate clustering, number of components estimated (DPM)
set.seed(1)
y <- c(rnorm(100, -3), rnorm(100, 3))
fit <- nimixClust(y, K_max = 8,
                  mcmcControl = list(niter = 2000, nburnin = 1000),
                  verbose = FALSE)
summary(fit)
plot(fit, type = "K")

## Mixture of linear regressions, fixed number of regimes (finite mixture)
x  <- runif(200, -3, 3)
g  <- rep(1:2, each = 100)
df <- data.frame(y = ifelse(g == 1, 2 * x, -2 * x) + rnorm(200, 0, 0.7), x = x)
reg <- nimixReg(y ~ x, df, K = 2, method = "fixedk",
                mcmcControl = list(niter = 2000, nburnin = 1000),
                verbose = FALSE)
summary(reg)
```

## Documentation

Function reference and vignettes are published at
<https://madsyair.github.io/nimix/>.

## License

GPL-3.
