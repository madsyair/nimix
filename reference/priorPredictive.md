# Prior predictive check for a mixture model

Simulates `nsim` datasets from the prior implied by `distribution`, `K`,
and the data-scaled default prior (optionally overridden), then compares
summary statistics of the observed data with their prior predictive
distributions. Run this *before* fitting: if the observed statistics
fall far outside the prior predictive range, the prior scale is off for
these data.

## Usage

``` r
priorPredictive(
  data,
  K,
  distribution = "normal",
  nsim = 200L,
  prior = NULL,
  conc = 1,
  seed = NULL
)
```

## Arguments

- data:

  Numeric vector of observations (univariate; the multivariate families
  are not yet covered).

- K:

  Number of mixture components to simulate under.

- distribution:

  Component family name, as in
  [`nimixClust`](https://madsyair.github.io/nimix/reference/nimixClust.md).

- nsim:

  Number of prior predictive datasets. Default 200.

- prior:

  Optional named list overriding entries of the data-scaled default
  prior (as in `nimixClust`).

- conc:

  Dirichlet concentration for the mixture weights. Default 1 (uniform on
  the simplex), matching the fixed-K engine default.

- seed:

  Optional RNG seed.

## Value

An object of class `nimixPriorPred`: a list with the observed statistics
(`obs`), the matrix of simulated statistics (`sim`, `nsim` rows), the
tail probability of each observed statistic under the prior predictive
(`pTail`, two-sided), and the simulated datasets' summary. Its `print`
method flags statistics with `pTail < 0.05`, and its `plot` method
overlays prior predictive densities on the observed data.

## References

Gelman, A., et al. (2020). Bayesian workflow. *arXiv:2011.01808*.

## Examples

``` r
if (FALSE) { # \dontrun{
y <- c(rnorm(80, -2), rnorm(120, 3))
pp <- priorPredictive(y, K = 2, distribution = "normal")
pp          # flags any statistic the prior cannot reach
plot(pp)    # observed density over prior predictive draws
} # }
```
