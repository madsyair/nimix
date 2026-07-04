# PSIS-LOO for a fitted mixture

Pareto-smoothed importance-sampling leave-one-out cross-validation
(Vehtari, Gelman & Gabry 2017) via the loo package, on the
label-invariant pointwise mixture log-likelihood. Requires loo; if it is
not installed, use
[`nimixWAIC`](https://madsyair.github.io/nimix/reference/nimixWAIC.md)
instead.

## Usage

``` r
nimixLOO(fit, maxDraws = 1000L)
```

## Arguments

- fit:

  A clustering
  [`FitResult`](https://madsyair.github.io/nimix/reference/FitResult-class.md).

- maxDraws:

  Cap on posterior draws used. Default 1000.

## Value

A `loo` object (see
[`loo::loo`](https://mc-stan.org/loo/reference/loo.html)); its
`estimates` carry `elpd_loo`, `p_loo`, and `looic`, and high Pareto-k
values flag observations where the approximation is unreliable.

## References

Vehtari, A., Gelman, A., & Gabry, J. (2017). Stat. Comput. 27,
1413–1432.
