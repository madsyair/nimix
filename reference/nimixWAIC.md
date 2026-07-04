# WAIC for a fitted mixture

Widely Applicable Information Criterion (Watanabe 2010; Gelman, Hwang &
Vehtari 2014) from the label-invariant pointwise mixture log-likelihood.
Lower WAIC (equivalently higher `elpd_waic`) indicates better expected
out-of-sample predictive fit. Useful for choosing K or the component
family.

## Usage

``` r
nimixWAIC(fit, maxDraws = 1000L)
```

## Arguments

- fit:

  A clustering
  [`FitResult`](https://madsyair.github.io/nimix/reference/FitResult-class.md).

- maxDraws:

  Cap on posterior draws used (thinned evenly). Default 1000.

## Value

A list with `waic`, `elpd_waic`, `p_waic`, and `se` (standard error of
elpd).

## References

Watanabe, S. (2010). JMLR 11, 3571–3594. Gelman, A., Hwang, J., &
Vehtari, A. (2014). Stat. Comput. 24, 997–1016.
