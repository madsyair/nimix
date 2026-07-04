# Ensemble several fitted mixtures

Combines several fitted mixtures into a single weighted predictive
model, rather than selecting one winner. Weights come from Bayesian
stacking or Pseudo-BMA+ (Yao et al. 2018) via loo when available, or
from WAIC (Akaike-style weights) natively. Stacking is the default and
is the most robust to model misspecification. Models must share the same
data.

## Usage

``` r
ensembleFit(..., method = c("stacking", "pseudobma", "waic"), maxDraws = 1000L)
```

## Arguments

- ...:

  Two or more clustering
  [`FitResult`](https://madsyair.github.io/nimix/reference/FitResult-class.md)
  objects, or a single named list of them.

- method:

  Weighting scheme: `"stacking"` or `"pseudobma"` (both need loo), or
  `"waic"` (native, no dependency).

- maxDraws:

  Cap on posterior draws used per fit. Default 1000.

## Value

A `nimixEnsemble` object carrying the fits and their weights, with a
`predict` method for the weighted predictive density.

## References

Yao, Y., Vehtari, A., Simpson, D., & Gelman, A. (2018). Bayesian
Analysis 13(3), 917–1007.
