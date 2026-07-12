# Posterior draws as an iterations x chains x parameters array

Returns a plain 3-D array in the layout `bayesplot`'s `mcmc_*` functions
accept natively (`iterations x chains x parameters`), so no extra
packages are required to hand a nimix fit to
[`bayesplot::mcmc_trace()`](https://mc-stan.org/bayesplot/reference/MCMC-traces.html),
`mcmc_rhat_hist()` and friends.

## Usage

``` r
drawsArray(fit, params = c("invariant", "components"))
```

## Arguments

- fit:

  A `FitResult`.

- params:

  `"invariant"` or `"components"`.

## Value

A numeric array `iterations x chains x parameters` with dimnames on the
parameter margin.

## Details

Two views are available, and the distinction is statistical, not
cosmetic:

- `"invariant"` (default):

  Label-invariant functionals – the number of occupied clusters, the
  allocation entropy, and (for DPM fits) the concentration parameter
  `alpha`. These are meaningful on raw draws and retain the per-chain
  structure, so cross-chain diagnostics like R-hat apply.

- `"components"`:

  Per-component parameters *after*
  [`relabel`](https://madsyair.github.io/nimix/reference/relabel.md).
  Refused if
  [`relabel()`](https://madsyair.github.io/nimix/reference/relabel.md)
  has not been run: under label switching, `muTilde[1]` names different
  components in different chains, and an R-hat computed on it looks
  valid while meaning nothing. Because relabelling conditions on the
  modal cluster count, chains lose equal lengths; the draws are
  therefore pooled into a single chain, suitable for posterior
  density/interval plots but not for cross-chain R-hat.

## See also

[`ppcData`](https://madsyair.github.io/nimix/reference/ppcData.md),
[`relabel`](https://madsyair.github.io/nimix/reference/relabel.md),
[`psm`](https://madsyair.github.io/nimix/reference/psm.md).
