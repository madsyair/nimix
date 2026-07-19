# Fitted nimix mixture result

Fitted nimix mixture result

## Usage

``` r
# S4 method for class 'FitResult'
show(object)

# S4 method for class 'FitResult'
summary(object, ...)

# S4 method for class 'FitResult,missing'
plot(
  x,
  y,
  type = c("K", "trace_raw", "trace_relabeled", "density", "cluster", "fitted", "series",
    "regime", "forecast"),
  ...
)

# S4 method for class 'FitResult'
predict(object, newdata, maxDraws = 500L, ...)
```

## Arguments

- object:

  A `FitResult`.

- ...:

  Passed to the underlying plot. For `type = "forecast"`, also `h`,
  `newdata`, `lags`, `draws` and `level`, as in
  [`nimixForecast`](https://madsyair.github.io/nimix/reference/nimixForecast.md).

- x:

  A `FitResult`.

- y:

  Ignored.

- type:

  Which view to draw. `"K"`, `"trace_raw"`, `"trace_relabeled"`,
  `"density"`, `"cluster"` and `"fitted"` are the mixture views.

  `"series"`, `"regime"` and `"forecast"` are the **time-series views**,
  for `method = "hmm"` only. A density plot of a regime-switching series
  is the wrong picture twice over: it discards the axis the model is
  about, and it shows a bimodal smear where the story is "the series sat
  in one regime, then moved". `"series"` draws the data with its decoded
  regimes shaded as blocks (blocks, because that is what the Markov
  chain models); `"regime"` draws the smoothed \\P(\mathrm{regime} \mid
  \mathrm{data})\\ through time, which is the honest companion to the
  Viterbi path – where the bands are mixed, the decode is a guess; and
  `"forecast"` draws the predictive fan past the end of the series (see
  [`nimixForecast`](https://madsyair.github.io/nimix/reference/nimixForecast.md)
  for what that fan can and cannot mean).

- newdata:

  Points at which to evaluate: a numeric vector (univariate clustering),
  a matrix with `d` columns (multivariate clustering), or a data frame
  of predictors (regression). Defaults to the training data.

- maxDraws:

  Integer cap on the number of posterior draws used (default 500); draws
  are thinned uniformly if exceeded.

## Value

Invisibly, a tidy data frame of exactly what was drawn (e.g.
`iteration`/`component`/`value` for traces, `x`/`density` for the
predictive density), so the plot can be reproduced with ggplot2 or any
other system without nimix depending on them.

## Methods (by generic)

- `show(FitResult)`: Compact display of a fitted result.

- `summary(FitResult)`: Posterior summary (relabelled component
  estimates, posterior of the number of clusters, and a mixing
  diagnostic).

- `plot(x = FitResult, y = missing)`: Diagnostic and result plots.

- `predict(FitResult)`: Posterior predictive density at new points.

  Returns the posterior predictive density averaged over MCMC draws,
  using the occupied-cluster mixture in each draw (weights = cluster
  sizes / n). This quantity is label-invariant, so no relabelling is
  required. For multivariate fits the density is evaluated at the
  supplied rows of `newdata`. For a regression fit (`nimixReg`) it
  instead returns the posterior predictive *mean* \\E\[y \mid x\]\\ per
  row of `newdata` (column `.fitted`). To keep evaluation tractable the
  draws are subsampled to at most `maxDraws`.

## Slots

- `mcmcSamples`:

  A matrix of monitored MCMC draws (iterations x parameters).

- `Kposterior`:

  Integer vector: number of occupied clusters per iteration.

- `clusterAllocation`:

  Integer matrix (iterations x n) of raw cluster indicators `xi`.

- `paramTrace`:

  A named list of raw cluster-parameter traces (each an iterations x
  Kmax matrix).

- `engineUsed`:

  Character scalar naming the engine.

- `distSpec`:

  The
  [`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md)
  used.

- `data`:

  The data the model was fit to.

- `Kmax`:

  Integer truncation level.

- `prior`:

  The prior list used.

- `relabeled`:

  A list cache populated by
  [`relabel`](https://madsyair.github.io/nimix/reference/relabel.md) (or
  empty).

- `mcmcControl`:

  The MCMC control list actually used.

- `diagnostics`:

  Multi-chain convergence diagnostics (Rhat, ESS) when more than one
  chain is run; otherwise a single-chain summary.

- `call`:

  The matched call.
