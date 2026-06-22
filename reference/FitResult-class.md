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
  type = c("K", "trace_raw", "trace_relabeled", "density", "cluster", "fitted"),
  ...
)

# S4 method for class 'FitResult'
predict(object, newdata, maxDraws = 500L, ...)
```

## Arguments

- object:

  A `FitResult`.

- ...:

  Passed to
  [`relabel`](https://madsyair.github.io/nimix/reference/relabel.md)
  when relabelling is triggered.

- x:

  A `FitResult`.

- y:

  Ignored.

- type:

  One of `"K"` (posterior of \#clusters), `"trace_raw"` (raw
  cluster-parameter traces; zig-zags reveal label switching),
  `"trace_relabeled"` (traces after relabelling), `"density"`
  (univariate clustering only: data histogram with posterior predictive
  overlay), `"cluster"` (multivariate clustering, \\d \ge 2\\: scatter
  coloured by MAP cluster), or `"fitted"` (regression only: observed
  response vs posterior predictive mean).

- newdata:

  Points at which to evaluate: a numeric vector (univariate clustering),
  a matrix with `d` columns (multivariate clustering), or a data frame
  of predictors (regression). Defaults to the training data.

- maxDraws:

  Integer cap on the number of posterior draws used (default 500); draws
  are thinned uniformly if exceeded.

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

- `call`:

  The matched call.
