# Bayesian mixture of linear regressions

Fit a mixture of Gaussian linear regressions. Each component is \\y \sim
N(x^\top \beta_k, \sigma^2_k)\\ with a conjugate Normal-Inverse-Gamma
cluster prior. The number of components can be inferred with a Dirichlet
Process Mixture (`method = "dpm"`, using `K_max`) or fixed
(`method = "fixedk"`, using `K`). Gating is constant: the mixing weights
do not depend on the covariates.

## Usage

``` r
nimixReg(
  formula,
  data,
  random = NULL,
  K = NULL,
  K_max = NULL,
  distribution = "normal",
  method = c("dpm", "fixedk", "mrf", "hmm"),
  gating = c("constant", "covariate"),
  prior = list(),
  mcmcControl = list(),
  initMethod = c("kmeans", "single", "spread"),
  seed = 1L,
  verbose = FALSE,
  spatialWeights = NULL
)
```

## Arguments

- formula:

  A model formula, e.g. `y ~ x1 + x2`.

- data:

  A data frame containing the formula variables.

- random:

  Optional one-sided formula for group-level random effects.
  `random = ~ region` gives a random intercept: every component's linear
  predictor gains a shared group offset \\b\_{g(i)} \sim N(0, \tau^2)\\.
  `random = ~ x | region` additionally gives a random slope for `x`
  (which must be a term of `formula`); as in `lme4`'s `(x|g)`, the
  intercept comes along. Both offsets use a sum-to-zero constraint – the
  parameterisation that mixes well – so the component coefficients
  absorb the group means and the reported `b` / `sRE` are centred. The
  two offsets are given independent priors by default. `tauRE` and
  `tauSlope` estimate the spread of *your* groups, which with few groups
  is itself variable. Currently supported with `method = "fixedk"` and
  `distribution = "normal"` or `"studentt"` (heavy-tailed residuals).

- K:

  Integer number of components for the finite mixture
  (`method = "fixedk"`). Required for that method.

- K_max:

  Integer truncation level for the DPM (`method = "dpm"`); it must sit
  comfortably above the expected number of clusters, since the dCRP
  sampler errors if the occupied-cluster count ever needs to exceed it.
  A generous data-aware default (with headroom) is used when missing.

- distribution:

  Component distribution. Gaussian (`"normal"`), heavy-tailed
  (`"studentt"`, `"normalgamma"`), GLM (`"poisson"`, `"binomial"`), and
  the nine neo-normal skew families (`"msnburr"`, `"msnburr2a"`,
  `"gmsnburr"`, `"fssn"`, `"fsst"`, `"sep"`, `"lep"`, `"fossep"`,
  `"jfst"`) are available; multivariate variants exist for the Gaussian
  and heavy-tailed families.

  **Non-Gaussian families cost effective sample size.** Only the
  Gaussian regression has a conjugate coefficient update; every other
  family (neo-normal, heavy-tailed, GLM) is sampled by NIMBLE's
  defaults, which mix more slowly. Measured on `"fixedk"` with two
  components, the Gaussian reaches about 10.7 ESS/s against 1.9 for
  `"msnburr"` and 1.4 for the three-shape `"gmsnburr"` – roughly a 5–7x
  factor. Budget proportionally more iterations for a skewed or
  heavy-tailed fit; the three-shape families (gmsnburr, fsst, fossep,
  jfst) are the slowest.

- method:

  Fitting method: `"dpm"`, `"fixedk"`, `"mrf"`, or `"hmm"`.

  `"hmm"` fits a **Markov-switching regression** (Hamilton 1989): the
  coefficients and error variance switch with a latent first-order
  Markov regime, so the rows of `data` are a *time series* rather than
  an exchangeable sample – their order carries meaning here as it does
  nowhere else in `nimixReg`. Give `K` for the number of regimes, as for
  `"fixedk"`; the regime path is marginalised out of the likelihood and
  decoded afterwards, so
  [`viterbiPath`](https://madsyair.github.io/nimix/reference/viterbiPath.md)
  gives the most probable regime sequence. Currently
  `distribution = "normal"` (Gaussian) or `"poisson"` (log-link counts),
  or `"studentt"` / `"normalgamma"` (heavy-tailed), or `"binomial"`
  (logit-link proportions, with the number of trials in
  `prior = list(size = )`).

  **Budget more iterations than for `"fixedk"`.** Marginalising the
  regime path rules out the conjugate Normal-Inverse-Gamma update, so
  the coefficients are sampled by NIMBLE's defaults instead. Measured on
  a two-regime series, that costs roughly four times the effective
  sample size per second (ESS/s 2.3 against 8.9 for the conjugate
  `"fixedk"` sampler), even though the marginalised chain has a shorter
  wall time. A light run gives usably-decoded regimes but wide
  coefficient intervals; raise `niter` until the intervals settle.

- gating:

  Mixing-weight model: `"constant"` (default; weights do not depend on
  covariates). `"covariate"` (concomitant gating) is a planned opt-in
  and currently errors.

- prior:

  A named list of prior overrides passed to
  [`defaultPrior`](https://madsyair.github.io/nimix/reference/defaultPrior.md)
  (e.g. `g` for the g-prior factor, `nu0` for the InvGamma shape, and
  `s2Guess` for its scale – see ‘Reading the error variance’ below; with
  a multivariate response, `sigmaGuess` plays the same role for the
  residual covariance) plus, for the DPM, optional
  `concPrior = c(shape, rate)`, or, for the finite mixture,
  `dirichletConc`.

- mcmcControl:

  A named list of MCMC controls: `niter`, `nburnin`, `thin`, and the
  optional `initRatio` – the fraction of the truncation / component cap
  (`K_max` or `K`) seeded by the dispersed cluster initialisation
  (default 0.8; must lie in (0, 1)). Lower it to leave more headroom
  below the truncation; raising it to 0.95 or above is allowed but
  warns, as it leaves little headroom.

- initMethod:

  Initialisation: `"kmeans"` (default) or `"single"`.

- seed:

  Integer RNG seed.

- verbose:

  Logical; print NIMBLE's configuration and progress output. Defaults to
  `FALSE` (quiet): NIMBLE's compilation notes and the benign dCRP
  truncation note are silenced, while nimix's own diagnostics (e.g. a
  censored-posterior warning) and any error still surface. Set `TRUE` to
  see NIMBLE's configuration and a progress bar.

- spatialWeights:

  Optional
  [`SpatialWeightSpec`](https://madsyair.github.io/nimix/reference/SpatialWeightSpec-class.md)
  (one region per observation). Required by, and only used with,
  `method = "mrf"`.

## Value

A
[`FitResult`](https://madsyair.github.io/nimix/reference/FitResult-class.md).
[`summary()`](https://rdrr.io/r/base/summary.html) reports relabelled
per-component regression coefficients and residual variances;
`predict(fit, newdata)` returns the posterior predictive mean;
`plot(fit, type = "fitted")` shows observed vs fitted.

## Reading the error variance

The prior on each component's error variance `s2` is centred on the
residual variance of a *global* OLS fit, which ignores the mixture. For
separated components that quantity measures the spread *between*
components as much as the spread within one, so the prior is
deliberately conservative and `s2` is biased upward. The bias is largest
exactly where mixtures are most useful – well-separated components at
moderate sample size – and it disappears as the per-component sample
size grows: with a prior/truth scale ratio of about 60, the measured
bias runs near 9x at 25 observations per component, 2.5x at 150, 1.2x at
1000, and 1.0x by 5000. It is a bias in the safe direction (wider
predictive intervals), and the coefficients are unaffected.

The conservatism is load-bearing rather than incidental: it also
regularises against splitting when `K` is over-specified, so it stays
the default.

If you know the within-component scale – from a pilot study, the
literature, or domain knowledge – say so directly:
`prior = list(s2Guess = 0.3)` makes 0.3 the prior mean of `s2`. (`s0`
sets the raw InvGamma scale instead, for callers who think in those
terms; give one or the other.) On a simulated benchmark with a true `s2`
of 0.25 that moved the estimate from 0.78 (3.1x) to 0.36 (1.4x), slopes
untouched. The override is deliberately absolute: a multiplier on the
automatic scale would only ever mean "a fraction of a quantity that
measures the wrong thing", which is a knob, not a statement.

**Where a tighter prior is safe.** A tenfold tighter scale left the
DPM's recovery of `K` untouched (modal `K` = 2 throughout a benchmark
with two true components), and is likewise safe for `fixedk` when `K` is
correct, where it halved the error in `s2`. It is **not** safe for
`fixedk` with an over-specified `K`: on the same data with `K = 4`
against two true components, the default occupied two components while a
tenfold tighter prior occupied three. Wide components make two suffice;
narrow ones make the model want more. If you are unsure whether `K` is
over-specified, leave the default alone or use the DPM.

**Recipe: two-stage empirical Bayes.** If you need the scale right and
have no external knowledge, fit once with the default, read the
within-component residual variance off that fit, and refit with it as
`s2Guess`. It costs a second compile-and-run, which is why it is not the
default, but it needs no guesswork:


    fit1 <- nimixReg(y ~ x, df, K = 2, method = "fixedk")
    fit1 <- relabel(fit1)
    # residuals of each point under its own MAP component
    z    <- binderPartition(fit1)$partition
    cf   <- fit1@relabeled$summary
    Xm   <- model.matrix(y ~ x, df)
    res  <- df$y - rowSums(Xm * as.matrix(cf[z, c("(Intercept)", "x")]))
    s2hat <- sum(res^2) / (nrow(df) - 2 * ncol(Xm))

    fit2 <- nimixReg(y ~ x, df, K = 2, method = "fixedk",
                     prior = list(s2Guess = s2hat))

The first fit's *allocation* is reliable even where its `s2` is not –
that is what makes the recipe work.

## References

Hurn, M., Justel, A., & Robert, C.P. (2003). Estimating mixtures of
regressions. *JCGS*, 12(1), 55–79.
[doi:10.1198/1061860031329](https://doi.org/10.1198/1061860031329)

Grün, B., & Leisch, F. (2008). FlexMix version 2. *JSS*, 28(4), 1–35.
[doi:10.18637/jss.v028.i04](https://doi.org/10.18637/jss.v028.i04)

## Examples

``` r
# \donttest{
set.seed(1)
x <- runif(200, -3, 3)
grp <- rep(1:2, each = 100)
y <- ifelse(grp == 1, 2 * x, -2 * x) + rnorm(200, 0, 0.7)
df <- data.frame(y = y, x = x)

## number of regimes estimated (DPM)
fit <- nimixReg(y ~ x, df, K_max = 8,
                mcmcControl = list(niter = 2000, nburnin = 1000),
                verbose = FALSE)
summary(fit)
#> Relabelling MCMC output before summarising (label switching)...
#> nimix mixture summary (engine: dpm, distribution: normal-reg)
#> Observations: 200 (dimension d = 1)
#> Relabelling: ECR-ITERATIVE-1 conditioned on modal K = 2 (834 draws)
#> 
#> Posterior of number of occupied clusters:
#> 
#>     2     3     4     5 
#> 0.834 0.141 0.024 0.001 
#> 
#> Relabelled component estimates (posterior mean; CIs for univariate):
#>  component weight (Intercept)     x s2_mean s2_med s2_lwr s2_upr
#>          1  0.498    -0.06331  1.97   0.993   0.97  0.749   1.33
#>          2  0.502     0.00544 -1.97   1.082   1.06  0.820   1.48
#> 
#> Mixing diagnostic (single chain): ESS(#clusters) = 421, ESS(alpha) = 832
#>   Set mcmcControl$nchains > 1 for cross-chain split-Rhat.
predict(fit, newdata = data.frame(x = c(-2, 0, 2)))
#>    x     .fitted
#> 1 -2 -0.01717332
#> 2  0 -0.03112922
#> 3  2 -0.04508512

## fixed number of regimes (finite mixture)
fit2 <- nimixReg(y ~ x, df, K = 2, method = "fixedk",
                 mcmcControl = list(niter = 2000, nburnin = 1000),
                 verbose = FALSE)
summary(fit2)
#> Relabelling MCMC output before summarising (label switching)...
#> Warning: Low effective sample size for the cluster count (ESS = 0 of 1000 draws): the chain may be mixing poorly across partitions. Consider a longer run or k-means initialisation.
#> nimix mixture summary (engine: fixedk, distribution: normal-reg)
#> Observations: 200 (dimension d = 1)
#> Relabelling: ECR-ITERATIVE-1 conditioned on modal K = 2 (1000 draws)
#> 
#> Posterior of number of occupied clusters:
#> 
#> 2 
#> 1 
#> 
#> Relabelled component estimates (posterior mean; CIs for univariate):
#>  component weight (Intercept)     x s2_mean s2_med s2_lwr s2_upr
#>          1  0.496    -0.07053  1.98    1.00  0.988  0.753   1.32
#>          2  0.504    -0.00105 -1.96    1.06  1.055  0.789   1.43
#> 
#> Mixing diagnostic (single chain): ESS(#clusters) = 0
#>   Set mcmcControl$nchains > 1 for cross-chain split-Rhat.
# }
```
