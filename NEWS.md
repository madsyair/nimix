# nimix 1.5.0

## ESS benchmark for regression across engines

* Documented the effective-sample-size cost of the non-Gaussian regressions.
  Only the Gaussian has a conjugate coefficient update; measured on fixed-K with
  two components, the Gaussian reaches ~10.7 ESS/s against 1.9 for MSNBurr and
  1.4 for the three-shape GMSNBurr (a 5-7x factor). HMM is slower still per
  second than fixed-K (regime path marginalised, coefficients more correlated),
  and the DPM mixes well on its concentration parameter. `?nimixReg` now states
  the factor so users can budget iterations; a full benchmark write-up ships
  with the project notes.

## Random effects extend to the DPM engine (and a note on HMM)

* Neo-normal (and Gaussian) regressions accept random effects under the DPM
  engine too: `nimixReg(y ~ x, method = "dpm", distribution = "msnburr", random = ~ grp)`.
  The random-cluster allocation and the fixed external grouping both act on iid
  exchangeable observations, so the offsets sit cleanly on the location, as they
  do under fixed-K. Measured: msnburr cor(b_hat, truth) = 0.96, sep 1.00, gmsnburr 1.00.
* Random effects remain refused under the HMM engine, now with a clear reason: a
  single time series has no exchangeable grouping, so the offset is unidentified
  and confounds with the regime transitions. Panel-HMM (multiple series) would be
  the prerequisite. See the design note for the full argument.

## Random effects for GLM regression -- the genuine GLMM case

* Poisson and Binomial regression now accept group-level random effects:
  `nimixReg(y ~ x, distribution = "poisson", random = ~ grp)`, and `~ x | grp`
  for a random slope. Unlike the location-scale families, a GLM random effect
  enters INSIDE the link (log or logit), so the group offset is multiplicative
  on the response scale -- this is the true GLMM the design study flagged. tauRE
  therefore lives on the link scale, so its prior bounds are fixed rather than
  scaled by sd(y) (a large count would otherwise give an absurd ceiling).
  Measured on six groups: Poisson cor(b_hat, truth) = 0.96, Binomial 0.99;
  Poisson random slope 0.93/0.99. Fixed-K engine.

## Random effects for MSNBurr regression (F10 prototype)

* MSNBurr regression accepts group-level random effects, like the Gaussian and
  Student-t regressions: `nimixReg(y ~ x, distribution = "msnburr", random = ~ grp)`
  for a random intercept, `random = ~ x | grp` to add a random slope. Because
  neo-normal families are location-scale, the random effect enters the location
  exactly as it does for a Gaussian (unlike a GLM, where it would act through
  the link), with the same sum-to-zero parameterisation. Measured on six groups:
  random-intercept recovery cor(b_hat, truth) = 0.999; with a random slope too,
  0.98 and 1.00. Fixed-K engine. All nine neo-normal families now support random
  effects through the same generic path (.neoRegREConstants + the re/reSlope
  arguments to .neoRegFixedKCode); measured recovery cor(b_hat, truth) >= 0.98
  across families.

## Neo-normal regression completes the engine trio with HMM

* All nine neo-normal regression families now run under the HMM engine too --
  `nimixReg(y ~ x, method = "hmm", distribution = "msnburr")` and the rest --
  giving Markov-switching skewed / heavy-tailed regression. The per-state
  location is a linear predictor, with scale and shape switching by regime.
  Measured on two regimes: MSNBurr intercepts -1.49/1.53 against -1.5/1.5.
* The forward FFBS kernel, its BUGSdist registration, and the three S4 methods
  (model code, allocation density, forecast RNG) are all generated from two
  facts -- the family's density kernel and its shape list -- via
  `.makeNeoHMMRegKernel` / `.makeNeoHMMRegMethods`. No hand-written
  nimbleFunction per family. Together with the fixed-K and DPM variants, every
  neo-normal family now spans the full engine trio.

## Neo-normal regression now runs under the DPM (nonparametric) engine

* Every neo-normal regression family gains a Dirichlet-process variant --
  `nimixReg(y ~ x, method = "dpm", distribution = "msnburr")` and the other
  eight -- with the number of components inferred rather than fixed. The DPM
  model code is generated generically from the family's density name and its
  shape-prior lines (the same ingredients the fixed-K code uses), so no
  per-family DPM code was written. No conjugate cluster wrapper (the likelihood
  is not Gaussian in the coefficients); NIMBLE samples the cluster parameters
  directly.

## Prediction for neo-normal regression is now complete

* `posteriorLinpred()` and `posteriorEpred()` already worked for the neo-normal
  regressions (identity link), but `posteriorPredictive()` fell back to a
  Gaussian draw and lost the skew. It now draws from each family via a generic
  `responseRng` the framework attaches from the family's own RNG -- so a
  MSNBurr predictive is skewed, an FSST one heavy-tailed, and so on. Works
  across all shape counts (2 or 3 parameters); the shape traces are pulled by
  name from a small registry keyed on the family.

## New: neo-normal regression (framework + MSNBurr and SEP)

* Skewed and heavy-tailed regression: each neo-normal component's location is
  a linear predictor `X beta_k`, with its shape parameters staying
  per-component. `nimixReg(y ~ x, K = 2, distribution = "msnburr")` and
  `"sep"`, , `"msnburr2a"` -- and now all nine neo-normal families: msnburr, msnburr2a, gmsnburr, fssn, fsst, sep, lep, fossep, jfst. Shape sets range from two (sigma, alpha or sigma, nu) to three (sigma, alpha, theta), with Gamma, log-normal, or truncated-Gamma shape priors; the framework absorbs all of these -- a new family is a short declaration (density, shape list, prior), and the traces, relabelling, and summary columns follow from the shape list automatically. Measured on two components: MSNBurr
  intercepts -1.44/1.64 against -1.5/1.5; SEP intercepts -1.26/1.10 against
  -1.2/1.2.
* Built on a generic helper (`.neoRegMethods`) so a new family needs only its
  density, its list of shape parameters, and a short prior -- the six
  boilerplate S4 methods (traces, relabelling, simulation, ...) are generated
  from the shape list. MSNBurr and MSNBurr-IIa use `(sigma, alpha)`, SEP uses `(sigma, nu)`;
  the summary columns follow automatically (`sigma_mean` + the family's own
  shape means). The remaining seven neo-normal families can be added the same
  way. No conjugate sampler (the likelihood is not Gaussian in the
  coefficients), so budget more iterations, as with the other non-Gaussian
  regressions.

## New: Markov-switching heavy-tail regression (Student-t and Normal-Gamma)

* `nimixReg(y ~ x, K = 2, method = "hmm", distribution = "studentt")` (and
  `"normalgamma"`) fits a regression whose location coefficients and scale
  switch with a latent Markov regime, with Student-t errors -- the heavy-tail
  sibling of the Gaussian, Poisson, and Binomial Markov-switching regressions.
  Own forward kernel (`dt_nonstandard`, df a fixed hyperparameter), not the
  inherited Gaussian one: the guard refuses that inheritance (a silent-wrong
  kernel otherwise) and this dedicated kernel runs instead. Measured on two
  regimes with df = 4 errors: intercepts -1.50/1.45 against -1.5/1.5, slopes
  -0.88/1.26 against -0.9/1.2, decoding 0.98.
* `"normalgamma"` shares the kernel (identical Student-t marginal) but has its
  OWN `buildModelCode` method, because the two heavy-tail specs are siblings
  under `NormalRegSpec`, not parent and child -- inheritance would have handed
  it the Gaussian default.
* `nimixForecast()` handles both: the predictive draws keep the heavy tails.
  This closes the long-standing backlog item; the HMM regression families are
  now normal, studentt, normalgamma, poisson, and binomial (18 HMM kernels).

## New: `initMethod = "spread"` for heterogeneous-variance data

* A third cluster-seeding option alongside `"kmeans"` (default) and
  `"single"`. For a univariate response, `"spread"` bands observations by
  `|y - median(y)|`, separating components by scale rather than location.
  This targets the one case measured to seed k-means poorly: components with
  very different variances but overlapping means (e.g. sd 0.3 vs 3, both
  centred near zero), where k-means is forced to draw a spatial boundary that
  cuts across the real, scale-based split. Initial-allocation accuracy on that
  case rises from 0.84 (k-means) to 0.90.
* It is a convenience, not a cure. As three initialisation studies found, a
  well-mixed Bayesian mixture recovers from a poor seed on its own -- the
  MCMC result is the same whether seeded by `"kmeans"` or `"spread"`. The
  value of `"spread"` is a shorter burn-in on the heterogeneous-variance case,
  not a different answer. For a multivariate response it has no natural
  analogue and falls back to k-means.
* Shared internally via a single `.initClusters()` helper, so future seeding
  strategies are written once rather than in each family's `componentInits`.

## New: Markov-switching GLM regression (Poisson and Binomial)

* `nimixReg(y ~ x, K = 2, method = "hmm", distribution = "poisson")` fits a
  count regression whose log-rate coefficients switch with a latent Markov
  regime -- the count sibling of the Gaussian Markov-switching regression.
  Own forward kernel (`dpois(exp(Xbeta))`, log link, no error variance), not
  an inherited one. Measured on two regimes: intercepts 0.93/2.23 against
  1.0/2.2, slopes 0.84/-0.51 against 0.8/-0.5, decoding 0.967.
* `distribution = "binomial"` does the same with a logit link and a known
  number of trials (`prior = list(size = )`): a proportion regression with
  regime-switching coefficients. Measured: intercepts -1.03/1.50 against
  -1.0/1.5, slopes 1.22/-0.81 against 1.2/-0.8, decoding 0.997.
* `nimixForecast()` handles both: the predictive draws are counts (Poisson)
  or proportions of `size` (Binomial). Fixed alongside: the forecaster tested
  `is(spec, "NormalRegSpec")` to decide whether a fit had covariates, which
  missed the sibling GLM specs; it now tests the shared regression trait.

## Fixed: `posteriorEpred()`/`posteriorPredictive()` now respect the link

* The first cut of the brms-style predictors was Gaussian-only in disguise:
  `posteriorEpred()` averaged the linear predictors `x'beta_k` directly, and
  `posteriorPredictive()` added Gaussian noise. For a Poisson (log link) or
  Binomial (logit link) regression that is simply wrong -- `epred` returned
  the log-mean instead of the mean, and the predictive drew Gaussian jitter
  instead of counts.
* Now matched to brms. `posteriorEpred()` is `E[Y|X]` on the **response
  scale**, applying each family's inverse link to every component *before*
  the mixture average: `sum_k w_k g^{-1}(x'beta_k)`. Verified on a Poisson
  fit -- epred came back at 0.87/1.79/3.68 against the truth exp(0.5 + 0.8x)
  = 0.74/1.65/3.67, where the linear predictor was 0.5 + 0.8x.
* `posteriorLinpred()` stays on the **linear-predictor scale** by design (the
  log-mean for a Poisson), and gains a `transform` argument -- the
  counterpart of brms's -- to apply the inverse link when you want each
  component's response mean instead. For a Normal fit the two coincide, so
  this is a no-op there.
* `posteriorPredictive()` now draws from the family via a new `responseRng()`
  generic, and this reaches every univariate regression family, not just the
  Gaussian one. A Poisson predictive draw is a count, a Binomial one a number
  of successes, and a **Student-t or Normal-Gamma one is heavy-tailed** -- the
  first cut drew Gaussian noise for those too, giving a predictive with
  measured kurtosis 3.14 where a t with df = 3 should be far heavier; it is
  now 5.0 (Student-t) and 7.1 (Normal-Gamma). Their `epred` was already right,
  since both use the identity link and E[Y|X] = Xbeta; only the predictive law
  was wrong. Multivariate-response fits are unaffected -- `posteriorLinpred()`
  and friends decline them by design, and `predict()` handles that case.
  Fixed alongside: a `K = 1` regression could not be predicted with `newdata`
  (it has no mixture weight node); it now takes weight one.

## brms-style prediction, time-series views, and Markov autoregression

* **Three prediction functions, because a mixture makes them three different
  questions.** `posteriorLinpred()` returns each component's own linear
  predictor (a `draws x n x K` array); `posteriorEpred()` averages the
  components into `E[Y|X]`; `posteriorPredictive()` adds the residual noise.
  For a mixture, reach for the first far more often than the second -- fit
  two crossing lines with slopes +1.5 and -1.5 and the expectation is a flat
  line through the middle (measured 0.028/0.069/0.111 at x = -1/0/1). Nothing
  is wrong with the number; it is a faithful summary of a distribution no
  component has. `?posteriorEpred` says so.
* The weights are what change meaning between calls, and the docs are
  explicit about it: in-sample, `posteriorEpred()` uses each row's posterior
  allocation and becomes useful (measured correlation 0.983 with `y` on the
  same fit whose `newdata` expectation was flat); with `newdata` it uses the
  mixture weights, since a new row's component is unknown. For
  `method = "hmm"` and `newdata` it is refused outright: a regime weight is a
  function of *when*, and a future row has no decoded regime --
  `nimixForecast()` is the function that projects it.
* **Fixed: a silent wrong answer in prediction with `newdata`.** If `newdata`
  lacked a predictor, `model.frame()` resolved the name against the formula's
  environment -- the fitted data -- and returned predictions for the original
  rows. Measured before the fix: one row of `newdata` in, a `10 x 60 x 2`
  array out. `predict()` surfaced this only as the downstream "replacement
  has 60 rows, data has 1"; both paths now check up front and name the
  missing predictors.
* **Time-series views for regime-switching fits**: `plot(fit, type =
  "series")` draws the data with its decoded regimes shaded as blocks,
  `type = "regime"` the smoothed regime probabilities through time (the
  honest companion to the Viterbi path -- where the bands are mixed, the
  decode is a guess), and `type = "forecast"` the predictive fan. A density
  plot of a switching series discards the axis the model is about; these put
  it back. Each returns the tidy frame it drew, as the other plot types do.
* **Markov-switching autoregression works today, with no new code**: lag the
  response into the design matrix. `nimixReg(y ~ ylag, df, K = 2, method =
  "hmm")` fits an MS-AR(1) whose intercept *and* AR coefficient switch with
  the regime -- measured, intercepts -1.42/1.48 against -1.5/1.5 and AR
  coefficients 0.76/0.32 against 0.75/0.30, decoding at 1.0. This is not a
  coincidence of the implementation: conditioning on lagged `y` is the
  standard conditional likelihood for an AR model, and the forward kernel
  takes an arbitrary design matrix.
* One practical note now that it is measured: the `"hmm"` regression path
  needs more iterations than `"fixedk"`. Marginalising the regime path rules
  out the conjugate sampler, and NIMBLE's defaults mix the coefficients about
  four times slower per second (ESS/s 2.3 vs 8.9 on a two-regime benchmark),
  even though the wall time is shorter. Light runs decode the regimes fine
  but leave the coefficient intervals wide; `?nimixReg` says so.
* **And forecasting one now works too**, via `nimixForecast(fit, h, lags =
  c(ylag = 1))`. `lags` marks the predictors that are the response in
  disguise; those columns are *generated by the forecast*, not supplied to
  it, with each posterior draw feeding its own trajectory back as its own
  lag. That is what makes the interval widen the way an autoregression's
  should -- measured on an MS-AR(1), from 2.15 at h = 1 to 5.12 at h = 10,
  with RMSE 0.54 against 3.97 for a constant-mean benchmark. Exogenous
  predictors go in `newdata` and lags in `lags`; give a predictor to one or
  the other, never both (measured on an MS-ARX: RMSE 0.83 against 5.29).
* **Do not fake a lag through `newdata`.** It would be taken as a known
  covariate, held fixed, and the forecast would inherit whatever you invented
  -- confidently and without complaint, because nothing in a design matrix
  says "this column is the past". `lags` exists so that you do not have to,
  and `?nimixForecast` says so where it will be read.

## New: forecasting from regime-switching fits (`nimixForecast()`)

* `nimixForecast(fit, h)` draws from the posterior predictive `h` steps past
  the end of a series, integrating over both parameter uncertainty (one path
  per posterior draw) and regime uncertainty (the regime is sampled, not
  fixed). It is the forward algorithm run one step further: filter the regime
  distribution to the last observation, push it through the transition matrix
  `h` times, sample a regime, draw from that regime's emission. All 13
  emission families are supported, plus Markov-switching regression (pass
  `newdata` with the future covariates -- the regime can be projected
  forward, the predictors cannot).
* Returns `$summary` (mean/median/interval), `$regime` (an `h` x `K` matrix
  of regime probabilities), and `$draws` (the raw predictive sample).
* **Calibration measured**, not asserted: interval coverage 0.917 against a
  nominal 0.90 over a 12-step horizon on a two-regime benchmark, with RMSE
  1.93 against 3.23 for a constant-mean forecast.
* **Three honest limits, all documented in `?nimixForecast` and all pinned by
  tests.** (1) The point forecast reverts to the stationary mixture as `h`
  grows -- measured, the median fell 1.8 to -1.1 over twelve steps while the
  regime probabilities went 0.11/0.89 to 0.59/0.41, converging on the
  stationary 0.67/0.33. That is the model working: an HMM forecasts the
  regime, and once the regime is unknowable the best guess is the long-run
  average. (2) Accuracy hinges on whether the regime persists: on a
  Markov-switching regression that held its regime across the boundary the
  forecast beat the benchmark fivefold (RMSE 0.50 vs 2.65); on one that
  switched at exactly the first forecast step -- a 5% event under the fitted
  persistence -- it lost (4.42 vs 2.45), while the interval still covered. A
  switch one step ahead is unpredictable by construction. (3) With separated
  regimes the predictive is **bimodal**, so its median sits in the trough
  between the modes where the model puts almost no mass; `$summary` is
  convenient, `$draws` is faithful, and `$regime` usually answers the
  question actually being asked.

## New: Markov-switching regression (`nimixReg(..., method = "hmm")`)

* `nimixReg` gains `method = "hmm"`: the regression coefficients and error
  variance switch with a latent first-order Markov regime -- Hamilton's
  (1989) model, and the one classic mixture-of-regressions variant nimix was
  missing. Give `K` regimes as for `"fixedk"`; the regime path is
  marginalised out by a forward kernel and decoded afterwards, so
  `viterbiPath()` returns the most probable regime sequence. Measured on two
  regimes with opposite slopes: intercepts 1.96/-1.96 against 2/-2, slopes
  1.50/-1.39 against 1.5/-1.5, decoding at 0.98.
* Unlike every other emission, this one's density depends on `t` through the
  design matrix rather than only through `y[t]`; the kernel carries `X`.
  Note that **the rows are a time series** here -- their order matters, which
  is true nowhere else in `nimixReg`.
* **No conjugate sampler here, by requirement.** The Normal-Inverse-Gamma
  update conditions on the allocations, and the HMM marginalises them away,
  so the coefficient full conditional is no longer NIG and NIMBLE's defaults
  are used instead. Installing the conjugate sampler anyway would be exactly
  the silent-wrong-sampler bug fixed for Student-t regression in 1.4.0; a
  test now asserts the mechanism, not just the outcome.
* `"studentt"` and `"normalgamma"` are refused for now. They are not merely
  unimplemented: both `contains` the Gaussian regression spec and so
  *inherit* its new HMM model code, which would have handed a user asking
  for heavy tails a plain Gaussian fit with no error. The guard matches on
  the exact class rather than `is()`, and a test pins that.
* The usual caveat about the default error-variance prior applies with force
  here: regimes with opposing slopes inflate the global OLS residual
  variance that the prior is centred on (measured `s2` ~2.2x on the
  benchmark above). `prior = list(s2Guess = )` is the escape hatch, as for
  the other regression paths.

## `sigmaGuess` reaches the multivariate regressions

* `prior = list(sigmaGuess = )` now works for a multivariate response too,
  closing a gap the v1.4.0 study left open: it had been added to the
  clustering prior only. `"studentt"` and `"normalgamma"` inherit it.
* Measuring it first turned out to unify the two earlier findings rather
  than repeat either. The **mechanism** is the univariate regression one:
  the multivariate prior centres Sigma on `cov()` of the *global* OLS
  residuals, which carry between-component variation as well as within. On
  two components differing only in their coefficients, with isotropic
  within-covariance, that prior was **22.4x** the truth in trace with
  condition number **40** where the truth is a circle -- wrong size *and*
  wrong shape. The **severity** is the multivariate one: the fitted
  covariance came back at ~1.5x with condition ~2, because
  `df0 = d + 2` keeps the InverseWishart prior worth exactly one observation
  where the univariate InvGamma is worth four. Same error, four times less
  weight, and it shows.
* The default is therefore untouched, as in the clustering case. Reach for
  `sigmaGuess = s` when you know the residual scale is isotropic with
  variance `s`, or pass a `d x d` matrix. Every number here is anchored by a
  test.

## HMM emissions complete: 13 families, and a silent-dispatch bug fixed

* Four more emission families for `method = "hmm"`: **`sep`**, **`lep`**
  (exponential-power, whose tail-shape `nu` proved *well* identified --
  measured 1.30/3.48 against a truth of 1.2/3.5), the four-parameter
  skew families **`fossep`** and **`jfst`**, and **`binomial`** (regime-switching
  proportions with known `size`, via `prior = list(size = )`; measured
  recovery prob 0.150/0.597 against 0.15/0.6, Viterbi 1.0). Measured decoding on simulated
  regimes: Viterbi 0.993 (sep), 1.0 (lep, fossep, jfst). JFST's `alpha` and
  `theta` jointly govern skew and tails and are weakly identified, like
  FSST's df -- documented and asserted loosely, per that precedent. This
  completes the planned univariate set at **12 families**; `normal-gamma`
  is excluded *by design*, since its augmented representation is exactly
  what the marginalised forward kernel exists to avoid, and direct
  `student-t` serves the heavy-tail case.
* **Fixed: a silent wrong-model dispatch.** `NormalGammaUvSpec` inherits
  from `NormalUvSpec`, so the engine's `is()`-based family guard accepted
  `distribution = "normal-gamma"` and S4 dispatch then fell through to the
  inherited Gaussian HMM method: the user asked for heavy tails and got a
  plain normal fit, with no error, ever since the HMM engine existed. The
  guard now matches the exact class. Found by the new *permanent* guard
  test the moment it moved onto `normal-gamma` -- after being relocated
  four times as the rollout advanced, it now sits on a family that is
  excluded forever, so it can never go stale again.

## Sparse adjacency: the graph is no longer the wall

* `SpatialWeightSpec` now stores its graph as an **edge list**; the dense
  matrix is derived on demand (`getAdjacency()`) and refused above 5000
  nodes with a message pointing at the new `getEdges()`. The motivation was
  measured, not assumed: a 10 000-node space-time graph (100 locations over
  100 times -- a *modest* spatio-temporal problem) was OOM-killed during
  graph construction, before any model existed, while the quantities the
  Potts engine actually needs came to ~0.6 MB, some **2700x smaller** than
  the dense path's transients. After the refactor that same graph builds in
  a quarter of a second within ~40 MB, and a 50 000-node graph in about a
  second.
* **Nothing changes numerically.** The edge list is kept in exactly the
  column-major order that `which(upper.tri(A) & A > 0)` produced on the old
  dense path, so the engine's structural constants -- and therefore fits
  under a fixed seed -- are bit-identical; the test suite asserts the
  equivalence for rook and queen grids and for space-time graphs.
  `spatialWeights(A)` with a dense matrix still works, with the same
  validation and the same error messages, and `getAdjacency()` still returns
  the same named matrix for small graphs.
* New sparse entry points: `spatialWeights(edges = , nNodes = )` builds a
  graph without ever touching a dense matrix (edge order and duplicates are
  normalised), and `getEdges()` returns the canonical two-column form.
  `gridAdjacency()` and `spacetimeAdjacency()` now build their edge lists
  directly.
* One incidental find while chasing a memory gate: `duplicated()` on a
  matrix coerces every row to a character string -- 55 MB of transient on a
  28 000-edge graph, and twice that again inside validity. Edges are now
  deduplicated on a numeric key `(i - 1) * n + j`, which is exact in a
  double for any realistic `n`.
* **The honest ledger.** The fit ceiling in a 4 GB container rose from
  ~5000 to ~7000 nodes (57 s, allocation accuracy 0.999 at 7000); at 10 000
  the fit still dies, but the binding constraint is now NIMBLE's per-node
  model memory (~1.6 GB at 5000, ~2.8 GB at 7000), not the graph. More RAM
  now buys more nodes roughly linearly, which was not true before.
* **Breaking change:** `SpatialWeightSpec` objects serialised (`saveRDS`)
  under earlier versions cannot be loaded; rebuild them from their adjacency
  or edges. Code that only used the exported constructors and accessors is
  unaffected.

# nimix 1.4.0

## New: spatio-temporal mixtures via `spacetimeAdjacency()`

* `spacetimeAdjacency(W, nTime)` expands a spatial adjacency over time, so
  that node (i, t) neighbours its spatial neighbours at t and itself at
  t +/- 1. Pass the result to `nimixClust(..., method = "mrf")` and the Potts
  prior couples allocations across space **and** time. There is no new
  engine: the architectural study behind this found the long-standing
  blocker -- "undirected spatial dependence breaks the HMM forward
  algorithm" -- to be true but beside the point, since the MRF engine never
  used the forward algorithm; it samples the allocations directly. A
  space-time graph is simply another adjacency.
* Measured on a 5x5 grid over 8 time points with deliberately overlapping
  components: a plain mixture failed outright (allocation accuracy 0.50, one
  component recovered), the spatial graph reached 0.95, and the space-time
  graph 0.995 -- the temporal edges earn their place rather than merely
  decorating.
* One honest limit: the coupling is isotropic, because the Potts prior reads
  the adjacency as unweighted, so a single `beta` governs spatial and
  temporal edges alike. Misspecifying it is mild rather than fatal -- on
  regimes random across space but perfectly persistent in time, imposing the
  spatial edges anyway cost 2.5 percentage points (0.900 against 0.925 for a
  temporal-only graph, which `spatial = FALSE` builds).
* For a pure time series, `method = "hmm"` remains the better tool: it
  marginalises the state path and offers `viterbiPath()`. The two are
  complementary -- `hmm` for time alone, a space-time Potts when space
  matters too.

## New: random slopes in the mixture of regressions

* `nimixReg(y ~ x, data, random = ~ x | group)` adds a group-varying slope for
  `x` on top of the random intercept (as in `lme4`'s `(x|g)`, the intercept
  comes along); `random = ~ group` still gives the intercept alone. A gate
  prototype tested the obvious worry -- that group-varying slopes would fight
  with component identification, since components in a mixture of regressions
  are distinguished *by* their slopes -- and it did not materialise: even with
  group slope deviations spanning -2.05 to 1.43 against a component
  separation of 4 (so one group's effective slope sat exactly between the two
  components), group slopes recovered at `cor = 0.997` and component
  allocation at 95.3% accuracy.
* Both offsets use the sum-to-zero parameterisation and independent priors.
  The gate measured why: free offsets produce two translation ridges
  (`cor(beta1, mean(s)) = -0.929`, `cor(beta0, mean(b)) = -0.953`) with min
  ESS 52, against 226 under the constraint -- and notably the free version
  received *conjugate* samplers from NIMBLE and still lost, i.e. posterior
  geometry beats sampler class. Independence was checked against correlated
  truth (`rho = 0.92`): the independent model recovered both sets of offsets
  (`cor = 0.97`) and the correlation itself reappeared empirically
  (`cor(b_hat, s_hat) = 0.876`), so a correlated prior buys little at
  moderate group sizes. With sparse groups it would matter more.
* Semantics, as for the intercept: the component slopes absorb `mean(s_g)`
  and the reported `sRE` are centred. `tauSlope` estimates the spread of
  *your* groups' slopes.
* Random effects (intercept and slope) now also work with heavy-tailed
  residuals: `distribution = "studentt"` accepts `random = ~ group` and
  `random = ~ x | group`. This needed no sampler work, because that family
  keeps NIMBLE's default samplers (see the dispatch fix below) -- the offsets
  simply enter the linear predictor. Measured on simulated t4 residuals:
  slopes -1.99/1.95 with `cor(b_hat, centred truth) = 0.992` for the
  intercept, and -1.82/2.09 against a sum-to-zero prediction of -1.88/2.12
  with `cor = 0.961` for the slope.

## The prior scale on the component error variance is now settable, and documented

* `defaultPrior` centres the InvGamma prior on `s2` at the residual variance
  of a **global** OLS fit, which ignores the mixture. For separated
  components that measures the spread *between* components as much as the
  spread within one, so `s2` is biased upward -- and the bias is largest
  exactly where mixtures are most useful. A study measured it: with a
  prior/truth scale ratio near 60 the bias runs ~9x at 25 observations per
  component, 2.5x at 150, 1.2x at 1000, and 1.0x by 5000; against component
  separation at n = 150 it runs 1.0x for overlapping components and 4.7x for
  well-separated ones. The prior's six nominal pseudo-observations carry a
  variance ~39x the truth, which makes them worth roughly 236 real ones.
* **The default does not change**, because the conservatism turned out to be
  load-bearing rather than incidental: shrinking the prior scale by 10x left
  the DPM's K recovery untouched (modal K = 2 throughout) but made an
  over-specified `fixedk` fit (K = 4 against 2 true components) occupy three
  components instead of two. Wide components make two suffice; narrow ones
  make the model want more. The bias is also in the safe direction (wider
  predictive intervals), and the coefficients are unaffected either way --
  which is why it stayed invisible.
* **What is new is the escape hatch.** `g` and `nu0` were overridable but the
  scale was not, so a user who knew the within-component variance had no way
  to say so. Now `prior = list(s2Guess = 0.3)` makes 0.3 the prior mean of
  `s2`; `distribution = "studentt"` inherits it, and `s0` still sets the raw
  InvGamma scale for callers who think in those terms. On the benchmark above
  (true `s2` = 0.25) that moved the estimate from 0.78 (3.1x) to 0.36 (1.4x),
  slopes unchanged. The override is deliberately absolute: a relative
  multiplier on the automatic scale was considered and dropped, since "a
  tenth of a quantity that measures the wrong thing" is a knob, not a
  statement a user can defend.
* **Where a tighter prior is safe** is documented and tested on both sides:
  safe for the DPM and for `fixedk` with a correct `K` (where it halved the
  error in `s2` at no cost), unsafe for `fixedk` with an over-specified `K`,
  where a tenfold tighter prior occupied three components against the
  default's two. `?nimixReg` gains a "Reading the error variance" section
  with the numbers, the boundary, and a two-stage empirical-Bayes recipe for
  callers who need the scale right without external knowledge: the first
  fit's allocation is reliable even where its `s2` is not, so it can supply
  `s2Guess` for a second fit. Every number is anchored by a test.

## The multivariate prior ellipse: measured, kept, and now overridable

* The multivariate analogue of the scale study above turned out to be a
  different failure. `cov(data)` is the **global** covariance, so for
  separated components it is inflated only *along* the direction that
  separates them: measured on isotropic components separated along a vector
  v, 37.6x along v against 0.9x across it -- a prior ellipse with condition
  number 42.8 where the truth is a circle. The prior has the wrong **shape**,
  not merely the wrong size.
* **And yet the default survives it**, unlike the univariate one: the fitted
  covariance came back at 1.2x with condition 1.4. The reason is a happy
  piece of existing design -- `df0 = d + 2` makes the InverseWishart prior
  worth exactly **one observation, for every d**, where the univariate
  `nu0 = 3` is worth four. The multivariate bias is therefore both milder and
  *constant in the dimension* (1.18x at 200 observations per component at any
  d), which is the opposite of what one might expect. The univariate cannot
  copy the trick: its weight is floored at two observations by the `nu0 > 2`
  requirement that keeps empty components finite.
* What the default cannot express is shape, so `prior = list(sigmaGuess = s)`
  now sets the prior mean of a component's covariance directly -- a positive
  scalar reads as isotropic, or pass a `d x d` positive-definite matrix. On
  the benchmark above `sigmaGuess = 0.25` restored the prior ellipse to a
  circle (condition 42.8 to 1). The default is untouched; the numbers here
  are anchored by tests.

## Fixed: Student-t regression was using a Gaussian-error sampler (affects earlier versions)

* `StudentTRegSpec` inherits from `NormalRegSpec`, so S4 dispatch handed it
  `NormalRegSpec`'s `customizeSamplers`, which installs the exact
  Normal-Inverse-Gamma Gibbs step on `(betaTilde, s2Tilde)`. That conditional
  is exact only under a Gaussian likelihood: with `distribution = "studentt"`
  it drew from the wrong conditional, with no accept/reject to correct it, so
  the chain targeted the wrong stationary distribution. It affected
  `method = "fixedk"` (and any engine whose allocation node is `z`); the DPM
  path was never affected, since its allocation node is `xi` and the method
  returned early.
* The symptom was easy to miss: the slopes were barely touched (symmetric
  errors), and only the scale moved. Measured against a correct `RW_block`
  reference on the same model and data, `s2` was biased ~17% at `df = 4`
  (0.99/1.08 vs 0.84/0.91, MCSE 0.003), and the gap shrank to ~1% at
  `df = 30` -- the bias vanishes as t approaches Normal, which pins the
  mechanism on the likelihood mismatch rather than on Monte Carlo error.
* `StudentTRegSpec` now keeps NIMBLE's default samplers, which are correct
  here (verified: `s2` 0.84/0.91, matching the reference; ESS 580 for the
  coefficients and 1010 for the scale). `normal-reg` and `normal-gamma-reg`
  are unaffected -- they use their own, correct samplers, and the test suite
  now asserts the dispatch for all three.

## Fixed: random-effect priors are now data-scaled (affects 1.3.0)

* The random-intercept prior shipped in 1.3.0 used fixed bounds
  `dunif(0.01, 5)` on `tauRE`, which silently broke the offsets whenever the
  response was on a large scale: with `y` multiplied by 1000 the needed
  `tauRE` was 771 against a hard ceiling of 5, and `cor(b_hat, truth)`
  collapsed from 0.992 to 0.091 -- no error, no warning, just wrong group
  effects. The bounds now scale with the data like the rest of nimix's
  priors (`tauRE` with `sd(y)`; `tauSlope` with `sd(y)/sd(x)`), and inits
  scale with them. Verified across a 1000x response rescale
  (`cor(b_hat, truth)` 0.991, `tauRE` 817 against a realized spread of 772)
  and a 1000x predictor rescale (`tauSlope` 0.000274 against a realized
  0.00029, where it previously pinned at its 0.01 floor).
* The existing scale-equivariance lock never caught this because it rescales
  *predictors* only. A response-rescaling test now closes that gap, and a new
  `test-prior-scale-invariance.R` locks the invariant per engine.
* An audit of the remaining priors found no other violation: the `fixedk`,
  `dpm` and `mrf` engines are *exactly* scale-equivariant in the response
  (their conjugate samplers reproduce the draws on data-scaled priors), and
  `hmm` agrees to within Monte Carlo error. The one caveat worth knowing for
  `hmm` on large-scale data: its adaptive random-walk samplers start from
  NIMBLE's default proposal scale of 1 regardless of the data, so short
  chains show a transient adaptation difference (0.08 at 2500 iterations,
  0.001 against an MCSE of 0.0013 at 8000) -- allow adequate burnin rather
  than reading a short chain as bias.

# nimix 1.3.0

This release adds a fourth inference engine -- hidden-Markov mixtures for
regime switching in time series, with six emission families -- and random
intercepts for grouped data in the mixture of regressions, both delivered
through gated prototypes whose measured findings are recorded below.

## New: random intercepts in the mixture of regressions

* `nimixReg(y ~ x, data, random = ~ group, method = "fixedk")` adds a shared
  group offset `b_g ~ N(0, tau^2)` to every component's linear predictor.
  Two design decisions came out of a measured gate prototype and are baked
  in: (i) **sum-to-zero parameterisation** -- with free `b`, the component
  intercepts and `mean(b)` form a pure translation ridge
  (`cor = -0.979`, min ESS 25 of 2500); the constraint restored min ESS to
  191-238 with recovery intact (`cor(b_hat, centred truth) = 0.992`,
  `tau_hat` 0.83 against a realized group spread of 0.77), and the reported
  `b` are centred with the component intercepts absorbing the group mean.
  Note `tauRE` estimates the spread of *your* groups: with few groups the
  realized spread is itself variable (for 12 groups drawn at sigma = 0.8 it
  spans roughly 0.5-1.1), so judge `tau_hat` against the groups you have,
  not the population sigma you may have simulated from. (ii) The **exact NIG
  Gibbs sampler gains a random-effect offset**: the gate found NIMBLE's
  conjugacy detection does handle dynamic indexing in additive *scalar*
  form, but not
  the `inprod` form production uses -- so the P1 sampler now conditions on
  the current `b`. The P1 scale-equivariance lock still holds with RE
  active (measured slope discrepancy 0.0026 under a 1000x rescale), and the
  test suite asserts it.
* Scope: `method = "fixedk"` with `distribution = "normal"`, one grouping
  factor, random intercept. Random slopes and further families follow the
  gated plan; other combinations are refused with a pointed message.


## New engine: hidden-Markov mixtures (regime switching), `method = "hmm"`

* Component labels can now follow a first-order Markov chain in time:
  `nimixClust(y, K = S, method = "hmm")` fits a regime-switching mixture in
  which the state path is **marginalised out of the likelihood by the forward
  algorithm** -- the MCMC only ever samples the continuous parameters. The
  gate prototype measured why: 4000 iterations in ~1 s at T = 300, min
  ESS/sec 456 vs 144 for the naive latent-state model (x3.2 on an easy
  setting, and the marginalised model removes T discrete nodes from the
  graph). The forward kernel is exact against a pure-R reference (difference
  0, compiled and uncompiled), asserted in the test suite.
* Allocation draws are recovered **post-hoc by forward-filter
  backward-sampling (FFBS)** per retained draw, so every existing tool works
  unchanged on HMM fits: `relabel()`, `psm()`, `binderPartition()`,
  `plot()`, and the bayesplot adaptors. Measured on simulated regimes:
  location recovery -1.99/2.14 (truth -2/2), self-transitions 0.965/0.895
  (truth 0.95/0.90), per-time MAP and Viterbi decoding both at accuracy 1.0.
* `viterbiPath(fit)` returns the jointly most probable state sequence at the
  posterior means -- complementary to `binderPartition()`, which summarises
  marginal co-clustering across the FFBS draws.
* `nimbleEcology` was evaluated and is not used: its `dHMM` family is
  categorical-emission only, while regime switching on continuous data needs
  continuous emissions. The gate showed nimix's own kernels compile exactly
  inside the forward pass, which is the path for extending the engine to the
  other emission families -- current scope is `"normal"`, `"student-t"`,
  `"poisson"` (count regimes), and the neo-normal skewed families
  `"msnburr"`, `"msnburr2a"`, `"gmsnburr"`, `"fssn"`, and `"fsst"` -- eight
  families, all univariate. Measured recovery on simulated regimes: FSSN
  mu -3.06/3.27 with skew alpha 0.63/1.54, and FSST (skewed *and*
  heavy-tailed at once) mu -4.05/4.18 with alpha 0.69/1.51, both decoding at
  Viterbi 1.0. FSST's degrees of freedom are deliberately not asserted
  tightly: nu is weakly identified once the tails are moderate (this run gave
  12.8/13.3 against a truth of 6 while everything else landed), which is a
  property of the family rather than a defect. New emission families
  implement one density method (`.hmmEmisDens`) and one forward kernel; the
  engine, FFBS, and `viterbiPath()` are family-generic, as the Poisson case
  (a non-location-scale emission) confirms.
* FFBS allocation decoding and Viterbi are now numerically hardened against
  emission underflow: for a thin-tailed family an outlying point can drive
  every state's density to 0 at some draw, which previously produced NaN
  weights and an "NA in probability vector" crash. The forward pass now falls
  back to a uniform over states in that degenerate case, and Viterbi floors
  zero densities before taking logs -- decoding stays correct on normal
  regimes and no longer crashes on outliers (asserted by a dedicated test).

* Over-parameterised fits (nStates above the true number of regimes) leave
  empty states without corruption, asserted by a dedicated test (the
  structural lesson of the 1.2.0 PPC bug applied to a new engine from day
  one).

# nimix 1.2.1

## New: internal cluster-validity indices

* `clusterValidity(fit)` computes silhouette width, the Dunn index, and
  Calinski-Harabasz for a clustering fit's point partition (default:
  `binderPartition()`, so every posterior draw contributes and no relabelling
  is needed), via the `cluster` and `fpc` packages -- both in Suggests.
  The documentation states, and the test suite asserts, the honest caveat:
  these indices reward *geometric* separation, while mixtures are
  *density*-based -- overlapping components can be exactly the right model
  and still score low (measured: silhouette 0.90/Dunn 1.6 for separated
  clusters vs 0.52/0.0 for a legitimate overlapping fit). They are a
  secondary comparison lens; `ppCheck()` and `psm()` remain the primary
  model-adequacy tools.

# nimix 1.2.0

Response wave to an external code review: two correctness fixes, five new
exported functions, sampler-default upgrades, and a test-harness overhaul.
Every claim below was measured in-session; see the file headers for the
numbers' provenance.

## Improved: plots return their data; PPC column lookup memoised

* Every `plot(fit, type = ...)` branch now returns, invisibly, the tidy data
  frame it drew (`iteration`/`component`/`value` for traces, `x`/`density`
  for the predictive density, `dim*`/`cluster` for the MAP scatter,
  `fitted`/`observed` for regression). Base `graphics` remains the only
  plotting dependency; users who want ggplot2/lattice/plotly replot from the
  returned data.
* `ppCheck()`/`posteriorPredict()` now resolve each monitored node's columns
  once per call instead of once per draw (memoised via an attribute cache).
  Measured cost of the old lookup was ~3e-05 s/call, so this is a clarity
  and large-monitor-set improvement, not a speedup claim.
* Considered and rejected: passing prior hyperparameters as data nodes to
  avoid recompilation on prior changes. It would touch every family's model
  code and risks breaking NIMBLE's conjugacy detection (which the FixedK
  regression fix showed is already fragile under dynamic indexing); the
  compile cache already forces a correct rebuild when priors change.

## Fixed: 65 tests silently erroring in installed-mode runs

* Running the suite against the *installed* package (as opposed to
  `pkgload::load_all()`) is this project's guard against a documented bug
  class -- kernels resolving under `load_all()` but not after
  `R CMD INSTALL`. That guard had a blind spot: many tests called internal
  helpers (`.rowPresence()`, `.nodeToArray()`, `.cacheKey()`,
  `buildConstants()`, ...) unqualified, which works under `load_all()` but
  errors when only the namespace exports are attached -- and testthat counts
  those as *errors*, not *failures*, so a summary reading "0 failed" hid
  them. All internal references in tests are now `nimix:::`-qualified, and
  both suite modes now report identically (588 passing, zero errors). Suite
  gates now check the error flag, not just the failure count.

## Improved: AF_slice defaults for correlated-parameter univariate families

* The nine 3-4 parameter univariate families (`msnburr`, `msnburr2a`,
  `gmsnburr`, `fssn`, `fossep`, `fsst`, `jfst`, `sep`, `lep`) now sample each
  component's parameters as a single automated-factor-slice block (Tibbits et
  al. 2014) instead of independent univariate samplers. Motivation was
  measured, not assumed: on `fssn`, `cor(mu, alpha)` reaches -0.94 and the
  default samplers delivered min ESS 12 of 1500 draws (0.28 ESS/sec); on
  `gmsnburr`, min ESS 46 (0.86 ESS/sec). The escalation ladder mattered --
  `RW_block` made `fssn` *worse* (min ESS 7) -- while AF_slice reached min ESS
  622 (`fssn`, x32 ESS/sec) and 417 (`gmsnburr`, x5.7), with parameter
  recovery unchanged. Verified under all three engines (FixedK, DPM, MRF).

## New: bayesplot interoperability -- `drawsArray()` and `ppcData()`

* `drawsArray(fit)` returns a plain `iterations x chains x parameters` array
  -- the layout `bayesplot::mcmc_*` functions accept natively -- and
  `ppcData(fit)` returns `list(y, yrep)` for `bayesplot::ppc_dens_overlay()`
  and friends (with a `margin` argument for multivariate fits). bayesplot
  enters `Suggests` only; the adaptors return base R objects and add no
  runtime dependency.
* The important part is the safety guard, not the plumbing:
  `drawsArray(fit, "components")` **refuses** to serve per-component draws
  before `relabel()`, explaining why -- under label switching, `muTilde[1]`
  names different components in different chains, so an R-hat computed on the
  raw trace looks valid and means nothing. The default `"invariant"` view
  (cluster count, allocation entropy, `alpha`) is safe on raw draws and keeps
  the true per-chain structure. After relabelling, the chain dimension is
  honestly collapsed to 1, because conditioning on the modal cluster count
  leaves chains with unequal lengths.

## New: interoperability foundations -- `yrep` and chain identity

* `posteriorPredict(fit, ndraws)` returns the posterior predictive replicates
  themselves (`ndraws x n`, or `ndraws x n x d` for multivariate fits), and
  `ppCheck(..., store_yrep = TRUE)` attaches `yrep`/`y`/`draws` attributes.
  Previously the replicates were computed and discarded, which made graphical
  PPC (e.g. `bayesplot::ppc_dens_overlay(y, yrep)`) impossible to drive from
  a nimix fit. Storage stays opt-in so the default result remains lean.
* Multi-chain fits now record `diagnostics$chainId`, marking which chain each
  pooled draw came from. Post-hoc per-chain diagnostics (R-hat on invariant
  functionals, per-chain traces, draws arrays) were previously impossible to
  reconstruct because chains were stacked without a marker. Stored in the
  diagnostics list, so the `FitResult` class is unchanged and existing objects
  remain valid.

## New: label-free partition summaries -- `psm()` and `binderPartition()`

* `psm(fit)` returns the posterior similarity matrix `P(z_i = z_j | y)`;
  `binderPartition(fit)` selects, among the partitions the chain actually
  visited, the one minimising the expected Binder loss (Dahl's least-squares
  criterion). Both are invariant to label permutations *and* to the number of
  occupied clusters, so **every draw contributes** -- unlike `relabel()`,
  which must condition on the modal cluster count to align component
  parameters (measured: 34% of DPM draws discarded on a two-cluster example).
  They are complements, not replacements: `relabel()` answers "what are the
  component parameters", `psm()`/`binderPartition()` answer "which
  observations belong together". On overlapping clusters the similarity
  matrix expresses genuine allocation uncertainty (mid-region pairs ~0.65)
  instead of forcing a hard 0/1 answer.

## Fixed: multivariate posterior predictive checks with empty components

* `ppCheck()` for `normal-mv` (and the inheriting `student-t-mv`,
  `normal-gamma-mv`) reconstructed each draw's covariance array with
  `dim = c(max(alloc), d, d)`. On draws where a component was empty,
  `max(alloc)` undercounts the monitored components, truncating and shifting
  the covariance entries -- usually a `chol()` error, but occasionally a
  shifted-yet-still-PD matrix returning wrong replicates silently. The
  dimension is now derived from the monitored component count, with a
  defensive length check. Regression tests deliberately over-fit K so empty
  components must occur; the six skew-mv families were audited and were never
  affected.

## Fixed: FixedK mixture regression is now scale-equivariant

* Under the FixedK engine, NIMBLE's conjugacy checker cannot see through the
  dynamically indexed `betaTilde[z[i], ]`, so `betaTilde`/`s2Tilde` fell back
  to adaptive random-walk samplers in raw units. With poorly scaled predictors
  this biased the fit visibly (slope 2.50 vs true 2.0 when X was multiplied by
  1000 -- systematic, not mixing). A conjugate Normal-Inverse-Gamma Gibbs
  sampler now replaces them on the FixedK path (the DPM path was already
  conjugate via `CRP_cluster_wrapper`). Measured: max slope discrepancy under
  a 1000x predictor rescale fell from 0.48 to 0.004, and minimum ESS/second on
  the coefficients rose ~35x (0.4 to 15.3). A scale-equivariance test locks
  the guarantee.

# nimix 1.1.0

## New: estimating the orthogonal factor O beyond two dimensions

All six skew multivariate families now run under the finite-K, DPM and MRF
engines. Budget more MCMC iterations for the general-m variants: each
component carries `m(m-1)/2` slice-sampled angles, which slows the Potts
sweep -- on an 8x8 grid with `m = 3` and heavy tails, 1200 iterations
recovered the spatial regions poorly and 3000 recovered them exactly.

* `distribution = "skewistudent-mv-o"` likewise accepts any `m >= 2`. The
  canonicalisation carries `nu` with the permutation but never inverts it: a
  sign flip inverts `gamma` and leaves `nu` alone, because the Student kernel
  is symmetric. Density invariance checked to 1e-15 for `m = 2, 3, 4`.

* `distribution = "skewnormal-mv-o"` now accepts any `m >= 2`. It routes on the
  data dimension: `m = 2` keeps its dedicated implementation, `m > 2` uses the
  general Householder parameterisation with `m(m-1)/2` angles (FS Lemma 2).
  New exported helpers `orthogonalFactor()` and `canonicaliseO()`.
* **Restriction (8) is a canonicalisation, not a sampling constraint.** FS write
  that confining the angles to their box `Theta^j` puts `O` in `O_m`; testing
  this directly, the fraction of box draws that satisfy (8) is 0.245 (`m = 2`),
  0.069 (`m = 3`) and 0.007 (`m = 4`). Constraining a sampler to a 0.7% slice of
  its own prior would mix badly. What *is* true, and what nimix uses: among the
  signed row permutations `P` of `A` with `|P| = +1`, exactly one `PO` satisfies
  (8) -- verified exhaustively for `m = 2, 3, 4`. So the angles are sampled
  unconstrained and each posterior draw is mapped to its unique representative,
  with `gamma` carried along (`gamma_i -> gamma_perm(i)` or its reciprocal, per
  the row sign). The density is invariant under the map, to 1e-14. The `m! 2^m`
  row ambiguity of `A` is label switching in the dimension index, and this
  package already prefers post-hoc relabelling to ordering constraints.
* Reassuringly, for `m = 2` the (8)-satisfying set is exactly
  `theta in (-pi/8, pi/8)` -- the prior support already shipped. The general
  treatment reduces to the bivariate one rather than replacing it.
* **Experimental, and read with care.** `gamma` and `O` are reported *after*
  canonicalisation, so comparing them with simulating values requires
  canonicalising those too. `O_mean` is an elementwise posterior mean and is not
  itself orthogonal. The mirror modes of the angle likelihood multiply with `m`:
  partition and location recovery are robust (accuracy 1.0 at `m = 3`), but
  individual angles are large-sample quantities.

## New: estimating the orthogonal factor O (bivariate)

All four skew multivariate families -- fixed-O and estimated-O, Gaussian and
heavy-tailed -- run under the finite-K, DPM and MRF engines.

* `distribution = "skewistudent-mv-o"`: the heavy-tailed counterpart, with `O`
  estimated the same way. **Its `theta` is better identified than the Gaussian
  one**, and the reason is instructive: at `gamma = 1` the skew-Normal density
  is theta-invariant, because spherical Normal errors carry no directional
  information (FS Lemma 1). Independent Student margins are *not* spherical, so
  the skew-IStudent density depends on `theta` even under symmetry -- verified
  numerically, with the profile likelihood recovering `theta` from symmetric
  data. Letting `nu -> Inf` restores sphericity and with it the invariance,
  exactly as the theory predicts.
* Grid-initialisation of `theta` uses each family's own density. Initialising a
  heavy-tailed family from a Gaussian profile picks the wrong angle, because
  outliers dominate the Gaussian fit -- this was observed, not assumed.

* `distribution = "skewnormal-mv-o"`: the FS skew multivariate Normal with the
  orthogonal factor of `A = OU` **estimated** rather than held fixed, via the
  Householder angle `theta` (FS 2007, Appendix A). This lifts the main scope
  limitation of `skewnormal-mv`: `O` is what determines the *directions* of
  asymmetry (FS Sec 3.3). Bivariate data only for now.
* `theta` has a uniform prior on `(-pi/8, pi/8)`, which is exactly FS's
  identifiability restriction (8) once written in the Householder
  parameterisation (`O11 = cos 2 theta`, `O21 = -sin 2 theta`, `|O| = -1`).
  Because `|O| = -1` always, `O = I` is *not* a member of FS's restricted set:
  `theta = 0` gives `O = diag(1, -1)`, which coincides with `skewnormal-mv`
  after replacing `gamma_2` by `1/gamma_2`. The fixed-`O` family is therefore
  nested here at `theta = 0`, up to that reflection.
* **Read `theta` with care.** At `gamma = 1` the density does not depend on
  `theta` at all, so `theta` is identified only through the skewness. Even with
  clear skewness the likelihood has a near-mirror secondary mode: at 150
  observations per component it sat within 1.65 log-likelihood units of the true
  mode (both above the log-likelihood at the true parameters) and chains settled
  on the wrong sign; at 500 per component the 95% intervals covered the
  simulating angles. `theta` is slice-sampled and grid-initialised, and should
  be treated as a large-sample quantity.

## New: Ferreira-Steel skew multivariate distributions

Both skew multivariate families run under all three engines: finite-K, DPM
and MRF (spatial Potts prior).
* `distribution = "skewistudent-mv"`: the FS skew multivariate
  independent-Student (Ferreira & Steel 2007, Sec 5.2) -- FS-skew Student-t
  margins with per-dimension degrees of freedom `nu` (stochastic, truncated
  below at 2), the same `A = chol(Sigma)` construction and harmonised `gamma`
  convention as `skewnormal-mv`. Closed-form (no lambda augmentation), which
  avoids the documented 25-38x partition-mixing penalty of the augmented
  skew-Student path and is the model most supported by the data in FS's own
  application. Validated: kernel equals the R reference to 1e-15, `nu -> Inf`
  recovers `skewnormal-mv`, `m = 1` equals the univariate `fsst`, the density
  integrates to one, and mixtures recover location, per-dimension skew
  direction and partition (accuracy 1.0; DPM modal K correct).


* `distribution = "skewnormal-mv"`: the skewed multivariate distribution of
  Ferreira & Steel (2007), `eta = A' eps + mu`, with independent FS-skew-Normal
  margins for `eps`, per-dimension skewness `gamma` (harmonised convention:
  `gamma_j > 1` skews dimension j right; `gamma = 1` recovers `dmnorm`
  exactly), and `A = chol(Sigma)` upper triangular with
  `Sigma ~ inverse-Wishart`. Available under the finite-K and DPM engines, with
  `dskewmvn()` / `rskewmvn()` and `ppCheck()` support. Validated: compiled
  kernel equals the R reference to 1e-14, the 2-D density integrates to one,
  the m = 1 case equals the univariate `fssn`, and mixtures recover location,
  per-dimension skew direction, and partition (accuracy 1.0; DPM modal K
  correct).
* **Scope, stated plainly**: the orthogonal factor `O` of `A = OU` (FS Lemma 1)
  is fixed at the identity. FS Sec 3.3 explains that under skewness `O`
  determines the *directions* of asymmetry, so this release ties those
  directions to the coordinate axes after the triangular transform. A
  Householder-parameterised `O` (FS Appendix A) is planned as a follow-up, as
  is the skew multivariate independent-Student family.


## Breaking change: one skewness convention across all Fernandez-Steel families

Previously `fssn` and `fsst` parameterised skewness as `alpha = 1/gamma`, while
`fossep` used `alpha = gamma`. The same `alpha = 2` therefore skewed *left* in
`fssn`/`fsst` and *right* in `fossep`. Fits and simulations were internally
consistent within each family, so this never produced wrong inference, but it
made the families incomparable and the parameter uninterpretable across them.

All Fernandez-Steel families (`fssn`, `fsst`, `fossep`) now share the convention
of Fernandez & Steel (1998, 2007): the exported `alpha` **is** the FS skewness
`gamma`, so

* `alpha = 1` is symmetric,
* `alpha > 1` skews right, `alpha < 1` skews left,
* `P(X > mu) = alpha^2 / (1 + alpha^2)` exactly, in every family.

**What this means for you.** `fssn` and `fsst` results computed with 1.0.1 or
earlier correspond to the reciprocal `alpha` under 1.1.0: an estimate of
`alpha = 2` then is `alpha = 0.5` now. Densities, CDFs, quantiles, RNG and the
compiled NIMBLE kernels were all updated together, and a regression test
(`test-skew-convention.R`) now pins the guarantee. `fossep` is unchanged.
The Jones-Faddy family (`jfst`) uses its own `alpha`/`theta` shapes and is not
affected; neither are the MSNBurr families.

# nimix 1.0.1

## Batch B: six new neo-normal component families

All six are univariate, non-conjugate, and available under the finite-K, DPM and
MRF engines, with `d/p/q/r` functions and `ppCheck()` support. Each NIMBLE
kernel is built and registered in the global environment (the pattern required
for scalar user-defined densities) and was validated against the R reference to
1e-8, with the density integrating to one and the family's known reduction
reproduced exactly.

* `distribution = "sep"` -- symmetric exponential power (`nu = 2` Normal,
  `nu = 1` Laplace).
* `distribution = "lep"` -- exponential power under the alternative
  parameterisation (`nu = 2` Normal).
* `distribution = "fssn"` -- Fernandez-Steel skew Normal (`alpha = 1` Normal),
  with a log-normal prior on `alpha` that treats left/right skew symmetrically.
* `distribution = "fossep"` -- Fernandez-Steel skew exponential power
  (`theta = 2` skew-Normal kernel).
* `distribution = "fsst"` -- Fernandez-Steel skew Student-t (`alpha = 1`
  symmetric-t). `nu` is a stochastic node truncated below at 2 so the variance
  exists; it is only weakly identified by the data. The t-kernel is inlined
  rather than calling NIMBLE's `dt`.
* `distribution = "jfst"` -- Jones-Faddy skew-t (`alpha = theta` symmetric).
  The density uses the branch-free identity
  `sign(z)/sqrt((a+th)/z^2 + 1) == z/sqrt(a + th + z^2)`, which is also finite
  at `z = 0`.

The `q*` functions for these families recycle vector parameters to the sample
length and subset them per branch, so per-observation parameter vectors (as used
by posterior predictive simulation) align without recycling warnings.

# nimix 1.0.0 (in development)

## Bug fixes (installed-package correctness)

* MRF engine: the Potts prior's `dPottsNimix` / `rPottsNimix` are now built and
  registered in the global environment (like the scalar neo-normal densities).
  Registering them from the package namespace made NIMBLE fail to find
  `rPottsNimix` during code generation for the latent label node once the
  package was installed (`library(nimix)`), so every MRF fit errored under a
  normal install while working under `load_all()`. Fixed.
* `ppCheck()` for MSNBurr / MSNBurr-IIa / GMSNBurr no longer emits recycling
  warnings: the quantile functions now recycle vector parameters to the sample
  length and index the interior subset consistently, so per-observation
  parameter vectors (as used by posterior predictive simulation) align.
* Multivariate cluster `summary()` columns are now `mu_<j>_mean` (with
  `_med` / `_lwr` / `_upr`) for consistency with the univariate summaries.

## New neo-normal family: GMSNBurr

* `distribution = "gmsnburr"`: the generalized MSNBurr (Iriawan 2000; Choir
  2020) with two shape parameters `alpha` and `theta`. `theta = 1` recovers
  MSNBurr, `alpha = 1` recovers MSNBurr-IIa, `alpha = theta` is symmetric, and
  `alpha = theta -> inf` converges to the Normal. Available under all three
  engines (finite-K, DPM, MRF), numerically stable (NIMBLE density matches the
  R reference to 1e-8; density integrates to one; exact reduction to MSNBurr /
  MSNBurr-IIa verified). `d/p/q/r` functions and `ppCheck()` support included.

## Parallel chains

* `mcmcControl` gains `parallel = TRUE` (with `nchains > 1`) to run chains in
  parallel via `parallel::mclapply`. Each worker builds and compiles its own
  model in a separate directory -- the only fork-safe way to parallelise
  NIMBLE, avoiding the shared-C++-object and temp-directory collisions that a
  naive `mclapply` would hit. Forking only (Unix/macOS; Windows falls back to
  sequential); `ncores` caps the worker count. Verified to recover the same
  label-invariant summary as the sequential path.

## Bug fix: robust source load order

* Added `@include` directives across the S4 class and spec files so the
  package's `Collate` order is derived topologically from the actual class
  dependencies. This fixes a load failure (`undefined slot classes in
  definition of "MRFEngine"` / `no definition found for superclass
  "NormalRegSpec"`) that could occur whenever files were sourced in
  alphabetical rather than `Collate` order.

## Cluster profiling, richer summaries, model selection and ensembling

* `summary()` of a relabelled clustering fit now reports the posterior
  **median** alongside the mean and 95% credible interval for every
  scalar component parameter (all univariate families) and per dimension
  for multivariate clustering.
* New `clusterProfile()`: assigns each observation to its MAP cluster and
  describes the observed data within each cluster (size, proportion, and
  per-variable mean / sd / median) -- the data-side complement to
  `summary()` (which reports fitted parameters) and `plot(fit, "cluster")`
  (which shows the partition).
* New predictive model-selection layer built on the label-invariant
  pointwise mixture log-likelihood: `nimixWAIC()` (native; Watanabe 2010),
  `nimixLOO()` (PSIS-LOO via the **loo** package; Vehtari, Gelman & Gabry
  2017), and `modelSelect()` to rank several fits (e.g. choosing K, or
  comparing Normal / Student-t / MSNBurr components on the same data).
* New `ensembleFit()`: combines several fits into one weighted predictive
  model via Bayesian stacking or Pseudo-BMA+ (Yao et al. 2018, needs
  **loo**) or Akaike-style WAIC weights (native). `predict()` returns the
  ensemble-weighted density. `loo` added to Suggests.


## Neo-normal components: MSNBurr and MSNBurr-IIa

* Two new univariate component families, `distribution = "msnburr"` and
  `"msnburr2a"` (Iriawan 2000; Choir 2020), with a location, scale, and a
  skewness shape `alpha` (`alpha = 1` is exactly the logistic distribution;
  MSNBurr accommodates left skew, MSNBurr-IIa its mirror). Available under all
  three engines: finite-K, DPM, and the spatial MRF.
* Numerically stable throughout, using the maintainer-contributed reference
  implementation: an asymptotic log-omega branch for `alpha -> 0`, a
  branch-free softplus in the compiled densities, and two-branch quantile
  inversion. Log-densities stay finite for standardized values in the
  hundreds; the density integrates to one from `alpha = 0.05` to `alpha = 100`;
  the NIMBLE and R densities agree to 1e-8; and `ppCheck()` supports both
  families.
* Verified recovery on synthetic skew-component mixtures (finite-K allocation
  accuracy 1.00 / 0.99, DPM recovers the true K, MRF accuracy 1.00).

## Bayesian workflow: modern mixture-appropriate convergence + predictive checks

## Bayesian workflow: modern mixture-appropriate convergence + predictive checks

* Full Vehtari et al. (2021) convergence suite -- rank-normalized split-Rhat
  (already present), **folded split-Rhat, bulk-ESS and tail-ESS** -- computed
  over **label-invariant functionals only** (occupied K, DPM alpha, MRF beta,
  and the new per-iteration allocation entropy): per-component traces are not
  identified under label switching, so chain diagnostics on them would be
  meaningless. `summary()` prints the functional table with the documented
  Rhat < 1.01 aim; helpers are unit-tested against theory (iid ESS ~ n,
  folded Rhat detects pure scale disagreement that location Rhat misses,
  AR(0.9) ESS matches the analytic value).
* New `ppCheck()`: posterior predictive checks (Gelman, Meng & Stern 1996;
  Gelman et al. 2020, Bayesian workflow) with replicates drawn conditionally
  on each draw's fitted allocation -- label-invariant by construction and
  valid for every engine including the spatial MRF. Built-in statistics
  (mean, sd, min, max, skew) or user functions; clustering families
  (Gaussian uv/mv, Student-t, Normal-Gamma, Poisson, Binomial) supported,
  regression checks planned. Verified to pass a well-specified Gaussian
  mixture and to flag a Poisson mixture on overdispersed counts
  (tail p = 0.013).

# nimix 0.9.0

## MRF engine across the distribution registry (batch 2) -- matrix complete

* `method = "mrf"` now covers **every registered family**: added multivariate
  Student-t clustering, multivariate-response regression (Gaussian, Student-t
  and Normal-Gamma responses), and the augmented Normal-Gamma routes
  (univariate/multivariate clustering and univariate regression), each with
  its own O(nK + edges) label sweep. Harness groups `mrfdist2` and
  `mrfdist2reg` pass on the synthetic lattice (7 combinations x 3 seeds).
* **Critical correctness fix.** The hand-written MRF kernels for the Gaussian
  specs shadowed the generic `.pottsify()` route for their heavy-tailed
  subclasses (S4 inheritance), so e.g. a Student-t regression under `mrf`
  silently used the Gaussian kernel for its parameter updates. The
  hand-written kernels were removed: every family's MRF kernel is now derived
  from its own fixed-K kernel, and a regression test pins kernel/family
  correspondence.
* An internal quiet-mode handler crashed on NIMBLE conditions with empty
  messages ("argument is of length zero"); it is now robust, and the benign
  "No samplers assigned" notice for the deliberately unsampled fixed `beta`
  node is muffled.
* **Mixing benchmark (documented commitment):** on identical heavy-tailed
  lattice data, the augmented Normal-Gamma route mixes the component means
  BETTER than the direct Student-t density under the MRF (ESS of the minimum
  component mean 1357 vs 550 at equal wall time) -- the reverse of the
  exchangeable-DPM finding, because with K fixed and labels pinned by the
  spatial field the omega-partition coupling penalty largely disappears while
  the conjugate Gibbs updates win. Guidance: prefer `normalgamma` for MRF
  parameter mixing; the direct-t route remains available and equally valid.

# nimix 0.8.0

## MRF engine across the distribution registry (batch 1: closed-form emissions)

* Per the maintainer-approved feasibility study, `method = "mrf"` now covers
  six more family/task combinations: **Poisson** and **Binomial** clustering,
  **Poisson-GLM** (log link) and **Binomial-GLM** (logit link) regression, and
  **direct Student-t** clustering and regression. Every DSL emission call was
  verified empirically against reference densities before implementation.
* New internal `.pottsify()` transformation derives any family's MRF kernel
  mechanically from its fixed-K kernel (drop the Dirichlet-categorical label
  layer, insert the joint Potts node), so a default
  `buildModelCode(<any spec>, MRFEngine)` now exists; the three original
  hand-written MRF kernels are unchanged. Label sweeps remain family-specific
  for O(nK + edges) performance.
* `method = "mrf"` now requires `K >= 2` (a one-state Potts field is
  degenerate). The augmented Normal-Gamma routes remain explicitly blocked
  pending batch 2 (with a planned mixing benchmark against the direct-t
  routes).
* Recovery harness gains `mrfdist1` and `mrfdist1reg` groups: all six new
  combinations pass on the synthetic two-block lattice across 3 seeds
  (parameter recovery + allocation accuracy 0.98-1.00 in verification runs).


## Spatial line hardening + official-statistics case study

* Two new packaged official-statistics datasets: `usStates2023` (SAIPE 2023
  median household income and poverty rates for the 48 contiguous states +
  DC) and `usStateAdj` (the state contiguity matrix derived from the Census
  Bureau's 2023 county adjacency file); full provenance in the help pages.
* New vignette `spatial-mixture.Rmd`: on the 2023 SAIPE poverty rates the MRF
  engine (with estimated interaction, posterior mean ~1.2) finds two regimes
  whose high-poverty regime is the spatially contiguous Southern belt
  (AL AR DC FL GA LA MS NC NM OK SC TN TX) -- structure an exchangeable
  mixture cannot represent; a plain `fixedk` fit agrees on only 82% of
  states. The spatial regression example honestly collapses to a single
  national income-poverty regime.
* The recovery harness gains an `mrfbeta` group (pseudo-likelihood beta
  estimation on the synthetic lattice, 3 seeds).

# nimix 0.7.0

## MRF engine: Bayesian estimation of the interaction beta

* `prior$estimateBeta = TRUE` (with optional `prior$betaMax`, default 2) puts a
  uniform prior on the Potts interaction and updates it by random-walk
  Metropolis against the Besag (1975) **pseudo-likelihood** -- the classical
  approximate route for hidden Potts fields, since the exact posterior of beta
  is doubly intractable. This is documented as an approximation; an
  exchange-algorithm refinement (Murray, Ghahramani & MacKay 2006) is a
  possible future upgrade. `beta` is now a model node in both modes; with
  `estimateBeta = FALSE` (default) no sampler touches it, reproducing the
  fixed-beta behaviour exactly.
* Correctness guards: NIMBLE's default sampler on `beta` (which would target
  the unnormalised Potts density -- wrong in beta) is removed unconditionally,
  and the compiled-model cache now distinguishes estimation from fixed mode
  (they compile different sampler sets).
* MRF-specific diagnostics: `summary()` reports split-Rhat / ESS and the
  posterior mean of `beta` when it is estimated, and the engine warns when the
  posterior piles up near `betaMax` (a near-saturated field).
* On the synthetic two-block lattice, the estimated interaction concentrates
  well above zero (posterior mean ~1.6, P(beta > 0.2) = 1) with perfect block
  recovery, while the fixed-beta path is unchanged.


## MRF engine: spatially clustered regressions

* `nimixReg()` gains `method = "mrf"` and a `spatialWeights` argument: a
  mixture of Gaussian linear regressions whose latent regime labels follow the
  Potts field on a neighbourhood graph -- regression coefficients that cluster
  spatially (e.g. growth patterns across adjacent regions). Same fixed
  `prior$beta` interaction as the clustering engine.
* On a synthetic two-block lattice with opposite slopes (+2 / -2) and heavy
  noise, the spatial smoothing lifts the regime-allocation accuracy from 0.867
  (`fixedk`) to 1.000 while recovering both slopes. The recovery harness gains
  an `mrfreg` group (3 seeds).
* Heavy-tailed regression responses are explicitly blocked under `mrf` with a
  clear message (they inherit from the Gaussian regression spec and would
  otherwise use the wrong emission density).


## MRF engine: multivariate Gaussian components

* `method = "mrf"` now also accepts multivariate data: multivariate Gaussian
  components under the Normal-Inverse-Wishart kernel, with the label sweep
  evaluating `dmnorm_chol` per component (Cholesky factors hoisted once per
  sweep). On a synthetic two-block lattice the spatial smoothing lifts the
  multivariate allocation accuracy from 0.942 (`fixedk`) to 0.992.
* Sampler selection is fully polymorphic (S4 dispatch on the component spec);
  heavy-tailed families that inherit from the Gaussian specs are explicitly
  blocked with a clear message instead of silently using the wrong emission
  density.

# nimix 0.6.0

## Spatially constrained mixtures: the MRF engine

* New `method = "mrf"` in `nimixClust()`: a finite mixture whose latent labels
  follow a Potts Markov random field on a `spatialWeights` neighbourhood
  (Potts 1952; Besag 1974; spatially variant finite mixtures, Blekas et al.
  2005), so neighbouring observations favour the same component. Requires a
  known `K` and a `SpatialWeightSpec`; the interaction strength is fixed at
  `prior$beta` (default 0.8; `beta = 0` removes the smoothing). Univariate
  Gaussian components in this release.
* Implementation: an (intentionally unnormalised) user-defined NIMBLE Potts
  distribution -- exact for MCMC because `beta` is fixed, so the intractable
  partition function is constant -- plus a custom single-site Gibbs sweep
  sampler over the labels; theta updates remain conjugate. Bayesian estimation
  of `beta` is planned for a later 1.x release.
* On a synthetic two-block lattice with overlapping components, spatial
  smoothing lifts the allocation accuracy from 0.883 (`fixedk`, no smoothing)
  to 0.967, and `beta = 0` reproduces the unsmoothed behaviour. The recovery
  harness gains an `mrf` group (3 seeds, all passing).


Opens the post-1.0 spatial line: spatially constrained mixtures in which the
latent component labels follow a Markov random field on a neighbourhood graph
(Besag 1974; spatially variant finite mixtures, Blekas et al. 2005).

* New S4 class `SpatialWeightSpec`: a validated neighbourhood structure
  (symmetric, zero-diagonal, non-negative weight matrix) deliberately
  orthogonal to `DistributionSpec`, so any registered component family can be
  paired with any graph. Constructors `spatialWeights()` (from a matrix) and
  `gridAdjacency()` (rook/queen contiguity on a regular lattice); accessors
  `nRegions()`, `getAdjacency()`, `neighborsOf()`.
* `nimixClust()` gains a `spatialWeights` argument (default `NULL`, fully
  backward-compatible). Supplying a structure validates it and points to the
  MRF engine planned for 0.6.0; the exchangeable mixture remains the default
  behaviour.


First stable release. The public API -- `nimixClust()`, `nimixReg()`,
`registerDistribution()`, `relabel()`, `summary()`, `plot()`, `predict()`,
`nimixClearCache()` -- is now considered stable; breaking changes will bump the
major version.

## Official-statistics data and vignettes

* New packaged dataset `wdi2022`: four development indicators for 207
  countries (World Bank World Development Indicators, 2022; CC BY 4.0; full
  provenance in `?wdi2022`).
* All four vignettes now run on this official-statistics data instead of
  simulations only: income-regime clustering (univariate), joint
  income-longevity clustering (multivariate), the Preston curve as a mixture
  of regressions (an instructive single-regime result: the DPM does not invent
  components), and the Student-t vs Normal-Gamma heavy-tail comparison.
  Reported numbers in the vignettes come from actual runs.

# nimix 0.5.0

This line opens the performance and hardening phase, built on the 0.4.3 feature
set. The two production engines are `method = "dpm"` (Dirichlet process /
Chinese restaurant process; the number of occupied components is estimated) and
`method = "fixedk"` (finite mixture with a known `K`).

## Implemented in 0.9.x

* **Compile-once / reuse.** A fit reuses a compiled NIMBLE model (and its MCMC)
  when a later fit has an identical structure -- same generated code, constants,
  monitors and component family -- resetting only data and initial values. This
  skips recompilation for repeated fits (multiple seeds, multiple chains) and is
  bit-for-bit identical to a fresh compile. Controlled by `mcmcControl$reuse`
  (default `TRUE`); `nimixClearCache()` releases the cached compiled models.
* **Reproducibility.** The dispersed (k-means) initialisation is now seeded by
  the fit's `seed`, so repeated fits with the same data and seed coincide
  exactly. The caller's global random stream is left untouched.
* **Multi-chain + convergence diagnostics.** `mcmcControl$nchains` (default 1)
  runs several chains from dispersed, separately seeded starts, reusing the
  compiled model so only the first chain compiles. `summary()` reports
  rank-normalized split-Rhat and effective sample size for label-invariant
  quantities (occupied-cluster count, and the DPM concentration), and warns when
  Rhat exceeds 1.1. Per-component parameters are deliberately excluded from Rhat
  to avoid label-switching artefacts.
* **Vectorised post-processing.** The per-draw occupied-cluster count, the
  allocation-trace parsing, and the relabelling recode/weight steps are now
  vectorised (single `tabulate`/matrix-index passes instead of per-row loops).
  Verified bit-for-bit identical to the previous implementation; the remaining
  relabelling cost is the external ECR algorithm itself.
* **Hardening harness.** `inst/harness/run_recovery_suite.R` is now a
  systematic distribution x engine recovery matrix: every released clustering
  family (Gaussian, Student-t, Normal-Gamma -- univariate and multivariate --
  Poisson, Binomial) and the regression path are fitted with BOTH engines on
  data with known truth, across three MCMC seeds (18 combinations, 54 fits).
  For the DPM on discrete counts the recovery criterion assesses the dominant
  components (weight >= 0.1), reflecting the known diffuseness of the DPM
  posterior on the number of components for such data (Miller & Harrison,
  2013). Previously untested engine pairings are additionally pinned in
  `tests/testthat/test-hardening-matrix.R`.

# nimix 0.4.3

## Robustness and ergonomics

* `verbose` now defaults to `FALSE` for `nimixClust()` and `nimixReg()`. Benign
  NIMBLE configuration chatter (e.g. the Chinese-restaurant-process truncation
  reminder) is muffled selectively; genuine warnings about potentially invalid
  MCMC draws and all errors always propagate, in both verbose modes.
* Dispersed k-means initialisation is now headroom-aware: the number of seeded
  clusters respects the truncation level (`K_max`) so early CRP transients no
  longer breach it on smaller `K_max`.
* The initialisation headroom is configurable through
  `mcmcControl$initRatio` (default 0.8). Values that leave little headroom
  (>= 0.95) warn; values outside `(0, 1)` error.

## Engines and distributions

* Two production engines: `method = "dpm"` and `method = "fixedk"`, available
  for univariate and multivariate clustering and for regression (including
  multivariate responses).
* Component distributions: Gaussian (univariate / multivariate), Student-t and
  Normal-Gamma (heavy-tailed, univariate / multivariate), and Poisson /
  Binomial counts.

# nimix 0.4.2

* `nimixReg()` gains multivariate responses (`cbind(y1, y2) ~ x`) for Normal,
  Student-t and Normal-Gamma components, with per-component coefficient
  matrices and error covariances.

# nimix 0.4.0

* Student-t and Normal-Gamma components (univariate and multivariate) and
  Poisson / Binomial counts; public `registerDistribution()`.

# nimix 0.3.0

* `nimixReg()` and the `RegressionMixModel` class: mixture-of-regressions with a
  Normal-Inverse-Gamma g-prior. `FixedKEngine` implemented across univariate,
  multivariate and regression models. Engine selection is polymorphic via
  `runEngine()`.

# nimix 0.2.0

* Multivariate Gaussian clustering (`NormalMvSpec`) with a Normal-Inverse-Wishart
  base measure. Engine generalised to be dimension-agnostic.

# nimix 0.1.0

* S4 foundation, univariate Gaussian clustering (`NormalUvSpec`),
  `nimixClust()` on the DPM and fixed-K engines.
