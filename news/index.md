# Changelog

## nimix 1.3.0

This release adds a fourth inference engine – hidden-Markov mixtures for
regime switching in time series, with six emission families – and random
intercepts for grouped data in the mixture of regressions, both
delivered through gated prototypes whose measured findings are recorded
below.

### New: random intercepts in the mixture of regressions

- `nimixReg(y ~ x, data, random = ~ group, method = "fixedk")` adds a
  shared group offset `b_g ~ N(0, tau^2)` to every component’s linear
  predictor. Two design decisions came out of a measured gate prototype
  and are baked in: (i) **sum-to-zero parameterisation** – with free
  `b`, the component intercepts and `mean(b)` form a pure translation
  ridge (`cor = -0.979`, min ESS 25 of 2500); the constraint restored
  min ESS to 205-238 with recovery intact
  (`cor(b_hat, centred truth) = 0.992`, `tau_hat` 0.83 vs 0.8), and the
  reported `b` are centred with the component intercepts absorbing the
  group mean. (ii) The **exact NIG Gibbs sampler gains a random-effect
  offset**: the gate found NIMBLE’s conjugacy detection does handle
  dynamic indexing in additive *scalar* form, but not the `inprod` form
  production uses – so the P1 sampler now conditions on the current `b`.
  The P1 scale-equivariance lock still holds with RE active (measured
  slope discrepancy 0.0026 under a 1000x rescale), and the test suite
  asserts it.
- Scope: `method = "fixedk"` with `distribution = "normal"`, one
  grouping factor, random intercept. Random slopes and further families
  follow the gated plan; other combinations are refused with a pointed
  message.

### New engine: hidden-Markov mixtures (regime switching), `method = "hmm"`

- Component labels can now follow a first-order Markov chain in time:
  `nimixClust(y, K = S, method = "hmm")` fits a regime-switching mixture
  in which the state path is **marginalised out of the likelihood by the
  forward algorithm** – the MCMC only ever samples the continuous
  parameters. The gate prototype measured why: 4000 iterations in ~1 s
  at T = 300, min ESS/sec 456 vs 144 for the naive latent-state model
  (x3.2 on an easy setting, and the marginalised model removes T
  discrete nodes from the graph). The forward kernel is exact against a
  pure-R reference (difference 0, compiled and uncompiled), asserted in
  the test suite.

- Allocation draws are recovered **post-hoc by forward-filter
  backward-sampling (FFBS)** per retained draw, so every existing tool
  works unchanged on HMM fits:
  [`relabel()`](https://madsyair.github.io/nimix/reference/relabel.md),
  [`psm()`](https://madsyair.github.io/nimix/reference/psm.md),
  [`binderPartition()`](https://madsyair.github.io/nimix/reference/binderPartition.md),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), and the
  bayesplot adaptors. Measured on simulated regimes: location recovery
  -1.99/2.14 (truth -2/2), self-transitions 0.965/0.895 (truth
  0.95/0.90), per-time MAP and Viterbi decoding both at accuracy 1.0.

- `viterbiPath(fit)` returns the jointly most probable state sequence at
  the posterior means – complementary to
  [`binderPartition()`](https://madsyair.github.io/nimix/reference/binderPartition.md),
  which summarises marginal co-clustering across the FFBS draws.

- `nimbleEcology` was evaluated and is not used: its `dHMM` family is
  categorical-emission only, while regime switching on continuous data
  needs continuous emissions. The gate showed nimix’s own kernels
  compile exactly inside the forward pass, which is the path for
  extending the engine to the other emission families – current scope is
  `"normal"`, `"student-t"`, `"poisson"` (count regimes; lambda
  2.94/14.98 on truth 3/15, Viterbi 0.993), and the neo-normal skewed
  families `"msnburr"`, `"msnburr2a"`, `"gmsnburr"` (measured recovery
  e.g. msnburr2a mu -3.01/2.93, gmsnburr mu -4.12/4.27, Viterbi
  0.99-1.0) – all univariate; the t family targets heavy-tailed regimes
  and recovered locations/transitions/decoding on simulated t4 regimes:
  mu -2.01/2.12 vs truth -2/2, Viterbi accuracy 1.0. further families
  follow the gated plan one at a time. New emission families implement
  one density method and one forward kernel; the engine, FFBS, and
  [`viterbiPath()`](https://madsyair.github.io/nimix/reference/viterbiPath.md)
  are family-generic.

- FFBS allocation decoding and Viterbi are now numerically hardened
  against emission underflow: for a thin-tailed family an outlying point
  can drive every state’s density to 0 at some draw, which previously
  produced NaN weights and an “NA in probability vector” crash. The
  forward pass now falls back to a uniform over states in that
  degenerate case, and Viterbi floors zero densities before taking logs
  – decoding stays correct on normal regimes and no longer crashes on
  outliers (asserted by a dedicated test).

- Over-parameterised fits (nStates above the true number of regimes)
  leave empty states without corruption, asserted by a dedicated test
  (the structural lesson of the 1.2.0 PPC bug applied to a new engine
  from day one).

## nimix 1.2.1

### New: internal cluster-validity indices

- `clusterValidity(fit)` computes silhouette width, the Dunn index, and
  Calinski-Harabasz for a clustering fit’s point partition (default:
  [`binderPartition()`](https://madsyair.github.io/nimix/reference/binderPartition.md),
  so every posterior draw contributes and no relabelling is needed), via
  the `cluster` and `fpc` packages – both in Suggests. The documentation
  states, and the test suite asserts, the honest caveat: these indices
  reward *geometric* separation, while mixtures are *density*-based –
  overlapping components can be exactly the right model and still score
  low (measured: silhouette 0.90/Dunn 1.6 for separated clusters vs
  0.52/0.0 for a legitimate overlapping fit). They are a secondary
  comparison lens;
  [`ppCheck()`](https://madsyair.github.io/nimix/reference/ppCheck.md)
  and [`psm()`](https://madsyair.github.io/nimix/reference/psm.md)
  remain the primary model-adequacy tools.

## nimix 1.2.0

Response wave to an external code review: two correctness fixes, five
new exported functions, sampler-default upgrades, and a test-harness
overhaul. Every claim below was measured in-session; see the file
headers for the numbers’ provenance.

### Improved: plots return their data; PPC column lookup memoised

- Every `plot(fit, type = ...)` branch now returns, invisibly, the tidy
  data frame it drew (`iteration`/`component`/`value` for traces,
  `x`/`density` for the predictive density, `dim*`/`cluster` for the MAP
  scatter, `fitted`/`observed` for regression). Base `graphics` remains
  the only plotting dependency; users who want ggplot2/lattice/plotly
  replot from the returned data.
- [`ppCheck()`](https://madsyair.github.io/nimix/reference/ppCheck.md)/[`posteriorPredict()`](https://madsyair.github.io/nimix/reference/posteriorPredict.md)
  now resolve each monitored node’s columns once per call instead of
  once per draw (memoised via an attribute cache). Measured cost of the
  old lookup was ~3e-05 s/call, so this is a clarity and
  large-monitor-set improvement, not a speedup claim.
- Considered and rejected: passing prior hyperparameters as data nodes
  to avoid recompilation on prior changes. It would touch every family’s
  model code and risks breaking NIMBLE’s conjugacy detection (which the
  FixedK regression fix showed is already fragile under dynamic
  indexing); the compile cache already forces a correct rebuild when
  priors change.

### Fixed: 65 tests silently erroring in installed-mode runs

- Running the suite against the *installed* package (as opposed to
  [`pkgload::load_all()`](https://pkgload.r-lib.org/reference/load_all.html))
  is this project’s guard against a documented bug class – kernels
  resolving under `load_all()` but not after `R CMD INSTALL`. That guard
  had a blind spot: many tests called internal helpers
  (`.rowPresence()`, `.nodeToArray()`, `.cacheKey()`,
  [`buildConstants()`](https://madsyair.github.io/nimix/reference/buildConstants.md),
  …) unqualified, which works under `load_all()` but errors when only
  the namespace exports are attached – and testthat counts those as
  *errors*, not *failures*, so a summary reading “0 failed” hid them.
  All internal references in tests are now `nimix:::`-qualified, and
  both suite modes now report identically (588 passing, zero errors).
  Suite gates now check the error flag, not just the failure count.

### Improved: AF_slice defaults for correlated-parameter univariate families

- The nine 3-4 parameter univariate families (`msnburr`, `msnburr2a`,
  `gmsnburr`, `fssn`, `fossep`, `fsst`, `jfst`, `sep`, `lep`) now sample
  each component’s parameters as a single automated-factor-slice block
  (Tibbits et al. 2014) instead of independent univariate samplers.
  Motivation was measured, not assumed: on `fssn`, `cor(mu, alpha)`
  reaches -0.94 and the default samplers delivered min ESS 12 of 1500
  draws (0.28 ESS/sec); on `gmsnburr`, min ESS 46 (0.86 ESS/sec). The
  escalation ladder mattered – `RW_block` made `fssn` *worse* (min
  ESS 7) – while AF_slice reached min ESS 622 (`fssn`, x32 ESS/sec) and
  417 (`gmsnburr`, x5.7), with parameter recovery unchanged. Verified
  under all three engines (FixedK, DPM, MRF).

### New: bayesplot interoperability – `drawsArray()` and `ppcData()`

- `drawsArray(fit)` returns a plain `iterations x chains x parameters`
  array – the layout `bayesplot::mcmc_*` functions accept natively – and
  `ppcData(fit)` returns `list(y, yrep)` for
  [`bayesplot::ppc_dens_overlay()`](https://mc-stan.org/bayesplot/reference/PPC-distributions.html)
  and friends (with a `margin` argument for multivariate fits).
  bayesplot enters `Suggests` only; the adaptors return base R objects
  and add no runtime dependency.
- The important part is the safety guard, not the plumbing:
  `drawsArray(fit, "components")` **refuses** to serve per-component
  draws before
  [`relabel()`](https://madsyair.github.io/nimix/reference/relabel.md),
  explaining why – under label switching, `muTilde[1]` names different
  components in different chains, so an R-hat computed on the raw trace
  looks valid and means nothing. The default `"invariant"` view (cluster
  count, allocation entropy, `alpha`) is safe on raw draws and keeps the
  true per-chain structure. After relabelling, the chain dimension is
  honestly collapsed to 1, because conditioning on the modal cluster
  count leaves chains with unequal lengths.

### New: interoperability foundations – `yrep` and chain identity

- `posteriorPredict(fit, ndraws)` returns the posterior predictive
  replicates themselves (`ndraws x n`, or `ndraws x n x d` for
  multivariate fits), and `ppCheck(..., store_yrep = TRUE)` attaches
  `yrep`/`y`/`draws` attributes. Previously the replicates were computed
  and discarded, which made graphical PPC
  (e.g. `bayesplot::ppc_dens_overlay(y, yrep)`) impossible to drive from
  a nimix fit. Storage stays opt-in so the default result remains lean.
- Multi-chain fits now record `diagnostics$chainId`, marking which chain
  each pooled draw came from. Post-hoc per-chain diagnostics (R-hat on
  invariant functionals, per-chain traces, draws arrays) were previously
  impossible to reconstruct because chains were stacked without a
  marker. Stored in the diagnostics list, so the `FitResult` class is
  unchanged and existing objects remain valid.

### New: label-free partition summaries – `psm()` and `binderPartition()`

- `psm(fit)` returns the posterior similarity matrix `P(z_i = z_j | y)`;
  `binderPartition(fit)` selects, among the partitions the chain
  actually visited, the one minimising the expected Binder loss (Dahl’s
  least-squares criterion). Both are invariant to label permutations
  *and* to the number of occupied clusters, so **every draw
  contributes** – unlike
  [`relabel()`](https://madsyair.github.io/nimix/reference/relabel.md),
  which must condition on the modal cluster count to align component
  parameters (measured: 34% of DPM draws discarded on a two-cluster
  example). They are complements, not replacements:
  [`relabel()`](https://madsyair.github.io/nimix/reference/relabel.md)
  answers “what are the component parameters”,
  [`psm()`](https://madsyair.github.io/nimix/reference/psm.md)/[`binderPartition()`](https://madsyair.github.io/nimix/reference/binderPartition.md)
  answer “which observations belong together”. On overlapping clusters
  the similarity matrix expresses genuine allocation uncertainty
  (mid-region pairs ~0.65) instead of forcing a hard 0/1 answer.

### Fixed: multivariate posterior predictive checks with empty components

- [`ppCheck()`](https://madsyair.github.io/nimix/reference/ppCheck.md)
  for `normal-mv` (and the inheriting `student-t-mv`, `normal-gamma-mv`)
  reconstructed each draw’s covariance array with
  `dim = c(max(alloc), d, d)`. On draws where a component was empty,
  `max(alloc)` undercounts the monitored components, truncating and
  shifting the covariance entries – usually a
  [`chol()`](https://rdrr.io/r/base/chol.html) error, but occasionally a
  shifted-yet-still-PD matrix returning wrong replicates silently. The
  dimension is now derived from the monitored component count, with a
  defensive length check. Regression tests deliberately over-fit K so
  empty components must occur; the six skew-mv families were audited and
  were never affected.

### Fixed: FixedK mixture regression is now scale-equivariant

- Under the FixedK engine, NIMBLE’s conjugacy checker cannot see through
  the dynamically indexed `betaTilde[z[i], ]`, so `betaTilde`/`s2Tilde`
  fell back to adaptive random-walk samplers in raw units. With poorly
  scaled predictors this biased the fit visibly (slope 2.50 vs true 2.0
  when X was multiplied by 1000 – systematic, not mixing). A conjugate
  Normal-Inverse-Gamma Gibbs sampler now replaces them on the FixedK
  path (the DPM path was already conjugate via `CRP_cluster_wrapper`).
  Measured: max slope discrepancy under a 1000x predictor rescale fell
  from 0.48 to 0.004, and minimum ESS/second on the coefficients rose
  ~35x (0.4 to 15.3). A scale-equivariance test locks the guarantee.

## nimix 1.1.0

### New: estimating the orthogonal factor O beyond two dimensions

All six skew multivariate families now run under the finite-K, DPM and
MRF engines. Budget more MCMC iterations for the general-m variants:
each component carries `m(m-1)/2` slice-sampled angles, which slows the
Potts sweep – on an 8x8 grid with `m = 3` and heavy tails, 1200
iterations recovered the spatial regions poorly and 3000 recovered them
exactly.

- `distribution = "skewistudent-mv-o"` likewise accepts any `m >= 2`.
  The canonicalisation carries `nu` with the permutation but never
  inverts it: a sign flip inverts `gamma` and leaves `nu` alone, because
  the Student kernel is symmetric. Density invariance checked to 1e-15
  for `m = 2, 3, 4`.

- `distribution = "skewnormal-mv-o"` now accepts any `m >= 2`. It routes
  on the data dimension: `m = 2` keeps its dedicated implementation,
  `m > 2` uses the general Householder parameterisation with `m(m-1)/2`
  angles (FS Lemma 2). New exported helpers
  [`orthogonalFactor()`](https://madsyair.github.io/nimix/reference/orthogonalFactor.md)
  and
  [`canonicaliseO()`](https://madsyair.github.io/nimix/reference/canonicaliseO.md).

- **Restriction (8) is a canonicalisation, not a sampling constraint.**
  FS write that confining the angles to their box `Theta^j` puts `O` in
  `O_m`; testing this directly, the fraction of box draws that
  satisfy (8) is 0.245 (`m = 2`), 0.069 (`m = 3`) and 0.007 (`m = 4`).
  Constraining a sampler to a 0.7% slice of its own prior would mix
  badly. What *is* true, and what nimix uses: among the signed row
  permutations `P` of `A` with `|P| = +1`, exactly one `PO` satisfies

  8.  – verified exhaustively for `m = 2, 3, 4`. So the angles are
      sampled unconstrained and each posterior draw is mapped to its
      unique representative, with `gamma` carried along
      (`gamma_i -> gamma_perm(i)` or its reciprocal, per the row sign).
      The density is invariant under the map, to 1e-14. The `m! 2^m` row
      ambiguity of `A` is label switching in the dimension index, and
      this package already prefers post-hoc relabelling to ordering
      constraints.

- Reassuringly, for `m = 2` the (8)-satisfying set is exactly
  `theta in (-pi/8, pi/8)` – the prior support already shipped. The
  general treatment reduces to the bivariate one rather than replacing
  it.

- **Experimental, and read with care.** `gamma` and `O` are reported
  *after* canonicalisation, so comparing them with simulating values
  requires canonicalising those too. `O_mean` is an elementwise
  posterior mean and is not itself orthogonal. The mirror modes of the
  angle likelihood multiply with `m`: partition and location recovery
  are robust (accuracy 1.0 at `m = 3`), but individual angles are
  large-sample quantities.

### New: estimating the orthogonal factor O (bivariate)

All four skew multivariate families – fixed-O and estimated-O, Gaussian
and heavy-tailed – run under the finite-K, DPM and MRF engines.

- `distribution = "skewistudent-mv-o"`: the heavy-tailed counterpart,
  with `O` estimated the same way. **Its `theta` is better identified
  than the Gaussian one**, and the reason is instructive: at `gamma = 1`
  the skew-Normal density is theta-invariant, because spherical Normal
  errors carry no directional information (FS Lemma 1). Independent
  Student margins are *not* spherical, so the skew-IStudent density
  depends on `theta` even under symmetry – verified numerically, with
  the profile likelihood recovering `theta` from symmetric data. Letting
  `nu -> Inf` restores sphericity and with it the invariance, exactly as
  the theory predicts.

- Grid-initialisation of `theta` uses each family’s own density.
  Initialising a heavy-tailed family from a Gaussian profile picks the
  wrong angle, because outliers dominate the Gaussian fit – this was
  observed, not assumed.

- `distribution = "skewnormal-mv-o"`: the FS skew multivariate Normal
  with the orthogonal factor of `A = OU` **estimated** rather than held
  fixed, via the Householder angle `theta` (FS 2007, Appendix A). This
  lifts the main scope limitation of `skewnormal-mv`: `O` is what
  determines the *directions* of asymmetry (FS Sec 3.3). Bivariate data
  only for now.

- `theta` has a uniform prior on `(-pi/8, pi/8)`, which is exactly FS’s
  identifiability restriction (8) once written in the Householder
  parameterisation (`O11 = cos 2 theta`, `O21 = -sin 2 theta`,
  `|O| = -1`). Because `|O| = -1` always, `O = I` is *not* a member of
  FS’s restricted set: `theta = 0` gives `O = diag(1, -1)`, which
  coincides with `skewnormal-mv` after replacing `gamma_2` by
  `1/gamma_2`. The fixed-`O` family is therefore nested here at
  `theta = 0`, up to that reflection.

- **Read `theta` with care.** At `gamma = 1` the density does not depend
  on `theta` at all, so `theta` is identified only through the skewness.
  Even with clear skewness the likelihood has a near-mirror secondary
  mode: at 150 observations per component it sat within 1.65
  log-likelihood units of the true mode (both above the log-likelihood
  at the true parameters) and chains settled on the wrong sign; at 500
  per component the 95% intervals covered the simulating angles. `theta`
  is slice-sampled and grid-initialised, and should be treated as a
  large-sample quantity.

### New: Ferreira-Steel skew multivariate distributions

Both skew multivariate families run under all three engines: finite-K,
DPM and MRF (spatial Potts prior). \*
`distribution = "skewistudent-mv"`: the FS skew multivariate
independent-Student (Ferreira & Steel 2007, Sec 5.2) – FS-skew Student-t
margins with per-dimension degrees of freedom `nu` (stochastic,
truncated below at 2), the same `A = chol(Sigma)` construction and
harmonised `gamma` convention as `skewnormal-mv`. Closed-form (no lambda
augmentation), which avoids the documented 25-38x partition-mixing
penalty of the augmented skew-Student path and is the model most
supported by the data in FS’s own application. Validated: kernel equals
the R reference to 1e-15, `nu -> Inf` recovers `skewnormal-mv`, `m = 1`
equals the univariate `fsst`, the density integrates to one, and
mixtures recover location, per-dimension skew direction and partition
(accuracy 1.0; DPM modal K correct).

- `distribution = "skewnormal-mv"`: the skewed multivariate distribution
  of Ferreira & Steel (2007), `eta = A' eps + mu`, with independent
  FS-skew-Normal margins for `eps`, per-dimension skewness `gamma`
  (harmonised convention: `gamma_j > 1` skews dimension j right;
  `gamma = 1` recovers `dmnorm` exactly), and `A = chol(Sigma)` upper
  triangular with `Sigma ~ inverse-Wishart`. Available under the
  finite-K and DPM engines, with
  [`dskewmvn()`](https://madsyair.github.io/nimix/reference/skewnormal-mv-distribution.md)
  /
  [`rskewmvn()`](https://madsyair.github.io/nimix/reference/skewnormal-mv-distribution.md)
  and
  [`ppCheck()`](https://madsyair.github.io/nimix/reference/ppCheck.md)
  support. Validated: compiled kernel equals the R reference to 1e-14,
  the 2-D density integrates to one, the m = 1 case equals the
  univariate `fssn`, and mixtures recover location, per-dimension skew
  direction, and partition (accuracy 1.0; DPM modal K correct).
- **Scope, stated plainly**: the orthogonal factor `O` of `A = OU` (FS
  Lemma 1) is fixed at the identity. FS Sec 3.3 explains that under
  skewness `O` determines the *directions* of asymmetry, so this release
  ties those directions to the coordinate axes after the triangular
  transform. A Householder-parameterised `O` (FS Appendix A) is planned
  as a follow-up, as is the skew multivariate independent-Student
  family.

### Breaking change: one skewness convention across all Fernandez-Steel families

Previously `fssn` and `fsst` parameterised skewness as
`alpha = 1/gamma`, while `fossep` used `alpha = gamma`. The same
`alpha = 2` therefore skewed *left* in `fssn`/`fsst` and *right* in
`fossep`. Fits and simulations were internally consistent within each
family, so this never produced wrong inference, but it made the families
incomparable and the parameter uninterpretable across them.

All Fernandez-Steel families (`fssn`, `fsst`, `fossep`) now share the
convention of Fernandez & Steel (1998, 2007): the exported `alpha`
**is** the FS skewness `gamma`, so

- `alpha = 1` is symmetric,
- `alpha > 1` skews right, `alpha < 1` skews left,
- `P(X > mu) = alpha^2 / (1 + alpha^2)` exactly, in every family.

**What this means for you.** `fssn` and `fsst` results computed with
1.0.1 or earlier correspond to the reciprocal `alpha` under 1.1.0: an
estimate of `alpha = 2` then is `alpha = 0.5` now. Densities, CDFs,
quantiles, RNG and the compiled NIMBLE kernels were all updated
together, and a regression test (`test-skew-convention.R`) now pins the
guarantee. `fossep` is unchanged. The Jones-Faddy family (`jfst`) uses
its own `alpha`/`theta` shapes and is not affected; neither are the
MSNBurr families.

## nimix 1.0.1

### Batch B: six new neo-normal component families

All six are univariate, non-conjugate, and available under the finite-K,
DPM and MRF engines, with `d/p/q/r` functions and
[`ppCheck()`](https://madsyair.github.io/nimix/reference/ppCheck.md)
support. Each NIMBLE kernel is built and registered in the global
environment (the pattern required for scalar user-defined densities) and
was validated against the R reference to 1e-8, with the density
integrating to one and the family’s known reduction reproduced exactly.

- `distribution = "sep"` – symmetric exponential power (`nu = 2` Normal,
  `nu = 1` Laplace).
- `distribution = "lep"` – exponential power under the alternative
  parameterisation (`nu = 2` Normal).
- `distribution = "fssn"` – Fernandez-Steel skew Normal (`alpha = 1`
  Normal), with a log-normal prior on `alpha` that treats left/right
  skew symmetrically.
- `distribution = "fossep"` – Fernandez-Steel skew exponential power
  (`theta = 2` skew-Normal kernel).
- `distribution = "fsst"` – Fernandez-Steel skew Student-t (`alpha = 1`
  symmetric-t). `nu` is a stochastic node truncated below at 2 so the
  variance exists; it is only weakly identified by the data. The
  t-kernel is inlined rather than calling NIMBLE’s `dt`.
- `distribution = "jfst"` – Jones-Faddy skew-t (`alpha = theta`
  symmetric). The density uses the branch-free identity
  `sign(z)/sqrt((a+th)/z^2 + 1) == z/sqrt(a + th + z^2)`, which is also
  finite at `z = 0`.

The `q*` functions for these families recycle vector parameters to the
sample length and subset them per branch, so per-observation parameter
vectors (as used by posterior predictive simulation) align without
recycling warnings.

## nimix 1.0.0 (in development)

### Bug fixes (installed-package correctness)

- MRF engine: the Potts prior’s `dPottsNimix` / `rPottsNimix` are now
  built and registered in the global environment (like the scalar
  neo-normal densities). Registering them from the package namespace
  made NIMBLE fail to find `rPottsNimix` during code generation for the
  latent label node once the package was installed
  ([`library(nimix)`](https://github.com/madsyair/nimix)), so every MRF
  fit errored under a normal install while working under `load_all()`.
  Fixed.
- [`ppCheck()`](https://madsyair.github.io/nimix/reference/ppCheck.md)
  for MSNBurr / MSNBurr-IIa / GMSNBurr no longer emits recycling
  warnings: the quantile functions now recycle vector parameters to the
  sample length and index the interior subset consistently, so
  per-observation parameter vectors (as used by posterior predictive
  simulation) align.
- Multivariate cluster
  [`summary()`](https://rdrr.io/r/base/summary.html) columns are now
  `mu_<j>_mean` (with `_med` / `_lwr` / `_upr`) for consistency with the
  univariate summaries.

### New neo-normal family: GMSNBurr

- `distribution = "gmsnburr"`: the generalized MSNBurr (Iriawan 2000;
  Choir
  2020. with two shape parameters `alpha` and `theta`. `theta = 1`
        recovers MSNBurr, `alpha = 1` recovers MSNBurr-IIa,
        `alpha = theta` is symmetric, and `alpha = theta -> inf`
        converges to the Normal. Available under all three engines
        (finite-K, DPM, MRF), numerically stable (NIMBLE density matches
        the R reference to 1e-8; density integrates to one; exact
        reduction to MSNBurr / MSNBurr-IIa verified). `d/p/q/r`
        functions and
        [`ppCheck()`](https://madsyair.github.io/nimix/reference/ppCheck.md)
        support included.

### Parallel chains

- `mcmcControl` gains `parallel = TRUE` (with `nchains > 1`) to run
  chains in parallel via
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html). Each
  worker builds and compiles its own model in a separate directory – the
  only fork-safe way to parallelise NIMBLE, avoiding the
  shared-C++-object and temp-directory collisions that a naive
  `mclapply` would hit. Forking only (Unix/macOS; Windows falls back to
  sequential); `ncores` caps the worker count. Verified to recover the
  same label-invariant summary as the sequential path.

### Bug fix: robust source load order

- Added `@include` directives across the S4 class and spec files so the
  package’s `Collate` order is derived topologically from the actual
  class dependencies. This fixes a load failure
  (`undefined slot classes in definition of "MRFEngine"` /
  `no definition found for superclass "NormalRegSpec"`) that could occur
  whenever files were sourced in alphabetical rather than `Collate`
  order.

### Cluster profiling, richer summaries, model selection and ensembling

- [`summary()`](https://rdrr.io/r/base/summary.html) of a relabelled
  clustering fit now reports the posterior **median** alongside the mean
  and 95% credible interval for every scalar component parameter (all
  univariate families) and per dimension for multivariate clustering.
- New
  [`clusterProfile()`](https://madsyair.github.io/nimix/reference/clusterProfile.md):
  assigns each observation to its MAP cluster and describes the observed
  data within each cluster (size, proportion, and per-variable mean / sd
  / median) – the data-side complement to
  [`summary()`](https://rdrr.io/r/base/summary.html) (which reports
  fitted parameters) and `plot(fit, "cluster")` (which shows the
  partition).
- New predictive model-selection layer built on the label-invariant
  pointwise mixture log-likelihood:
  [`nimixWAIC()`](https://madsyair.github.io/nimix/reference/nimixWAIC.md)
  (native; Watanabe 2010),
  [`nimixLOO()`](https://madsyair.github.io/nimix/reference/nimixLOO.md)
  (PSIS-LOO via the **loo** package; Vehtari, Gelman & Gabry 2017), and
  [`modelSelect()`](https://madsyair.github.io/nimix/reference/modelSelect.md)
  to rank several fits (e.g. choosing K, or comparing Normal / Student-t
  / MSNBurr components on the same data).
- New
  [`ensembleFit()`](https://madsyair.github.io/nimix/reference/ensembleFit.md):
  combines several fits into one weighted predictive model via Bayesian
  stacking or Pseudo-BMA+ (Yao et al. 2018, needs **loo**) or
  Akaike-style WAIC weights (native).
  [`predict()`](https://rdrr.io/r/stats/predict.html) returns the
  ensemble-weighted density. `loo` added to Suggests.

### Neo-normal components: MSNBurr and MSNBurr-IIa

- Two new univariate component families, `distribution = "msnburr"` and
  `"msnburr2a"` (Iriawan 2000; Choir 2020), with a location, scale, and
  a skewness shape `alpha` (`alpha = 1` is exactly the logistic
  distribution; MSNBurr accommodates left skew, MSNBurr-IIa its mirror).
  Available under all three engines: finite-K, DPM, and the spatial MRF.
- Numerically stable throughout, using the maintainer-contributed
  reference implementation: an asymptotic log-omega branch for
  `alpha -> 0`, a branch-free softplus in the compiled densities, and
  two-branch quantile inversion. Log-densities stay finite for
  standardized values in the hundreds; the density integrates to one
  from `alpha = 0.05` to `alpha = 100`; the NIMBLE and R densities agree
  to 1e-8; and
  [`ppCheck()`](https://madsyair.github.io/nimix/reference/ppCheck.md)
  supports both families.
- Verified recovery on synthetic skew-component mixtures (finite-K
  allocation accuracy 1.00 / 0.99, DPM recovers the true K, MRF accuracy
  1.00).

### Bayesian workflow: modern mixture-appropriate convergence + predictive checks

### Bayesian workflow: modern mixture-appropriate convergence + predictive checks

- Full Vehtari et al. (2021) convergence suite – rank-normalized
  split-Rhat (already present), **folded split-Rhat, bulk-ESS and
  tail-ESS** – computed over **label-invariant functionals only**
  (occupied K, DPM alpha, MRF beta, and the new per-iteration allocation
  entropy): per-component traces are not identified under label
  switching, so chain diagnostics on them would be meaningless.
  [`summary()`](https://rdrr.io/r/base/summary.html) prints the
  functional table with the documented Rhat \< 1.01 aim; helpers are
  unit-tested against theory (iid ESS ~ n, folded Rhat detects pure
  scale disagreement that location Rhat misses, AR(0.9) ESS matches the
  analytic value).
- New
  [`ppCheck()`](https://madsyair.github.io/nimix/reference/ppCheck.md):
  posterior predictive checks (Gelman, Meng & Stern 1996; Gelman et
  al. 2020, Bayesian workflow) with replicates drawn conditionally on
  each draw’s fitted allocation – label-invariant by construction and
  valid for every engine including the spatial MRF. Built-in statistics
  (mean, sd, min, max, skew) or user functions; clustering families
  (Gaussian uv/mv, Student-t, Normal-Gamma, Poisson, Binomial)
  supported, regression checks planned. Verified to pass a
  well-specified Gaussian mixture and to flag a Poisson mixture on
  overdispersed counts (tail p = 0.013).

## nimix 0.9.0

### MRF engine across the distribution registry (batch 2) – matrix complete

- `method = "mrf"` now covers **every registered family**: added
  multivariate Student-t clustering, multivariate-response regression
  (Gaussian, Student-t and Normal-Gamma responses), and the augmented
  Normal-Gamma routes (univariate/multivariate clustering and univariate
  regression), each with its own O(nK + edges) label sweep. Harness
  groups `mrfdist2` and `mrfdist2reg` pass on the synthetic lattice (7
  combinations x 3 seeds).
- **Critical correctness fix.** The hand-written MRF kernels for the
  Gaussian specs shadowed the generic `.pottsify()` route for their
  heavy-tailed subclasses (S4 inheritance), so e.g. a Student-t
  regression under `mrf` silently used the Gaussian kernel for its
  parameter updates. The hand-written kernels were removed: every
  family’s MRF kernel is now derived from its own fixed-K kernel, and a
  regression test pins kernel/family correspondence.
- An internal quiet-mode handler crashed on NIMBLE conditions with empty
  messages (“argument is of length zero”); it is now robust, and the
  benign “No samplers assigned” notice for the deliberately unsampled
  fixed `beta` node is muffled.
- **Mixing benchmark (documented commitment):** on identical
  heavy-tailed lattice data, the augmented Normal-Gamma route mixes the
  component means BETTER than the direct Student-t density under the MRF
  (ESS of the minimum component mean 1357 vs 550 at equal wall time) –
  the reverse of the exchangeable-DPM finding, because with K fixed and
  labels pinned by the spatial field the omega-partition coupling
  penalty largely disappears while the conjugate Gibbs updates win.
  Guidance: prefer `normalgamma` for MRF parameter mixing; the direct-t
  route remains available and equally valid.

## nimix 0.8.0

### MRF engine across the distribution registry (batch 1: closed-form emissions)

- Per the maintainer-approved feasibility study, `method = "mrf"` now
  covers six more family/task combinations: **Poisson** and **Binomial**
  clustering, **Poisson-GLM** (log link) and **Binomial-GLM** (logit
  link) regression, and **direct Student-t** clustering and regression.
  Every DSL emission call was verified empirically against reference
  densities before implementation.
- New internal `.pottsify()` transformation derives any family’s MRF
  kernel mechanically from its fixed-K kernel (drop the
  Dirichlet-categorical label layer, insert the joint Potts node), so a
  default `buildModelCode(<any spec>, MRFEngine)` now exists; the three
  original hand-written MRF kernels are unchanged. Label sweeps remain
  family-specific for O(nK + edges) performance.
- `method = "mrf"` now requires `K >= 2` (a one-state Potts field is
  degenerate). The augmented Normal-Gamma routes remain explicitly
  blocked pending batch 2 (with a planned mixing benchmark against the
  direct-t routes).
- Recovery harness gains `mrfdist1` and `mrfdist1reg` groups: all six
  new combinations pass on the synthetic two-block lattice across 3
  seeds (parameter recovery + allocation accuracy 0.98-1.00 in
  verification runs).

### Spatial line hardening + official-statistics case study

- Two new packaged official-statistics datasets: `usStates2023` (SAIPE
  2023 median household income and poverty rates for the 48 contiguous
  states +
  600. and `usStateAdj` (the state contiguity matrix derived from the
       Census Bureau’s 2023 county adjacency file); full provenance in
       the help pages.
- New vignette `spatial-mixture.Rmd`: on the 2023 SAIPE poverty rates
  the MRF engine (with estimated interaction, posterior mean ~1.2) finds
  two regimes whose high-poverty regime is the spatially contiguous
  Southern belt (AL AR DC FL GA LA MS NC NM OK SC TN TX) – structure an
  exchangeable mixture cannot represent; a plain `fixedk` fit agrees on
  only 82% of states. The spatial regression example honestly collapses
  to a single national income-poverty regime.
- The recovery harness gains an `mrfbeta` group (pseudo-likelihood beta
  estimation on the synthetic lattice, 3 seeds).

## nimix 0.7.0

### MRF engine: Bayesian estimation of the interaction beta

- `prior$estimateBeta = TRUE` (with optional `prior$betaMax`, default 2)
  puts a uniform prior on the Potts interaction and updates it by
  random-walk Metropolis against the Besag (1975) **pseudo-likelihood**
  – the classical approximate route for hidden Potts fields, since the
  exact posterior of beta is doubly intractable. This is documented as
  an approximation; an exchange-algorithm refinement (Murray, Ghahramani
  & MacKay 2006) is a possible future upgrade. `beta` is now a model
  node in both modes; with `estimateBeta = FALSE` (default) no sampler
  touches it, reproducing the fixed-beta behaviour exactly.
- Correctness guards: NIMBLE’s default sampler on `beta` (which would
  target the unnormalised Potts density – wrong in beta) is removed
  unconditionally, and the compiled-model cache now distinguishes
  estimation from fixed mode (they compile different sampler sets).
- MRF-specific diagnostics:
  [`summary()`](https://rdrr.io/r/base/summary.html) reports split-Rhat
  / ESS and the posterior mean of `beta` when it is estimated, and the
  engine warns when the posterior piles up near `betaMax` (a
  near-saturated field).
- On the synthetic two-block lattice, the estimated interaction
  concentrates well above zero (posterior mean ~1.6, P(beta \> 0.2) = 1)
  with perfect block recovery, while the fixed-beta path is unchanged.

### MRF engine: spatially clustered regressions

- [`nimixReg()`](https://madsyair.github.io/nimix/reference/nimixReg.md)
  gains `method = "mrf"` and a `spatialWeights` argument: a mixture of
  Gaussian linear regressions whose latent regime labels follow the
  Potts field on a neighbourhood graph – regression coefficients that
  cluster spatially (e.g. growth patterns across adjacent regions). Same
  fixed `prior$beta` interaction as the clustering engine.
- On a synthetic two-block lattice with opposite slopes (+2 / -2) and
  heavy noise, the spatial smoothing lifts the regime-allocation
  accuracy from 0.867 (`fixedk`) to 1.000 while recovering both slopes.
  The recovery harness gains an `mrfreg` group (3 seeds).
- Heavy-tailed regression responses are explicitly blocked under `mrf`
  with a clear message (they inherit from the Gaussian regression spec
  and would otherwise use the wrong emission density).

### MRF engine: multivariate Gaussian components

- `method = "mrf"` now also accepts multivariate data: multivariate
  Gaussian components under the Normal-Inverse-Wishart kernel, with the
  label sweep evaluating `dmnorm_chol` per component (Cholesky factors
  hoisted once per sweep). On a synthetic two-block lattice the spatial
  smoothing lifts the multivariate allocation accuracy from 0.942
  (`fixedk`) to 0.992.
- Sampler selection is fully polymorphic (S4 dispatch on the component
  spec); heavy-tailed families that inherit from the Gaussian specs are
  explicitly blocked with a clear message instead of silently using the
  wrong emission density.

## nimix 0.6.0

### Spatially constrained mixtures: the MRF engine

- New `method = "mrf"` in
  [`nimixClust()`](https://madsyair.github.io/nimix/reference/nimixClust.md):
  a finite mixture whose latent labels follow a Potts Markov random
  field on a `spatialWeights` neighbourhood (Potts 1952; Besag 1974;
  spatially variant finite mixtures, Blekas et al. 2005), so
  neighbouring observations favour the same component. Requires a known
  `K` and a `SpatialWeightSpec`; the interaction strength is fixed at
  `prior$beta` (default 0.8; `beta = 0` removes the smoothing).
  Univariate Gaussian components in this release.
- Implementation: an (intentionally unnormalised) user-defined NIMBLE
  Potts distribution – exact for MCMC because `beta` is fixed, so the
  intractable partition function is constant – plus a custom single-site
  Gibbs sweep sampler over the labels; theta updates remain conjugate.
  Bayesian estimation of `beta` is planned for a later 1.x release.
- On a synthetic two-block lattice with overlapping components, spatial
  smoothing lifts the allocation accuracy from 0.883 (`fixedk`, no
  smoothing) to 0.967, and `beta = 0` reproduces the unsmoothed
  behaviour. The recovery harness gains an `mrf` group (3 seeds, all
  passing).

Opens the post-1.0 spatial line: spatially constrained mixtures in which
the latent component labels follow a Markov random field on a
neighbourhood graph (Besag 1974; spatially variant finite mixtures,
Blekas et al. 2005).

- New S4 class `SpatialWeightSpec`: a validated neighbourhood structure
  (symmetric, zero-diagonal, non-negative weight matrix) deliberately
  orthogonal to `DistributionSpec`, so any registered component family
  can be paired with any graph. Constructors
  [`spatialWeights()`](https://madsyair.github.io/nimix/reference/spatialWeights.md)
  (from a matrix) and
  [`gridAdjacency()`](https://madsyair.github.io/nimix/reference/gridAdjacency.md)
  (rook/queen contiguity on a regular lattice); accessors
  [`nRegions()`](https://madsyair.github.io/nimix/reference/nRegions.md),
  [`getAdjacency()`](https://madsyair.github.io/nimix/reference/getAdjacency.md),
  [`neighborsOf()`](https://madsyair.github.io/nimix/reference/neighborsOf.md).
- [`nimixClust()`](https://madsyair.github.io/nimix/reference/nimixClust.md)
  gains a `spatialWeights` argument (default `NULL`, fully
  backward-compatible). Supplying a structure validates it and points to
  the MRF engine planned for 0.6.0; the exchangeable mixture remains the
  default behaviour.

First stable release. The public API –
[`nimixClust()`](https://madsyair.github.io/nimix/reference/nimixClust.md),
[`nimixReg()`](https://madsyair.github.io/nimix/reference/nimixReg.md),
[`registerDistribution()`](https://madsyair.github.io/nimix/reference/registerDistribution.md),
[`relabel()`](https://madsyair.github.io/nimix/reference/relabel.md),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`nimixClearCache()`](https://madsyair.github.io/nimix/reference/nimixClearCache.md)
– is now considered stable; breaking changes will bump the major
version.

### Official-statistics data and vignettes

- New packaged dataset `wdi2022`: four development indicators for 207
  countries (World Bank World Development Indicators, 2022; CC BY 4.0;
  full provenance in
  [`?wdi2022`](https://madsyair.github.io/nimix/reference/wdi2022.md)).
- All four vignettes now run on this official-statistics data instead of
  simulations only: income-regime clustering (univariate), joint
  income-longevity clustering (multivariate), the Preston curve as a
  mixture of regressions (an instructive single-regime result: the DPM
  does not invent components), and the Student-t vs Normal-Gamma
  heavy-tail comparison. Reported numbers in the vignettes come from
  actual runs.

## nimix 0.5.0

This line opens the performance and hardening phase, built on the 0.4.3
feature set. The two production engines are `method = "dpm"` (Dirichlet
process / Chinese restaurant process; the number of occupied components
is estimated) and `method = "fixedk"` (finite mixture with a known `K`).

### Implemented in 0.9.x

- **Compile-once / reuse.** A fit reuses a compiled NIMBLE model (and
  its MCMC) when a later fit has an identical structure – same generated
  code, constants, monitors and component family – resetting only data
  and initial values. This skips recompilation for repeated fits
  (multiple seeds, multiple chains) and is bit-for-bit identical to a
  fresh compile. Controlled by `mcmcControl$reuse` (default `TRUE`);
  [`nimixClearCache()`](https://madsyair.github.io/nimix/reference/nimixClearCache.md)
  releases the cached compiled models.
- **Reproducibility.** The dispersed (k-means) initialisation is now
  seeded by the fit’s `seed`, so repeated fits with the same data and
  seed coincide exactly. The caller’s global random stream is left
  untouched.
- **Multi-chain + convergence diagnostics.** `mcmcControl$nchains`
  (default 1) runs several chains from dispersed, separately seeded
  starts, reusing the compiled model so only the first chain compiles.
  [`summary()`](https://rdrr.io/r/base/summary.html) reports
  rank-normalized split-Rhat and effective sample size for
  label-invariant quantities (occupied-cluster count, and the DPM
  concentration), and warns when Rhat exceeds 1.1. Per-component
  parameters are deliberately excluded from Rhat to avoid
  label-switching artefacts.
- **Vectorised post-processing.** The per-draw occupied-cluster count,
  the allocation-trace parsing, and the relabelling recode/weight steps
  are now vectorised (single `tabulate`/matrix-index passes instead of
  per-row loops). Verified bit-for-bit identical to the previous
  implementation; the remaining relabelling cost is the external ECR
  algorithm itself.
- **Hardening harness.** `inst/harness/run_recovery_suite.R` is now a
  systematic distribution x engine recovery matrix: every released
  clustering family (Gaussian, Student-t, Normal-Gamma – univariate and
  multivariate – Poisson, Binomial) and the regression path are fitted
  with BOTH engines on data with known truth, across three MCMC seeds
  (18 combinations, 54 fits). For the DPM on discrete counts the
  recovery criterion assesses the dominant components (weight \>= 0.1),
  reflecting the known diffuseness of the DPM posterior on the number of
  components for such data (Miller & Harrison, 2013). Previously
  untested engine pairings are additionally pinned in
  `tests/testthat/test-hardening-matrix.R`.

## nimix 0.4.3

### Robustness and ergonomics

- `verbose` now defaults to `FALSE` for
  [`nimixClust()`](https://madsyair.github.io/nimix/reference/nimixClust.md)
  and
  [`nimixReg()`](https://madsyair.github.io/nimix/reference/nimixReg.md).
  Benign NIMBLE configuration chatter (e.g. the
  Chinese-restaurant-process truncation reminder) is muffled
  selectively; genuine warnings about potentially invalid MCMC draws and
  all errors always propagate, in both verbose modes.
- Dispersed k-means initialisation is now headroom-aware: the number of
  seeded clusters respects the truncation level (`K_max`) so early CRP
  transients no longer breach it on smaller `K_max`.
- The initialisation headroom is configurable through
  `mcmcControl$initRatio` (default 0.8). Values that leave little
  headroom (\>= 0.95) warn; values outside `(0, 1)` error.

### Engines and distributions

- Two production engines: `method = "dpm"` and `method = "fixedk"`,
  available for univariate and multivariate clustering and for
  regression (including multivariate responses).
- Component distributions: Gaussian (univariate / multivariate),
  Student-t and Normal-Gamma (heavy-tailed, univariate / multivariate),
  and Poisson / Binomial counts.

## nimix 0.4.2

- [`nimixReg()`](https://madsyair.github.io/nimix/reference/nimixReg.md)
  gains multivariate responses (`cbind(y1, y2) ~ x`) for Normal,
  Student-t and Normal-Gamma components, with per-component coefficient
  matrices and error covariances.

## nimix 0.4.0

- Student-t and Normal-Gamma components (univariate and multivariate)
  and Poisson / Binomial counts; public
  [`registerDistribution()`](https://madsyair.github.io/nimix/reference/registerDistribution.md).

## nimix 0.3.0

- [`nimixReg()`](https://madsyair.github.io/nimix/reference/nimixReg.md)
  and the `RegressionMixModel` class: mixture-of-regressions with a
  Normal-Inverse-Gamma g-prior. `FixedKEngine` implemented across
  univariate, multivariate and regression models. Engine selection is
  polymorphic via
  [`runEngine()`](https://madsyair.github.io/nimix/reference/runEngine.md).

## nimix 0.2.0

- Multivariate Gaussian clustering (`NormalMvSpec`) with a
  Normal-Inverse-Wishart base measure. Engine generalised to be
  dimension-agnostic.

## nimix 0.1.0

- S4 foundation, univariate Gaussian clustering (`NormalUvSpec`),
  [`nimixClust()`](https://madsyair.github.io/nimix/reference/nimixClust.md)
  on the DPM and fixed-K engines.
