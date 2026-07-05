# Changelog

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
