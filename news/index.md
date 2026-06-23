# Changelog

## nimix 0.4.3

Bug-fix and packaging/infrastructure release. No new modelling features;
this clears `R CMD check` and adds continuous integration and a
documentation site ahead of the v0.5.x reversible-jump work.

### Fixed

- **Namespace imports under `R_CHECK_DEPENDS_ONLY`.** Added roxygen
  `@import` / `@importFrom` tags (notably `@import methods`) on the
  package doc so a regenerated `NAMESPACE` keeps its import directives.
  Previously a `document()`-regenerated `NAMESPACE` dropped them, so
  `.onLoad()`’s `new("NormalUvSpec")` failed with *could not find
  function “new”* when the package was loaded with only its stated
  dependencies.

- **Poisson / Binomial regression recovery.** Replaced the global-GLM
  initialisation (one pooled fit copied to every cluster, k-means on the
  raw response) with a dispersed-slope start plus one hard
  classification step. Components that differ only by the sign of their
  slope (crossing regression lines) are now separated at initialisation,
  so the finite-mixture sampler recovers the per-component slopes.

- **Single-component finite mixture (`K = 1`).**
  `nimixReg(..., K = 1, method = "fixedk")` (and the clustering
  analogue) failed to build with *inconsistent dimensionality provided
  for node ‘alphaVec’*: a length-1 Dirichlet parameter was demoted to a
  scalar. The weight-vector dimensions are now declared explicitly, so
  the `K = 1` baseline (useful for WAIC model comparison) builds.

- **Count-response validation.**
  [`nimixReg()`](https://madsyair.github.io/nimix/reference/nimixReg.md)
  /
  [`nimixClust()`](https://madsyair.github.io/nimix/reference/nimixClust.md)
  now reject a continuous response for `distribution = "poisson"` /
  `"binomial"` early with a clear message instead of surfacing a
  numerical underflow inside the sampler.

### Infrastructure

- Added GitHub Actions workflows for `R CMD check` (multi-OS) and a
  `pkgdown` documentation site, a `_pkgdown.yml` reference layout, and
  the corresponding `.Rbuildignore` entries.

### Changed

- **DPM truncation handling.** `K_max` is the dCRP truncation level: the
  sampler errors if the occupied-cluster count ever needs to exceed it.
  Three improvements make this robust and transparent:

  - The default `K_max` (when unspecified) now scales with `n` and
    carries generous headroom above the expected cluster count, instead
    of the old `min(10, floor(n / 5))` which could sit right at the
    operating point.
  - If the sampler still exceeds the truncation, NIMBLE’s raw *“not a
    proper model”* error is translated into an actionable message naming
    `K_max` and suggesting a concrete larger value.
  - After a DPM fit, a warning is issued when the posterior is
    *sustainedly* censored at the truncation (a meaningful share of
    iterations occupy every slot), so a too-tight `K_max` is flagged
    rather than silently distorting the posterior on the number of
    clusters.

- **`verbose` now defaults to `FALSE`** in
  [`nimixClust()`](https://madsyair.github.io/nimix/reference/nimixClust.md)
  and
  [`nimixReg()`](https://madsyair.github.io/nimix/reference/nimixReg.md),
  matching the documented design. Quiet mode is now **selective, not a
  blanket warning suppressor**: only NIMBLE’s known-benign configuration
  chatter (the *“number of clusters … is less than the number of
  potential clusters”* reminder and the *“model is not fully
  initialized”* note) is muffled. Any other warning – in particular
  anything that could indicate the sampler produced invalid draws –
  propagates to the user even when `verbose = FALSE`, as do nimix’s own
  diagnostics (the censored-posterior warning) and any error. Pass
  `verbose = TRUE` for NIMBLE’s full configuration output and a progress
  bar.

- **Headroom-aware cluster initialisation.** The k-means / hard-E-step
  starts used to seed up to `count - 1` clusters (`count` = the
  truncation level `K_max` for the DPM), which sat right against the
  ceiling for a modest `K_max` and left no room for the early CRP
  transient (which briefly occupies more clusters than the modal K
  before merging). Initialisation now seeds at most `floor(0.8 * count)`
  clusters, guaranteeing at least `0.2 * K_max` slots of headroom above
  the dispersed start; the GLM-regression E-step keeps a floor of 2 so
  small fixed-K starts still separate. For large `K_max` the
  `ceil(sqrt(n))` cap still binds, so dispersed-start behaviour there is
  unchanged. This lowers the chance a modest explicit `K_max` breaches
  the truncation; the headroom default remains the primary safeguard.

## nimix 0.4.2

Multivariate-response mixture regression.

### Added – multivariate-response regression (matrix coefficients)

- [`nimixReg()`](https://madsyair.github.io/nimix/reference/nimixReg.md)
  now accepts a matrix response (e.g. `cbind(y1, y2) ~ x`): each mixture
  component is a multivariate linear regression with its own coefficient
  matrix `B_k` (p x d), mean `X B_k`, and a d x d error covariance. The
  variant is chosen automatically from the response shape.
  - `NormalMvRegSpec`: `y ~ N_d(X B_k, Sigma_k)`.
  - `StudentTMvRegSpec`: direct multivariate-t errors via the
    user-defined `dmvt_nimix` density (non-conjugate).
  - `NormalGammaMvRegSpec`: conjugate scale-mixture
    (`omega ~ Gamma(df/2, df/2)`, `y ~ N_d(X B_k, Sigma_k / omega)`),
    with `omega` slice-sampled.
- Cluster prior: the conjugate matrix-Normal-Inverse-Wishart of Backlund
  & Hobert (2020) – inverse-Wishart on the error covariance and, for
  each coefficient row, prior covariance `v0[l] * Sigma` with a g-prior
  among-row scale `v0 = g * diag((X'X)^{-1})`. Conditional on `Sigma`
  the coefficients are Gaussian, so the collapsed CRP cluster sampler is
  kept. [`predict()`](https://rdrr.io/r/stats/predict.html) returns one
  fitted column per response.
- Verified: all three recover the component coefficient matrices under
  DPM and fixed-K; the `dmvt_nimix` kernel and the scale mixture both
  run inside the CRP with dynamically indexed coefficient matrices.

### Benchmark

- `inst/harness/run_benchmark_heavytail.R` and the
  `normal-gamma-vs-studentt` vignette now report measured
  effective-sample-size-per-second. The direct Student-t route mixes the
  partition far better for **clustering** (about 7x), while for
  **regression** the two routes are comparable (the conjugate
  coefficient update offsets the augmentation penalty) – choose by
  measurement, not assumption.

### References

- Zellner (1976); Fernandez & Steel (1999); Backlund & Hobert (2020) for
  multivariate regression with heavy-tailed / scale-mixture errors.

## nimix 0.4.1

Heavy-tailed mixture regression and a regression-routing fix.

### Added – univariate heavy-tail mixture regression

- [`nimixReg()`](https://madsyair.github.io/nimix/reference/nimixReg.md)
  now accepts `distribution = "studentt"` and `"normalgamma"`: linear
  mixture regression with Student-t errors, in the two equivalent
  parameterisations already offered for clustering.
  - `StudentTRegSpec`: direct t kernel `y ~ t(X beta_k, tau_k, df)`
    (non-conjugate).
  - `NormalGammaRegSpec`: conjugate scale mixture
    `omega ~ Gamma(df/2, df/2)`, `y ~ N(X beta_k, s2_k / omega)` – a
    pure Gibbs sampler keeping the conjugate NIG g-prior on the
    coefficients (Geweke 1993), with `omega` slice-sampled and
    interpretable as robustness weights. Both inherit the g-prior and
    trace/relabel handling of `NormalRegSpec`; `df` is a fixed
    hyperparameter (\> 2). Verified: both recover two regression regimes
    under DPM and fixed-K on heavy-tailed, outlier-contaminated data.

### Fixed

- [`isRegressionSpec()`](https://madsyair.github.io/nimix/reference/isRegressionSpec.md)
  now returns `TRUE` for `NormalRegSpec`, so
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) route
  Gaussian linear regression through the regression branch (regression
  introduced in 0.4.0 via the GLM hooks; this restores the routing for
  the Normal-linear case).

### Examined and deferred (theory + feasibility confirmed)

- Multivariate-response heavy-tail regression (multivariate-t /
  Normal-Gamma errors, matrix coefficients) is theoretically supported
  (Zellner 1976; Fernandez & Steel 1999; Backlund & Hobert 2020) and a
  working prototype recovers matrix coefficients under the CRP. It needs
  a matrix-response `RegressionMixModel` (a real architecture step) and
  is proposed as v0.4.2 rather than bundled here, to keep verification
  rigorous.

## nimix 0.4.0

Fourth development release. Roadmap scope: **finalise the public
[`registerDistribution()`](https://madsyair.github.io/nimix/reference/registerDistribution.md),
and add Student-t and Normal-Gamma (univariate AND multivariate) plus
Poisson and Binomial components** – all on the existing DPM and fixed-K
engines, adding only `DistributionSpec` subclasses (no engine changes),
which is the extensibility claim made concrete.

### Added – heavy-tail components (a matched pair, identical marginal)

- `StudentTUvSpec` / `StudentTMvSpec`: Student-t components evaluated as
  a **direct** t density. NIMBLE has no multivariate-t, so the
  multivariate kernel is supplied as a **user-defined NIMBLE
  distribution** (`dmvt_nimix`, registered at load) – the same mechanism
  [`registerDistribution()`](https://madsyair.github.io/nimix/reference/registerDistribution.md)
  offers users. Non-conjugate cluster updates.
- `NormalGammaUvSpec` / `NormalGammaMvSpec`: the **conjugate
  scale-mixture** route to the *same* Student-t / multivariate-t
  marginal, via a per-observation latent precision multiplier `omega`
  (`y ~ N(mu, Sigma/omega)`, `omega ~ Gamma(df/2, df/2)`). Conditional
  on `omega` the kernel is Gaussian, so the conjugate
  Normal-Inverse-Gamma / Normal-Inverse-Wishart cluster updates are kept
  (cheaper per iteration than the direct t density). The latent `omega`
  double as robustness weights. These inherit from `NormalUvSpec` /
  `NormalMvSpec`.
- The marginal equivalence (scale-mixture == analytic t /
  multivariate-t) is verified numerically to ~1e-10, both univariate and
  multivariate.
- `df` is a fixed hyperparameter (default 4 univariate, 5 multivariate;
  must exceed 2). These are the SAME marginal reached two ways, offered
  as distinct choices for their different sampling cost – not different
  distributions, and unrelated to the Normal-Gamma shrinkage prior.

### Added – discrete count components

- `PoissonSpec` (Gamma-Poisson conjugate) and `BinomialSpec`
  (Beta-Binomial conjugate; the number of trials `size` is supplied in
  `prior`). These show the `DistributionSpec` contract is not tied to
  continuous data.

### Distribution selection

- `distribution = "studentt"` / `"normalgamma"` auto-selects the
  univariate or multivariate variant from the data shape (as `"normal"`
  already did); `"poisson"` / `"binomial"` select the count components.
  All work with both `method = "dpm"` and `method = "fixedk"`.
- [`registerDistribution()`](https://madsyair.github.io/nimix/reference/registerDistribution.md)
  (public) registers user `DistributionSpec` objects; `dmvt_nimix`
  demonstrates pairing a custom spec with a custom NIMBLE density.

### Added – GLM mixture-of-regression (link functions)

- [`nimixReg()`](https://madsyair.github.io/nimix/reference/nimixReg.md)
  now accepts `distribution = "poisson"` (**log link**) and `"binomial"`
  (**logit link**, trials in `prior = list(size = ...)`), in addition to
  `"normal"` (identity link). Each component is a GLM with its own
  coefficient vector: Poisson `y ~ Poisson(exp(X beta_k))`, Binomial
  `y ~ Binomial(size, plogis(X beta_k))`. This is the “extend to GLM”
  item for the regression module.
- New `PoissonRegSpec` / `BinomialRegSpec` (g-prior on coefficients).
  The GLM likelihood is non-conjugate, so cluster coefficients are
  updated non-conjugately (no Polya-Gamma augmentation, keeping
  dependencies unchanged and avoiding the augmentation mixing penalty).
  Verified: both recover two regression regimes under DPM and fixed-K;
  [`predict()`](https://rdrr.io/r/stats/predict.html) applies the
  inverse link.
- New polymorphic hooks
  [`isRegressionSpec()`](https://madsyair.github.io/nimix/reference/isRegressionSpec.md)
  and
  [`linkInv()`](https://madsyair.github.io/nimix/reference/linkInv.md)
  route the predict path and apply the component link without class-name
  branching.

### Sampler strategy (Normal-Gamma)

- The latent precision multipliers `omega` are now **slice-sampled**
  (via a new polymorphic
  [`customizeSamplers()`](https://madsyair.github.io/nimix/reference/customizeSamplers.md)
  hook), not random-walk sampled. On a controlled run this raised the
  effective sample size of the number of clusters from ~41 to ~314
  (about 7x). Even so, an empirical comparison found the **direct
  Student-t route mixes the partition markedly better than the
  Normal-Gamma route at equal or lower cost**; Normal-Gamma is preferred
  when the per-observation robustness weights `omega` are wanted, not
  for speed. This corrects the earlier expectation; see the knowledge
  patch.

### Verified (R 4.3.3, NIMBLE 1.4.2)

- Student-t / Normal-Gamma, univariate and multivariate, recover two
  heavy-tailed clusters under DPM and fixed-K; the custom `dmvt_nimix`
  kernel compiles and runs inside the CRP model. Poisson (rates 2/12)
  and Binomial (probs 0.2/0.75, size 20) recover their components.

### Deferred (still planned for the 0.4.x line / later)

- Heavy-tail residuals in
  [`nimixReg()`](https://madsyair.github.io/nimix/reference/nimixReg.md)
  (regression with t / Normal-Gamma errors) and a benchmark vignette
  comparing the two routes are the remaining 0.4.x items. RJMCMC remains
  v0.5.0+ (experimental, conditional).

## nimix 0.3.0

Third development release. Scope follows the official roadmap (project
knowledge, Section 4): **v0.3.0 = `RegressionMixModel` +
[`nimixReg()`](https://madsyair.github.io/nimix/reference/nimixReg.md),
a DPM mixture of linear regressions.**

### Added

- **`FixedKEngine` — finite mixture with a fixed, known number of
  components** (`method = "fixedk"`, argument `K`). Planned since the
  start of the roadmap as the simplest engine; now implemented. Uses a
  symmetric Dirichlet prior on the mixing weights and a categorical
  allocation (no Chinese Restaurant Process), so NIMBLE assigns
  conjugate samplers to the weights and component parameters and the DPM
  truncation reminder does not arise. Works for univariate,
  multivariate, and regression components. Useful as a fast baseline
  when `K` is known, and for classical model selection by comparing fits
  across values of `K`.
- [`nimixClust()`](https://madsyair.github.io/nimix/reference/nimixClust.md)
  /
  [`nimixReg()`](https://madsyair.github.io/nimix/reference/nimixReg.md)
  gained the argument `K` (for `method = "fixedk"`) alongside `K_max`
  (for `method = "dpm"`), with validation that errors clearly if the two
  are swapped.
- `nimixReg(formula, data, ...)`: Bayesian Dirichlet Process Mixture of
  Gaussian **linear regressions**. Each component is ; the number of
  components is inferred via the Chinese Restaurant Process (or fixed
  with `method = "fixedk"`). **Gating is constant** (mixing weights do
  not depend on covariates); covariate-dependent (concomitant) gating is
  a planned opt-in and currently errors with a pointer.
- `NormalRegSpec`: Normal-linear regression component with a conjugate
  **Normal-Inverse-Gamma g-prior**
  (`beta | s2 ~ N(b0, s2 * g * (X'X)^{-1})`, unit-information `g = n`;
  `s2 ~ InvGamma(nu0, s0)` data-scaled to the global OLS residual
  variance, `nu0 >= 3` so the residual variance cannot collapse on small
  components). Verified: NIMBLE assigns the conjugate collapsed CRP
  samplers (`CRP_cluster_wrapper`) to both `betaTilde` and `s2Tilde`.
- `RegressionMixModel` (extends `MixtureModel`) carrying the `formula`
  and the design matrix `X`.
- [`summary()`](https://rdrr.io/r/base/summary.html) reports relabelled
  per-component regression coefficients (named after the design columns)
  and residual variances; `predict(fit, newdata)` returns the posterior
  predictive mean ; `plot(fit, type = "fitted")` shows observed vs
  fitted.

### Changed

- Engines are now selected polymorphically through a `runEngine` generic
  that dispatches on the engine (`DPMEngine` vs `FixedKEngine`); a
  shared helper does the build / compile / run / extract, so the two
  engines share one code path.
- Component initialisation is engine-agnostic (`componentInits`),
  returning an allocation plus component-parameter starts that both
  engines reuse.
- When `verbose = FALSE`, NIMBLE’s compilation notes and the CRP
  truncation reminder are silenced, so quiet runs and examples produce
  clean output.
- Code comments and roxygen no longer cite internal design-document
  section numbers; the explanatory content is kept inline.
- The deterministic-node patterns established for the multivariate model
  are reused for regression: the `s2`-scaled coefficient covariance is
  bound to a node `covBeta`, and dynamically indexed cluster
  coefficients are resolved via `betaObs`.

### Deferred (per roadmap — intentionally NOT in 0.3.0)

- Covariate-dependent (concomitant) gating — planned opt-in (Section
  9.8).
- Student-t / Normal-Gamma / Poisson / Binomial components are
  **v0.4.0**.
- RJMCMC engine is **v0.5.0+** (experimental, conditional — Section
  4.0).

### Verification status

- Pure-R unit tests (g-prior scaling, NIG invariants, `simulateParams`)
  and an end-to-end **unbalanced 80/20** recovery (two regression
  regimes with slopes +2 / -2) pass in R 4.3.3 with NIMBLE 1.4.2: modal
  K = 2, recovered slopes 2.00 / -1.93 and weights 0.78 / 0.22. The full
  CRAN harness should still be run before any release is considered
  CRAN-ready (Section 8, 0.3.4).

## nimix 0.2.0

Second development release. Scope follows the official roadmap (project
knowledge, Section 4): **v0.2.0 = `NormalMvSpec` (multivariate
Gaussian) + multivariate
[`nimixClust()`](https://madsyair.github.io/nimix/reference/nimixClust.md)
on the DPM (dCRP) engine + label switching.**

### Added

- `NormalMvSpec`: multivariate Gaussian component with a conjugate
  **Normal-Inverse-Wishart** cluster base measure (`dmnorm(cov=)` +
  `dinvwish`), data-scaled by default (`mu0 = colMeans`, prior mean of
  the cluster covariance `= cov(data)`, mean dispersion `cLoc`-scaled).
  The inverse-Wishart degrees of freedom are validated to exceed `d + 1`
  so prior draws on empty components are finite and non-singular
  (project knowledge Section 9.3).
- [`nimixClust()`](https://madsyair.github.io/nimix/reference/nimixClust.md)
  now accepts a numeric **matrix** (one row per observation) and routes
  it to `NormalMvSpec`; a numeric vector still routes to `NormalUvSpec`.
  `distribution = "normal"` chooses univariate/multivariate by data
  shape; `"normal-uv"` / `"normal-mv"` force a specific family.
- Dimension-agnostic DPM engine: the engine now delegates all
  dimension-specific work (constants, data, inits, trace parsing,
  relabelled summaries) to new `DistributionSpec` generics —
  [`buildConstants()`](https://madsyair.github.io/nimix/reference/buildConstants.md),
  [`buildDataList()`](https://madsyair.github.io/nimix/reference/buildDataList.md),
  `dpmInits()`,
  [`extractParamTraces()`](https://madsyair.github.io/nimix/reference/extractParamTraces.md),
  [`relabelComponents()`](https://madsyair.github.io/nimix/reference/relabelComponents.md).
  Adding a distribution no longer touches the engine (project knowledge
  Section 3, extensibility seam).
- [`relabel()`](https://madsyair.github.io/nimix/reference/relabel.md)
  generalised: the label-permutation derivation (wrapping
  `label.switching` ECR) depends only on the allocations, so it is
  identical for univariate and multivariate; parameter permutation and
  component summaries are delegated to the spec.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) gains
  `type = "cluster"` (multivariate : scatter of the first two dimensions
  coloured by MAP cluster);
  [`summary()`](https://rdrr.io/r/base/summary.html)/[`predict()`](https://rdrr.io/r/stats/predict.html)
  are dimension-aware (multivariate predictive density via the spec’s
  [`componentDensity()`](https://madsyair.github.io/nimix/reference/componentDensity.md),
  with draw subsampling for tractability).
- Multivariate recovery scenario added to
  `inst/harness/run_recovery_suite.R` and
  `tests/testthat/test-dist-normal-mv.R` (verified to recover two
  well-separated 2-D clusters).

### Changed / fixed

- **`nimble` moved from `Imports` to `Depends`.** NIMBLE’s BUGS-language
  distributions (`dCRP`, `dinvwish`, …) are only resolvable during model
  building when `nimble` is attached, not merely imported (this is the
  convention used by other NIMBLE-extension packages). This fixes a
  latent model-build error (`R function 'dCRP' ... does not exist`).
- Multivariate parameterisation note: the conjugate base measure is
  written with `cov=` + `dinvwish` rather than `prec=` + `dwish`. The
  two are mathematically equivalent, but the precision/`dwish` path
  triggers a Cholesky-lifting failure under dynamic CRP indexing in
  NIMBLE 1.4.x, whereas `cov=` + `dinvwish` compiles, runs, and is
  assigned NIMBLE’s conjugate collapsed CRP samplers (verified
  empirically).
- `ClusterModel` data slot widened to accept matrices while preserving
  their dimensions; univariate behaviour is unchanged.
- The univariate engine/relabel logic was relocated into `NormalUvSpec`
  methods (no behavioural change; required for the dimension-agnostic
  engine).

### Deferred (per roadmap — intentionally NOT in 0.2.0)

- [`nimixReg()`](https://madsyair.github.io/nimix/reference/nimixReg.md)
  (mixture-of-regressions) is **v0.3.0** (stub still errors).
- Student-t / Normal-Gamma / Poisson / Binomial are **v0.4.0**.
- RJMCMC engine is **v0.5.0+** (experimental, conditional — project
  knowledge Section 4.0).

### Verification status

- Pure-R unit tests (spec routing, data-scaled prior, NIW invariants,
  multivariate-normal density, `simulateParams`) and an end-to-end
  multivariate DPM recovery (two 2-D clusters) pass in R 4.3.3 with
  NIMBLE 1.4.2. The full CRAN harness (`R CMD check --as-cran`,
  multi-seed recovery suite, vignette build) should still be run before
  any release is considered CRAN-ready (project knowledge Section 8,
  0.3.4).

## nimix 0.1.0

First development release. Scope follows the official roadmap (project
knowledge, Section 4): **v0.1.0 = S4 foundation + `NormalUvSpec` +
univariate
[`nimixClust()`](https://madsyair.github.io/nimix/reference/nimixClust.md)
on the DPM (dCRP) engine.**

### Added

- S4 architecture foundation:
  - `DistributionSpec` (virtual) with generics
    [`validateParams()`](https://madsyair.github.io/nimix/reference/validateParams.md),
    [`defaultPrior()`](https://madsyair.github.io/nimix/reference/defaultPrior.md),
    [`buildModelCode()`](https://madsyair.github.io/nimix/reference/buildModelCode.md).
  - `NormalUvSpec`: univariate Gaussian component with a conjugate
    Normal-Inverse-Gamma cluster prior (data-scaled by default; see
    project knowledge Section 9.2).
  - `MixtureModel` (virtual) and `ClusterModel`.
  - `EngineConfig` (virtual) and `DPMEngine` wrapping NIMBLE’s native
    `dCRP` / CRP samplers (Neal 2000).
  - `FitResult` with [`summary()`](https://rdrr.io/r/base/summary.html),
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
    [`predict()`](https://rdrr.io/r/stats/predict.html),
    [`relabel()`](https://madsyair.github.io/nimix/reference/relabel.md).
- [`nimixClust()`](https://madsyair.github.io/nimix/reference/nimixClust.md)
  user-facing function (univariate, `method = "dpm"`).
- [`relabel()`](https://madsyair.github.io/nimix/reference/relabel.md)
  wraps `label.switching` (ECR-ITERATIVE-1 default, Papastamoulis &
  Iliopoulos 2010), conditioning on the modal number of occupied
  clusters.
- [`registerDistribution()`](https://madsyair.github.io/nimix/reference/registerDistribution.md)
  /
  [`getDistribution()`](https://madsyair.github.io/nimix/reference/getDistribution.md)
  registry skeleton.
- Recovery-test harness (`inst/harness/run_recovery_suite.R`) and CRAN
  check harness (`inst/harness/run_cran_check.R`).

### Deferred (per roadmap — intentionally NOT in 0.1.0)

- [`nimixReg()`](https://madsyair.github.io/nimix/reference/nimixReg.md)
  (mixture-of-regressions) is **v0.3.0**. A stub is exported that errors
  with an informative message rather than silently reordering the
  roadmap (project knowledge Section 0.3.1).
- Multivariate clustering (`NormalMvSpec`) is **v0.2.0**.
- RJMCMC engine is **v0.5.0+**.

### Naming

- Package and user-facing functions are consistently named `nimix`
  ([`nimixClust()`](https://madsyair.github.io/nimix/reference/nimixClust.md),
  [`nimixReg()`](https://madsyair.github.io/nimix/reference/nimixReg.md)),
  matching the interface signatures in the project-knowledge document
  (Section 3). An earlier draft used the mixed-case spelling `nimMix*`;
  this has been renamed throughout for consistency, so no
  package/document naming discrepancy remains.
