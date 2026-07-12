<!-- README.md is generated from README.Rmd if you keep one; this file is the
     pkgdown home page and the GitHub landing page. -->

# nimix <img src="man/figures/logo.png" align="right" height="139" alt="nimix logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/madsyair/nimix/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/madsyair/nimix/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/madsyair/nimix/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/madsyair/nimix/actions/workflows/pkgdown.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**nimix** is a Bayesian mixture-modelling package built on top of
[NIMBLE](https://r-nimble.org/). It provides mixture *clustering* (univariate
and multivariate), mixtures *of regressions* (univariate and multivariate
response), and *spatially coupled* mixtures, through an extensible S4
`DistributionSpec` contract: new component families plug into every inference
engine without touching the engine code. The registry currently holds **32
component distributions**.

## Component families

- **Classical** — Normal (univariate and multivariate), Student-t (both),
  Normal-Gamma heavy-tailed (both), Poisson, Binomial.
- **Neo-normal univariate** — MSNBurr, MSNBurr-IIa, GMSNBurr, the exponential
  power pair `sep`/`lep`, and the Fernandez–Steel / Jones–Faddy skew families
  `fssn`, `fossep`, `fsst`, `jfst`. All follow one skewness convention:
  `alpha` *is* the FS `gamma`, so `alpha > 1` skews right and
  `P(X > mu) = alpha^2 / (1 + alpha^2)`, uniformly across families.
- **Ferreira–Steel skew multivariate** — skew-Normal and skew-independent-
  Student components (`"skewnormal-mv"`, `"skewistudent-mv"`), with the
  orthogonal factor of `A = OU` either fixed or *estimated*
  (`"skewnormal-mv-o"`, `"skewistudent-mv-o"`). The estimated-O families
  accept any dimension `m >= 2`: `m = 2` uses a dedicated implementation,
  `m > 2` a general Householder parameterisation with `m(m-1)/2` angles and
  FS's identifiability restriction applied as a post-hoc canonicalisation.
  The general-`m` variants are **experimental** — partition and location
  recovery are robust, but individual angles are large-sample quantities and
  the reported `gamma`/`O` are canonicalised (compare against canonicalised
  truth; see `?canonicaliseO`).

## Inference engines

Selected with `method`:

- `"dpm"` — Dirichlet Process Mixture that *infers* the number of components
  (NIMBLE's CRP), bounded by `K_max`.
- `"fixedk"` — finite mixture with known `K` (including the `K = 1` baseline
  for model comparison).
- `"mrf"` — spatially coupled allocations via a Potts prior over a
  user-supplied adjacency (`spatialWeights`; `gridAdjacency()` builds regular
  grids).
- `"hmm"` — regime switching in time: labels follow a first-order Markov
  chain, the state path is marginalised out by the forward algorithm
  (measured: min ESS/sec 456 vs 144 for naive latent-state sampling), and
  exact allocation draws are recovered post hoc by FFBS -- so `relabel()`,
  `psm()`, `binderPartition()` and the plots work unchanged.
  `viterbiPath(fit)` decodes the jointly most probable state sequence.
  Current emissions (univariate): Gaussian, Student-t (heavy tails),
  Poisson (count regimes), and the neo-normal skewed families MSNBurr,
  MSNBurr-IIa and GMSNBurr; further families follow the gated roadmap.

All 32 families run under the dpm, fixedk and mrf engines.

## Sampler defaults that were measured, not assumed

- FixedK mixture regression uses an **exact Normal-Inverse-Gamma Gibbs**
  update for coefficients and residual variances (NIMBLE's conjugacy checker
  cannot see through the mixture's dynamic indexing). This makes the fit
  scale-equivariant in the predictors — measured slope discrepancy under a
  1000× rescale fell from 0.48 to 0.004 — and raised coefficient ESS/second
  ~35×. The DPM path was already conjugate via NIMBLE's CRP machinery.
- The nine 3–4 parameter univariate families sample each component's
  parameters as an **automated-factor-slice block** (`AF_slice`; Tibbits et
  al. 2014). On `fssn`, where `cor(mu, alpha)` reaches −0.94, minimum ESS
  went from 12 to 622 per 1500 draws; a naive `RW_block` made it *worse*.

## Post-processing and interoperability

- **Component parameters:** `relabel()` (post-hoc label-switching correction
  via the [`label.switching`](https://CRAN.R-project.org/package=label.switching)
  package, conditioning on the modal number of clusters), then `summary()`,
  `plot()`, `predict()`.
- **Partitions, label-free:** `psm()` (posterior similarity matrix) and
  `binderPartition()` (Dahl's least-squares / expected-Binder-loss point
  partition) use *every* draw — no relabelling, no fixed `K` required. They
  complement `relabel()`: it answers "what are the component parameters",
  these answer "which observations belong together".
- **Internal validity, with an honest caveat:** `clusterValidity(fit)`
  computes silhouette, Dunn and Calinski-Harabasz on the Binder partition
  (via `cluster`/`fpc`, Suggests-only). These indices reward *geometric*
  separation while mixtures are *density*-based: overlapping components can
  be exactly the right model and still score low -- use them to compare
  partitions, not to judge model adequacy.
- **Predictive checking:** `ppCheck()` for tail-probability summaries;
  `posteriorPredict()` (or `ppCheck(..., store_yrep = TRUE)`) for the
  replicates themselves.
- **bayesplot, without the dependency:** `drawsArray(fit)` returns a plain
  `iterations × chains × parameters` array that `bayesplot::mcmc_*` accept
  natively, and `ppcData(fit)` returns `list(y, yrep)` for
  `bayesplot::ppc_dens_overlay()` and friends. bayesplot sits in `Suggests`
  only. The adaptors enforce one statistical guard: per-component draws are
  refused before `relabel()`, because under label switching an R-hat on a raw
  `muTilde` trace looks valid and means nothing — the default `"invariant"`
  view (cluster count, allocation entropy, `alpha`) is what is safe on raw
  draws.
- Every `plot(fit, type = ...)` invisibly returns the tidy data frame it
  drew, so you can replot with ggplot2/lattice/plotly without nimix carrying
  those packages.

## Installation

`nimix` depends on `nimble`, which compiles model code, so you need a working
C/C++ toolchain (Rtools on Windows, Xcode CLT on macOS, build-essential on
Linux).

``` r
# install.packages("pak")
pak::pak("madsyair/nimix")

# or
# install.packages("remotes")
remotes::install_github("madsyair/nimix")
```

## Quick start

``` r
library(nimix)

## Univariate clustering, number of components estimated (DPM)
set.seed(1)
y <- c(rnorm(100, -3), rnorm(100, 3))
fit <- nimixClust(y, K_max = 8,
                  mcmcControl = list(niter = 2000, nburnin = 1000),
                  verbose = FALSE)
summary(fit)
plot(fit, type = "K")

## Label-free partition summary using every draw
bp <- binderPartition(fit)
table(bp$partition)

## Skew multivariate clustering with an estimated orthogonal factor (m = 3)
Y <- rbind(matrix(rnorm(300, -4), 100, 3), matrix(rnorm(300, 4), 100, 3))
sk <- relabel(nimixClust(Y, K = 2, method = "fixedk",
                         distribution = "skewnormal-mv-o",
                         mcmcControl = list(niter = 2500, nburnin = 1000)))

## Regime-switching time series (hidden-Markov engine); data order = time order
zt  <- cumsum(c(1, runif(299) < 0.07)) %% 2 + 1     # two persistent regimes
yts <- rnorm(300, c(-2, 2)[zt], 0.7)
hm  <- nimixClust(yts, K = 2, method = "hmm",
                  mcmcControl = list(niter = 4000, nburnin = 1500))
viterbiPath(hm)                        # decoded regime per time point

## Mixture of linear regressions, fixed number of regimes (finite mixture)
x  <- runif(200, -3, 3)
g  <- rep(1:2, each = 100)
df <- data.frame(y = ifelse(g == 1, 2 * x, -2 * x) + rnorm(200, 0, 0.7), x = x)
reg <- nimixReg(y ~ x, df, K = 2, method = "fixedk",
                mcmcControl = list(niter = 2000, nburnin = 1000),
                verbose = FALSE)
summary(reg)

## ... with a random intercept for grouped data (fixedk + normal)
# nimixReg(y ~ x, df, random = ~ region, K = 2, method = "fixedk")

## Hand the fit to bayesplot (if installed)
# bayesplot::mcmc_trace(drawsArray(fit))               # invariant functionals
# pd <- ppcData(fit); bayesplot::ppc_dens_overlay(pd$y, pd$yrep)
```

## Documentation

Function reference and vignettes are published at
<https://madsyair.github.io/nimix/>. `NEWS.md` documents each release,
including the measured numbers behind sampler-default decisions and two
breaking-change notes (the v1.1.0 skewness-convention harmonisation, with its
reciprocal-`alpha` migration).

## License

GPL-3.
