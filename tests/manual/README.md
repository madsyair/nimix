# Manual test suite

The tests in this directory are **not** run by `R CMD check` or by
`devtools::test()`. They fit real MCMC models, so each file compiles NIMBLE
code and runs chains; the whole suite takes minutes to tens of minutes.

## Why the split

`tests/testthat/` holds the automatic suite: distribution mathematics
(densities integrate to one, `d`/`p`/`q`/`r` agree), S4 contracts, prior
scaling, validation and error messages, diagnostics and helper functions.
These are cheap, and they are what silently breaks when internals are
refactored, so they run on every check.

This directory holds what is expensive but only needs running deliberately:

* parameter **recovery** for each engine and component family,
* **engine integration** (fixed-K, DPM, MRF/Potts, HMM) end to end,
* **regression** paths, random effects, and GLM/GLMM fits,
* **prediction, forecasting and posterior predictive** workflows,
* model selection, convergence workflow, compile-cache reuse, parallel runs.

Splitting them this way keeps `R CMD check` fast enough for CRAN without
throwing away the tests that actually verify the statistics.

## Running them

From the package root:

```sh
Rscript tests/manual/run-manual-tests.R            # everything
Rscript tests/manual/run-manual-tests.R hmm        # files matching "hmm"
Rscript tests/manual/run-manual-tests.R engine-mrf reg-hmm
```

The runner loads the working tree with `pkgload::load_all()` and sets
`NOT_CRAN=true`, so `skip_on_cran()` tests execute — that is the point.
It exits non-zero if anything fails, so it can be wired into CI.

## When to run them

* before a release,
* after changing an engine, a sampler, or a component distribution,
* whenever a recovery or performance claim needs re-checking.
