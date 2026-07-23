# cran-comments

## Test environments

* local Ubuntu 24.04, R 4.3.3
* (add win-builder / R-hub / macOS results before submitting)

## R CMD check results

0 errors | 0 warnings | 1 note

### Note: possible error in `integer(0, default = 0)`

```
.nimixEnsureHMM : <anonymous>: possible error in integer(0, default = 0):
  unused argument (default = 0)
```

This is a false positive from `codetools`. The expression appears inside the
argument list of a `nimbleFunction`:

```r
nimble::nimbleFunction(run = function(x = double(1), mu = double(1),
                                      log = integer(0, default = 0)) { ... })
```

In NIMBLE, `double(1)` and `integer(0, default = 0)` are **type declarations**,
not function calls. NIMBLE parses the argument list symbolically to build the
C++ signature; the expressions are never evaluated as R code, so
`base::integer()` is never called with a `default` argument. `integer(0,
default = 0)` is the documented NIMBLE idiom for an optional `log` argument in
a user-defined distribution (see the NIMBLE user manual, "User-defined
distributions"). Static analysis cannot see this, so the note is unavoidable
without abandoning the documented NIMBLE interface.

## Notes seen only in restricted environments

If suggested packages are unavailable, an additional note lists them
(`rmarkdown`, `bayesplot`, `loo`, `cluster`, `fpc`, `sf`). All are used
conditionally, guarded by `requireNamespace()`, and all corresponding tests
and vignette chunks are skipped when the package is absent.

## Tests

`tests/testthat/` holds the automatic suite (distribution mathematics, S4
contracts, prior scaling, validation, diagnostics) and completes in well under
a minute.

MCMC-heavy tests — parameter recovery, engine integration, prediction and
workflow checks — live in `tests/manual/` and are **not** run by `R CMD check`,
because each fits real models and the suite takes tens of minutes. They are run
deliberately before releases via `Rscript tests/manual/run-manual-tests.R`.
