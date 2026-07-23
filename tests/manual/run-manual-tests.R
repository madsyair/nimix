#!/usr/bin/env Rscript
#
# Runner for the manual test suite.
#
# These tests fit real MCMC models -- parameter recovery, engine integration,
# prediction and workflow checks. Each one compiles NIMBLE code and runs
# chains, so the whole suite takes minutes to tens of minutes, far too long
# for `R CMD check`. They are kept out of tests/testthat/ so that the
# automatic suite stays fast, and are run deliberately instead: before a
# release, after touching an engine or a sampler, or when a recovery claim
# needs re-checking.
#
# Usage, from the package root:
#
#   Rscript tests/manual/run-manual-tests.R              # everything
#   Rscript tests/manual/run-manual-tests.R hmm          # files matching "hmm"
#   Rscript tests/manual/run-manual-tests.R engine-mrf reg-hmm
#
# The package is loaded with pkgload::load_all() so the working tree is
# tested, not the installed copy. NOT_CRAN is set to "true" so that
# skip_on_cran() tests actually run -- that is the point of this suite.

suppressMessages({
  if (!requireNamespace("pkgload", quietly = TRUE))
    stop("pkgload is required to run the manual tests.", call. = FALSE)
  if (!requireNamespace("testthat", quietly = TRUE))
    stop("testthat is required to run the manual tests.", call. = FALSE)
})

Sys.setenv(NOT_CRAN = "true")

root <- normalizePath(getwd(), mustWork = FALSE)
if (!dir.exists(file.path(root, "R")) &&
    dir.exists(file.path(root, "..", "..", "R")))
  root <- normalizePath(file.path(root, "..", ".."), mustWork = FALSE)
if (!dir.exists(file.path(root, "R")))
  stop("Run this from the package root (the directory containing R/ and ",
       "tests/), e.g. Rscript tests/manual/run-manual-tests.R", call. = FALSE)
manualDir <- file.path(root, "tests", "manual")
if (!dir.exists(manualDir))
  stop("Cannot locate tests/manual under ", root, ".", call. = FALSE)

suppressMessages(pkgload::load_all(root, quiet = TRUE))

patterns <- commandArgs(trailingOnly = TRUE)
files <- list.files(manualDir, pattern = "^test-.*\\.R$", full.names = TRUE)
if (length(patterns))
  files <- files[vapply(basename(files), function(f)
    any(vapply(patterns, function(p) grepl(p, f, fixed = TRUE), logical(1))),
    logical(1))]

if (!length(files)) {
  message("No manual test files matched.")
  quit(save = "no", status = 0L)
}

message(sprintf("Running %d manual test file(s). This is slow by design.\n",
                length(files)))

t0 <- Sys.time()
failed <- 0L
for (f in files) {
  message("--- ", basename(f))
  ft <- Sys.time()
  res <- testthat::test_file(f, reporter = "summary")
  df <- as.data.frame(res)
  nf <- sum(df$failed) + sum(df$error)
  failed <- failed + nf
  message(sprintf("    %s: %d passed, %d failed, %.0fs\n", basename(f),
                  sum(df$passed), nf,
                  as.numeric(difftime(Sys.time(), ft, units = "secs"))))
}

message(sprintf("Manual suite finished in %.1f min; %d failure(s).",
                as.numeric(difftime(Sys.time(), t0, units = "mins")), failed))
quit(save = "no", status = if (failed > 0L) 1L else 0L)
