## inst/harness/run_cran_check.R
## Layer 1 of the engineering harness (project knowledge Section 8.1).
## Run from the package root in a clean R session with devtools installed.
##
##   Rscript inst/harness/run_cran_check.R
##
## Pass criteria: 0 ERROR, 0 WARNING. NOTEs minimised and each remaining NOTE
## documented (e.g. compiled-code size from NIMBLE).

if (!requireNamespace("devtools", quietly = TRUE))
  stop("Install 'devtools' to run the CRAN-check harness.")

devtools::document()
res <- devtools::check(cran = TRUE, error_on = "warning")

cat("\n--- check summary ---\n")
print(res)
cat("\nReminder: do NOT claim CRAN-readiness until Layers 1-6 pass per the\n",
    "gating policy in project knowledge Section 8.2.\n", sep = "")
