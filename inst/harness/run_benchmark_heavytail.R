## run_benchmark_heavytail.R
## Effective-sample-size-per-second comparison of the two heavy-tail routes:
## direct Student-t vs the Normal-Gamma scale mixture (slice-sampled omega).
## ESS is measured on the number of occupied clusters K (label-invariant).
## Run with: Rscript inst/harness/run_benchmark_heavytail.R

suppressMessages({library(nimix); library(coda)})
essK <- function(fit) as.numeric(effectiveSize(as.numeric(fit@Kposterior)))
bench <- function(label, expr) {
  t0 <- proc.time()[3]; fit <- force(expr); el <- proc.time()[3] - t0
  e <- essK(fit)
  cat(sprintf("%-22s %6.1fs  ESS(K)=%5.0f  ESS(K)/s=%6.2f  modalK=%d\n",
              label, el, e, e / el,
              as.integer(names(sort(table(fit@Kposterior), decreasing = TRUE))[1])))
  invisible(fit)
}
ctrl <- list(niter = 6000, nburnin = 2000)

cat("== Univariate clustering (two heavy-tailed clusters + outliers) ==\n")
set.seed(7); y <- c(rt(80, 4) - 5, rt(80, 4) + 5, 18, -20)
bench("studentt (direct)",
      nimixClust(y, K_max = 8, distribution = "studentt", prior = list(df = 4),
                 mcmcControl = ctrl, seed = 7, verbose = FALSE))
bench("normalgamma (slice)",
      nimixClust(y, K_max = 8, distribution = "normalgamma", prior = list(df = 4),
                 mcmcControl = ctrl, seed = 7, verbose = FALSE))

cat("\n== Univariate regression (two regimes, heavy-tailed errors) ==\n")
set.seed(1); x <- rnorm(200)
dr <- data.frame(x = x, y = c(2 * x[1:100], -2 * x[101:200]) +
                            c(rt(100, 4), rt(100, 4)) * 0.8)
bench("studentt-reg (direct)",
      nimixReg(y ~ x, dr, K_max = 8, distribution = "studentt", prior = list(df = 4),
               mcmcControl = ctrl, seed = 1, verbose = FALSE))
bench("normalgamma-reg (slice)",
      nimixReg(y ~ x, dr, K_max = 8, distribution = "normalgamma", prior = list(df = 4),
               mcmcControl = ctrl, seed = 1, verbose = FALSE))
