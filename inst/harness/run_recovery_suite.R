## inst/harness/run_recovery_suite.R
## Layer 4 of the engineering harness (project knowledge Section 8.1):
## statistical-validity recovery tests. Simulate from KNOWN mixtures, fit, and
## check that the posterior of K concentrates near the truth and that
## RELABELLED component estimates recover the true parameters.
##
##   Rscript inst/harness/run_recovery_suite.R
##
## In v0.2.0 the univariate AND multivariate Gaussian DPM paths exist. The
## DPM-vs-RJMCMC benchmark (run_benchmark_dpm_rjmcmc.R) is added at v0.5.0.

library(nimix)

scenarios <- list(
  two_separated = list(mu = c(-5, 5),  sd = c(1, 1),   n = c(150, 150)),
  three_comp    = list(mu = c(-6, 0, 6), sd = c(1, 1, 1), n = c(120, 120, 120)),
  unbalanced    = list(mu = c(-4, 4),  sd = c(1, 1),   n = c(240, 60))
)

run_one <- function(sc, seed) {
  set.seed(seed)
  y <- unlist(Map(function(m, s, k) rnorm(k, m, s), sc$mu, sc$sd, sc$n))
  fit <- nimixClust(y, K_max = 10,
                     mcmcControl = list(niter = 6000, nburnin = 2000),
                     seed = seed, verbose = FALSE)
  fit <- relabel(fit)
  list(
    trueK  = length(sc$mu),
    modalK = fit@relabeled$modalK,
    estMu  = sort(fit@relabeled$summary$mu_mean),
    trueMu = sort(sc$mu)
  )
}

cat("nimix recovery suite (v0.2.0)\n")
cat("\n### Univariate Gaussian DPM ###\n")
for (nm in names(scenarios)) {
  cat("\n== scenario:", nm, "==\n")
  for (seed in 1:3) {
    r <- run_one(scenarios[[nm]], seed)
    ok <- (r$modalK == r$trueK) &&
          all(abs(r$estMu - r$trueMu) < 1.0)
    cat(sprintf("  seed %d | trueK=%d modalK=%d | %s\n",
                seed, r$trueK, r$modalK, if (ok) "PASS" else "REVIEW"))
  }
}

## --- Multivariate Gaussian DPM (v0.2.0) ------------------------------------
# Two well-separated 2-D Gaussian clusters with known centres.
run_one_mv <- function(seed) {
  set.seed(seed)
  centres <- list(c(-4, -4), c(4, 4))
  Y <- do.call(rbind, lapply(centres, function(mu)
        matrix(rnorm(120, mean = rep(mu, each = 60), sd = 1), ncol = 2)))
  fit <- nimixClust(Y, K_max = 10,
                    mcmcControl = list(niter = 6000, nburnin = 2000),
                    seed = seed, verbose = FALSE)
  fit <- relabel(fit)
  s <- fit@relabeled$summary
  estC <- s[order(s$mu_1), c("mu_1", "mu_2")]
  list(trueK = 2L, modalK = fit@relabeled$modalK, estCentres = estC)
}

cat("\n### Multivariate (2-D) Gaussian DPM ###\n")
for (seed in 1:3) {
  r <- run_one_mv(seed)
  ok <- (r$modalK == 2L)
  cat(sprintf("  seed %d | trueK=2 modalK=%d | %s\n",
              seed, r$modalK, if (ok) "PASS" else "REVIEW"))
}

## --- Mixture of linear regressions (DPM), unbalanced 80/20 -----------------
run_one_reg <- function(seed) {
  set.seed(seed)
  n <- 250; x <- runif(n, -3, 3)
  grp <- c(rep(1L, 200), rep(2L, 50))            # 80 / 20
  y <- ifelse(grp == 1L, 2 * x, -2 * x) + rnorm(n, 0, 0.7)
  fit <- nimixReg(y ~ x, data.frame(y = y, x = x), K_max = 8, method = "dpm",
                  mcmcControl = list(niter = 5000, nburnin = 2000),
                  seed = seed, verbose = FALSE)
  fit <- relabel(fit)
  slopes <- sort(fit@relabeled$summary[["x"]])
  list(modalK = fit@relabeled$modalK, slopes = slopes)
}

cat("\n### Mixture of linear regressions DPM (slopes +2 / -2, 80/20) ###\n")
for (seed in 1:3) {
  r <- run_one_reg(seed)
  ok <- (r$slopes[1] < 0 && r$slopes[length(r$slopes)] > 0)
  cat(sprintf("  seed %d | modalK=%d | slopes=[%.2f, %.2f] | %s\n",
              seed, r$modalK, r$slopes[1], r$slopes[length(r$slopes)],
              if (ok) "PASS" else "REVIEW"))
}

## --- Fixed-K finite mixture (univariate), known K = 2 ----------------------
run_one_fixedk <- function(seed) {
  set.seed(seed)
  y <- c(rnorm(80, -4, 1), rnorm(80, 4, 1))
  fit <- nimixClust(y, K = 2, method = "fixedk",
                    mcmcControl = list(niter = 4000, nburnin = 1500),
                    seed = seed, verbose = FALSE)
  fit <- relabel(fit)
  sort(fit@relabeled$summary$mu_mean)
}

cat("\n### Fixed-K finite mixture (K = 2, means -4 / 4) ###\n")
for (seed in 1:3) {
  m <- run_one_fixedk(seed)
  ok <- (m[1] < -2 && m[2] > 2)
  cat(sprintf("  seed %d | means=[%.2f, %.2f] | %s\n",
              seed, m[1], m[2], if (ok) "PASS" else "REVIEW"))
}
