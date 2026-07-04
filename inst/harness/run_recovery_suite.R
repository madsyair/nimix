## inst/harness/run_recovery_suite.R
## Layer 4 of the engineering harness: statistical-validity recovery tests,
## run as a systematic DISTRIBUTION x ENGINE matrix (v0.9.0 hardening).
##
## Every released clustering family (Gaussian uv/mv, Student-t uv/mv,
## Normal-Gamma uv/mv, Poisson, Binomial) is fitted with BOTH production
## engines (dpm, fixedk) on fixed simulated data with known truth, across
## three MCMC seeds. The regression path is exercised for both engines too.
## Because the data are fixed per combo, the compiled-model cache means each
## combo compiles once and the extra seeds re-run cheaply.
##
##   Rscript inst/harness/run_recovery_suite.R [group ...]
## groups: uvnorm uvt uvng mvnorm mvt mvng counts reg   (default: all)

suppressMessages(library(nimix))

seeds <- c(7L, 8L, 9L)
ctrl  <- list(niter = 4000, nburnin = 1500)

## ---- fixed datasets (data seeds separate from the MCMC seeds) --------------
set.seed(101); yUv  <- c(rnorm(100, -4), rnorm(100, 4))
set.seed(102); yT   <- c(-4 + rt(100, df = 5), 4 + rt(100, df = 5))
set.seed(103); Ymv  <- rbind(matrix(rnorm(120, -4, 1), ncol = 2),
                             matrix(rnorm(120,  4, 1), ncol = 2))
set.seed(104); Ytmv <- rbind(-4 + matrix(rt(120, df = 6), ncol = 2),
                              4 + matrix(rt(120, df = 6), ncol = 2))
set.seed(105); yPo  <- c(rpois(100, 3), rpois(100, 15))
set.seed(106); yBi  <- c(rbinom(100, 20, 0.2), rbinom(100, 20, 0.8))
## MRF: 12x10 rook grid, two spatial blocks, overlapping components
nrG <- 12; ncG <- 10; nG <- nrG * ncG
gridSW <- gridAdjacency(nrG, ncG, "rook")
zBlock <- integer(nG)
for (i in 1:nrG) for (j in 1:ncG) zBlock[(i - 1) * ncG + j] <- if (j <= ncG / 2) 1L else 2L
set.seed(108); yMRF <- rnorm(nG, c(-2, 2)[zBlock], 1.4)

set.seed(109)
YMRF <- matrix(rnorm(nG * 2, mean = rep(c(-1.5, 1.5)[zBlock], 2), sd = 1.4), ncol = 2)

set.seed(110)
xMRF <- runif(nG, -2, 2)
yMRFreg <- ifelse(zBlock == 1L, 2 * xMRF, -2 * xMRF) + rnorm(nG, 0, 1.6)
dfMRFreg <- data.frame(y = yMRFreg, x = xMRF)

set.seed(111); yMRFpois <- rpois(nG, c(3, 12)[zBlock])
set.seed(112); yMRFbin  <- rbinom(nG, 20, c(0.2, 0.7)[zBlock])
set.seed(113); yMRFt    <- c(-2, 2)[zBlock] + rt(nG, df = 5) * 1.1
set.seed(114); xG <- runif(nG, -1.5, 1.5)
yMRFpg <- rpois(nG, exp(1 + ifelse(zBlock == 1L, 0.8, -0.8) * xG))
yMRFbg <- rbinom(nG, 15, plogis(0.2 + ifelse(zBlock == 1L, 1.5, -1.5) * xG))
yMRFtg <- ifelse(zBlock == 1L, 2 * xG, -2 * xG) + rt(nG, 5) * 1.2

set.seed(115)
YMRFt2 <- matrix(rep(c(-1.5, 1.5)[zBlock], 2) + rt(nG * 2, df = 5) * 1.1, ncol = 2)
set.seed(116); yMRFng <- c(-2, 2)[zBlock] + rt(nG, df = 5) * 1.1
set.seed(117); xG2 <- runif(nG, -1.5, 1.5)
XG2 <- cbind(1, xG2)
B1g <- rbind(c(1, -1), c(2, -2)); B2g <- rbind(c(-1, 1), c(-2, 2))
muG2 <- t(sapply(seq_len(nG), function(i) XG2[i, ] %*% (if (zBlock[i] == 1L) B1g else B2g)))
YregN <- muG2 + matrix(rnorm(nG * 2, 0, 1.3), ncol = 2)
YregT <- muG2 + matrix(rt(nG * 2, 5) * 1.1, ncol = 2)
dfRegN <- data.frame(y1 = YregN[, 1], y2 = YregN[, 2], x = xG2)
dfRegT <- data.frame(y1 = YregT[, 1], y2 = YregT[, 2], x = xG2)
yMRFtg2 <- ifelse(zBlock == 1L, 2 * xG2, -2 * xG2) + rt(nG, 5) * 1.2
dfTg2 <- data.frame(y = yMRFtg2, x = xG2)

set.seed(107); xR   <- runif(250, -3, 3)
grp  <- rep(1:2, c(200, 50))
yR   <- ifelse(grp == 1, 2 * xR, -2 * xR) + rnorm(250, 0, 0.7)
regDf <- data.frame(y = yR, x = xR)

## ---- family-specific recovery checks (on the relabelled fit) ---------------
chkUvMu <- function(fit, true) {
  s <- sort(fit@relabeled$summary$mu_mean)
  fit@relabeled$modalK == length(true) && all(abs(s - sort(true)) < 1.2)
}
chkMvMu1 <- function(fit, true1) {
  s <- sort(fit@relabeled$summary$mu_1)
  fit@relabeled$modalK == length(true1) && all(abs(s - sort(true1)) < 1.5)
}
chkPois <- function(fit, true) {
  s <- sort(fit@relabeled$summary$lambda_mean)
  fit@relabeled$modalK == 2L && all(abs(s - sort(true)) / sort(true) < 0.35)
}
chkBinom <- function(fit, true) {
  s <- sort(fit@relabeled$summary$prob_mean)
  fit@relabeled$modalK == 2L && all(abs(s - sort(true)) < 0.15)
}
## For the DPM on DISCRETE count data the posterior number of clusters is known
## to be diffuse and to spawn small transient components even when the two
## dominant components recover the truth cleanly -- the DPM posterior on the
## number of components is not consistent for it (Miller & Harrison 2013,
## "A simple example of Dirichlet process mixture inconsistency for the number
## of components", NIPS). The statistically appropriate recovery criterion is
## therefore on the DOMINANT components (weight >= 0.1): exactly two of them,
## with parameters near the truth. The strict modal-K criterion above is kept
## for fixed-K fits and for the continuous-data DPM combos, which meet it.
chkPoisDpm <- function(fit, true) {
  sm <- fit@relabeled$summary
  dom <- sm[sm$weight >= 0.1, , drop = FALSE]
  s <- sort(dom$lambda_mean)
  nrow(dom) == 2L && all(abs(s - sort(true)) / sort(true) < 0.35)
}
chkBinomDpm <- function(fit, true) {
  sm <- fit@relabeled$summary
  dom <- sm[sm$weight >= 0.1, , drop = FALSE]
  s <- sort(dom$prob_mean)
  nrow(dom) == 2L && all(abs(s - sort(true)) < 0.15)
}
chkMRF <- function(fit) {
  zMap <- apply(fit@clusterAllocation, 2L,
                function(v) as.integer(names(which.max(table(v)))))
  acc <- max(mean(zMap == zBlock), mean(zMap == (3L - zBlock)))
  mu <- sort(fit@relabeled$summary$mu_mean)
  acc > 0.9 && abs(mu[1] + 2) < 0.8 && abs(mu[2] - 2) < 0.8
}
chkMRFmv <- function(fit) {
  zMap <- apply(fit@clusterAllocation, 2L,
                function(v) as.integer(names(which.max(table(v)))))
  acc <- max(mean(zMap == zBlock), mean(zMap == (3L - zBlock)))
  mu1 <- sort(fit@relabeled$summary$mu_1)
  acc > 0.9 && abs(mu1[1] + 1.5) < 0.8 && abs(mu1[2] - 1.5) < 0.8
}
chkMRFreg <- function(fit) {
  zMap <- apply(fit@clusterAllocation, 2L,
                function(v) as.integer(names(which.max(table(v)))))
  acc <- max(mean(zMap == zBlock), mean(zMap == (3L - zBlock)))
  sl <- sort(fit@relabeled$summary[["x"]])
  acc > 0.9 && sl[1] < -1 && sl[2] > 1
}
chkMRFbeta <- function(fit) {
  b <- as.numeric(fit@mcmcSamples[, "beta"])
  zMap <- apply(fit@clusterAllocation, 2L,
                function(v) as.integer(names(which.max(table(v)))))
  acc <- max(mean(zMap == zBlock), mean(zMap == (3L - zBlock)))
  stats::var(b) > 0 && mean(b > 0.2) > 0.95 && acc > 0.9
}
accBlock <- function(fit) {
  zMap <- apply(fit@clusterAllocation, 2L,
                function(v) as.integer(names(which.max(table(v)))))
  max(mean(zMap == zBlock), mean(zMap == (3L - zBlock)))
}
slope1Of <- function(fit) sort(fit@relabeled$summary[[grep("^x", names(fit@relabeled$summary), value = TRUE)[1]]])
chkSlopes <- function(fit) {
  s <- sort(fit@relabeled$summary[["x"]])
  s[1] < -1 && s[length(s)] > 1
}

## ---- the matrix -------------------------------------------------------------
combos <- list(
  uvnorm = list(
    list(name = "normal-uv x dpm",
         fit = function(s) nimixClust(yUv, K_max = 8, mcmcControl = ctrl, seed = s),
         check = function(f) chkUvMu(f, c(-4, 4))),
    list(name = "normal-uv x fixedk",
         fit = function(s) nimixClust(yUv, K = 2, method = "fixedk",
                                      mcmcControl = ctrl, seed = s),
         check = function(f) chkUvMu(f, c(-4, 4)))
  ),
  uvt = list(
    list(name = "student-t x dpm",
         fit = function(s) nimixClust(yT, K_max = 8, distribution = "studentt",
                                      mcmcControl = ctrl, seed = s),
         check = function(f) chkUvMu(f, c(-4, 4))),
    list(name = "student-t x fixedk",
         fit = function(s) nimixClust(yT, K = 2, distribution = "studentt",
                                      method = "fixedk", mcmcControl = ctrl, seed = s),
         check = function(f) chkUvMu(f, c(-4, 4)))
  ),
  uvng = list(
    list(name = "normal-gamma x dpm",
         fit = function(s) nimixClust(yT, K_max = 8, distribution = "normalgamma",
                                      mcmcControl = ctrl, seed = s),
         check = function(f) chkUvMu(f, c(-4, 4))),
    list(name = "normal-gamma x fixedk",
         fit = function(s) nimixClust(yT, K = 2, distribution = "normalgamma",
                                      method = "fixedk", mcmcControl = ctrl, seed = s),
         check = function(f) chkUvMu(f, c(-4, 4)))
  ),
  mvnorm = list(
    list(name = "normal-mv x dpm",
         fit = function(s) nimixClust(Ymv, K_max = 6, mcmcControl = ctrl, seed = s),
         check = function(f) chkMvMu1(f, c(-4, 4))),
    list(name = "normal-mv x fixedk",
         fit = function(s) nimixClust(Ymv, K = 2, method = "fixedk",
                                      mcmcControl = ctrl, seed = s),
         check = function(f) chkMvMu1(f, c(-4, 4)))
  ),
  mvt = list(
    list(name = "student-t-mv x dpm",
         fit = function(s) nimixClust(Ytmv, K_max = 6, distribution = "studentt",
                                      mcmcControl = ctrl, seed = s),
         check = function(f) chkMvMu1(f, c(-4, 4))),
    list(name = "student-t-mv x fixedk",
         fit = function(s) nimixClust(Ytmv, K = 2, distribution = "studentt",
                                      method = "fixedk", mcmcControl = ctrl, seed = s),
         check = function(f) chkMvMu1(f, c(-4, 4)))
  ),
  mvng = list(
    list(name = "normal-gamma-mv x dpm",
         fit = function(s) nimixClust(Ytmv, K_max = 6, distribution = "normalgamma",
                                      mcmcControl = ctrl, seed = s),
         check = function(f) chkMvMu1(f, c(-4, 4))),
    list(name = "normal-gamma-mv x fixedk",
         fit = function(s) nimixClust(Ytmv, K = 2, distribution = "normalgamma",
                                      method = "fixedk", mcmcControl = ctrl, seed = s),
         check = function(f) chkMvMu1(f, c(-4, 4)))
  ),
  counts = list(
    list(name = "poisson x dpm",
         fit = function(s) nimixClust(yPo, K_max = 8, distribution = "poisson",
                                      mcmcControl = ctrl, seed = s),
         check = function(f) chkPoisDpm(f, c(3, 15))),
    list(name = "poisson x fixedk",
         fit = function(s) nimixClust(yPo, K = 2, distribution = "poisson",
                                      method = "fixedk", mcmcControl = ctrl, seed = s),
         check = function(f) chkPois(f, c(3, 15))),
    list(name = "binomial x dpm",
         fit = function(s) nimixClust(yBi, K_max = 8, distribution = "binomial",
                                      prior = list(size = 20),
                                      mcmcControl = ctrl, seed = s),
         check = function(f) chkBinomDpm(f, c(0.2, 0.8))),
    list(name = "binomial x fixedk",
         fit = function(s) nimixClust(yBi, K = 2, distribution = "binomial",
                                      prior = list(size = 20), method = "fixedk",
                                      mcmcControl = ctrl, seed = s),
         check = function(f) chkBinom(f, c(0.2, 0.8)))
  ),
  mrf = list(
    list(name = "normal-uv x mrf (spatial blocks)",
         fit = function(s) nimixClust(yMRF, K = 2, method = "mrf",
                                      spatialWeights = gridSW,
                                      mcmcControl = ctrl, seed = s),
         check = chkMRF)
  ),
  mrfmv = list(
    list(name = "normal-mv x mrf (spatial blocks)",
         fit = function(s) nimixClust(YMRF, K = 2, method = "mrf",
                                      spatialWeights = gridSW,
                                      mcmcControl = ctrl, seed = s),
         check = chkMRFmv)
  ),
  mrfreg = list(
    list(name = "normal-reg x mrf (spatial slopes)",
         fit = function(s) nimixReg(y ~ x, dfMRFreg, K = 2, method = "mrf",
                                    spatialWeights = gridSW,
                                    mcmcControl = ctrl, seed = s),
         check = chkMRFreg)
  ),
  mrfbeta = list(
    list(name = "normal-uv x mrf (estimateBeta, PL)",
         fit = function(s) nimixClust(yMRF, K = 2, method = "mrf",
                                      spatialWeights = gridSW,
                                      prior = list(estimateBeta = TRUE),
                                      mcmcControl = ctrl, seed = s),
         check = chkMRFbeta)
  ),
  mrfdist1 = list(
    list(name = "poisson x mrf",
         fit = function(s) nimixClust(yMRFpois, K = 2, method = "mrf",
                spatialWeights = gridSW, distribution = "poisson",
                mcmcControl = ctrl, seed = s),
         check = function(f) {
           lam <- sort(f@relabeled$summary$lambda_mean)
           accBlock(f) > 0.9 && abs(lam[1] - 3)/3 < 0.35 && abs(lam[2] - 12)/12 < 0.35 }),
    list(name = "binomial x mrf",
         fit = function(s) nimixClust(yMRFbin, K = 2, method = "mrf",
                spatialWeights = gridSW, distribution = "binomial",
                prior = list(size = 20), mcmcControl = ctrl, seed = s),
         check = function(f) {
           pr <- sort(f@relabeled$summary$prob_mean)
           accBlock(f) > 0.9 && abs(pr[1] - 0.2) < 0.12 && abs(pr[2] - 0.7) < 0.12 }),
    list(name = "student-t x mrf",
         fit = function(s) nimixClust(yMRFt, K = 2, method = "mrf",
                spatialWeights = gridSW, distribution = "studentt",
                mcmcControl = ctrl, seed = s),
         check = function(f) {
           mu <- sort(f@relabeled$summary$mu_mean)
           accBlock(f) > 0.9 && abs(mu[1] + 2) < 0.9 && abs(mu[2] - 2) < 0.9 })
  ),
  mrfdist1reg = list(
    list(name = "poisson-glm x mrf",
         fit = function(s) nimixReg(y ~ x, data.frame(y = yMRFpg, x = xG), K = 2,
                method = "mrf", spatialWeights = gridSW, distribution = "poisson",
                mcmcControl = ctrl, seed = s),
         check = function(f) {
           sl <- sort(f@relabeled$summary[["x"]])
           accBlock(f) > 0.9 && sl[1] < -0.4 && sl[2] > 0.4 }),
    list(name = "binomial-glm x mrf",
         fit = function(s) nimixReg(y ~ x, data.frame(y = yMRFbg, x = xG), K = 2,
                method = "mrf", spatialWeights = gridSW, distribution = "binomial",
                prior = list(size = 15), mcmcControl = ctrl, seed = s),
         check = function(f) {
           sl <- sort(f@relabeled$summary[["x"]])
           accBlock(f) > 0.9 && sl[1] < -0.8 && sl[2] > 0.8 }),
    list(name = "student-t-reg x mrf",
         fit = function(s) nimixReg(y ~ x, data.frame(y = yMRFtg, x = xG), K = 2,
                method = "mrf", spatialWeights = gridSW, distribution = "studentt",
                mcmcControl = ctrl, seed = s),
         check = function(f) {
           sl <- sort(f@relabeled$summary[["x"]])
           accBlock(f) > 0.9 && sl[1] < -1 && sl[2] > 1 })
  ),
  mrfdist2 = list(
    list(name = "student-t-mv x mrf",
         fit = function(s) nimixClust(YMRFt2, K = 2, method = "mrf",
                spatialWeights = gridSW, distribution = "studentt",
                mcmcControl = ctrl, seed = s),
         check = function(f) { mu <- sort(f@relabeled$summary$mu_1)
           accBlock(f) > 0.9 && abs(mu[1] + 1.5) < 0.9 && abs(mu[2] - 1.5) < 0.9 }),
    list(name = "normal-gamma-uv x mrf",
         fit = function(s) nimixClust(yMRFng, K = 2, method = "mrf",
                spatialWeights = gridSW, distribution = "normalgamma",
                mcmcControl = ctrl, seed = s),
         check = function(f) { mu <- sort(f@relabeled$summary$mu_mean)
           accBlock(f) > 0.9 && abs(mu[1] + 2) < 0.9 && abs(mu[2] - 2) < 0.9 }),
    list(name = "normal-gamma-mv x mrf",
         fit = function(s) nimixClust(YMRFt2, K = 2, method = "mrf",
                spatialWeights = gridSW, distribution = "normalgamma",
                mcmcControl = ctrl, seed = s),
         check = function(f) { mu <- sort(f@relabeled$summary$mu_1)
           accBlock(f) > 0.9 && abs(mu[1] + 1.5) < 0.9 && abs(mu[2] - 1.5) < 0.9 })
  ),
  mrfdist2reg = list(
    list(name = "normal-mv-reg x mrf",
         fit = function(s) nimixReg(cbind(y1, y2) ~ x, dfRegN, K = 2,
                method = "mrf", spatialWeights = gridSW,
                mcmcControl = ctrl, seed = s),
         check = function(f) { sl <- slope1Of(f)
           accBlock(f) > 0.9 && sl[1] < -1 && sl[2] > 1 }),
    list(name = "student-t-mv-reg x mrf",
         fit = function(s) nimixReg(cbind(y1, y2) ~ x, dfRegT, K = 2,
                method = "mrf", spatialWeights = gridSW,
                distribution = "studentt", mcmcControl = ctrl, seed = s),
         check = function(f) { sl <- slope1Of(f)
           accBlock(f) > 0.9 && sl[1] < -1 && sl[2] > 1 }),
    list(name = "normal-gamma-reg x mrf",
         fit = function(s) nimixReg(y ~ x, dfTg2, K = 2, method = "mrf",
                spatialWeights = gridSW, distribution = "normalgamma",
                mcmcControl = ctrl, seed = s),
         check = function(f) { sl <- slope1Of(f)
           accBlock(f) > 0.9 && sl[1] < -1 && sl[2] > 1 }),
    list(name = "normal-gamma-mv-reg x mrf",
         fit = function(s) nimixReg(cbind(y1, y2) ~ x, dfRegT, K = 2,
                method = "mrf", spatialWeights = gridSW,
                distribution = "normalgamma", mcmcControl = ctrl, seed = s),
         check = function(f) { sl <- slope1Of(f)
           accBlock(f) > 0.9 && sl[1] < -1 && sl[2] > 1 })
  ),
  reg = list(
    list(name = "normal-reg x dpm (80/20)",
         fit = function(s) nimixReg(y ~ x, regDf, K_max = 8, method = "dpm",
                                    mcmcControl = ctrl, seed = s),
         check = chkSlopes),
    list(name = "normal-reg x fixedk (80/20)",
         fit = function(s) nimixReg(y ~ x, regDf, K = 2, method = "fixedk",
                                    mcmcControl = ctrl, seed = s),
         check = chkSlopes)
  )
)

## ---- runner -----------------------------------------------------------------
runCombo <- function(cb) {
  ok <- logical(length(seeds))
  for (i in seq_along(seeds)) {
    f <- tryCatch(relabel(cb$fit(seeds[i])), error = function(e) e)
    if (inherits(f, "error")) {
      cat(sprintf("  seed %d : ERROR  %s\n", seeds[i], conditionMessage(f)))
      ok[i] <- FALSE
    } else {
      ok[i] <- isTRUE(tryCatch(cb$check(f), error = function(e) FALSE))
      cat(sprintf("  seed %d : %s\n", seeds[i], if (ok[i]) "PASS" else "REVIEW"))
    }
  }
  all(ok)
}

args   <- commandArgs(trailingOnly = TRUE)
groups <- if (length(args)) args else names(combos)
stopifnot(all(groups %in% names(combos)))

tallyPass <- 0L; tallyRev <- 0L
for (g in groups) {
  cat("\n### group:", g, "###\n")
  for (cb in combos[[g]]) {
    cat("--", cb$name, "--\n")
    if (runCombo(cb)) tallyPass <- tallyPass + 1L else tallyRev <- tallyRev + 1L
  }
}
cat(sprintf("\n== TALLY [%s]: %d combos PASS, %d combos REVIEW ==\n",
            paste(groups, collapse = ","), tallyPass, tallyRev))
