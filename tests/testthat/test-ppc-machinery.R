# Posterior-predictive machinery: the empty-component reconstruction fix,
# yrep exposure, and chain identity. Split from test-convergence-workflow.R
# so each file stays within the CI time budget (each NIMBLE compile ~35s).

test_that("mv ppCheck survives empty components (K > true clusters)", {
  # Regression for a silent-corruption bug: .ppcSimulate reshaped the
  # covariance trace with dim = c(max(alloc), d, d). On draws where a
  # component was empty, max(alloc) < K_monitor, so the Kmon*d*d-length trace
  # was truncated and shifted -- producing a non-symmetric Sigma (chol error)
  # or, worse, a shifted-but-still-PD Sigma returning wrong numbers silently.
  # The fix derives the array dimension from the monitored count. This test
  # deliberately over-fits K so empty components MUST occur; the old harness
  # only ever used K = true clusters, which is why the bug survived it.
  skip_on_cran()
  set.seed(4)
  Y <- rbind(matrix(rnorm(300, -3), 100, 3), matrix(rnorm(300, 3), 100, 3))
  for (dist in c("normal-mv", "student-t-mv")) {
    f <- nimixClust(Y, K = 3, method = "fixedk", distribution = dist,
                    mcmcControl = list(niter = 500, nburnin = 200), seed = 1)
    # the scenario is only exercised if some draws really have empty components
    kOcc <- apply(f@clusterAllocation, 1L, function(r) length(unique(r)))
    expect_true(any(kOcc < 3L))
    pc <- ppCheck(f, nrep = 50)
    expect_true(all(is.finite(pc$ppp)))
  }
})

test_that("mv ppc reconstruction matches column-name indexing draw by draw", {
  # Stronger than "no error": the reconstructed Sigma must equal the values
  # addressed by explicit covTilde[k, r, c] column names, and be symmetric PD,
  # on every draw -- including those with empty components.
  skip_on_cran()
  set.seed(4)
  Y <- rbind(matrix(rnorm(300, -3), 100, 3), matrix(rnorm(300, 3), 100, 3))
  f <- nimixClust(Y, K = 3, method = "fixedk", distribution = "normal-mv",
                  mcmcControl = list(niter = 400, nburnin = 150), seed = 1)
  S <- f@mcmcSamples; d <- 3L
  cs <- nimix:::.ppcCols(S, "covTilde")
  cm <- nimix:::.ppcCols(S, "muTilde")
  Kmon <- length(cm) / d
  drs <- unique(round(seq(1L, nrow(S), length.out = 40L)))
  for (dr in drs) {
    covA <- array(unname(S[dr, cs]), dim = c(Kmon, d, d))
    for (k in seq_len(Kmon)) {
      Sk <- matrix(covA[k, , ], d, d)
      nmMat <- matrix(0, d, d)
      for (r in seq_len(d)) for (c2 in seq_len(d))
        nmMat[r, c2] <- unname(S[dr, sprintf("covTilde[%d, %d, %d]", k, r, c2)])
      expect_equal(Sk, nmMat)
      expect_gt(min(eigen((Sk + t(Sk)) / 2, only.values = TRUE)$values), 0)
    }
  }
})

test_that("chainId preserves per-chain structure of pooled draws", {
  # Chains used to be rbind-ed with no marker, making post-hoc per-chain
  # R-hat and bayesplot's iter x chain x param array impossible to build.
  skip_on_cran()
  set.seed(2); y <- c(rnorm(50, -3), rnorm(50, 3))
  f <- nimixClust(y, K = 2, method = "fixedk",
                  mcmcControl = list(niter = 400, nburnin = 150, nchains = 3),
                  seed = 1)
  cid <- f@diagnostics$chainId
  expect_equal(length(cid), nrow(f@mcmcSamples))
  expect_equal(sort(unique(cid)), 1:3)
  expect_true(all(table(cid) == nrow(f@mcmcSamples) / 3))
  # single chain still tagged
  f1 <- nimixClust(y, K = 2, method = "fixedk",
                   mcmcControl = list(niter = 300, nburnin = 100), seed = 1)
  expect_equal(unique(f1@diagnostics$chainId), 1L)
})

test_that("ppCheck can retain yrep, and posteriorPredict exposes replicates", {
  skip_on_cran()
  set.seed(2); y <- c(rnorm(50, -3), rnorm(50, 3))
  f <- nimixClust(y, K = 2, method = "fixedk",
                  mcmcControl = list(niter = 400, nburnin = 150), seed = 1)
  # default stays lean
  pc0 <- ppCheck(f, nrep = 20)
  expect_null(attr(pc0, "yrep"))
  # opt-in attaches y, yrep, draws -- the graphical-PPC inputs
  pc <- ppCheck(f, nrep = 20, store_yrep = TRUE)
  yr <- attr(pc, "yrep")
  expect_equal(dim(yr), c(length(attr(pc, "draws")), length(y)))
  expect_identical(attr(pc, "y"), f@data)
  # standalone replicate generator, univariate
  pp <- posteriorPredict(f, ndraws = 15)
  expect_equal(nrow(pp), length(attr(pp, "draws")))
  expect_equal(ncol(pp), length(y))
  # multivariate returns ndraws x n x d -- and exercises empty components
  # (K = 3 on 2 true clusters), which the P0 fix made safe
  set.seed(4)
  Y <- rbind(matrix(rnorm(150, -3), 50, 3), matrix(rnorm(150, 3), 50, 3))
  fm <- nimixClust(Y, K = 3, method = "fixedk", distribution = "normal-mv",
                   mcmcControl = list(niter = 400, nburnin = 150), seed = 1)
  ppm <- posteriorPredict(fm, ndraws = 10)
  expect_equal(length(dim(ppm)), 3L)
  expect_equal(dim(ppm)[2:3], c(100L, 3L))
  # replicates statistically mirror the data (guards against silent shifts)
  expect_lt(max(abs(apply(ppm, 3, mean) - colMeans(Y))), 0.5)
})
