# Regression guard for the covariance-reconstruction bug in .ppcSimulate.
#
# When a mixture is fit with more components than the data support, some
# components are empty on some MCMC draws, so max(alloc) < the number of
# monitored components. The covariance trace always holds Kmon * d * d entries;
# reshaping it to c(max(alloc), d, d) truncates and shifts those entries,
# producing a non-symmetric (or, worse, still-PD but wrong) matrix. chol() then
# either errors or -- the dangerous case -- silently returns garbage.
#
# Every mv ppCheck test elsewhere uses K equal to the true cluster count, so no
# component is ever empty and this class of bug slips through. These tests use
# K strictly greater than the truth on purpose.

test_that(".ppcSimulate reconstructs each Sigma symmetric and PD when a component is empty", {
  skip_on_cran()
  ns <- asNamespace("nimix")
  set.seed(4)
  d <- 3
  Y <- rbind(matrix(rnorm(100 * d, -3), 100, d),
             matrix(rnorm(100 * d,  3), 100, d))          # two true clusters
  f <- nimixClust(Y, K = 4, method = "fixedk", distribution = "normal-mv",
                  mcmcControl = list(niter = 500, nburnin = 200), seed = 1)
  S <- f@mcmcSamples
  A <- f@clusterAllocation
  cs <- nimix:::.ppcCols(S, "covTilde")
  cm <- nimix:::.ppcCols(S, "muTilde")
  Kmon <- length(cm) / d

  # there must actually be draws with an empty component, or the test is vacuous
  emptyDraws <- which(apply(A, 1L, max) < Kmon)
  expect_gt(length(emptyDraws), 0L)

  # on every draw with an empty component, each occupied Sigma must be
  # symmetric and positive definite after correct reconstruction
  for (dr in emptyDraws[seq_len(min(30L, length(emptyDraws)))]) {
    covA <- array(S[dr, cs], dim = c(Kmon, d, d))          # correct: Kmon, not max(alloc)
    for (k in seq_len(max(A[dr, ]))) {
      if (!any(A[dr, ] == k)) next
      Sk <- matrix(covA[k, , ], d, d)
      expect_equal(unname(Sk), unname(t(Sk)), tolerance = 1e-10)
      expect_gt(min(eigen(Sk, symmetric = TRUE, only.values = TRUE)$values), 0)
    }
  }
})

test_that("ppCheck runs on over-fitted K for every multivariate family", {
  skip_on_cran()
  set.seed(4)
  d <- 3
  Y <- rbind(matrix(rnorm(100 * d, -3), 100, d),
             matrix(rnorm(100 * d,  3), 100, d))
  for (dist in c("normal-mv", "student-t-mv", "normal-gamma-mv")) {
    f <- nimixClust(Y, K = 4, method = "fixedk", distribution = dist,
                    mcmcControl = list(niter = 400, nburnin = 150), seed = 1)
    pc <- ppCheck(f, nrep = 40)
    expect_s3_class(pc, "data.frame")
    expect_true(all(is.finite(pc$observed)))
  }
})

test_that(".ppcSimulate errors loudly if the covariance trace length is inconsistent", {
  skip_on_cran()
  ns <- asNamespace("nimix")
  set.seed(4)
  d <- 3
  Y <- rbind(matrix(rnorm(60 * d, -3), 60, d),
             matrix(rnorm(60 * d,  3), 60, d))
  f <- nimixClust(Y, K = 3, method = "fixedk", distribution = "normal-mv",
                  mcmcControl = list(niter = 300, nburnin = 100), seed = 1)
  S <- f@mcmcSamples
  # drop one covariance column to simulate an inconsistent trace
  cs <- nimix:::.ppcCols(S, "covTilde")
  Sbad <- S[, -cs[length(cs)], drop = FALSE]
  spec <- f@distSpec
  expect_error(
    nimix:::.ppcSimulate(spec, Sbad, 1L, f@clusterAllocation[1, ], list(d = d)),
    "not Kmon"
  )
})
