# HMMEngine: regime-switching mixture with the state path marginalised by a
# forward-algorithm likelihood and recovered post-hoc by FFBS. Design and
# measured motivation are in the F2 gate report; the numbers asserted here
# (recovery, decoding accuracy 1.0 on well-separated regimes) reproduce it.

.simRegime <- function(T = 300L, seed = 11L) {
  set.seed(seed)
  P <- rbind(c(.95, .05), c(.10, .90)); mu <- c(-2, 2); sg <- c(.6, .8)
  z <- integer(T); z[1] <- 1L
  for (t in 2:T) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  list(y = rnorm(T, mu[z], sg[z]), z = z, P = P, mu = mu)
}

test_that("hmm engine recovers states, parameters, and transitions", {
  skip_on_cran()
  d <- .simRegime()
  f <- nimixClust(d$y, K = 2, method = "hmm",
                  mcmcControl = list(niter = 3000, nburnin = 1200), seed = 1)
  expect_identical(f@engineUsed, "hmm")
  # FFBS allocation draws exist with the right shape and range
  expect_identical(dim(f@clusterAllocation)[2], length(d$y))
  expect_true(all(f@clusterAllocation %in% 1:2))
  # location recovery through the standard relabel path
  fr <- relabel(f)
  s <- fr@relabeled$summary
  expect_lt(max(abs(sort(s$mu_mean) - c(-2, 2))), 0.5)
  # decoded path: MAP-per-time from FFBS and Viterbi both match the truth
  zmap <- apply(f@clusterAllocation, 2L,
                function(v) as.integer(names(which.max(table(v)))))
  accM <- max(mean(zmap == d$z), mean((3L - zmap) == d$z))
  expect_gt(accM, 0.95)
  zv <- viterbiPath(f)
  accV <- max(mean(zv == d$z), mean((3L - zv) == d$z))
  expect_gt(accV, 0.95)
  # self-transitions recovered (order-free: compare sorted diagonals)
  pd <- sort(c(mean(f@mcmcSamples[, "P[1, 1]"]),
               mean(f@mcmcSamples[, "P[2, 2]"])))
  expect_lt(max(abs(pd - c(0.90, 0.95))), 0.08)
})

test_that("label-free partition tools work unchanged on hmm fits", {
  skip_on_cran()
  d <- .simRegime(T = 200L)
  f <- nimixClust(d$y, K = 2, method = "hmm",
                  mcmcControl = list(niter = 2000, nburnin = 800), seed = 1)
  S <- psm(f)
  expect_identical(dim(S), c(200L, 200L))
  bp <- binderPartition(f, S)
  expect_identical(bp$K, 2L)
})

test_that("over-parameterised nStates leaves empty states without corruption", {
  skip_on_cran()
  # K = 4 states against 2 true regimes: the 9.16 lesson made structural --
  # empty states MUST occur and everything must stay consistent.
  d <- .simRegime()
  f <- nimixClust(d$y, K = 4, method = "hmm",
                  mcmcControl = list(niter = 2500, nburnin = 1000), seed = 1)
  expect_true(any(f@Kposterior < 4L))          # scenario is non-vacuous
  expect_true(all(f@clusterAllocation %in% 1:4))
  zmap <- apply(f@clusterAllocation, 2L,
                function(v) as.integer(names(which.max(table(v)))))
  agree <- function(a, b) mean(outer(a, a, "==") == outer(b, b, "=="))
  expect_gt(agree(zmap, d$z), 0.95)            # partition structure survives
})

test_that("hmm guards: multivariate, non-normal, and K_max are refused", {
  skip_on_cran()
  y <- rnorm(50)
  expect_error(nimixClust(matrix(rnorm(60), ncol = 2), K = 2, method = "hmm"),
               "univariate|normal")
  # Emission families join the engine incrementally, so this guard must
  # track the rollout (the 9.19 stale-test class -- caught twice already):
  # use a family still OUTSIDE the current set {normal, student-t, poisson,
  # msnburr, msnburr2a, gmsnburr}.
  expect_error(nimixClust(y, K = 2, method = "hmm", distribution = "fssn"),
               "gated plan")
  expect_error(nimixClust(y, K_max = 5, method = "hmm"),
               "needs the number of components K")
  expect_error(viterbiPath(structure(list(), class = "logical")))
})

test_that("forward kernel is exact against a pure-R reference", {
  skip_on_cran()
  nimix:::.nimixEnsureHMM()
  fwdR <- function(x, mu, sg, P, init) {
    a <- init * dnorm(x[1], mu, sg); ll <- log(sum(a)); a <- a / sum(a)
    for (t in 2:length(x)) {
      a <- as.vector(t(P) %*% a) * dnorm(x[t], mu, sg)
      ll <- ll + log(sum(a)); a <- a / sum(a)
    }
    ll
  }
  set.seed(3)
  y <- rnorm(100)
  P <- rbind(c(.9, .1), c(.2, .8))
  llR <- fwdR(y, c(-1, 1), c(1, 0.5), P, c(.5, .5))
  llK <- get("dRegimeHMMNorm_k", envir = globalenv())(
    y, c(-1, 1), c(1, 0.25), P, c(.5, .5), log = 1)
  expect_equal(llK, llR, tolerance = 1e-12)
})

test_that("student-t emissions recover heavy-tailed regimes", {
  skip_on_cran()
  set.seed(21)
  P <- rbind(c(.93, .07), c(.08, .92)); mu <- c(-2, 2)
  z <- integer(350); z[1] <- 1L
  for (t in 2:350) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  y <- mu[z] + 0.6 * rt(350, df = 4)
  f <- nimixClust(y, K = 2, method = "hmm", distribution = "student-t",
                  mcmcControl = list(niter = 3000, nburnin = 1200), seed = 1)
  fr <- relabel(f)
  s <- fr@relabeled$summary
  expect_lt(max(abs(sort(s$mu_mean) - c(-2, 2))), 0.5)
  zv <- viterbiPath(f)
  expect_gt(max(mean(zv == z), mean((3L - zv) == z)), 0.95)
  pd <- sort(c(mean(f@mcmcSamples[, "P[1, 1]"]),
               mean(f@mcmcSamples[, "P[2, 2]"])))
  expect_lt(max(abs(pd - c(0.92, 0.93))), 0.08)
})

test_that("msnburr emissions recover skewed regimes", {
  skip_on_cran()
  set.seed(41)
  P <- rbind(c(.94, .06), c(.09, .91))
  z <- integer(320); z[1] <- 1L
  for (t in 2:320) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  y <- ifelse(z == 1, rmsnburr(320, -3, .8, .6), rmsnburr(320, 3, .8, 2.2))
  f <- nimixClust(y, K = 2, method = "hmm", distribution = "msnburr",
                  mcmcControl = list(niter = 3000, nburnin = 1200), seed = 1)
  s <- relabel(f)@relabeled$summary
  expect_lt(max(abs(sort(s$mu_mean) - c(-3, 3))), 0.6)
  zv <- viterbiPath(f)
  expect_gt(max(mean(zv == z), mean((3L - zv) == z)), 0.95)
})

test_that("poisson emissions recover count regimes", {
  skip_on_cran()
  set.seed(43)
  P <- rbind(c(.95, .05), c(.10, .90)); lam <- c(3, 15)
  z <- integer(300); z[1] <- 1L
  for (t in 2:300) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  y <- rpois(300, lam[z])
  f <- nimixClust(y, K = 2, method = "hmm", distribution = "poisson",
                  mcmcControl = list(niter = 2500, nburnin = 1000), seed = 1)
  lh <- sort(colMeans(nimix:::.nodeToArray(f@mcmcSamples, "lambda", 2)))
  expect_lt(max(abs(lh - c(3, 15))), 1.5)
  zv <- viterbiPath(f)
  expect_gt(max(mean(zv == z), mean((3L - zv) == z)), 0.95)
  expect_identical(binderPartition(f)$K, 2L)
})

test_that("msnburr2a emissions recover skewed regimes", {
  skip_on_cran()
  set.seed(61)
  P <- rbind(c(.92, .08), c(.09, .91)); mu <- c(-3, 3); al <- c(0.6, 2.2)
  z <- integer(300); z[1] <- 1L
  for (t in 2:300) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  y <- numeric(300)
  for (t in 1:300) y[t] <- rmsnburr2a(1, mu[z[t]], 0.8, al[z[t]])
  f <- nimixClust(y, K = 2, method = "hmm", distribution = "msnburr2a",
                  mcmcControl = list(niter = 2500, nburnin = 1000), seed = 1)
  s <- relabel(f)@relabeled$summary
  expect_lt(max(abs(sort(s$mu_mean) - c(-3, 3))), 0.6)
  zv <- viterbiPath(f)
  expect_gt(max(mean(zv == z), mean((3L - zv) == z)), 0.95)
})

test_that("gmsnburr (four-parameter) emissions recover skewed regimes", {
  skip_on_cran()
  set.seed(71)
  P <- rbind(c(.93, .07), c(.08, .92)); mu <- c(-4, 4)
  al <- c(0.6, 2.0); th <- c(1.5, 1.0)
  z <- integer(320); z[1] <- 1L
  for (t in 2:320) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  y <- numeric(320)
  for (t in 1:320) y[t] <- rgmsnburr(1, mu[z[t]], 1, al[z[t]], th[z[t]])
  f <- nimixClust(y, K = 2, method = "hmm", distribution = "gmsnburr",
                  mcmcControl = list(niter = 2500, nburnin = 1000), seed = 1)
  s <- relabel(f)@relabeled$summary
  expect_lt(max(abs(sort(s$mu_mean) - c(-4, 4))), 0.6)
  zv <- viterbiPath(f)
  expect_gt(max(mean(zv == z), mean((3L - zv) == z)), 0.95)
})

test_that("FFBS is numerically stable under emission underflow", {
  # A thin-tailed emission can underflow to 0 across all states for an
  # outlying point at some draw; the safeNorm guard must keep FFBS from
  # producing NaN weights (which crashed sample.int before). Construct a
  # draw matrix by hand and decode.
  skip_on_cran()
  set.seed(5)
  P <- rbind(c(.9, .1), c(.1, .9)); mu <- c(-3, 3)
  z <- integer(120); z[1] <- 1L
  for (t in 2:120) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  y <- rnorm(120, mu[z], 0.4)
  y[60] <- 500          # wild outlier: underflows under every state
  f <- nimixClust(y, K = 2, method = "hmm", distribution = "normal",
                  mcmcControl = list(niter = 1500, nburnin = 600), seed = 1)
  expect_true(all(f@clusterAllocation %in% 1:2))   # no NaN crash
  expect_false(anyNA(f@clusterAllocation))
})
