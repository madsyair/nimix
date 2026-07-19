# Sampler dispatch for regression families.
#
# StudentTRegSpec contains NormalRegSpec, so S4 dispatch used to hand it
# NormalRegSpec's customizeSamplers, which installs the exact
# Normal-Inverse-Gamma Gibbs step on (betaTilde, s2Tilde). That conditional is
# exact ONLY under a Gaussian likelihood; under a t likelihood it is a Gibbs
# step drawing from the wrong conditional with no accept/reject to correct it,
# so the chain targets the wrong stationary distribution.
#
# It hid well: the slopes were barely affected (symmetric errors), and only
# the scale moved -- s2 biased ~17% at df = 4 against a correct RW_block
# reference, shrinking to ~1% at df = 30 as t approaches Normal. That df
# gradient is what pinned the mechanism on the likelihood mismatch.

test_that("student-t regression does not inherit the Gaussian NIG sampler", {
  spec <- getDistribution("student-t-reg")
  # the guard is the dispatch itself: an own method, not NormalRegSpec's
  m <- selectMethod("customizeSamplers", class(spec))
  expect_identical(as.character(m@defined[[1]]), "StudentTRegSpec")

  # and the families that legitimately use their own samplers keep them
  expect_identical(
    as.character(selectMethod("customizeSamplers", "NormalRegSpec")@defined[[1]]),
    "NormalRegSpec")
  expect_identical(
    as.character(selectMethod("customizeSamplers",
                              "NormalGammaRegSpec")@defined[[1]]),
    "NormalGammaRegSpec")
})

test_that("student-t regression installs no conjugate sampler on betaTilde", {
  skip_on_cran()
  set.seed(5)
  n <- 120L
  x <- runif(n, -2, 2)
  zc <- rep(1:2, length.out = n)
  y <- c(3, -3)[zc] + c(2, -2)[zc] * x + rt(n, 4) * 0.5
  X <- stats::model.matrix(~ x)
  spec <- getDistribution("student-t-reg")
  pr <- nimix:::defaultPrior(spec, y, control = list(X = X))
  mc <- nimix:::buildModelCode(spec, FixedKEngine(), n = n, L = 2, d = 1)
  cn <- c(nimix:::buildConstants(spec, pr, n), list(K = 2, alphaVec = rep(1, 2)))
  m <- nimble::nimbleModel(
    mc$code, constants = cn, data = list(y = y),
    inits = list(z = zc, weights = c(.5, .5),
                 betaTilde = rbind(c(3, 2), c(-3, -2)), s2Tilde = c(1, 1)),
    calculate = TRUE)
  conf <- nimble::configureMCMC(m, monitors = c("betaTilde", "s2Tilde"),
                                print = FALSE)
  nimix:::customizeSamplers(spec, conf, m)
  nm <- vapply(conf$getSamplers(), function(u) u$name, character(1))
  tg <- vapply(conf$getSamplers(), function(u) u$target[1], character(1))
  beta1 <- nm[grep("betaTilde\\[1", tg)]
  expect_true(length(beta1) > 0L)
  expect_false(any(grepl("ConjSampler", beta1)))
})
