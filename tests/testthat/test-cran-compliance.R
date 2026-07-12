# Compliance tests for the v0.1.0 scope guards. These are fast (no NIMBLE
# compilation) and assert that out-of-roadmap requests fail with clear errors
# rather than silently doing something (project knowledge 0.3.1 / 0.3.7).

test_that("nimixReg scope guards: covariate gating, count response", {
  df <- data.frame(y = rnorm(20), x = rnorm(20))
  expect_error(nimixReg(y ~ x, df, gating = "covariate"), "9.8|gating")
  # Poisson regression is supported (v0.4.x) but a continuous response must be
  # rejected early with a clear message, not crash inside the sampler.
  expect_error(nimixReg(y ~ x, df, distribution = "poisson"), "count|integer")
  # a genuinely unknown family still points the user at the available ones
  expect_error(nimixReg(y ~ x, df, distribution = "gamma"), "not available")
})

test_that("multivariate data routes to NormalMvSpec, univariate to NormalUvSpec", {
  expect_s4_class(nimix:::.selectClusterSpec("normal", TRUE,  2L), "NormalMvSpec")
  expect_s4_class(nimix:::.selectClusterSpec("normal", FALSE, 1L), "NormalUvSpec")
  # forcing the wrong family for the data shape errors clearly
  expect_error(nimix:::.selectClusterSpec("normal-mv", FALSE, 1L), "matrix")
  expect_error(nimix:::.selectClusterSpec("normal-uv", TRUE,  2L), "vector")
  # Student-t (v0.4.0) now routes to the univariate / multivariate spec.
  expect_s4_class(nimix:::.selectClusterSpec("student-t", FALSE, 1L), "StudentTUvSpec")
  expect_s4_class(nimix:::.selectClusterSpec("studentt",  TRUE,  2L), "StudentTMvSpec")
  # a family that is genuinely not provided still errors clearly
  expect_error(nimix:::.selectClusterSpec("gamma", FALSE, 1L), "not available")
})

test_that("unsupported distributions error clearly", {
  expect_error(nimixClust(rnorm(20), distribution = "gamma"), "not available")
})

test_that("input validation: NAs, tiny n, K_max bounds", {
  expect_error(nimixClust(c(1, NA, 3)))
  expect_error(nimixClust(1))                       # n < 2
  expect_error(nimixClust(rnorm(20), K_max = 1))    # K_max < 2
  expect_error(nimixClust(rnorm(5), K_max = 50))    # K_max > n
})

test_that("DPMEngine validates its concentration prior", {
  expect_error(DPMEngine(concPrior = c(-1, 2)))
  expect_error(DPMEngine(concPrior = 2))
  expect_s4_class(DPMEngine(), "DPMEngine")
})
