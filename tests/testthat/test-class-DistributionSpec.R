test_that("DistributionSpec is virtual and cannot be instantiated", {
  expect_error(new("DistributionSpec"))
})

test_that("NormalUvSpec constructs with expected slots", {
  spec <- NormalUvSpec()
  expect_s4_class(spec, "NormalUvSpec")
  expect_s4_class(spec, "DistributionSpec")
  expect_identical(spec@name, "normal-uv")
  expect_identical(spec@paramNames, c("mu", "s2"))
  expect_identical(spec@dataDim, 1L)
})

test_that("registry registers, retrieves, and lists built-ins", {
  expect_true(all(c("normal", "normal-uv") %in% listDistributions()))
  expect_s4_class(getDistribution("normal"), "NormalUvSpec")
  expect_error(getDistribution("does-not-exist"))
})

test_that("registerDistribution guards duplicates and types", {
  expect_error(registerDistribution(42))
  expect_error(registerDistribution(NormalUvSpec()))           # already there
  expect_silent(registerDistribution(NormalUvSpec(), overwrite = TRUE))
})
