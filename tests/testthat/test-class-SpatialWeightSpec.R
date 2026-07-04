# Tests for SpatialWeightSpec (v1.1.0): the neighbourhood-structure class that
# the MRF engine (planned v1.2.0) will pair with any DistributionSpec.

test_that("validity rejects malformed adjacency structures", {
  ok <- matrix(0, 3, 3); ok[1, 2] <- ok[2, 1] <- 1
  expect_s4_class(spatialWeights(ok), "SpatialWeightSpec")

  asym <- ok; asym[1, 3] <- 1                       # not symmetric
  expect_error(spatialWeights(asym), "symmetric")

  selfd <- ok; diag(selfd) <- 1                     # self-neighbours
  expect_error(spatialWeights(selfd), "diagonal")

  neg <- ok; neg[2, 3] <- neg[3, 2] <- -1           # negative weight
  expect_error(spatialWeights(neg), "non-negative")

  expect_error(spatialWeights(matrix(0, 2, 3)), "square")
  expect_error(spatialWeights(ok, regionIds = c("a", "a", "b")), "unique")
  expect_error(spatialWeights(ok, regionIds = c("a", "b")), "one entry per")
})

test_that("grid contiguity has the textbook neighbour counts", {
  rook <- gridAdjacency(3, 3, "rook")
  expect_equal(nRegions(rook), 9L)
  expect_length(neighborsOf(rook, "r1c1"), 2L)      # corner, rook
  expect_length(neighborsOf(rook, "r1c2"), 3L)      # edge, rook
  expect_length(neighborsOf(rook, "r2c2"), 4L)      # interior, rook

  queen <- gridAdjacency(3, 3, "queen")
  expect_length(neighborsOf(queen, "r1c1"), 3L)     # corner, queen
  expect_length(neighborsOf(queen, "r1c2"), 5L)     # edge, queen
  expect_length(neighborsOf(queen, "r2c2"), 8L)     # interior, queen

  # symmetry of the generated graph and named accessors
  expect_true(isSymmetric(unname(getAdjacency(queen))))
  expect_true("r2c2" %in% neighborsOf(queen, "r1c1"))
})

test_that("neighborsOf accepts index or id and rejects unknown regions", {
  g <- gridAdjacency(2, 2, "rook")
  expect_identical(neighborsOf(g, 1L), neighborsOf(g, "r1c1"))
  expect_error(neighborsOf(g, "nowhere"), "Unknown region")
  expect_error(neighborsOf(g, 99L), "Unknown region")
})

test_that("nimixClust validates spatialWeights and ties it to method = 'mrf'", {
  y <- rnorm(20)
  expect_error(nimixClust(y, spatialWeights = "not-a-spec"),
               "SpatialWeightSpec")
  g <- gridAdjacency(4, 5)
  expect_error(nimixClust(y, spatialWeights = g), "only used by method = 'mrf'")
})
