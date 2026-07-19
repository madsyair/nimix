# spacetimeAdjacency(): the whole of the MRF-HMM fusion, as it turned out.
#
# The roadmap had this gated as "RED, needs research", on the reasoning that
# undirected spatial dependence breaks the HMM forward algorithm. True, and
# beside the point: the MRF engine never used the forward algorithm -- it
# samples z directly. So a space-time graph is just another adjacency, and the
# existing Potts engine handles it with no changes at all. This function only
# builds the graph.
#
# Measured on a 5x5 grid over 8 times with deliberately overlapping components
# (mu +-1.2, sd 1.0): a plain mixture failed outright (accuracy 0.50, K = 1);
# the spatial graph reached 0.95; the space-time graph reached 0.995.

test_that("spacetimeAdjacency builds the expected graph", {
  W <- gridAdjacency(5, 5)                 # 25 locations, 320 spatial edges
  st <- getAdjacency(spacetimeAdjacency(W, nTime = 8))
  sp <- getAdjacency(spacetimeAdjacency(W, nTime = 8, temporal = FALSE))
  tm <- getAdjacency(spacetimeAdjacency(W, nTime = 8, spatial = FALSE))

  expect_identical(nrow(st), 200L)                  # 25 x 8
  expect_identical(sum(tm) / 2, 25 * 7)             # one chain per location
  expect_identical((sum(st) - sum(sp)) / 2, 25 * 7) # space-time = space + time
  expect_true(all(diag(st) == 0))
  expect_true(isSymmetric(st))

  expect_error(spacetimeAdjacency(W, nTime = 1), "integer >= 2")
  expect_error(spacetimeAdjacency(W, nTime = 4, spatial = FALSE,
                                  temporal = FALSE), "no edges")
})

test_that("temporal edges earn their place on spatio-temporal data", {
  skip_on_cran()
  # Components overlap badly on purpose, so structure has to do the work.
  set.seed(21)
  nr <- 5L; nc <- 5L; nLoc <- nr * nc; nT <- 8L
  gridcol <- rep(seq_len(nc), times = nr)
  ztrue <- matrix(0L, nLoc, nT)
  for (t in seq_len(nT)) ztrue[, t] <- ifelse(gridcol <= 2L, 1L, 2L)
  ztrue[gridcol == 3L, 5:nT] <- 1L        # one column switches regime at t = 5
  z <- as.vector(ztrue)                   # time varies slowest: the required order
  y <- rnorm(nLoc * nT, c(-1.2, 1.2)[z], 1.0)

  W <- gridAdjacency(nr, nc)
  mc <- list(niter = 2500, nburnin = 1000)
  acc <- function(f) {
    zh <- binderPartition(f)$partition
    max(mean(zh == z), mean((3L - zh) == z))
  }

  fSp <- nimixClust(y, K = 2, method = "mrf",
                    spatialWeights = spacetimeAdjacency(W, nT, temporal = FALSE),
                    mcmcControl = mc, seed = 1, verbose = FALSE)
  fSt <- nimixClust(y, K = 2, method = "mrf",
                    spatialWeights = spacetimeAdjacency(W, nT),
                    mcmcControl = mc, seed = 1, verbose = FALSE)

  expect_gt(acc(fSp), 0.85)          # spatial structure already rescues a lot
  expect_gt(acc(fSt), acc(fSp))      # and the temporal edges add on top
  expect_gt(acc(fSt), 0.95)
  expect_identical(binderPartition(fSt)$K, 2L)
})
