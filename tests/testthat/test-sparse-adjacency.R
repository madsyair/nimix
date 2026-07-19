# F7: the canonical representation of SpatialWeightSpec is the edge list.
#
# The dense matrix was what capped the usable problem size, and it was the
# FIRST binding constraint: a 10 000-node space-time graph was OOM-killed
# during graph construction, before any model existed, while the quantities
# the engine actually needs came to ~0.6 MB (~2700x smaller than the dense
# path's transients). After the refactor the same graph builds in a quarter
# of a second within ~40 MB, a 50 000-node graph in ~1.3 s, and the fit
# ceiling in this 4 GB container rose from ~5000 to ~7000 nodes -- the wall
# now sits in NIMBLE's per-node model memory, which is the next problem, not
# this one.

test_that("the sparse path reproduces the dense path's constants exactly", {
  # Gate 1. The engine derives e1/e2/deg/nbrs from getEdges(); with the old
  # code they came from which(upper.tri(A) & A > 0) and rowSums/which on the
  # dense matrix. Bit-identical constants => identical fits under a seed.
  buildDense <- function(nr, nc, cont) {
    n <- nr * nc; A <- matrix(0, n, n)
    cell <- function(i, j) (i - 1L) * nc + j
    offs <- if (cont == "rook")
      list(c(-1, 0), c(1, 0), c(0, -1), c(0, 1))
    else list(c(-1, 0), c(1, 0), c(0, -1), c(0, 1),
              c(-1, -1), c(-1, 1), c(1, -1), c(1, 1))
    for (i in seq_len(nr)) for (j in seq_len(nc)) for (o in offs) {
      ii <- i + o[1]; jj <- j + o[2]
      if (ii >= 1 && ii <= nr && jj >= 1 && jj <= nc)
        A[cell(i, j), cell(ii, jj)] <- 1
    }
    A
  }
  for (cont in c("rook", "queen")) {
    A <- buildDense(5, 5, cont)
    ut <- which(upper.tri(A) & A > 0, arr.ind = TRUE)
    sw <- gridAdjacency(5, 5, cont)
    E <- getEdges(sw)
    expect_identical(as.integer(ut[, 1]), E[, 1L])
    expect_identical(as.integer(ut[, 2]), E[, 2L])
    # neighbour rows ascending, as which(A[i, ] > 0) gave them
    deg <- tabulate(c(E[, 1L], E[, 2L]), nbins = 25L)
    expect_identical(deg, as.integer(rowSums(A > 0)))
    # dense round-trip for small graphs is unchanged
    expect_identical(unname(getAdjacency(sw)), A)
  }
})

test_that("spacetimeAdjacency matches its dense definition", {
  g <- gridAdjacency(5, 5)
  As <- unname(getAdjacency(g))
  manual <- matrix(0, 100, 100)
  for (t in 1:4) {
    o <- (t - 1) * 25
    manual[(o + 1):(o + 25), (o + 1):(o + 25)] <- As
  }
  for (t in 1:3) for (i in 1:25) {
    a <- (t - 1) * 25 + i; b <- t * 25 + i
    manual[a, b] <- 1; manual[b, a] <- 1
  }
  expect_identical(unname(getAdjacency(spacetimeAdjacency(g, 4))), manual)
})

test_that("the edge constructor validates and normalises", {
  # unordered + duplicate edges are normalised to canonical unique i < j
  sw <- spatialWeights(edges = rbind(c(2L, 1L), c(1L, 2L), c(3L, 1L)),
                       nNodes = 3)
  expect_identical(getEdges(sw), rbind(c(1L, 2L), c(1L, 3L)))
  expect_identical(nRegions(sw), 3L)
  expect_identical(neighborsOf(sw, 1L), c("region2", "region3"))

  expect_error(spatialWeights(edges = rbind(c(1L, 1L)), nNodes = 3),
               "no self-neighbours")
  expect_error(spatialWeights(edges = rbind(c(1L, 4L)), nNodes = 3),
               "1..nNodes")
  expect_error(spatialWeights(edges = rbind(c(1L, 2L))), "nNodes")
  expect_error(spatialWeights(), "exactly one")
  expect_error(spatialWeights(matrix(0, 2, 2), edges = rbind(c(1L, 2L)),
                              nNodes = 2), "exactly one")
})

test_that("getAdjacency refuses to materialise large graphs, pointing at getEdges", {
  # Gate 4. Silently allocating hundreds of MB is how the old wall was hit;
  # the refusal must tell the user where the sparse form lives.
  big <- spatialWeights(edges = cbind(seq_len(5999L), seq_len(5999L) + 1L),
                        nNodes = 6000)
  expect_error(getAdjacency(big), "getEdges")
  expect_identical(nrow(getEdges(big)), 5999L)   # the sparse form still works
})

test_that("formerly impossible graph sizes now build", {
  skip_on_cran()
  # Gate 2. This exact construction was OOM-killed before the refactor.
  t0 <- system.time(ST <- spacetimeAdjacency(gridAdjacency(10, 10), 100))
  expect_lt(t0[["elapsed"]], 10)
  expect_identical(nRegions(ST), 10000L)
  expect_identical(nrow(getEdges(ST)), 100L * 180L + 100L * 99L)
  # and an order of magnitude beyond
  ST2 <- spacetimeAdjacency(gridAdjacency(20, 25), 100)
  expect_identical(nRegions(ST2), 50000L)
})
