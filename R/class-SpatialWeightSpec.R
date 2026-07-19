## class-SpatialWeightSpec.R ---------------------------------------------------
## Spatial neighbourhood structure for spatially constrained mixtures (v0.6.0).
##
## This class is deliberately ORTHOGONAL to DistributionSpec: it describes the
## graph structure on the OBSERVATIONS/REGIONS (who neighbours whom), not the
## shape of a component density. A future MRF engine (planned for v0.6.0) will
## pair any registered DistributionSpec with a SpatialWeightSpec, so no
## distribution class needs rewriting to gain spatial awareness.
##
## The stored object is a plain dense adjacency/weight matrix: symmetric,
## zero-diagonal, non-negative. Binary contiguity (0/1) is the common case for
## administrative regions (queen/rook contiguity as used throughout spatial
## econometrics; Anselin 1988, Ch. 3); general non-negative weights are
## accepted for weighted graphs. No new package dependency is introduced.

#' Spatial neighbourhood structure for spatially constrained mixtures
#'
#' Represents the neighbourhood graph of \code{n} regions (or observations) as
#' a symmetric, zero-diagonal, non-negative weight matrix. It is the structural
#' ingredient of the spatially constrained mixture models planned for the
#' 1.x series, in which the latent component labels follow a Markov random
#' field on this graph rather than being independent across observations
#' (Besag 1974; spatially variant finite mixtures, Blekas et al. 2005).
#'
#' A \code{SpatialWeightSpec} is intentionally independent of
#' \code{\linkS4class{DistributionSpec}}: the same neighbourhood structure can
#' be paired with any registered component distribution.
#'
#' @slot edges Integer matrix, one row per undirected edge with endpoints
#'   \code{i < j}, in canonical column-major order (see
#'   \code{\link{getEdges}}). This is the canonical representation since
#'   v1.5.0; the dense matrix is derived on demand and refused for large
#'   graphs.
#' @slot edgeWeights Positive numeric weights, one per edge (1 for binary
#'   contiguity).
#' @slot nNodes Integer number of regions.
#' @slot regionIds Character vector of length \code{nNodes} naming the
#'   regions (defaults to row names of the matrix, or \code{"region1"}, ...
#'   when absent).
#'
#' @references
#' Besag, J. (1974). Spatial interaction and the statistical analysis of
#' lattice systems. \emph{Journal of the Royal Statistical Society B}, 36(2),
#' 192--236.
#'
#' Blekas, K., Likas, A., Galatsanos, N.P., & Lagaris, I.E. (2005). A spatially
#' constrained mixture model for image segmentation. \emph{IEEE Transactions on
#' Neural Networks}, 16(2), 494--498. \doi{10.1109/TNN.2004.841773}
#'
#' Anselin, L. (1988). \emph{Spatial Econometrics: Methods and Models}.
#' Kluwer, Dordrecht. (Queen/rook contiguity conventions, Ch. 3.)
#'
#' @seealso \code{\link{spatialWeights}}, \code{\link{gridAdjacency}},
#'   \code{\link{neighborsOf}}
#' @export
# Since F7 the canonical representation is the EDGE LIST, not the dense
# matrix. The dense form is derived on demand (getAdjacency) and refused for
# large graphs: the Potts engine only ever needs the edge list and a
# neighbour table, and the dense intermediate was what capped the usable
# problem size -- measured, a 10 000-node space-time graph died during graph
# CONSTRUCTION (before any model existed) while the quantities actually
# needed came to ~0.6 MB, some 2700x smaller than the dense path's footprint.
#
# Edge order is canonical: sorted by second endpoint then first, i.e. exactly
# the column-major order of which(upper.tri(A) & A > 0) on the old dense
# path. The MRF engine's e1/e2/deg/nbrs are bit-identical to what the dense
# extraction produced, so fits with the same seed are unchanged.
setClass(
  "SpatialWeightSpec",
  representation(edges = "matrix", edgeWeights = "numeric",
                 nNodes = "integer", regionIds = "character")
)

setValidity("SpatialWeightSpec", function(object) {
  E <- object@edges
  n <- object@nNodes
  msgs <- character(0)
  if (length(n) != 1L || is.na(n) || n < 1L)
    msgs <- c(msgs, "nNodes must be a single positive integer.")
  if (!is.integer(E) || (nrow(E) > 0L && ncol(E) != 2L))
    msgs <- c(msgs, "edges must be an integer matrix with two columns.")
  else if (nrow(E) > 0L) {
    if (anyNA(E) || any(E < 1L) || any(E > n))
      msgs <- c(msgs, "edge endpoints must lie in 1..nNodes.")
    else {
      if (any(E[, 1L] >= E[, 2L]))
        msgs <- c(msgs, "edges must satisfy i < j (undirected, no ",
                  "self-neighbours).")
      # duplicated() on a matrix coerces every row to a string -- measured 55
      # MB of transient on a 28k-edge graph, and twice that through new().
      # An integer pair with i < j <= n keys exactly into one double for any
      # realistic n (n^2 < 2^53), so dedup on the key instead.
      else if (anyDuplicated((as.numeric(E[, 1L]) - 1) * n + E[, 2L]))
        msgs <- c(msgs, "edges must be unique.")
    }
  }
  if (length(object@edgeWeights) != nrow(E))
    msgs <- c(msgs, "edgeWeights must have one entry per edge.")
  else if (nrow(E) > 0L &&
           (any(!is.finite(object@edgeWeights)) ||
            any(object@edgeWeights <= 0)))
    msgs <- c(msgs, "edgeWeights must be finite and positive.")
  if (length(object@regionIds) != n)
    msgs <- c(msgs, "regionIds must have one entry per region.")
  if (anyDuplicated(object@regionIds))
    msgs <- c(msgs, "regionIds must be unique.")
  if (length(msgs)) paste(msgs, collapse = " ") else TRUE
})

# canonical ordering: by second endpoint, then first (column-major upper.tri)
.canonEdgeOrder <- function(E) order(E[, 2L], E[, 1L])

#' Construct a SpatialWeightSpec from an adjacency/weight matrix
#'
#' @param adjacency Numeric n x n matrix: symmetric, zero diagonal,
#'   non-negative. Binary 0/1 contiguity is the typical case. Give either
#'   this \emph{or} \code{edges}, not both.
#' @param regionIds Optional character vector of n unique region names.
#'   Defaults to the matrix row names, or \code{"region1"}, ....
#' @param edges Alternative sparse input: a two-column matrix of node
#'   indices, one row per undirected edge (order and duplicates are
#'   normalised). Never allocates an \eqn{n \times n} matrix, so it is the
#'   route for large graphs. Requires \code{nNodes}.
#' @param nNodes Number of nodes when constructing from \code{edges}.
#' @param edgeWeights Optional positive weights, one per row of
#'   \code{edges}; defaults to 1 (binary contiguity).
#'
#' @return A validated \code{\linkS4class{SpatialWeightSpec}}.
#' @examples
#' A <- matrix(0, 3, 3); A[1, 2] <- A[2, 1] <- 1; A[2, 3] <- A[3, 2] <- 1
#' sw <- spatialWeights(A, regionIds = c("west", "centre", "east"))
#' nRegions(sw)
#' neighborsOf(sw, "centre")
#' @export
spatialWeights <- function(adjacency = NULL, regionIds = NULL, edges = NULL,
                           nNodes = NULL, edgeWeights = NULL) {
  if (is.null(adjacency) == is.null(edges))
    stop("Give exactly one of `adjacency` (dense matrix) or `edges` ",
         "(two-column edge list with `nNodes`).", call. = FALSE)

  if (!is.null(adjacency)) {
    # Dense path, backward compatible: same checks, same messages as before.
    if (!is.matrix(adjacency))
      stop("adjacency must be a matrix.", call. = FALSE)
    if (nrow(adjacency) != ncol(adjacency))
      stop("adjacency must be square.", call. = FALSE)
    if (!is.numeric(adjacency))
      stop("adjacency must be a numeric matrix.", call. = FALSE)
    if (any(!is.finite(adjacency)))
      stop("adjacency must not contain NA/NaN/Inf.", call. = FALSE)
    if (any(adjacency < 0))
      stop("adjacency weights must be non-negative.", call. = FALSE)
    if (any(diag(adjacency) != 0))
      stop("adjacency diagonal must be zero (no self-neighbours).",
           call. = FALSE)
    if (!isSymmetric(unname(adjacency)))
      stop("adjacency must be symmetric (undirected graph).", call. = FALSE)
    if (!is.null(regionIds) && length(regionIds) != nrow(adjacency))
      stop("regionIds must have one entry per region.", call. = FALSE)
    n <- nrow(adjacency)
    if (is.null(regionIds)) {
      regionIds <- rownames(adjacency)
      if (is.null(regionIds)) regionIds <- paste0("region", seq_len(n))
    }
    # which() is column-major, so this IS the canonical order already.
    ut <- which(upper.tri(adjacency) & adjacency > 0, arr.ind = TRUE)
    E <- matrix(as.integer(ut), ncol = 2L)
    w <- as.numeric(adjacency[ut])
  } else {
    # Edge path: never touches a dense matrix. This is what makes large
    # space-time graphs possible at all.
    if (is.null(nNodes))
      stop("`nNodes` is required with `edges`.", call. = FALSE)
    n <- as.integer(nNodes)
    E <- edges
    if (is.data.frame(E)) E <- as.matrix(E)
    if (!is.matrix(E) || (nrow(E) > 0L && ncol(E) != 2L))
      stop("edges must be a two-column matrix of node indices.",
           call. = FALSE)
    storage.mode(E) <- "integer"
    if (nrow(E) > 0L) {
      if (anyNA(E) || any(E < 1L) || any(E > n))
        stop("edge endpoints must lie in 1..nNodes.", call. = FALSE)
      if (any(E[, 1L] == E[, 2L]))
        stop("adjacency diagonal must be zero (no self-neighbours).",
             call. = FALSE)
      # normalise to i < j, drop duplicates, canonical order
      flip <- E[, 1L] > E[, 2L]
      if (any(flip)) E[flip, ] <- E[flip, c(2L, 1L)]
      if (is.null(edgeWeights)) edgeWeights <- rep(1, nrow(E))
      if (length(edgeWeights) != nrow(E))
        stop("edgeWeights must have one entry per edge.", call. = FALSE)
      keep <- !duplicated((as.numeric(E[, 1L]) - 1) * n + E[, 2L])
      E <- E[keep, , drop = FALSE]
      w <- as.numeric(edgeWeights[keep])
    } else {
      E <- matrix(integer(0), ncol = 2L)
      w <- numeric(0)
    }
    if (is.null(regionIds)) regionIds <- paste0("region", seq_len(n))
    if (length(regionIds) != n)
      stop("regionIds must have one entry per region.", call. = FALSE)
  }
  o <- .canonEdgeOrder(E)
  methods::new("SpatialWeightSpec",
               edges = E[o, , drop = FALSE], edgeWeights = w[o],
               nNodes = as.integer(n), regionIds = as.character(regionIds))
}

#' Rook/queen contiguity on a regular grid
#'
#' Builds the \code{\linkS4class{SpatialWeightSpec}} of a \code{nrow} x
#' \code{ncol} regular lattice under rook (shared edge) or queen (shared edge
#' or corner) contiguity -- the standard conventions of spatial econometrics
#' (Anselin 1988, Ch. 3). Regular grids with known block structure are also the
#' prescribed synthetic-graph setting for the spatial recovery tests of the
#' MRF engine planned for v0.6.0.
#'
#' @param nrow,ncol Grid dimensions (each >= 1).
#' @param contiguity \code{"rook"} (default) or \code{"queen"}.
#' @return A \code{\linkS4class{SpatialWeightSpec}} with regions named
#'   \code{"r<i>c<j>"} in row-major order.
#' @examples
#' g <- gridAdjacency(3, 3, "queen")
#' neighborsOf(g, "r2c2")   # interior cell: 8 queen neighbours
#' @export
gridAdjacency <- function(nrow, ncol, contiguity = c("rook", "queen")) {
  contiguity <- match.arg(contiguity)
  if (nrow < 1L || ncol < 1L) stop("nrow and ncol must be >= 1.", call. = FALSE)
  nrow <- as.integer(nrow); ncol <- as.integer(ncol)
  n <- nrow * ncol
  ids <- as.vector(t(outer(seq_len(nrow), seq_len(ncol),
                           function(i, j) paste0("r", i, "c", j))))
  # Build the edge list directly -- no n x n matrix. Each undirected edge is
  # emitted once via "forward" offsets only.
  cell <- function(i, j) (i - 1L) * ncol + j
  ii <- rep(seq_len(nrow), each = ncol)
  jj <- rep(seq_len(ncol), times = nrow)
  offs <- list(c(0L, 1L), c(1L, 0L))
  if (contiguity == "queen") offs <- c(offs, list(c(1L, 1L), c(1L, -1L)))
  es <- vector("list", length(offs))
  for (k in seq_along(offs)) {
    o <- offs[[k]]
    ok <- (ii + o[1L]) >= 1L & (ii + o[1L]) <= nrow &
          (jj + o[2L]) >= 1L & (jj + o[2L]) <= ncol
    es[[k]] <- cbind(cell(ii[ok], jj[ok]),
                     cell(ii[ok] + o[1L], jj[ok] + o[2L]))
  }
  E <- do.call(rbind, es)
  spatialWeights(edges = E, nNodes = n, regionIds = ids)
}
#' Number of regions in a spatial weight structure
#' @param spec A \code{\linkS4class{SpatialWeightSpec}}.
#' @return Integer number of regions.
#' @export
setGeneric("nRegions", function(spec) standardGeneric("nRegions"))

#' @rdname nRegions
#' @export
setMethod("nRegions", "SpatialWeightSpec",
          function(spec) spec@nNodes)

#' Adjacency matrix of a spatial weight structure
#' @param spec A \code{\linkS4class{SpatialWeightSpec}}.
#' @return The (named) numeric adjacency matrix.
#' @export
setGeneric("getAdjacency", function(spec) standardGeneric("getAdjacency"))

#' @rdname getAdjacency
#' @export
setMethod("getAdjacency", "SpatialWeightSpec",
  function(spec) {
    n <- spec@nNodes
    # Materialising is O(n^2) memory, which is exactly what F7 removed from
    # the pipeline; refuse rather than silently allocate gigabytes. 5000
    # nodes is ~200 MB of intermediates -- workable, but the ceiling.
    if (n > 5000L)
      stop("Refusing to materialise a dense ", n, " x ", n, " adjacency (",
           round(n^2 * 8 / 1024^2), " MB). The graph itself is available ",
           "sparsely via getEdges(); the mixture engines use that directly.",
           call. = FALSE)
    A <- matrix(0, n, n)
    E <- spec@edges
    if (nrow(E) > 0L) {
      A[E] <- spec@edgeWeights
      A[E[, c(2L, 1L), drop = FALSE]] <- spec@edgeWeights
    }
    dimnames(A) <- list(spec@regionIds, spec@regionIds)
    A
  })

#' Edge list of a spatial weight structure
#'
#' The canonical sparse form: one row per undirected edge, endpoints as
#' integer node indices with \code{i < j}, ordered by second endpoint then
#' first (the column-major order of \code{which(upper.tri(A) & A > 0)} on
#' the dense form, so downstream constants are identical either way). Unlike
#' \code{\link{getAdjacency}} this never allocates an \eqn{n \times n}
#' matrix, so it works at any graph size.
#'
#' @param spec A \code{\linkS4class{SpatialWeightSpec}}.
#' @return Integer matrix with two columns.
#' @export
setGeneric("getEdges", function(spec) standardGeneric("getEdges"))

#' @rdname getEdges
#' @export
setMethod("getEdges", "SpatialWeightSpec", function(spec) spec@edges)

#' Neighbours of one region
#' @param spec A \code{\linkS4class{SpatialWeightSpec}}.
#' @param region A region id (character) or index (integer).
#' @return Character vector of neighbouring region ids.
#' @export
setGeneric("neighborsOf", function(spec, region) standardGeneric("neighborsOf"))

#' @rdname neighborsOf
#' @export
setMethod("neighborsOf", "SpatialWeightSpec", function(spec, region) {
  ids <- spec@regionIds
  i <- if (is.character(region)) match(region, ids) else as.integer(region)
  if (is.na(i) || i < 1L || i > length(ids))
    stop("Unknown region: ", region, call. = FALSE)
  E <- spec@edges
  nb <- sort(c(E[E[, 1L] == i, 2L], E[E[, 2L] == i, 1L]))
  ids[nb]
})

#' @describeIn SpatialWeightSpec-class Compact printout: number of regions,
#'   edges, degree range, and weight type.
#' @param object A \code{SpatialWeightSpec}.
#' @export
setMethod("show", "SpatialWeightSpec", function(object) {
  E <- object@edges
  deg <- tabulate(c(E[, 1L], E[, 2L]), nbins = object@nNodes)
  cat("SpatialWeightSpec:", object@nNodes, "regions,",
      nrow(E), "undirected edges\n")
  cat("  degree range:", min(deg), "-", max(deg),
      "| weights:", if (all(object@edgeWeights == 1)) "binary (contiguity)"
                    else "general non-negative", "\n")
})

#' Build a space-time adjacency for spatio-temporal mixtures
#'
#' Expands a spatial adjacency over \code{nTime} time points, so that node
#' \eqn{(i, t)} neighbours its spatial neighbours at the same time and itself
#' at \eqn{t \pm 1}. The result is an ordinary
#' \code{\linkS4class{SpatialWeightSpec}}: pass it to
#' \code{nimixClust(..., method = "mrf", spatialWeights = )} and the Potts
#' prior couples allocations across space \emph{and} time, with no other
#' change.
#'
#' Observations must be ordered with time varying slowest -- node \eqn{(i, t)}
#' is row \eqn{(t - 1) n_{loc} + i} -- which is what \code{as.vector()} of a
#' \eqn{n_{loc} \times n_{time}} matrix gives.
#'
#' On a 5x5 grid over 8 time points with deliberately overlapping components,
#' a plain mixture failed outright (allocation accuracy 0.50, one component
#' recovered); the spatial graph reached 0.95; the space-time graph reached
#' 0.995. The temporal edges earn their place.
#'
#' Note the coupling is isotropic: one \code{beta} governs spatial and
#' temporal edges alike, because the Potts prior reads the adjacency as
#' unweighted. Misspecifying it is mild rather than fatal -- on regimes that
#' were random across space but perfectly persistent in time, imposing the
#' spatial edges anyway cost 2.5 percentage points (0.900 against 0.925 for a
#' temporal-only graph). If your structure is purely temporal, build the
#' temporal-only graph with \code{spatial = FALSE}.
#'
#' For a pure time series with no spatial component, \code{method = "hmm"} is
#' the better tool: it marginalises the state path (better mixing) and offers
#' \code{\link{viterbiPath}}. This function is for when space matters too.
#'
#' @section Scale limit:
#' The graph itself no longer limits the problem size:
#' \code{SpatialWeightSpec} stores an edge list, so a 10 000-node space-time
#' graph builds in a fraction of a second within a few tens of MB (it was
#' OOM-killed outright before v1.5.0, during construction), and 50 000 nodes
#' take about a second. The binding constraint is now NIMBLE's per-node model
#' memory during \code{nimbleModel}/compilation: measured in a 4 GB
#' container, fits ran at 5000 nodes (~1.6 GB) and 7000 nodes (~2.8 GB) but
#' died near 10 000. If you need more, more RAM buys it roughly linearly --
#' the graph is no longer the wall.
#'
#' @param spaceWeights A \code{\linkS4class{SpatialWeightSpec}} over
#'   \eqn{n_{loc}} locations, e.g. from \code{\link{gridAdjacency}}.
#' @param nTime Integer, number of time points (>= 2).
#' @param spatial Logical; keep the spatial edges within each time point.
#'   \code{FALSE} gives a temporal-only graph (independent chains per
#'   location).
#' @param temporal Logical; keep the edges linking consecutive time points.
#' @return A \code{\linkS4class{SpatialWeightSpec}} over
#'   \eqn{n_{loc} \times n_{time}} nodes.
#' @examples
#' W  <- gridAdjacency(3, 3)
#' ST <- spacetimeAdjacency(W, nTime = 4)
#' nrow(getAdjacency(ST))   # 36 = 9 locations x 4 times
#' @seealso \code{\link{gridAdjacency}}, \code{\link{viterbiPath}}
#' @export
spacetimeAdjacency <- function(spaceWeights, nTime, spatial = TRUE,
                               temporal = TRUE) {
  Es <- getEdges(spaceWeights)
  nLoc <- nRegions(spaceWeights)
  nTime <- as.integer(nTime)
  if (length(nTime) != 1L || is.na(nTime) || nTime < 2L)
    stop("nTime must be a single integer >= 2.", call. = FALSE)
  if (!spatial && !temporal)
    stop("A graph with neither spatial nor temporal edges has no edges; ",
         "use method = 'fixedk' for an unstructured mixture.", call. = FALSE)
  N <- nLoc * nTime
  # Built directly as an edge list: the dense version of this function was
  # what died first at scale -- graph construction alone was OOM-killed at
  # 10 000 nodes, before any model existed.
  parts <- list()
  if (spatial && nrow(Es) > 0L) {
    offs <- (seq_len(nTime) - 1L) * nLoc
    parts$space <- cbind(rep(Es[, 1L], nTime) + rep(offs, each = nrow(Es)),
                         rep(Es[, 2L], nTime) + rep(offs, each = nrow(Es)))
  }
  if (temporal && nTime >= 2L) {
    a <- rep(seq_len(nLoc), nTime - 1L) +
         rep((seq_len(nTime - 1L) - 1L) * nLoc, each = nLoc)
    parts$time <- cbind(a, a + nLoc)
  }
  E <- do.call(rbind, parts)
  if (is.null(E) || nrow(E) == 0L)
    stop("The resulting graph has no edges: check spaceWeights and nTime.",
         call. = FALSE)
  spatialWeights(edges = E, nNodes = N)
}
