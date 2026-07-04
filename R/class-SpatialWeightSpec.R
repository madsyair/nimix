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
#' @slot adjacency Numeric matrix (n x n): symmetric, zero diagonal,
#'   non-negative entries. Entry (i, j) > 0 means regions i and j are
#'   neighbours (with that weight).
#' @slot regionIds Character vector of length n naming the regions (defaults
#'   to row names of the matrix, or \code{"region1"}, ... when absent).
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
setClass(
  "SpatialWeightSpec",
  representation(adjacency = "matrix", regionIds = "character")
)

setValidity("SpatialWeightSpec", function(object) {
  A <- object@adjacency
  n <- nrow(A)
  msgs <- character(0)
  if (!is.numeric(A)) msgs <- c(msgs, "adjacency must be a numeric matrix.")
  if (n != ncol(A))   msgs <- c(msgs, "adjacency must be square.")
  else {
    if (any(!is.finite(A)))
      msgs <- c(msgs, "adjacency must not contain NA/NaN/Inf.")
    else {
      if (any(A < 0))
        msgs <- c(msgs, "adjacency weights must be non-negative.")
      if (any(diag(A) != 0))
        msgs <- c(msgs, "adjacency diagonal must be zero (no self-neighbours).")
      if (!isSymmetric(unname(A)))
        msgs <- c(msgs, "adjacency must be symmetric (undirected graph).")
    }
  }
  if (length(object@regionIds) != n)
    msgs <- c(msgs, "regionIds must have one entry per region.")
  if (anyDuplicated(object@regionIds))
    msgs <- c(msgs, "regionIds must be unique.")
  if (length(msgs)) msgs else TRUE
})

#' Construct a SpatialWeightSpec from an adjacency/weight matrix
#'
#' @param adjacency Numeric n x n matrix: symmetric, zero diagonal,
#'   non-negative. Binary 0/1 contiguity is the typical case.
#' @param regionIds Optional character vector of n unique region names.
#'   Defaults to the matrix row names, or \code{"region1"}, ....
#'
#' @return A validated \code{\linkS4class{SpatialWeightSpec}}.
#' @examples
#' A <- matrix(0, 3, 3); A[1, 2] <- A[2, 1] <- 1; A[2, 3] <- A[3, 2] <- 1
#' sw <- spatialWeights(A, regionIds = c("west", "centre", "east"))
#' nRegions(sw)
#' neighborsOf(sw, "centre")
#' @export
spatialWeights <- function(adjacency, regionIds = NULL) {
  if (!is.matrix(adjacency))
    stop("adjacency must be a matrix.", call. = FALSE)
  if (nrow(adjacency) != ncol(adjacency))
    stop("adjacency must be square.", call. = FALSE)
  if (!is.null(regionIds) && length(regionIds) != nrow(adjacency))
    stop("regionIds must have one entry per region.", call. = FALSE)
  storage.mode(adjacency) <- "double"
  if (is.null(regionIds)) {
    regionIds <- rownames(adjacency)
    if (is.null(regionIds)) regionIds <- paste0("region", seq_len(nrow(adjacency)))
  }
  dimnames(adjacency) <- list(regionIds, regionIds)
  methods::new("SpatialWeightSpec", adjacency = adjacency,
               regionIds = as.character(regionIds))
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
  n <- nrow * ncol
  ids <- as.vector(t(outer(seq_len(nrow), seq_len(ncol),
                           function(i, j) paste0("r", i, "c", j))))
  A <- matrix(0, n, n)
  cell <- function(i, j) (i - 1L) * ncol + j
  offs <- if (contiguity == "rook")
    list(c(-1L, 0L), c(1L, 0L), c(0L, -1L), c(0L, 1L))
  else
    list(c(-1L, 0L), c(1L, 0L), c(0L, -1L), c(0L, 1L),
         c(-1L, -1L), c(-1L, 1L), c(1L, -1L), c(1L, 1L))
  for (i in seq_len(nrow)) for (j in seq_len(ncol)) for (o in offs) {
    ii <- i + o[1L]; jj <- j + o[2L]
    if (ii >= 1L && ii <= nrow && jj >= 1L && jj <= ncol)
      A[cell(i, j), cell(ii, jj)] <- 1
  }
  spatialWeights(A, regionIds = ids)
}

#' Number of regions in a spatial weight structure
#' @param spec A \code{\linkS4class{SpatialWeightSpec}}.
#' @return Integer number of regions.
#' @export
setGeneric("nRegions", function(spec) standardGeneric("nRegions"))

#' @rdname nRegions
#' @export
setMethod("nRegions", "SpatialWeightSpec",
          function(spec) nrow(spec@adjacency))

#' Adjacency matrix of a spatial weight structure
#' @param spec A \code{\linkS4class{SpatialWeightSpec}}.
#' @return The (named) numeric adjacency matrix.
#' @export
setGeneric("getAdjacency", function(spec) standardGeneric("getAdjacency"))

#' @rdname getAdjacency
#' @export
setMethod("getAdjacency", "SpatialWeightSpec",
          function(spec) spec@adjacency)

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
  ids[spec@adjacency[i, ] > 0]
})

#' @describeIn SpatialWeightSpec-class Compact printout: number of regions,
#'   edges, degree range, and weight type.
#' @param object A \code{SpatialWeightSpec}.
#' @export
setMethod("show", "SpatialWeightSpec", function(object) {
  A <- object@adjacency
  deg <- rowSums(A > 0)
  cat("SpatialWeightSpec:", nrow(A), "regions,",
      sum(A > 0) / 2, "undirected edges\n")
  cat("  degree range:", min(deg), "-", max(deg),
      "| weights:", if (all(A %in% c(0, 1))) "binary (contiguity)"
                    else "general non-negative", "\n")
})
