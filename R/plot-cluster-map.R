# plot-cluster-map.R
#
# Map visualisation of a fitted clustering when the observations are regions
# with known geometry (a shapefile / sf object). Natural companion to the MRF
# engine -- whose whole point is spatial coherence -- but works for any
# clustering fit whose observations correspond one-to-one with map features.
#
# sf is Suggests-only: everything here is guarded by requireNamespace("sf").

#' Map the clusters of a spatial mixture fit
#'
#' Draws the fitted partition as a choropleth: each map feature (region) is
#' coloured by its cluster label. Optionally shades each region by allocation
#' uncertainty, so that regions whose membership the posterior is unsure about
#' stand out.
#'
#' @param fit A clustering \code{\linkS4class{FitResult}} (any engine; the MRF
#'   engine is the usual source).
#' @param shp The map: either a path to a shapefile (\code{.shp}, read with
#'   \code{sf::st_read}) or an \code{sf} object already in memory. Features
#'   must correspond to observations.
#' @param partition How to summarise the posterior into one label per region:
#'   \code{"binder"} (default) uses \code{\link{binderPartition}}, the draw
#'   minimising Binder loss against the posterior similarity matrix -- a
#'   label-invariant summary; \code{"modal"} uses the per-region posterior mode
#'   of the allocation trace.
#' @param idCol Optional name of a column in \code{shp} giving the observation
#'   index (1-based) of each feature. Default \code{NULL} assumes the features
#'   are already in observation order.
#' @param uncertainty Logical; if \code{TRUE}, regions are shaded towards white
#'   in proportion to their allocation entropy (0 = certain, log K =
#'   maximally uncertain), so pale regions are the ones the posterior cannot
#'   place.
#' @param palette Colours, one per cluster. Default is
#'   \code{grDevices::hcl.colors(K, "Dark 3")}.
#' @param main Plot title.
#' @param legendPos Legend position keyword (see \code{\link[graphics]{legend}}),
#'   or \code{NA} to suppress the legend.
#' @param ... Passed on to \code{plot} of the sf geometry (e.g. \code{border},
#'   \code{lwd}).
#' @return Invisibly, a data.frame with one row per region: \code{cluster} and
#'   (if requested) \code{entropy}. Called for its side effect, the plot.
#' @examples
#' \dontrun{
#' sw  <- spatialWeights(nb)                      # neighbourhood used in the fit
#' fit <- nimixClust(y, K = 3, method = "mrf", spatialWeights = sw)
#' plotClusterMap(fit, "regions.shp")             # from a shapefile on disk
#' plotClusterMap(fit, mysf, uncertainty = TRUE)  # from an sf object
#' }
#' @export
plotClusterMap <- function(fit, shp, partition = c("binder", "modal"),
                           idCol = NULL, uncertainty = FALSE, palette = NULL,
                           main = NULL, legendPos = "topright", ...) {
  if (!methods::is(fit, "FitResult"))
    stop("plotClusterMap() expects a FitResult.", call. = FALSE)
  if (!requireNamespace("sf", quietly = TRUE))
    stop("plotClusterMap() needs the 'sf' package for map geometry. ",
         "Install it with install.packages(\"sf\").", call. = FALSE)
  partition <- match.arg(partition)

  geom <- if (inherits(shp, "sf")) shp
  else if (is.character(shp) && length(shp) == 1L) {
    if (!file.exists(shp))
      stop("Shapefile not found: ", shp, call. = FALSE)
    sf::st_read(shp, quiet = TRUE)
  } else
    stop("`shp` must be a path to a shapefile or an sf object.", call. = FALSE)

  A <- fit@clusterAllocation
  if (is.null(A) || !nrow(A))
    stop("The fit carries no allocation trace; was it run with monitors ",
         "including the allocation node?", call. = FALSE)
  n <- ncol(A)
  if (nrow(geom) != n)
    stop("The map has ", nrow(geom), " features but the fit has ", n,
         " observations; they must correspond one-to-one.", call. = FALSE)

  # Feature -> observation mapping.
  ord <- seq_len(n)
  if (!is.null(idCol)) {
    if (!idCol %in% names(geom))
      stop("Column '", idCol, "' not found in the map attributes.",
           call. = FALSE)
    ord <- as.integer(geom[[idCol]])
    if (anyNA(ord) || !setequal(ord, seq_len(n)))
      stop("`idCol` must hold a permutation of 1..", n, ".", call. = FALSE)
  }

  # One label per region, label-invariantly.
  z <- if (partition == "binder") binderPartition(fit)$partition
       else apply(A, 2L, function(col) {
         tb <- tabulate(col)
         which.max(tb)
       })
  z <- as.integer(factor(z))          # compact labels 1..K for colouring
  K <- max(z)

  if (is.null(palette)) palette <- grDevices::hcl.colors(K, "Dark 3")
  if (length(palette) < K)
    stop("palette has ", length(palette), " colours but the partition has ",
         K, " clusters.", call. = FALSE)
  cols <- palette[z][ord]

  ent <- NULL
  if (isTRUE(uncertainty)) {
    # Per-region allocation entropy over the posterior draws, normalised by
    # log K so it lies in [0, 1]; used to fade uncertain regions towards white.
    ent <- apply(A, 2L, function(col) {
      p <- tabulate(col) / length(col)
      p <- p[p > 0]
      -sum(p * log(p))
    })
    if (K > 1L) ent <- ent / log(K)
    mixWhite <- function(col, w) {
      rgb0 <- grDevices::col2rgb(col) / 255
      grDevices::rgb(rgb0[1] * (1 - w) + w, rgb0[2] * (1 - w) + w,
                     rgb0[3] * (1 - w) + w)
    }
    cols <- vapply(seq_len(n), function(i) mixWhite(cols[i], 0.75 * ent[ord][i]),
                   character(1))
  }

  if (is.null(main))
    main <- sprintf("Cluster map (%s partition, K = %d)", partition, K)
  plot(sf::st_geometry(geom), col = cols, main = main, ...)
  if (!is.na(legendPos) && length(legendPos) == 1L)
    graphics::legend(legendPos, legend = paste("Cluster", seq_len(K)),
                     fill = palette[seq_len(K)], bty = "n", cex = 0.9)

  invisible(data.frame(cluster = z,
                       entropy = if (is.null(ent)) NA_real_ else ent))
}
