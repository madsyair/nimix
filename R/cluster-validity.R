## cluster-validity.R -----------------------------------------------------------
## Internal cluster-validity indices for a fitted clustering, computed on the
## Binder/Dahl point partition (or any partition the user supplies) via the
## established backends 'cluster' (silhouette) and 'fpc' (Dunn,
## Calinski-Harabasz). Both live in Suggests: nothing here is needed to fit or
## summarise a model, and the package must pass checks without them.
##
## STATISTICAL CAVEAT (deliberate, documented, and asserted in the tests):
## these indices assume geometric clusters -- compact, well-separated in
## distance. nimix mixtures are *density-based*: two overlapping Gaussian
## components can be exactly the right model and still earn a low silhouette.
## So these numbers are a secondary, model-free lens for comparing partitions,
## NOT a verdict on model correctness; posterior predictive checks (ppCheck)
## and partition uncertainty (psm) remain the primary tools.

#' Internal validity indices for a fitted clustering
#'
#' Computes standard internal cluster-validity indices -- silhouette width,
#' Dunn index, Calinski-Harabasz -- for the point partition of a clustering
#' fit. By default the partition is [binderPartition()]'s least-squares
#' partition, so every posterior draw informs it and no relabelling is
#' required.
#'
#' @param fit A [FitResult][FitResult-class] from a *clustering* fit
#'   ([nimixClust()]). Regression fits are refused: distance-based indices
#'   have no meaning for mixtures of regressions.
#' @param metrics Character vector, any of `"silhouette"`, `"dunn"`, `"ch"`.
#' @param partition Optional integer vector of cluster labels (length `n`).
#'   Default: `binderPartition(fit)$partition`.
#' @param dist Optional `stats::dist` object. Default: Euclidean distance on
#'   the fitted data. Supply your own for scaled or non-Euclidean analyses.
#'
#' @return Named numeric vector with one entry per requested metric.
#'
#' @section Interpretation caveat:
#' Internal indices reward *geometric* separation. A mixture model is
#' *density*-based, and a fit with genuinely overlapping components -- often
#' the scientifically correct model -- will score a low silhouette even when
#' [ppCheck()] says the model reproduces the data perfectly. Treat these
#' indices as a secondary lens (useful, e.g., for comparing partitions from
#' different `K` on equal footing), never as a model-adequacy verdict; that
#' job belongs to [ppCheck()] and [psm()].
#'
#' For the full battery of indices beyond these three, call the backends
#' directly, e.g. `fpc::cluster.stats(d, part)` or
#' `clusterCrit::intCriteria(X, part, "all")` -- both accept exactly the
#' `(partition, distance/data)` pair this function assembles.
#'
#' @examples
#' \donttest{
#' y <- c(rnorm(60, -3), rnorm(60, 3))
#' fit <- nimixClust(y, K_max = 6,
#'                   mcmcControl = list(niter = 800, nburnin = 300),
#'                   verbose = FALSE)
#' clusterValidity(fit)
#' }
#' @seealso [binderPartition()], [psm()], [ppCheck()]
#' @export
clusterValidity <- function(fit,
                            metrics = c("silhouette", "dunn", "ch"),
                            partition = NULL, dist = NULL) {
  stopifnot(is(fit, "FitResult"))
  metrics <- match.arg(metrics, several.ok = TRUE)
  if (isRegressionSpec(fit@distSpec))
    stop("clusterValidity() is for clustering fits; distance-based indices ",
         "have no meaning for a mixture of regressions.", call. = FALSE)

  X <- as.matrix(fit@data)
  if (is.null(partition)) partition <- binderPartition(fit)$partition
  partition <- as.integer(partition)
  if (length(partition) != nrow(X))
    stop("`partition` has length ", length(partition), " but the fit has ",
         nrow(X), " observations.", call. = FALSE)
  K <- length(unique(partition))
  if (K < 2L)
    stop("The partition has a single cluster; internal validity indices ",
         "need at least two.", call. = FALSE)

  if (is.null(dist)) dist <- stats::dist(X)
  if (!inherits(dist, "dist"))
    stop("`dist` must be a stats::dist object.", call. = FALSE)

  out <- c()

  if ("silhouette" %in% metrics) {
    if (!requireNamespace("cluster", quietly = TRUE))
      stop("metric \"silhouette\" needs the 'cluster' package; ",
           "install.packages(\"cluster\").", call. = FALSE)
    sil <- cluster::silhouette(partition, dist)
    out["silhouette"] <- mean(sil[, "sil_width"])
  }

  if (any(c("dunn", "ch") %in% metrics)) {
    if (!requireNamespace("fpc", quietly = TRUE))
      stop("metrics \"dunn\"/\"ch\" need the 'fpc' package; ",
           "install.packages(\"fpc\").", call. = FALSE)
    cs <- fpc::cluster.stats(dist, partition, silhouette = FALSE)
    if ("dunn" %in% metrics) out["dunn"] <- cs$dunn
    if ("ch"   %in% metrics) out["ch"]   <- cs$ch
  }

  out[metrics]
}
