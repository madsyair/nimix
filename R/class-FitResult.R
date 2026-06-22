## ---------------------------------------------------------------------------
## class-FitResult.R
##
## FitResult is the object returned by nimixClust(). It stores raw MCMC
## output plus derived quantities (per-iteration number of occupied clusters,
## cluster allocations) and a cache for relabelled summaries.
##
## IMPORTANT: raw per-component
## posterior summaries are meaningless under label switching. summary()/plot()
## therefore route through relabel() before reporting component parameters.
## ---------------------------------------------------------------------------

#' Fitted nimix mixture result
#'
#' @slot mcmcSamples A matrix of monitored MCMC draws (iterations x parameters).
#' @slot Kposterior Integer vector: number of occupied clusters per iteration.
#' @slot clusterAllocation Integer matrix (iterations x n) of raw cluster
#'   indicators \code{xi}.
#' @slot paramTrace A named list of raw cluster-parameter traces
#'   (each an iterations x Kmax matrix).
#' @slot engineUsed Character scalar naming the engine.
#' @slot distSpec The \code{\linkS4class{DistributionSpec}} used.
#' @slot data The data the model was fit to.
#' @slot Kmax Integer truncation level.
#' @slot prior The prior list used.
#' @slot relabeled A list cache populated by \code{\link{relabel}} (or empty).
#' @slot mcmcControl The MCMC control list actually used.
#' @slot call The matched call.
#' @export
setClass(
  "FitResult",
  representation(
    mcmcSamples       = "matrix",
    Kposterior        = "integer",
    clusterAllocation = "matrix",
    paramTrace        = "list",
    engineUsed        = "character",
    distSpec          = "DistributionSpec",
    data              = "ANY",
    Kmax              = "integer",
    prior             = "list",
    relabeled         = "list",
    mcmcControl       = "list",
    call              = "ANY"
  ),
  prototype = prototype(relabeled = list())
)

#' @describeIn FitResult Compact display of a fitted result.
#' @param object A \code{FitResult}.
#' @export
setMethod("show", "FitResult", function(object) {
  kt <- sort(table(object@Kposterior), decreasing = TRUE)
  modalK <- as.integer(names(kt)[1])
  cat("<nimix FitResult>\n")
  cat(sprintf("  engine        : %s\n", object@engineUsed))
  cat(sprintf("  distribution  : %s\n", object@distSpec@name))
  cat(sprintf("  observations  : %d (dimension d = %d)\n",
              .nObs(object@data), .dataDimOf(object@data)))
  cat(sprintf("  K_max (trunc) : %d\n", object@Kmax))
  cat(sprintf("  posterior draws: %d\n", nrow(object@mcmcSamples)))
  cat(sprintf("  modal #clusters: %d (%.1f%% of draws)\n",
              modalK, 100 * kt[1] / sum(kt)))
  cat("  Use summary(fit) for relabelled component estimates,\n")
  cat("  plot(fit, type=...) for diagnostics, predict(fit, newdata=...).\n")
  invisible(object)
})
