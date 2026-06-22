## ---------------------------------------------------------------------------
## class-MixtureModel.R
##
## MixtureModel bundles the data, the component DistributionSpec, the
## truncation level K_max, and the engine into one object that the fitting
## pipeline consumes. In v0.1.0 only ClusterModel (for mixClust) is concrete;
## RegressionMixModel is deferred to v0.3.0.
## ---------------------------------------------------------------------------

#' Virtual base class for mixture models
#'
#' @slot data Numeric data (vector for univariate clustering in v0.1.0).
#' @slot distSpec A \code{\linkS4class{DistributionSpec}}.
#' @slot engine An \code{\linkS4class{EngineConfig}}.
#' @slot Kmax Integer truncation level for the number of components.
#' @slot prior A named list of prior hyperparameters.
#' @export
setClass(
  "MixtureModel",
  representation(
    "VIRTUAL",
    data     = "ANY",
    distSpec = "DistributionSpec",
    engine   = "EngineConfig",
    Kmax     = "integer",
    prior    = "list"
  )
)

#' Mixture clustering model (for nimixClust)
#'
#' @slot data Numeric vector (univariate) or numeric matrix with one row per
#'   observation and one column per dimension (multivariate, v0.2.0).
#' @export
setClass(
  "ClusterModel",
  contains = "MixtureModel",
  representation(data = "ANY")
)

#' Construct a ClusterModel
#'
#' @param data Numeric vector (univariate) or numeric matrix (multivariate,
#'   one row per observation).
#' @param distSpec A \code{\linkS4class{DistributionSpec}}.
#' @param engine An \code{\linkS4class{EngineConfig}}.
#' @param Kmax Integer truncation level.
#' @param prior Named list of prior hyperparameters.
#' @return A \code{ClusterModel}.
#' @keywords internal
ClusterModel <- function(data, distSpec, engine, Kmax, prior) {
  # Preserve matrix structure for multivariate data; coerce only true vectors.
  if (!is.matrix(data)) data <- as.numeric(data)
  new("ClusterModel",
      data = data, distSpec = distSpec, engine = engine,
      Kmax = as.integer(Kmax), prior = prior)
}

#' Mixture-of-regressions model (for nimixReg)
#'
#' @slot data Numeric response vector.
#' @slot formula The model formula.
#' @slot X The numeric design matrix (one row per observation, including the
#'   intercept column when present).
#' @export
setClass(
  "RegressionMixModel",
  contains = "MixtureModel",
  representation(data = "ANY", formula = "ANY", X = "matrix")
)

#' Construct a RegressionMixModel
#'
#' @param data Numeric response vector.
#' @param X Numeric design matrix.
#' @param formula The model formula.
#' @param distSpec A \code{\linkS4class{DistributionSpec}}.
#' @param engine An \code{\linkS4class{EngineConfig}}.
#' @param Kmax Integer truncation level.
#' @param prior Named list of prior hyperparameters.
#' @return A \code{RegressionMixModel}.
#' @keywords internal
RegressionMixModel <- function(data, X, formula, distSpec, engine, Kmax,
                               prior) {
  new("RegressionMixModel",
      data = if (is.matrix(data)) data else as.numeric(data),
      X = as.matrix(X), formula = formula,
      distSpec = distSpec, engine = engine, Kmax = as.integer(Kmax),
      prior = prior)
}
