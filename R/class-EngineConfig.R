## ---------------------------------------------------------------------------
## class-EngineConfig.R
##
## EngineConfig abstracts how the number of components is handled.
##   FixedKEngine : finite mixture with a fixed, known K (Dirichlet weights +
##                  categorical allocation; the simplest, fastest baseline).
##   DPMEngine    : Dirichlet Process Mixture; K is inferred via NIMBLE's native
##                  Chinese Restaurant Process. These CRP samplers are optimised
##                  by the NIMBLE core team and must not be reimplemented by hand.
##
## A reversible-jump engine (custom trans-dimensional sampler) is planned for a
## later release and is intentionally absent here.
## ---------------------------------------------------------------------------

#' Virtual base class for mixture sampling engines
#'
#' @slot name Character scalar identifying the engine.
#' @export
setClass(
  "EngineConfig",
  representation("VIRTUAL", name = "character")
)

#' Run a mixture engine on a model (internal generic)
#'
#' Dispatches on the engine: \code{\linkS4class{DPMEngine}} builds a CRP model,
#' \code{\linkS4class{FixedKEngine}} builds a finite-mixture model. Returns the
#' raw pieces used to construct a \code{\linkS4class{FitResult}}.
#'
#' @param engine An \code{\linkS4class{EngineConfig}}.
#' @param model A \code{\linkS4class{MixtureModel}}.
#' @param mcmcControl,initMethod,seed,verbose Sampler controls.
#' @param ... Reserved.
#' @return A named list with the MCMC samples, the posterior of the number of
#'   occupied components, the allocation matrix, the parsed parameter traces,
#'   and the resolved MCMC control list.
#' @keywords internal
setGeneric("runEngine",
  function(engine, model, mcmcControl = list(), initMethod = "kmeans",
           seed = 1L, verbose = TRUE, ...) {
    standardGeneric("runEngine")
  })

#' Dirichlet Process Mixture engine (native NIMBLE dCRP)
#'
#' Wraps NIMBLE's Chinese Restaurant Process distribution \code{dCRP} and the
#' specialised \code{CRP}, \code{CRP_cluster_wrapper}, and
#' \code{CRP_concentration} samplers. Empty components are handled natively:
#' the collapsed sampler only updates parameters of occupied clusters.
#'
#' @slot name Engine identifier, fixed to \code{"dpm"}.
#' @slot concPrior A length-2 numeric \code{c(shape, rate)} for the Gamma
#'   hyperprior on the concentration parameter \eqn{\alpha}. A Gamma prior is
#'   required for NIMBLE to assign the \code{CRP_concentration} sampler, so the
#'   data can inform the level of concentration rather than fixing it.
#'
#' @references
#' Neal, R.M. (2000). Markov chain sampling methods for Dirichlet process
#' mixture models. \emph{JCGS}, 9(2), 249--265.
#' \doi{10.1080/10618600.2000.10474879}
#'
#' Ferguson, T.S. (1973). A Bayesian analysis of some nonparametric problems.
#' \emph{The Annals of Statistics}, 1(2), 209--230.
#' \doi{10.1214/aos/1176342360}
#'
#' @export
setClass(
  "DPMEngine",
  contains = "EngineConfig",
  representation(concPrior = "numeric"),
  prototype = prototype(name = "dpm", concPrior = c(2, 4))
)

#' Construct a DPM engine configuration
#'
#' @param concPrior Length-2 numeric \code{c(shape, rate)} for the Gamma prior
#'   on the DP concentration \eqn{\alpha}. Defaults to \code{c(2, 4)}
#'   (weakly informative, prior mean 0.5).
#' @return A \code{\linkS4class{DPMEngine}} object.
#' @examples
#' eng <- DPMEngine(concPrior = c(2, 4))
#' eng
#' @export
DPMEngine <- function(concPrior = c(2, 4)) {
  if (length(concPrior) != 2L || any(concPrior <= 0))
    stop("concPrior must be a length-2 positive numeric c(shape, rate).",
         call. = FALSE)
  new("DPMEngine", name = "dpm", concPrior = as.numeric(concPrior))
}

#' Finite-mixture engine with fixed, known K
#'
#' The simplest engine: the number of components K is fixed (not inferred).
#' Mixing weights have a symmetric Dirichlet prior and each observation has a
#' categorical allocation. Because there is no Chinese Restaurant Process, the
#' truncation considerations of the DPM do not apply, and NIMBLE assigns
#' conjugate samplers to the weights and component parameters. Useful as a fast
#' baseline when K is known or assumed, and for classical model selection by
#' comparing fits across several values of K.
#'
#' @slot name Engine identifier, fixed to \code{"fixedk"}.
#' @slot dirichletConc Positive scalar concentration of the symmetric Dirichlet
#'   prior on the mixing weights (\eqn{1} is uniform on the simplex; values
#'   below 1 favour sparser weight vectors).
#'
#' @references
#' McLachlan, G.J., & Peel, D. (2000). \emph{Finite Mixture Models}. Wiley.
#' \doi{10.1002/0471721182}
#'
#' Frühwirth-Schnatter, S. (2006). \emph{Finite Mixture and Markov Switching
#' Models}. Springer. \doi{10.1007/978-0-387-35768-3}
#'
#' @export
setClass(
  "FixedKEngine",
  contains = "EngineConfig",
  representation(dirichletConc = "numeric"),
  prototype = prototype(name = "fixedk", dirichletConc = 1)
)

#' Construct a fixed-K finite-mixture engine configuration
#'
#' @param dirichletConc Positive scalar concentration of the symmetric Dirichlet
#'   prior on the mixing weights. Defaults to \code{1} (uniform on the simplex).
#' @return A \code{\linkS4class{FixedKEngine}} object.
#' @examples
#' eng <- FixedKEngine(dirichletConc = 1)
#' eng
#' @export
FixedKEngine <- function(dirichletConc = 1) {
  if (length(dirichletConc) != 1L || dirichletConc <= 0)
    stop("dirichletConc must be a positive scalar.", call. = FALSE)
  new("FixedKEngine", name = "fixedk",
      dirichletConc = as.numeric(dirichletConc))
}
