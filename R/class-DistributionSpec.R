## ---------------------------------------------------------------------------
## class-DistributionSpec.R
##
## The DistributionSpec is the unit of extensibility in nimix (project
## knowledge ). A concrete spec describes ONE family of mixture
## component distributions: its parameter names, a sensible data-scaled default
## prior, a parameter validator, and -- through the buildModelCode() generic --
## the NIMBLE code fragment for that component under a given engine.
##
## Distributions registered for the DPM engine in v0.1.0-v0.4.0 are reused
## verbatim across engines (DPM and fixed-K share the same component contract),
## never a rewrite of DistributionSpec.
## ---------------------------------------------------------------------------

#' Virtual base class for mixture component distributions
#'
#' \code{DistributionSpec} is the abstract S4 class that every component
#' distribution in nimix extends. It is never instantiated directly; use a
#' concrete subclass such as \code{\linkS4class{NormalUvSpec}}.
#'
#' @slot name Character scalar, a short identifier (e.g. \code{"normal-uv"}).
#' @slot paramNames Character vector of component parameter names.
#' @slot priorSpec A named list of prior hyperparameters. Empty until
#'   \code{\link{defaultPrior}} (or the user) fills it.
#' @slot dataDim Integer, the data dimension the spec is meant for (1 for
#'   univariate). Used for early validation in \code{\link{nimixClust}}.
#'
#' @references
#' Frühwirth-Schnatter, S. (2006). \emph{Finite Mixture and Markov Switching
#' Models}. Springer. \doi{10.1007/978-0-387-35768-3}
#'
#' @seealso \code{\linkS4class{NormalUvSpec}}, \code{\link{registerDistribution}}
#' @export
setClass(
  "DistributionSpec",
  representation(
    "VIRTUAL",
    name       = "character",
    paramNames = "character",
    priorSpec  = "list",
    dataDim    = "integer"
  ),
  prototype = prototype(priorSpec = list(), dataDim = 1L)
)

# --- Generics --------------------------------------------------------------

#' Validate component parameters or a prior specification
#'
#' Checks that a candidate parameter / prior list is internally consistent for
#' a given \code{\linkS4class{DistributionSpec}}. Each concrete distribution
#' must implement this.
#'
#' @param spec A \code{\linkS4class{DistributionSpec}}.
#' @param params A named list of parameters or hyperparameters to validate.
#' @param ... Reserved for methods.
#' @return Invisibly \code{TRUE} if valid; otherwise an error is raised.
#' @export
setGeneric("validateParams", function(spec, params, ...) {
  standardGeneric("validateParams")
})

#' Build a data-scaled default prior for a distribution
#'
#' Returns a named list of prior hyperparameters scaled to the observed data,
#' following the weakly-informative, data-scaled philosophy in project
#' knowledge (priors for location parameters must not be made
#' arbitrarily vague when \code{K_max} is large).
#'
#' @param spec A \code{\linkS4class{DistributionSpec}}.
#' @param data Numeric data used to scale the prior.
#' @param control A named list of user overrides (merged over the defaults).
#' @param ... Reserved for methods.
#' @return A named list of prior hyperparameters.
#' @export
setGeneric("defaultPrior", function(spec, data, control = list(), ...) {
  standardGeneric("defaultPrior")
})

#' Build the NIMBLE model code for a (distribution, engine) pair
#'
#' Dispatches on BOTH the component distribution and the engine. This is the
#' extensibility seam described in adding a new
#' distribution to the DPM engine means adding one method here, not editing the
#' engine or the spec base class.
#'
#' @param spec A \code{\linkS4class{DistributionSpec}}.
#' @param engine An \code{\linkS4class{EngineConfig}}.
#' @param n Integer, number of observations.
#' @param L Integer, cluster-parameter truncation length (\code{= K_max}).
#' @param ... Reserved for methods.
#' @return A list with elements \code{code} (a \code{nimbleCode} object),
#'   \code{monitors} (character vector), and \code{paramNodes}
#'   (named character vector mapping logical parameter names to model nodes).
#' @export
setGeneric("buildModelCode", function(spec, engine, n, L, ...) {
  standardGeneric("buildModelCode")
})

#' Simulate component parameters from a prior (for inits / recovery tests)
#'
#' @param spec A \code{\linkS4class{DistributionSpec}}.
#' @param prior A prior list (typically from \code{\link{defaultPrior}}).
#' @param nClust Integer number of components to simulate.
#' @param ... Reserved for methods.
#' @return A named list of simulated parameter vectors/matrices.
#' @keywords internal
setGeneric("simulateParams", function(spec, prior, nClust, ...) {
  standardGeneric("simulateParams")
})

#' Component density evaluator (R-level, for posterior predictive checks)
#'
#' Returns a function \code{f(x, params)} giving the component density at
#' \code{x}. Used by \code{predict()} to build posterior predictive densities.
#' For univariate specs \code{x} is a scalar and \code{params} a list with
#' scalar entries; for multivariate specs \code{x} is a length-\code{d} vector
#' and \code{params$mu} a vector, \code{params$Sigma} a \code{d x d} matrix.
#'
#' @param spec A \code{\linkS4class{DistributionSpec}}.
#' @param ... Reserved for methods.
#' @return A function of \code{(x, params)}.
#' @keywords internal
setGeneric("componentDensity", function(spec, ...) {
  standardGeneric("componentDensity")
})

# --- Engine-facing generics (dimension-agnostic DPM orchestration) ---------
# These let engine-dpm.R stay free of any univariate/multivariate branching:
# every dimension-specific decision (how priors map to NIMBLE constants, how
# data and inits are shaped, how raw parameter traces are parsed, and how
# relabelled component summaries are built) lives on the spec. Adding a new
# distribution therefore means adding methods here, never editing the engine
# (DistributionSpec is the extensibility unit).

#' Assemble distribution-specific NIMBLE constants for the DPM engine
#'
#' @param spec A \code{\linkS4class{DistributionSpec}}.
#' @param prior The prior list (from \code{\link{defaultPrior}}).
#' @param n Integer number of observations.
#' @param ... Reserved for methods.
#' @return A named list of constants (excluding the concentration hyperprior,
#'   which the engine appends).
#' @keywords internal
setGeneric("buildConstants", function(spec, prior, n, ...) {
  standardGeneric("buildConstants")
})

#' Shape the observed data into the NIMBLE \code{data} list
#'
#' @param spec A \code{\linkS4class{DistributionSpec}}.
#' @param data The observed data (vector or matrix).
#' @param ... Reserved for methods.
#' @return A named list, typically \code{list(y = ...)}.
#' @keywords internal
setGeneric("buildDataList", function(spec, data, ...) {
  standardGeneric("buildDataList")
})

#' Build dispersed initial values (engine-agnostic)
#'
#' @param spec A \code{\linkS4class{DistributionSpec}}.
#' @param prior The prior list.
#' @param data The observed data.
#' @param count Integer number of component slots (K_max for the DPM, K for
#'   the finite mixture).
#' @param initMethod \code{"kmeans"} (default) or \code{"single"}.
#' @param ... Reserved for methods.
#' @return A list with \code{alloc} (an integer allocation vector) and
#'   \code{params} (a named list of component-parameter initial values).
#' @keywords internal
setGeneric("componentInits", function(spec, prior, data, count, initMethod = "kmeans", ...) {
  standardGeneric("componentInits")
})

#' Parse raw cluster-parameter traces from the MCMC sample matrix
#'
#' @param spec A \code{\linkS4class{DistributionSpec}}.
#' @param samples The MCMC sample matrix (iterations x monitored nodes).
#' @param L Integer truncation length.
#' @param ... Reserved for methods.
#' @return A named list of raw parameter traces (matrices/arrays) keyed by the
#'   logical parameter names of the distribution.
#' @keywords internal
setGeneric("extractParamTraces", function(spec, samples, L, ...) {
  standardGeneric("extractParamTraces")
})

#' Permute cluster parameters and build the relabelled component summary
#'
#' Called by \code{\link{relabel}} after the label-permutation matrix has been
#' derived from the allocation vectors (which is distribution-independent). The
#' spec is responsible only for permuting its own parameters and producing a
#' tidy per-component summary, so multivariate covariance handling stays inside
#' \code{\linkS4class{NormalMvSpec}}.
#'
#' @param spec A \code{\linkS4class{DistributionSpec}}.
#' @param paramTrace The raw parameter trace list (from
#'   \code{extractParamTraces}).
#' @param idx Integer indices of retained (modal-K) iterations.
#' @param occList A list (length = number of retained iterations) of the sorted
#'   occupied cluster labels at each retained iteration.
#' @param perms The permutation matrix (retained iterations x modalK).
#' @param modalK Integer modal number of occupied clusters.
#' @param weights A retained-iterations x modalK matrix of already-permuted
#'   mixing weights.
#' @param ... Reserved for methods.
#' @return A named list with at least \code{summary} (a data.frame) plus
#'   permuted parameter arrays cached on the \code{FitResult}.
#' @keywords internal
setGeneric("relabelComponents",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...) {
    standardGeneric("relabelComponents")
  })

#' Customise MCMC samplers for a component spec
#'
#' Hook called by the engine after \code{configureMCMC()} and before
#' \code{buildMCMC()}, letting a spec swap NIMBLE's default sampler on its own
#' nodes. The default is a no-op; the scale-mixture specs override it to put a
#' slice sampler on the latent precision multipliers, which mixes the partition
#' markedly better than the default random walk.
#'
#' @param spec A \code{\linkS4class{DistributionSpec}}.
#' @param conf An MCMC configuration object.
#' @param model The (uncompiled) NIMBLE model.
#' @param ... Unused.
#' @return The (possibly modified) \code{conf}, invisibly.
#' @export
setGeneric("customizeSamplers", function(spec, conf, model, ...) {
  standardGeneric("customizeSamplers")
})

#' @describeIn customizeSamplers Default: leave NIMBLE's samplers unchanged.
#' @export
setMethod("customizeSamplers", "DistributionSpec",
  function(spec, conf, model, ...) invisible(conf))

#' Is this a regression component spec?
#'
#' Predicate the predict path uses to route to the regression branch. Default
#' \code{FALSE}; regression specs override it to \code{TRUE}.
#' @param spec A \code{\linkS4class{DistributionSpec}}.
#' @param ... Unused.
#' @return Logical scalar.
#' @export
setGeneric("isRegressionSpec", function(spec, ...) {
  standardGeneric("isRegressionSpec")
})

#' @describeIn isRegressionSpec Default: not a regression spec.
#' @export
setMethod("isRegressionSpec", "DistributionSpec", function(spec, ...) FALSE)

#' Inverse link for a regression component
#'
#' Maps the linear predictor to the response mean. Default is the identity
#' link (Normal-linear); GLM specs override it (log, logit).
#' @param spec A \code{\linkS4class{DistributionSpec}}.
#' @param eta Linear predictor value(s).
#' @param prior Optional prior list (for e.g. the Binomial \code{size}).
#' @param ... Unused.
#' @return The response mean.
#' @export
setGeneric("linkInv", function(spec, eta, prior = NULL, ...) {
  standardGeneric("linkInv")
})

#' @describeIn linkInv Identity link (Normal-linear).
#' @export
setMethod("linkInv", "DistributionSpec",
  function(spec, eta, prior = NULL, ...) eta)
