## dist-lep.R --------------------------------------------------------------------
## LEP: exponential power under an alternative (rate-nu) parameterisation. Shares
## the symmetric-EP machinery with SEP (dist-sep.R); differs only in the kernel
## dLEP_k, built and registered in globalenv by .nimixEnsureMSNBurr().

#' @include class-DistributionSpec.R
#' @include dist-sep.R
#' @include dist-lep-core.R
NULL

#' LEP mixture components (exponential power, alternative parameterisation)
#'
#' Symmetric exponential-power component with location \code{mu}, scale
#' \code{sigma}, and shape \code{nu}; an alternative parameterisation to
#' \code{\link{SEPUvSpec}}. Non-conjugate.
#'
#' @references
#' Choir, A. S. (2020). The New Neo-Normal Distributions and their Properties.
#' Dissertation. Institut Teknologi Sepuluh Nopember.
#' @keywords internal
#' @rdname LEPUvSpec-class
#' @export
setClass("LEPUvSpec", contains = "DistributionSpec",
  prototype = prototype(name = "lep",
                        paramNames = c("mu", "sigma", "nu"), dataDim = 1L))

#' @rdname LEPUvSpec-class
#' @export
LEPUvSpec <- function() methods::new("LEPUvSpec")

#' @describeIn defaultPrior Data-scaled LEP prior.
setMethod("defaultPrior", "LEPUvSpec",
  function(spec, data, control = list(), ...) .epDefaultPrior(data, control))

#' @describeIn validateParams LEP hyperparameter checks.
setMethod("validateParams", "LEPUvSpec", function(spec, params, ...) {
  stopifnot(is.finite(params$mu0), params$muSd > 0, params$aSig > 1,
            params$bSig > 0, params$aNu > 0, params$bNu > 0)
  invisible(TRUE)
})

#' @describeIn simulateParams Draw LEP component parameters from the prior.
setMethod("simulateParams", "LEPUvSpec", function(spec, prior, nClust, ...) {
  list(mu = stats::rnorm(nClust, prior$mu0, prior$muSd),
       sigma = 1 / stats::rgamma(nClust, prior$aSig, rate = prior$bSig),
       nu = stats::rgamma(nClust, prior$aNu, rate = prior$bNu))
})

#' @describeIn componentDensity LEP density closure (stable reference form).
setMethod("componentDensity", "LEPUvSpec", function(spec, ...) {
  function(x, params) dlep(x, params[["mu"]], params[["sigma"]], params[["nu"]])
})

#' @describeIn buildModelCode LEP finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("LEPUvSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) .epMc(.epFixedKCode("dLEP_k"), "z"))

#' @describeIn buildModelCode LEP DPM mixture.
#' @export
setMethod("buildModelCode", signature("LEPUvSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) .epMc(.epDPMCode("dLEP_k"), "xi"))

#' @describeIn buildConstants LEP constants.
setMethod("buildConstants", "LEPUvSpec", function(spec, prior, n, ...)
  list(n = n, mu0 = prior$mu0, muSd = prior$muSd, aSig = prior$aSig,
       bSig = prior$bSig, aNu = prior$aNu, bNu = prior$bNu))

#' @describeIn buildDataList Univariate response.
setMethod("buildDataList", "LEPUvSpec",
  function(spec, data, ...) list(y = as.numeric(data)))

#' @describeIn componentInits Dispersed k-means start for LEP.
setMethod("componentInits", "LEPUvSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...)
    .epComponentInits(prior, data, count, initMethod, ...))

#' @describeIn extractParamTraces Parse mu / sigma / nu traces.
setMethod("extractParamTraces", "LEPUvSpec", function(spec, samples, L, ...)
  list(mu = .nodeToArray(samples, "muTilde", L),
       sigma = .nodeToArray(samples, "sigmaTilde", L),
       nu = .nodeToArray(samples, "nuTilde", L)))

#' @describeIn relabelComponents Permute mu / sigma / nu, summarise.
setMethod("relabelComponents", "LEPUvSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...)
    .epRelabel(paramTrace, idx, occList, perms, modalK, weights))
