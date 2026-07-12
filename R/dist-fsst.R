## dist-fsst.R -------------------------------------------------------------------
## FSST: Fernandez-Steel skew Student-t. Skewness alpha (alpha = 1 symmetric-t)
## and degrees of freedom nu. The t-kernel is inlined in dFSST_k (NIMBLE's dt is
## avoided), built/registered in globalenv by .nimixEnsureMSNBurr().
##
## nuTilde is a *stochastic* node truncated below at 2 (so the variance exists).
## It must not be written as a deterministic transform of another node: dCRP
## refuses to cluster deterministic nodes, which breaks the DPM engine.

#' @include class-DistributionSpec.R
#' @include dist-fossep.R
#' @include dist-fsst-core.R
NULL

#' FSST mixture components (Fernandez-Steel skew Student-t)
#'
#' Four-parameter component with location \code{mu}, scale \code{sigma},
#' skewness \code{alpha} (\code{alpha = 1} symmetric) and degrees of freedom
#' \code{nu}, truncated below at 2 so the component variance exists. Heavy tails
#' for small \code{nu}. Non-conjugate. Note the skew-\eqn{t} pitfalls of
#' Fernandez & Steel (1999): very small \code{nu} weakens identifiability, and
#' \code{nu} is typically only weakly identified by the data.
#'
#' @references
#' Fernandez, C. & Steel, M. F. J. (1999). Multivariate Student-t regression
#' models: Pitfalls and inference. Biometrika 86, 153--167.
#' Choir, A. S. (2020). The New Neo-Normal Distributions and their Properties.
#' @keywords internal
#' @rdname FSSTUvSpec-class
#' @export
setClass("FSSTUvSpec", contains = "DistributionSpec",
  prototype = prototype(name = "fsst",
                        paramNames = c("mu", "sigma", "alpha", "nu"),
                        dataDim = 1L))

#' @rdname FSSTUvSpec-class
#' @export
FSSTUvSpec <- function() methods::new("FSSTUvSpec")

.fsstDefaultPrior <- function(data, control = list()) {
  p <- .fsSkewPrior(data, control)               # mu, sigma, alpha (log-normal)
  p$aNu <- if (!is.null(control$aNuShape)) control$aNuShape else 2
  p$bNu <- if (!is.null(control$bNuRate))  control$bNuRate  else 0.15
  p
}

#' @describeIn defaultPrior Data-scaled FSST prior.
setMethod("defaultPrior", "FSSTUvSpec",
  function(spec, data, control = list(), ...) .fsstDefaultPrior(data, control))

#' @describeIn validateParams FSST hyperparameter checks.
setMethod("validateParams", "FSSTUvSpec", function(spec, params, ...) {
  stopifnot(is.finite(params$mu0), params$muSd > 0, params$aSig > 1,
            params$bSig > 0, params$aScale > 0, params$aNu > 0, params$bNu > 0)
  invisible(TRUE)
})

#' @describeIn simulateParams Draw FSST component parameters from the prior.
setMethod("simulateParams", "FSSTUvSpec", function(spec, prior, K, ...) {
  list(mu = stats::rnorm(K, prior$mu0, prior$muSd),
       sigma = 1 / stats::rgamma(K, prior$aSig, rate = prior$bSig),
       alpha = stats::rlnorm(K, 0, prior$aScale),
       nu = 2 + stats::rgamma(K, prior$aNu, rate = prior$bNu))
})

#' @describeIn componentDensity FSST density closure (stable reference form).
setMethod("componentDensity", "FSSTUvSpec", function(spec, ...) {
  function(x, params) dfsst(x, params[["mu"]], params[["sigma"]],
                            params[["alpha"]], params[["nu"]])
})

.fsstFixedKCode <- function() str2lang("{
  for (i in 1:n) {
    z[i] ~ dcat(weights[1:K])
    y[i] ~ dFSST_k(muTilde[z[i]], sigmaTilde[z[i]], alphaTilde[z[i]],
                   nuTilde[z[i]])
  }
  weights[1:K] ~ ddirch(alphaVec[1:K])
  for (j in 1:K) {
    muTilde[j]    ~ dnorm(mu0, sd = muSd)
    sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
    alphaTilde[j] ~ dlnorm(0, sd = aScale)
    nuTilde[j]    ~ T(dgamma(shape = aNu, rate = bNu), 2, )
  }
}")

.fsstDPMCode <- function() str2lang("{
  for (i in 1:n) {
    y[i] ~ dFSST_k(muTilde[xi[i]], sigmaTilde[xi[i]], alphaTilde[xi[i]],
                   nuTilde[xi[i]])
  }
  xi[1:n] ~ dCRP(alpha, size = n)
  for (j in 1:L) {
    muTilde[j]    ~ dnorm(mu0, sd = muSd)
    sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
    alphaTilde[j] ~ dlnorm(0, sd = aScale)
    nuTilde[j]    ~ T(dgamma(shape = aNu, rate = bNu), 2, )
  }
  alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
}")

#' @describeIn buildModelCode FSST finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("FSSTUvSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...)
    .skew4Mc(.fsstFixedKCode(), "z", "alphaTilde", "nuTilde"))

#' @describeIn buildModelCode FSST DPM mixture.
#' @export
setMethod("buildModelCode", signature("FSSTUvSpec", "DPMEngine"),
  function(spec, engine, n, L, ...)
    .skew4Mc(.fsstDPMCode(), "xi", "alphaTilde", "nuTilde"))

#' @describeIn buildConstants FSST constants.
setMethod("buildConstants", "FSSTUvSpec", function(spec, prior, n, ...)
  list(n = n, mu0 = prior$mu0, muSd = prior$muSd, aSig = prior$aSig,
       bSig = prior$bSig, aScale = prior$aScale, aNu = prior$aNu,
       bNu = prior$bNu))

#' @describeIn buildDataList Univariate response.
setMethod("buildDataList", "FSSTUvSpec",
  function(spec, data, ...) list(y = as.numeric(data)))

#' @describeIn componentInits Dispersed k-means start for FSST.
setMethod("componentInits", "FSSTUvSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...)
    .skew4ComponentInits(prior, data, count, initMethod,
                         "alphaTilde", 1, "nuTilde", 8, ...))

#' @describeIn extractParamTraces Parse mu / sigma / alpha / nu traces.
setMethod("extractParamTraces", "FSSTUvSpec", function(spec, samples, L, ...)
  list(mu = .nodeToArray(samples, "muTilde", L),
       sigma = .nodeToArray(samples, "sigmaTilde", L),
       alpha = .nodeToArray(samples, "alphaTilde", L),
       nu = .nodeToArray(samples, "nuTilde", L)))

#' @describeIn relabelComponents Permute mu / sigma / alpha / nu, summarise.
setMethod("relabelComponents", "FSSTUvSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...)
    .skew4Relabel(paramTrace, idx, occList, perms, modalK, weights,
                  "alpha", "nu"))
