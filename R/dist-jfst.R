## dist-jfst.R -------------------------------------------------------------------
## JFST: Jones-Faddy skew-t. Two shape parameters alpha, theta: alpha = theta is
## symmetric (a scaled t with 2*alpha df), alpha > theta right-skews and
## alpha < theta left-skews. Kernel dJFST_k built/registered in globalenv by
## .nimixEnsureMSNBurr(); it uses the branch-free identity
## sign(z)/sqrt((a+th)/z^2 + 1) == z/sqrt(a + th + z^2), which is also finite
## at z = 0.

#' @include class-DistributionSpec.R
#' @include dist-fossep.R
#' @include dist-jfst-core.R
NULL

#' JFST mixture components (Jones-Faddy skew-t)
#'
#' Four-parameter component with location \code{mu}, scale \code{sigma}, and two
#' shape parameters \code{alpha}, \code{theta}. \code{alpha = theta} is
#' symmetric; \code{alpha > theta} skews right and \code{alpha < theta} skews
#' left. Both control tail weight. Non-conjugate.
#'
#' @references
#' Jones, M. C. & Faddy, M. J. (2003). A skew extension of the t-distribution,
#' with applications. JRSS-B 65, 159--174.
#' Choir, A. S. (2020). The New Neo-Normal Distributions and their Properties.
#' @keywords internal
#' @rdname JFSTUvSpec-class
#' @export
setClass("JFSTUvSpec", contains = "DistributionSpec",
  prototype = prototype(name = "jfst",
                        paramNames = c("mu", "sigma", "alpha", "theta"),
                        dataDim = 1L))

#' @rdname JFSTUvSpec-class
#' @export
JFSTUvSpec <- function() methods::new("JFSTUvSpec")

.jfstDefaultPrior <- function(data, control = list()) {
  y <- as.numeric(data); sy <- stats::sd(y)
  if (!is.finite(sy) || sy <= 0) sy <- 1
  cLoc <- if (!is.null(control$cLoc)) control$cLoc else 2
  aSig <- if (!is.null(control$aSig)) control$aSig else 3
  bSig <- if (!is.null(control$bSig)) control$bSig else (aSig - 1) * sy
  # Both shapes get the same weakly-informative Gamma (mode 3), keeping the
  # prior symmetric in alpha <-> theta, i.e. centred on symmetry.
  aSh <- if (!is.null(control$aShapeShape)) control$aShapeShape else 4
  bSh <- if (!is.null(control$bShapeRate))  control$bShapeRate  else 1
  list(mu0 = mean(y), muSd = cLoc * sy, aSig = aSig, bSig = bSig,
       aSh = aSh, bSh = bSh)
}

#' @describeIn defaultPrior Data-scaled JFST prior (symmetric in alpha/theta).
setMethod("defaultPrior", "JFSTUvSpec",
  function(spec, data, control = list(), ...) .jfstDefaultPrior(data, control))

#' @describeIn validateParams JFST hyperparameter checks.
setMethod("validateParams", "JFSTUvSpec", function(spec, params, ...) {
  stopifnot(is.finite(params$mu0), params$muSd > 0, params$aSig > 1,
            params$bSig > 0, params$aSh > 0, params$bSh > 0)
  invisible(TRUE)
})

#' @describeIn simulateParams Draw JFST component parameters from the prior.
setMethod("simulateParams", "JFSTUvSpec", function(spec, prior, K, ...) {
  list(mu = stats::rnorm(K, prior$mu0, prior$muSd),
       sigma = 1 / stats::rgamma(K, prior$aSig, rate = prior$bSig),
       alpha = stats::rgamma(K, prior$aSh, rate = prior$bSh),
       theta = stats::rgamma(K, prior$aSh, rate = prior$bSh))
})

#' @describeIn componentDensity JFST density closure (stable reference form).
setMethod("componentDensity", "JFSTUvSpec", function(spec, ...) {
  function(x, params) djfst(x, params[["mu"]], params[["sigma"]],
                            params[["alpha"]], params[["theta"]])
})

.jfstFixedKCode <- function() str2lang("{
  for (i in 1:n) {
    z[i] ~ dcat(weights[1:K])
    y[i] ~ dJFST_k(muTilde[z[i]], sigmaTilde[z[i]], alphaTilde[z[i]],
                   thetaTilde[z[i]])
  }
  weights[1:K] ~ ddirch(alphaVec[1:K])
  for (j in 1:K) {
    muTilde[j]    ~ dnorm(mu0, sd = muSd)
    sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
    alphaTilde[j] ~ dgamma(shape = aSh, rate = bSh)
    thetaTilde[j] ~ dgamma(shape = aSh, rate = bSh)
  }
}")

.jfstDPMCode <- function() str2lang("{
  for (i in 1:n) {
    y[i] ~ dJFST_k(muTilde[xi[i]], sigmaTilde[xi[i]], alphaTilde[xi[i]],
                   thetaTilde[xi[i]])
  }
  xi[1:n] ~ dCRP(alpha, size = n)
  for (j in 1:L) {
    muTilde[j]    ~ dnorm(mu0, sd = muSd)
    sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
    alphaTilde[j] ~ dgamma(shape = aSh, rate = bSh)
    thetaTilde[j] ~ dgamma(shape = aSh, rate = bSh)
  }
  alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
}")

#' @describeIn buildModelCode JFST finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("JFSTUvSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) .skew4Mc(.jfstFixedKCode(), "z"))

#' @describeIn buildModelCode JFST DPM mixture.
#' @export
setMethod("buildModelCode", signature("JFSTUvSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) .skew4Mc(.jfstDPMCode(), "xi"))

#' @describeIn buildConstants JFST constants.
setMethod("buildConstants", "JFSTUvSpec", function(spec, prior, n, ...)
  list(n = n, mu0 = prior$mu0, muSd = prior$muSd, aSig = prior$aSig,
       bSig = prior$bSig, aSh = prior$aSh, bSh = prior$bSh))

#' @describeIn buildDataList Univariate response.
setMethod("buildDataList", "JFSTUvSpec",
  function(spec, data, ...) list(y = as.numeric(data)))

#' @describeIn componentInits Dispersed k-means start for JFST (symmetric start).
setMethod("componentInits", "JFSTUvSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...)
    .skew4ComponentInits(prior, data, count, initMethod,
                         "alphaTilde", 3, "thetaTilde", 3, ...))

#' @describeIn extractParamTraces Parse mu / sigma / alpha / theta traces.
setMethod("extractParamTraces", "JFSTUvSpec", function(spec, samples, L, ...)
  list(mu = .nodeToArray(samples, "muTilde", L),
       sigma = .nodeToArray(samples, "sigmaTilde", L),
       alpha = .nodeToArray(samples, "alphaTilde", L),
       theta = .nodeToArray(samples, "thetaTilde", L)))

#' @describeIn relabelComponents Permute mu / sigma / alpha / theta, summarise.
setMethod("relabelComponents", "JFSTUvSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...)
    .skew4Relabel(paramTrace, idx, occList, perms, modalK, weights,
                  "alpha", "theta"))
