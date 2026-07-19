## dist-gmsnburr.R ---------------------------------------------------------------
## GMSNBurr neo-normal component (Iriawan 2000; Choir 2020): two shape
## parameters alpha, theta. theta = 1 recovers MSNBurr, alpha = 1 recovers
## MSNBurr-IIa, alpha = theta -> inf converges to the Normal. The compiled
## NIMBLE density dGMSNBurr_k is built and registered lazily in globalenv by
## .nimixEnsureMSNBurr() (see registerDistribution.R), the fork-safe pattern
## required for scalar user-defined densities.

#' @include class-DistributionSpec.R
#' @include dist-gmsnburr-core.R
NULL

#' GMSNBurr mixture components (generalized neo-normal)
#'
#' Univariate neo-normal component family with location \code{mu}, scale
#' \code{sigma}, and two shape parameters \code{alpha}, \code{theta} governing
#' skewness (Iriawan 2000; Choir 2020). \code{theta = 1} is MSNBurr,
#' \code{alpha = 1} is MSNBurr-IIa, and \code{alpha = theta} is symmetric.
#' Non-conjugate; NIMBLE assigns adaptive samplers to the component parameters.
#'
#' @references
#' Choir, A. S. (2020). The New Neo-Normal Distributions and their Properties.
#' Dissertation. Institut Teknologi Sepuluh Nopember.
#' @keywords internal
#' @export
setClass(
  "GMSNBurrUvSpec",
  contains = "DistributionSpec",
  prototype = prototype(name = "gmsnburr",
                        paramNames = c("mu", "sigma", "alpha", "theta"),
                        dataDim = 1L)
)

#' @rdname GMSNBurrUvSpec-class
#' @export
GMSNBurrUvSpec <- function() methods::new("GMSNBurrUvSpec")

.gmsnburrDefaultPrior <- function(spec, data, control = list()) {
  y <- as.numeric(data)
  sy <- stats::sd(y); if (!is.finite(sy) || sy <= 0) sy <- 1
  cLoc <- if (!is.null(control$cLoc)) control$cLoc else 2
  aSig <- if (!is.null(control$aSig)) control$aSig else 3
  bSig <- if (!is.null(control$bSig)) control$bSig else (aSig - 1) * sy
  aA <- if (!is.null(control$aAlphaShape)) control$aAlphaShape else 2
  bA <- if (!is.null(control$bAlphaRate))  control$bAlphaRate  else 2
  aT <- if (!is.null(control$aThetaShape)) control$aThetaShape else 2
  bT <- if (!is.null(control$bThetaRate))  control$bThetaRate  else 2
  list(mu0 = mean(y), muSd = cLoc * sy,
       aSig = aSig, bSig = bSig, aA = aA, bA = bA, aT = aT, bT = bT,
       cLoc = cLoc)
}

#' @describeIn defaultPrior Data-scaled GMSNBurr prior.
setMethod("defaultPrior", "GMSNBurrUvSpec",
  function(spec, data, control = list(), ...)
    .gmsnburrDefaultPrior(spec, data, control))

#' @describeIn validateParams GMSNBurr hyperparameter checks.
setMethod("validateParams", "GMSNBurrUvSpec",
  function(spec, params, ...) {
    stopifnot(is.finite(params$mu0), params$muSd > 0,
              params$aSig > 1, params$bSig > 0,
              params$aA > 0, params$bA > 0, params$aT > 0, params$bT > 0)
    invisible(TRUE)
  })

#' @describeIn simulateParams Draw GMSNBurr component parameters from the prior.
setMethod("simulateParams", "GMSNBurrUvSpec",
  function(spec, prior, K, ...) {
    list(mu = stats::rnorm(K, prior$mu0, prior$muSd),
         sigma = 1 / stats::rgamma(K, prior$aSig, rate = prior$bSig),
         alpha = stats::rgamma(K, prior$aA, rate = prior$bA),
         theta = stats::rgamma(K, prior$aT, rate = prior$bT))
  })

#' @describeIn componentDensity GMSNBurr density closure (stable reference form).
setMethod("componentDensity", "GMSNBurrUvSpec",
  function(spec, ...) {
    function(x, params) dgmsnburr(x, params[["mu"]], params[["sigma"]],
                                  params[["alpha"]], params[["theta"]])
  })

.gmsnburrFixedKCode <- function() str2lang("{
  for (i in 1:n) {
    z[i] ~ dcat(weights[1:K])
    y[i] ~ dGMSNBurr_k(muTilde[z[i]], sigmaTilde[z[i]],
                       alphaTilde[z[i]], thetaTilde[z[i]])
  }
  weights[1:K] ~ ddirch(alphaVec[1:K])
  for (j in 1:K) {
    muTilde[j]    ~ dnorm(mu0, sd = muSd)
    sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
    alphaTilde[j] ~ dgamma(shape = aA, rate = bA)
    thetaTilde[j] ~ dgamma(shape = aT, rate = bT)
  }
}")

.gmsnburrDPMCode <- function() str2lang("{
  for (i in 1:n) {
    y[i] ~ dGMSNBurr_k(muTilde[xi[i]], sigmaTilde[xi[i]],
                       alphaTilde[xi[i]], thetaTilde[xi[i]])
  }
  xi[1:n] ~ dCRP(alpha, size = n)
  for (j in 1:L) {
    muTilde[j]    ~ dnorm(mu0, sd = muSd)
    sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
    alphaTilde[j] ~ dgamma(shape = aA, rate = bA)
    thetaTilde[j] ~ dgamma(shape = aT, rate = bT)
  }
  alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
}")

.gmsnburrMc <- function(code, alloc) {
  list(code = code,
       monitors = c(alloc, "muTilde", "sigmaTilde", "alphaTilde", "thetaTilde",
                    if (alloc == "z") "weights" else "alpha"),
       paramNodes = c(mu = "muTilde", sigma = "sigmaTilde",
                      alpha = "alphaTilde", theta = "thetaTilde"),
       allocNode = alloc)
}

#' @describeIn buildModelCode GMSNBurr finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("GMSNBurrUvSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) .gmsnburrMc(.gmsnburrFixedKCode(), "z"))

#' @describeIn buildModelCode GMSNBurr DPM mixture.
#' @export
setMethod("buildModelCode", signature("GMSNBurrUvSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) .gmsnburrMc(.gmsnburrDPMCode(), "xi"))

#' @describeIn buildConstants GMSNBurr constants.
setMethod("buildConstants", "GMSNBurrUvSpec",
  function(spec, prior, n, ...)
    list(n = n, mu0 = prior$mu0, muSd = prior$muSd,
         aSig = prior$aSig, bSig = prior$bSig, aA = prior$aA, bA = prior$bA,
         aT = prior$aT, bT = prior$bT))

#' @describeIn buildDataList Univariate response.
setMethod("buildDataList", "GMSNBurrUvSpec",
  function(spec, data, ...) list(y = as.numeric(data)))

#' @describeIn componentInits Dispersed k-means start for GMSNBurr.
setMethod("componentInits", "GMSNBurrUvSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    y <- as.numeric(data); n <- length(y)
    nUnique <- length(unique(y))
    initRatio <- .initRatioArg(...)
    k0 <- max(1L, min(as.integer(floor(initRatio * count)),
                      as.integer(ceiling(sqrt(n)))))
    k0 <- min(k0, max(1L, nUnique))
    alloc <- rep(1L, n); centers <- mean(y)
    sds <- stats::sd(y); if (!is.finite(sds) || sds <= 0) sds <- 1
    if (!identical(initMethod, "single") && k0 >= 2L && nUnique >= k0) {
      cl <- .initClusters(y, k0, initMethod)
      if (!is.null(cl)) {
        alloc <- cl; centers <- vapply(sort(unique(cl)), function(j) mean(y[cl == j]), numeric(1))
        sds <- vapply(seq_len(k0), function(j) {
          s <- stats::sd(y[alloc == j]); if (!is.finite(s) || s <= 0) 1 else s
        }, numeric(1))
      }
    }
    muInit  <- rep(prior$mu0, count)
    sigInit <- rep(prior$bSig / (prior$aSig - 1), count)
    alInit  <- rep(1, count); thInit <- rep(1, count)   # symmetric start
    occ <- sort(unique(alloc))
    for (idx in seq_along(occ)) {
      j <- occ[idx]
      if (length(centers) >= idx) muInit[j] <- centers[idx]
      if (length(sds) >= idx && is.finite(sds[idx]) && sds[idx] > 0)
        sigInit[j] <- sds[idx]
    }
    list(alloc = alloc,
         params = list(muTilde = muInit, sigmaTilde = sigInit,
                       alphaTilde = alInit, thetaTilde = thInit))
  })

#' @describeIn extractParamTraces Parse mu / sigma / alpha / theta traces.
setMethod("extractParamTraces", "GMSNBurrUvSpec",
  function(spec, samples, L, ...) {
    list(mu    = .nodeToArray(samples, "muTilde",    L),
         sigma = .nodeToArray(samples, "sigmaTilde", L),
         alpha = .nodeToArray(samples, "alphaTilde", L),
         theta = .nodeToArray(samples, "thetaTilde", L))
  })

#' @describeIn relabelComponents Permute mu / sigma / alpha / theta, summarise.
setMethod("relabelComponents", "GMSNBurrUvSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...) {
    muTr <- paramTrace$mu[idx, , drop = FALSE]
    sgTr <- paramTrace$sigma[idx, , drop = FALSE]
    alTr <- paramTrace$alpha[idx, , drop = FALSE]
    thTr <- paramTrace$theta[idx, , drop = FALSE]
    m <- length(idx)
    muRe <- matrix(NA_real_, m, modalK); sgRe <- matrix(NA_real_, m, modalK)
    alRe <- matrix(NA_real_, m, modalK); thRe <- matrix(NA_real_, m, modalK)
    for (r in seq_len(m)) {
      occ <- occList[[r]]; pr <- perms[r, ]
      muRe[r, ] <- muTr[r, occ][pr]; sgRe[r, ] <- sgTr[r, occ][pr]
      alRe[r, ] <- alTr[r, occ][pr]; thRe[r, ] <- thTr[r, occ][pr]
    }
    q <- function(M, p) apply(M, 2L, stats::quantile, probs = p, names = FALSE)
    summ <- data.frame(
      component = seq_len(modalK), weight = colMeans(weights),
      mu_mean = colMeans(muRe), mu_lwr = q(muRe, 0.025), mu_upr = q(muRe, 0.975),
      sigma_mean = colMeans(sgRe), sigma_lwr = q(sgRe, 0.025),
      sigma_upr = q(sgRe, 0.975),
      alpha_mean = colMeans(alRe), alpha_lwr = q(alRe, 0.025),
      alpha_upr = q(alRe, 0.975),
      theta_mean = colMeans(thRe), theta_lwr = q(thRe, 0.025),
      theta_upr = q(thRe, 0.975))
    list(mu = muRe, sigma = sgRe, alpha = alRe, theta = thRe, summary = summ)
  })
