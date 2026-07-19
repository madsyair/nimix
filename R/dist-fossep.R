## dist-fossep.R -----------------------------------------------------------------
## FOSSEP: Fernandez-Steel skew exponential power. Skewness alpha (alpha = 1
## symmetric) and shape theta (theta = 2 gives the skew-Normal kernel). Kernel
## dFOSSEP_k built/registered in globalenv by .nimixEnsureMSNBurr().

#' @include class-DistributionSpec.R
#' @include dist-fssn.R
#' @include dist-fossep-core.R
NULL

#' FOSSEP mixture components (Fernandez-Steel skew exponential power)
#'
#' Four-parameter neo-normal component with location \code{mu}, scale
#' \code{sigma}, skewness \code{alpha} (\code{alpha = 1} symmetric) and shape
#' \code{theta} (\code{theta = 2} skew-Normal kernel). Non-conjugate.
#'
#' @references
#' Choir, A. S. (2020). The New Neo-Normal Distributions and their Properties.
#' Dissertation. Institut Teknologi Sepuluh Nopember.
#' @keywords internal
#' @rdname FOSSEPUvSpec-class
#' @export
setClass("FOSSEPUvSpec", contains = "DistributionSpec",
  prototype = prototype(name = "fossep",
                        paramNames = c("mu", "sigma", "alpha", "theta"),
                        dataDim = 1L))

#' @rdname FOSSEPUvSpec-class
#' @export
FOSSEPUvSpec <- function() methods::new("FOSSEPUvSpec")

.fossepDefaultPrior <- function(data, control = list()) {
  p <- .fsSkewPrior(data, control)               # mu, sigma, alpha (log-normal)
  p$aTheta <- if (!is.null(control$aThetaShape)) control$aThetaShape else 4
  p$bTheta <- if (!is.null(control$bThetaRate))  control$bThetaRate  else 2
  p
}

#' @describeIn defaultPrior Data-scaled FOSSEP prior.
setMethod("defaultPrior", "FOSSEPUvSpec",
  function(spec, data, control = list(), ...) .fossepDefaultPrior(data, control))

#' @describeIn validateParams FOSSEP hyperparameter checks.
setMethod("validateParams", "FOSSEPUvSpec", function(spec, params, ...) {
  stopifnot(is.finite(params$mu0), params$muSd > 0, params$aSig > 1,
            params$bSig > 0, params$aScale > 0, params$aTheta > 0,
            params$bTheta > 0)
  invisible(TRUE)
})

#' @describeIn simulateParams Draw FOSSEP component parameters from the prior.
setMethod("simulateParams", "FOSSEPUvSpec", function(spec, prior, K, ...) {
  list(mu = stats::rnorm(K, prior$mu0, prior$muSd),
       sigma = 1 / stats::rgamma(K, prior$aSig, rate = prior$bSig),
       alpha = stats::rlnorm(K, 0, prior$aScale),
       theta = stats::rgamma(K, prior$aTheta, rate = prior$bTheta))
})

#' @describeIn componentDensity FOSSEP density closure (stable reference form).
setMethod("componentDensity", "FOSSEPUvSpec", function(spec, ...) {
  function(x, params) dfossep(x, params[["mu"]], params[["sigma"]],
                              params[["alpha"]], params[["theta"]])
})

.fossepFixedKCode <- function() str2lang("{
  for (i in 1:n) {
    z[i] ~ dcat(weights[1:K])
    y[i] ~ dFOSSEP_k(muTilde[z[i]], sigmaTilde[z[i]], alphaTilde[z[i]],
                     thetaTilde[z[i]])
  }
  weights[1:K] ~ ddirch(alphaVec[1:K])
  for (j in 1:K) {
    muTilde[j]    ~ dnorm(mu0, sd = muSd)
    sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
    alphaTilde[j] ~ dlnorm(0, sd = aScale)
    thetaTilde[j] ~ dgamma(shape = aTheta, rate = bTheta)
  }
}")

.fossepDPMCode <- function() str2lang("{
  for (i in 1:n) {
    y[i] ~ dFOSSEP_k(muTilde[xi[i]], sigmaTilde[xi[i]], alphaTilde[xi[i]],
                     thetaTilde[xi[i]])
  }
  xi[1:n] ~ dCRP(alpha, size = n)
  for (j in 1:L) {
    muTilde[j]    ~ dnorm(mu0, sd = muSd)
    sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
    alphaTilde[j] ~ dlnorm(0, sd = aScale)
    thetaTilde[j] ~ dgamma(shape = aTheta, rate = bTheta)
  }
  alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
}")

#' @describeIn buildModelCode FOSSEP finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("FOSSEPUvSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) .skew4Mc(.fossepFixedKCode(), "z"))

#' @describeIn buildModelCode FOSSEP DPM mixture.
#' @export
setMethod("buildModelCode", signature("FOSSEPUvSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) .skew4Mc(.fossepDPMCode(), "xi"))

#' @describeIn buildConstants FOSSEP constants.
setMethod("buildConstants", "FOSSEPUvSpec", function(spec, prior, n, ...)
  list(n = n, mu0 = prior$mu0, muSd = prior$muSd, aSig = prior$aSig,
       bSig = prior$bSig, aScale = prior$aScale, aTheta = prior$aTheta,
       bTheta = prior$bTheta))

#' @describeIn buildDataList Univariate response.
setMethod("buildDataList", "FOSSEPUvSpec",
  function(spec, data, ...) list(y = as.numeric(data)))

#' @describeIn componentInits Dispersed k-means start for FOSSEP.
setMethod("componentInits", "FOSSEPUvSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...)
    .skew4ComponentInits(prior, data, count, initMethod,
                         "alphaTilde", 1, "thetaTilde", 2, ...))

#' @describeIn extractParamTraces Parse mu / sigma / alpha / theta traces.
setMethod("extractParamTraces", "FOSSEPUvSpec", function(spec, samples, L, ...)
  list(mu = .nodeToArray(samples, "muTilde", L),
       sigma = .nodeToArray(samples, "sigmaTilde", L),
       alpha = .nodeToArray(samples, "alphaTilde", L),
       theta = .nodeToArray(samples, "thetaTilde", L)))

#' @describeIn relabelComponents Permute mu / sigma / alpha / theta, summarise.
setMethod("relabelComponents", "FOSSEPUvSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...)
    .skew4Relabel(paramTrace, idx, occList, perms, modalK, weights,
                  "alpha", "theta"))

## shared helpers for 4-parameter skew families (FOSSEP, FSST, JFST) ------------
## The fourth node is named per family (thetaTilde, or nuTilde for FSST), so the
## helpers take the node names explicitly rather than assuming them.

.skew4Mc <- function(code, alloc, s1 = "alphaTilde", s2 = "thetaTilde") list(
  code = code,
  monitors = c(alloc, "muTilde", "sigmaTilde", s1, s2,
               if (alloc == "z") "weights" else "alpha"),
  paramNodes = stats::setNames(c("muTilde", "sigmaTilde", s1, s2),
                               c("mu", "sigma", "alpha", "theta")),
  allocNode = alloc)

.skew4ComponentInits <- function(prior, data, count, initMethod,
                                 s1Node, s1Init, s2Node, s2Init, ...) {
  y <- as.numeric(data); n <- length(y); nUnique <- length(unique(y))
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
  muInit <- rep(prior$mu0, count)
  sigInit <- rep(prior$bSig / (prior$aSig - 1), count)
  occ <- sort(unique(alloc))
  for (idx in seq_along(occ)) {
    j <- occ[idx]
    if (length(centers) >= idx) muInit[j] <- centers[idx]
    if (length(sds) >= idx && is.finite(sds[idx]) && sds[idx] > 0)
      sigInit[j] <- sds[idx]
  }
  params <- list(muTilde = muInit, sigmaTilde = sigInit)
  params[[s1Node]] <- rep(s1Init, count)
  params[[s2Node]] <- rep(s2Init, count)
  list(alloc = alloc, params = params)
}

.skew4Relabel <- function(paramTrace, idx, occList, perms, modalK, weights,
                          n1, n2) {
  muTr <- paramTrace$mu[idx, , drop = FALSE]
  sgTr <- paramTrace$sigma[idx, , drop = FALSE]
  s1Tr <- paramTrace[[n1]][idx, , drop = FALSE]
  s2Tr <- paramTrace[[n2]][idx, , drop = FALSE]
  m <- length(idx)
  muRe <- matrix(NA_real_, m, modalK); sgRe <- matrix(NA_real_, m, modalK)
  s1Re <- matrix(NA_real_, m, modalK); s2Re <- matrix(NA_real_, m, modalK)
  for (r in seq_len(m)) {
    occ <- occList[[r]]; pr <- perms[r, ]
    muRe[r, ] <- muTr[r, occ][pr]; sgRe[r, ] <- sgTr[r, occ][pr]
    s1Re[r, ] <- s1Tr[r, occ][pr]; s2Re[r, ] <- s2Tr[r, occ][pr]
  }
  q <- function(M, p) apply(M, 2L, stats::quantile, probs = p, names = FALSE)
  summ <- data.frame(
    component = seq_len(modalK), weight = colMeans(weights),
    mu_mean = colMeans(muRe), mu_lwr = q(muRe, 0.025), mu_upr = q(muRe, 0.975),
    sigma_mean = colMeans(sgRe), sigma_lwr = q(sgRe, 0.025),
    sigma_upr = q(sgRe, 0.975))
  summ[[paste0(n1, "_mean")]] <- colMeans(s1Re)
  summ[[paste0(n1, "_lwr")]]  <- q(s1Re, 0.025)
  summ[[paste0(n1, "_upr")]]  <- q(s1Re, 0.975)
  summ[[paste0(n2, "_mean")]] <- colMeans(s2Re)
  summ[[paste0(n2, "_lwr")]]  <- q(s2Re, 0.025)
  summ[[paste0(n2, "_upr")]]  <- q(s2Re, 0.975)
  out <- list(mu = muRe, sigma = sgRe, summary = summ)
  out[[n1]] <- s1Re; out[[n2]] <- s2Re
  out
}
