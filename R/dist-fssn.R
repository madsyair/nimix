## dist-fssn.R -------------------------------------------------------------------
## FSSN: Fernandez-Steel skew Normal. Skewness alpha > 0 (alpha = 1 Normal);
## alpha and 1/alpha are treated symmetrically by a log-normal prior on alpha.
## Kernel dFSSN_k built/registered in globalenv by .nimixEnsureMSNBurr().

#' @include class-DistributionSpec.R
#' @include dist-fssn-core.R
NULL

#' FSSN mixture components (Fernandez-Steel skew Normal)
#'
#' Fernandez-Steel skew-Normal component with location \code{mu}, scale
#' \code{sigma}, and skewness \code{alpha} (\code{alpha = 1} Normal). A
#' log-normal prior on \code{alpha} treats left/right skew symmetrically.
#' Non-conjugate.
#'
#' @references
#' Fernandez, C. & Steel, M. F. J. (1998). On Bayesian modeling of fat tails
#' and skewness. JASA 93, 359--371.
#' Choir, A. S. (2020). The New Neo-Normal Distributions and their Properties.
#' @keywords internal
#' @rdname FSSNUvSpec-class
#' @export
setClass("FSSNUvSpec", contains = "DistributionSpec",
  prototype = prototype(name = "fssn",
                        paramNames = c("mu", "sigma", "alpha"), dataDim = 1L))

#' @rdname FSSNUvSpec-class
#' @export
FSSNUvSpec <- function() methods::new("FSSNUvSpec")

.fsSkewPrior <- function(data, control = list()) {
  y <- as.numeric(data); sy <- stats::sd(y)
  if (!is.finite(sy) || sy <= 0) sy <- 1
  cLoc <- if (!is.null(control$cLoc)) control$cLoc else 2
  aSig <- if (!is.null(control$aSig)) control$aSig else 3
  bSig <- if (!is.null(control$bSig)) control$bSig else (aSig - 1) * sy
  aScale <- if (!is.null(control$alphaLogSd)) control$alphaLogSd else 1
  list(mu0 = mean(y), muSd = cLoc * sy, aSig = aSig, bSig = bSig,
       aScale = aScale)
}

#' @describeIn defaultPrior Data-scaled FSSN prior (log-normal skewness).
setMethod("defaultPrior", "FSSNUvSpec",
  function(spec, data, control = list(), ...) .fsSkewPrior(data, control))

#' @describeIn validateParams FSSN hyperparameter checks.
setMethod("validateParams", "FSSNUvSpec", function(spec, params, ...) {
  stopifnot(is.finite(params$mu0), params$muSd > 0, params$aSig > 1,
            params$bSig > 0, params$aScale > 0)
  invisible(TRUE)
})

#' @describeIn simulateParams Draw FSSN component parameters from the prior.
setMethod("simulateParams", "FSSNUvSpec", function(spec, prior, K, ...) {
  list(mu = stats::rnorm(K, prior$mu0, prior$muSd),
       sigma = 1 / stats::rgamma(K, prior$aSig, rate = prior$bSig),
       alpha = stats::rlnorm(K, 0, prior$aScale))
})

#' @describeIn componentDensity FSSN density closure (stable reference form).
setMethod("componentDensity", "FSSNUvSpec", function(spec, ...) {
  function(x, params) dfssn(x, params[["mu"]], params[["sigma"]],
                            params[["alpha"]])
})

.fssnFixedKCode <- function() str2lang("{
  for (i in 1:n) {
    z[i] ~ dcat(weights[1:K])
    y[i] ~ dFSSN_k(muTilde[z[i]], sigmaTilde[z[i]], alphaTilde[z[i]])
  }
  weights[1:K] ~ ddirch(alphaVec[1:K])
  for (j in 1:K) {
    muTilde[j]    ~ dnorm(mu0, sd = muSd)
    sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
    alphaTilde[j] ~ dlnorm(0, sd = aScale)
  }
}")

.fssnDPMCode <- function() str2lang("{
  for (i in 1:n) {
    y[i] ~ dFSSN_k(muTilde[xi[i]], sigmaTilde[xi[i]], alphaTilde[xi[i]])
  }
  xi[1:n] ~ dCRP(alpha, size = n)
  for (j in 1:L) {
    muTilde[j]    ~ dnorm(mu0, sd = muSd)
    sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
    alphaTilde[j] ~ dlnorm(0, sd = aScale)
  }
  alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
}")

.skewMc <- function(code, alloc) list(
  code = code,
  monitors = c(alloc, "muTilde", "sigmaTilde", "alphaTilde",
               if (alloc == "z") "weights" else "alpha"),
  paramNodes = c(mu = "muTilde", sigma = "sigmaTilde", alpha = "alphaTilde"),
  allocNode = alloc)

#' @describeIn buildModelCode FSSN finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("FSSNUvSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) .skewMc(.fssnFixedKCode(), "z"))

#' @describeIn buildModelCode FSSN DPM mixture.
#' @export
setMethod("buildModelCode", signature("FSSNUvSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) .skewMc(.fssnDPMCode(), "xi"))

#' @describeIn buildConstants FSSN constants.
setMethod("buildConstants", "FSSNUvSpec", function(spec, prior, n, ...)
  list(n = n, mu0 = prior$mu0, muSd = prior$muSd, aSig = prior$aSig,
       bSig = prior$bSig, aScale = prior$aScale))

#' @describeIn buildDataList Univariate response.
setMethod("buildDataList", "FSSNUvSpec",
  function(spec, data, ...) list(y = as.numeric(data)))

#' @describeIn componentInits Dispersed k-means start for FSSN (alpha at 1).
setMethod("componentInits", "FSSNUvSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...)
    .skewComponentInits(prior, data, count, initMethod, "alphaTilde", 1, ...))

#' @describeIn extractParamTraces Parse mu / sigma / alpha traces.
setMethod("extractParamTraces", "FSSNUvSpec", function(spec, samples, L, ...)
  list(mu = .nodeToArray(samples, "muTilde", L),
       sigma = .nodeToArray(samples, "sigmaTilde", L),
       alpha = .nodeToArray(samples, "alphaTilde", L)))

#' @describeIn relabelComponents Permute mu / sigma / alpha, summarise.
setMethod("relabelComponents", "FSSNUvSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...)
    .skew3Relabel(paramTrace, idx, occList, perms, modalK, weights, "alpha"))

## shared helpers for 3-parameter skew families (FSSN; reused by FSST) ----------

.skewComponentInits <- function(prior, data, count, initMethod, shapeNode,
                                shapeInit, ...) {
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
  params[[shapeNode]] <- rep(shapeInit, count)
  list(alloc = alloc, params = params)
}

.skew3Relabel <- function(paramTrace, idx, occList, perms, modalK, weights,
                          shapeName) {
  muTr <- paramTrace$mu[idx, , drop = FALSE]
  sgTr <- paramTrace$sigma[idx, , drop = FALSE]
  shTr <- paramTrace[[shapeName]][idx, , drop = FALSE]
  m <- length(idx)
  muRe <- matrix(NA_real_, m, modalK); sgRe <- matrix(NA_real_, m, modalK)
  shRe <- matrix(NA_real_, m, modalK)
  for (r in seq_len(m)) {
    occ <- occList[[r]]; pr <- perms[r, ]
    muRe[r, ] <- muTr[r, occ][pr]; sgRe[r, ] <- sgTr[r, occ][pr]
    shRe[r, ] <- shTr[r, occ][pr]
  }
  q <- function(M, p) apply(M, 2L, stats::quantile, probs = p, names = FALSE)
  summ <- data.frame(
    component = seq_len(modalK), weight = colMeans(weights),
    mu_mean = colMeans(muRe), mu_lwr = q(muRe, 0.025), mu_upr = q(muRe, 0.975),
    sigma_mean = colMeans(sgRe), sigma_lwr = q(sgRe, 0.025),
    sigma_upr = q(sgRe, 0.975))
  summ[[paste0(shapeName, "_mean")]] <- colMeans(shRe)
  summ[[paste0(shapeName, "_lwr")]]  <- q(shRe, 0.025)
  summ[[paste0(shapeName, "_upr")]]  <- q(shRe, 0.975)
  out <- list(mu = muRe, sigma = sgRe, summary = summ)
  out[[shapeName]] <- shRe
  out
}
