## dist-sep.R --------------------------------------------------------------------
## SEP: symmetric exponential power (generalised error) component. Shape nu
## controls tail weight: nu = 2 is Normal, nu = 1 is Laplace, nu -> inf is
## uniform-like. Non-conjugate; the compiled NIMBLE density dSEP_k is built and
## registered in globalenv by .nimixEnsureMSNBurr() (registerDistribution.R).

#' @include class-DistributionSpec.R
#' @include dist-sep-core.R
NULL

#' SEP mixture components (symmetric exponential power)
#'
#' Symmetric exponential-power (generalised error) component family with
#' location \code{mu}, scale \code{sigma}, and shape \code{nu} governing tail
#' weight (\code{nu = 2} Normal, \code{nu = 1} Laplace). Non-conjugate; NIMBLE
#' assigns adaptive samplers to the component parameters.
#'
#' @references
#' Choir, A. S. (2020). The New Neo-Normal Distributions and their Properties.
#' Dissertation. Institut Teknologi Sepuluh Nopember.
#' @keywords internal
#' @rdname SEPUvSpec-class
#' @export
setClass("SEPUvSpec", contains = "DistributionSpec",
  prototype = prototype(name = "sep",
                        paramNames = c("mu", "sigma", "nu"), dataDim = 1L))

#' @rdname SEPUvSpec-class
#' @export
SEPUvSpec <- function() methods::new("SEPUvSpec")

.epDefaultPrior <- function(data, control = list()) {
  y <- as.numeric(data); sy <- stats::sd(y)
  if (!is.finite(sy) || sy <= 0) sy <- 1
  cLoc <- if (!is.null(control$cLoc)) control$cLoc else 2
  aSig <- if (!is.null(control$aSig)) control$aSig else 3
  bSig <- if (!is.null(control$bSig)) control$bSig else (aSig - 1) * sy
  aNu <- if (!is.null(control$aNuShape)) control$aNuShape else 4
  bNu <- if (!is.null(control$bNuRate))  control$bNuRate  else 2   # mode = 1.5
  list(mu0 = mean(y), muSd = cLoc * sy, aSig = aSig, bSig = bSig,
       aNu = aNu, bNu = bNu)
}

#' @describeIn defaultPrior Data-scaled SEP prior.
setMethod("defaultPrior", "SEPUvSpec",
  function(spec, data, control = list(), ...) .epDefaultPrior(data, control))

#' @describeIn validateParams SEP hyperparameter checks.
setMethod("validateParams", "SEPUvSpec", function(spec, params, ...) {
  stopifnot(is.finite(params$mu0), params$muSd > 0, params$aSig > 1,
            params$bSig > 0, params$aNu > 0, params$bNu > 0)
  invisible(TRUE)
})

#' @describeIn simulateParams Draw SEP component parameters from the prior.
setMethod("simulateParams", "SEPUvSpec", function(spec, prior, K, ...) {
  list(mu = stats::rnorm(K, prior$mu0, prior$muSd),
       sigma = 1 / stats::rgamma(K, prior$aSig, rate = prior$bSig),
       nu = stats::rgamma(K, prior$aNu, rate = prior$bNu))
})

#' @describeIn componentDensity SEP density closure (stable reference form).
setMethod("componentDensity", "SEPUvSpec", function(spec, ...) {
  function(x, params) dsep(x, params[["mu"]], params[["sigma"]], params[["nu"]])
})

.epFixedKCode <- function(dfun) str2lang(sprintf("{
  for (i in 1:n) {
    z[i] ~ dcat(weights[1:K])
    y[i] ~ %s(muTilde[z[i]], sigmaTilde[z[i]], nuTilde[z[i]])
  }
  weights[1:K] ~ ddirch(alphaVec[1:K])
  for (j in 1:K) {
    muTilde[j]    ~ dnorm(mu0, sd = muSd)
    sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
    nuTilde[j]    ~ dgamma(shape = aNu, rate = bNu)
  }
}", dfun))

.epDPMCode <- function(dfun) str2lang(sprintf("{
  for (i in 1:n) {
    y[i] ~ %s(muTilde[xi[i]], sigmaTilde[xi[i]], nuTilde[xi[i]])
  }
  xi[1:n] ~ dCRP(alpha, size = n)
  for (j in 1:L) {
    muTilde[j]    ~ dnorm(mu0, sd = muSd)
    sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
    nuTilde[j]    ~ dgamma(shape = aNu, rate = bNu)
  }
  alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
}", dfun))

.epMc <- function(code, alloc) list(
  code = code,
  monitors = c(alloc, "muTilde", "sigmaTilde", "nuTilde",
               if (alloc == "z") "weights" else "alpha"),
  paramNodes = c(mu = "muTilde", sigma = "sigmaTilde", nu = "nuTilde"),
  allocNode = alloc)

#' @describeIn buildModelCode SEP finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("SEPUvSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) .epMc(.epFixedKCode("dSEP_k"), "z"))

#' @describeIn buildModelCode SEP DPM mixture.
#' @export
setMethod("buildModelCode", signature("SEPUvSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) .epMc(.epDPMCode("dSEP_k"), "xi"))

#' @describeIn buildConstants SEP constants.
setMethod("buildConstants", "SEPUvSpec", function(spec, prior, n, ...)
  list(n = n, mu0 = prior$mu0, muSd = prior$muSd, aSig = prior$aSig,
       bSig = prior$bSig, aNu = prior$aNu, bNu = prior$bNu))

#' @describeIn buildDataList Univariate response.
setMethod("buildDataList", "SEPUvSpec",
  function(spec, data, ...) list(y = as.numeric(data)))

#' @describeIn componentInits Dispersed k-means start for SEP.
setMethod("componentInits", "SEPUvSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...)
    .epComponentInits(prior, data, count, initMethod, ...))

#' @describeIn extractParamTraces Parse mu / sigma / nu traces.
setMethod("extractParamTraces", "SEPUvSpec", function(spec, samples, L, ...)
  list(mu = .nodeToArray(samples, "muTilde", L),
       sigma = .nodeToArray(samples, "sigmaTilde", L),
       nu = .nodeToArray(samples, "nuTilde", L)))

#' @describeIn relabelComponents Permute mu / sigma / nu, summarise.
setMethod("relabelComponents", "SEPUvSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...)
    .epRelabel(paramTrace, idx, occList, perms, modalK, weights))

## shared helpers for the two symmetric EP families (SEP, LEP) ------------------

.epComponentInits <- function(prior, data, count, initMethod = "kmeans", ...) {
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
  nuInit <- rep(2, count)                       # start at Normal
  occ <- sort(unique(alloc))
  for (idx in seq_along(occ)) {
    j <- occ[idx]
    if (length(centers) >= idx) muInit[j] <- centers[idx]
    if (length(sds) >= idx && is.finite(sds[idx]) && sds[idx] > 0)
      sigInit[j] <- sds[idx]
  }
  list(alloc = alloc,
       params = list(muTilde = muInit, sigmaTilde = sigInit, nuTilde = nuInit))
}

.epRelabel <- function(paramTrace, idx, occList, perms, modalK, weights) {
  muTr <- paramTrace$mu[idx, , drop = FALSE]
  sgTr <- paramTrace$sigma[idx, , drop = FALSE]
  nuTr <- paramTrace$nu[idx, , drop = FALSE]
  m <- length(idx)
  muRe <- matrix(NA_real_, m, modalK); sgRe <- matrix(NA_real_, m, modalK)
  nuRe <- matrix(NA_real_, m, modalK)
  for (r in seq_len(m)) {
    occ <- occList[[r]]; pr <- perms[r, ]
    muRe[r, ] <- muTr[r, occ][pr]; sgRe[r, ] <- sgTr[r, occ][pr]
    nuRe[r, ] <- nuTr[r, occ][pr]
  }
  q <- function(M, p) apply(M, 2L, stats::quantile, probs = p, names = FALSE)
  summ <- data.frame(
    component = seq_len(modalK), weight = colMeans(weights),
    mu_mean = colMeans(muRe), mu_lwr = q(muRe, 0.025), mu_upr = q(muRe, 0.975),
    sigma_mean = colMeans(sgRe), sigma_lwr = q(sgRe, 0.025),
    sigma_upr = q(sgRe, 0.975),
    nu_mean = colMeans(nuRe), nu_lwr = q(nuRe, 0.025), nu_upr = q(nuRe, 0.975))
  list(mu = muRe, sigma = sgRe, nu = nuRe, summary = summ)
}
