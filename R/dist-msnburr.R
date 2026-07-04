#' @include class-DistributionSpec.R
#' @include dist-msnburr-core.R
NULL

## dist-msnburr.R ----------------------------------------------------------------
## MSNBurr and MSNBurr-IIa neo-normal mixture components (Iriawan 2000;
## Choir 2020). alpha controls skewness: MSNBurr accommodates left skew,
## MSNBurr-IIa mirrors it for right skew; alpha = 1 is the logistic
## distribution (exact identity, pinned in the tests). The NIMBLE densities
## below reproduce the maintainer-contributed stable reference implementation
## (dist-msnburr-core.R) exactly: asymptotic log-omega branch for alpha -> 0,
## Maechler-style log1pexp thresholds, and tail-stable algebraic arrangements
## -- so log-densities remain finite for |z| in the hundreds.

# --- S4 specs -------------------------------------------------------------------

#' MSNBurr mixture components (neo-normal, left-skew capable)
#'
#' Univariate neo-normal component family with location \code{mu}, scale
#' \code{sigma} and shape \code{alpha} controlling skewness (Iriawan 2000;
#' Choir 2020). \code{alpha = 1} is the logistic distribution. Non-conjugate;
#' NIMBLE assigns adaptive samplers to the component parameters.
#'
#' @references
#' Iriawan, N. (2000). Computationally Intensive Approaches to Inference in
#' Neo-Normal Linear Models. PhD Thesis, Curtin University of Technology.
#'
#' Choir, A. S. (2020). The New Neo-Normal Distributions and their
#' Properties. Dissertation. Institut Teknologi Sepuluh Nopember.
#' @keywords internal
#' @export
setClass(
  "MSNBurrUvSpec",
  contains = "DistributionSpec",
  prototype = prototype(name = "msnburr",
                        paramNames = c("mu", "sigma", "alpha"), dataDim = 1L)
)

#' @rdname MSNBurrUvSpec-class
#' @export
MSNBurrUvSpec <- function() methods::new("MSNBurrUvSpec")

#' MSNBurr-IIa mixture components (neo-normal, right-skew capable)
#'
#' Mirror image of \code{\link{MSNBurrUvSpec}}: the MSNBurr-IIa family
#' (Iriawan 2000; Choir 2020) accommodates right skew, with \code{alpha = 1}
#' again the logistic distribution.
#'
#' @references
#' Choir, A. S. (2020). The New Neo-Normal Distributions and their Properties.
#' Dissertation. Institut Teknologi Sepuluh Nopember.
#' @keywords internal
#' @rdname MSNBurr2aUvSpec-class
#' @export
setClass(
  "MSNBurr2aUvSpec",
  contains = "DistributionSpec",
  prototype = prototype(name = "msnburr2a",
                        paramNames = c("mu", "sigma", "alpha"), dataDim = 1L)
)

#' @rdname MSNBurr2aUvSpec-class
#' @export
MSNBurr2aUvSpec <- function() methods::new("MSNBurr2aUvSpec")

# --- contract methods (shared implementations dispatched per class) ------------

.msnburrDefaultPrior <- function(spec, data, control = list()) {
  y <- as.numeric(data)
  sy <- stats::sd(y)
  if (!is.finite(sy) || sy <= 0) sy <- 1
  cLoc <- if (!is.null(control$cLoc)) control$cLoc else 2
  ## E[sigma] = bSig / (aSig - 1) = sd(y); weak shape
  aSig <- if (!is.null(control$aSig)) control$aSig else 3
  bSig <- if (!is.null(control$bSig)) control$bSig else (aSig - 1) * sy
  ## alpha ~ Gamma(aA, bA), E = 1 around the logistic/symmetric-ish region
  aA <- if (!is.null(control$aAlphaShape)) control$aAlphaShape else 2
  bA <- if (!is.null(control$bAlphaRate))  control$bAlphaRate  else 2
  list(mu0 = mean(y), muSd = cLoc * sy,
       aSig = aSig, bSig = bSig, aA = aA, bA = bA, cLoc = cLoc)
}

#' @describeIn defaultPrior Data-scaled MSNBurr prior (location/scale/shape).
setMethod("defaultPrior", "MSNBurrUvSpec",
  function(spec, data, control = list(), ...)
    .msnburrDefaultPrior(spec, data, control))

#' @describeIn defaultPrior Data-scaled MSNBurr-IIa prior.
setMethod("defaultPrior", "MSNBurr2aUvSpec",
  function(spec, data, control = list(), ...)
    .msnburrDefaultPrior(spec, data, control))

.msnburrValidate <- function(prior) {
  stopifnot(is.finite(prior$mu0), prior$muSd > 0,
            prior$aSig > 1, prior$bSig > 0, prior$aA > 0, prior$bA > 0)
  invisible(TRUE)
}

#' @describeIn validateParams MSNBurr hyperparameter checks.
setMethod("validateParams", "MSNBurrUvSpec",
  function(spec, params, ...) .msnburrValidate(params))

#' @describeIn validateParams MSNBurr-IIa hyperparameter checks.
setMethod("validateParams", "MSNBurr2aUvSpec",
  function(spec, params, ...) .msnburrValidate(params))

#' @describeIn simulateParams Draw component parameters from the prior.
setMethod("simulateParams", "MSNBurrUvSpec",
  function(spec, prior, K, ...) {
    list(mu = stats::rnorm(K, prior$mu0, prior$muSd),
         sigma = 1 / stats::rgamma(K, prior$aSig, rate = prior$bSig),
         alpha = stats::rgamma(K, prior$aA, rate = prior$bA))
  })

#' @describeIn simulateParams Draw component parameters from the prior.
setMethod("simulateParams", "MSNBurr2aUvSpec",
  function(spec, prior, K, ...) {
    list(mu = stats::rnorm(K, prior$mu0, prior$muSd),
         sigma = 1 / stats::rgamma(K, prior$aSig, rate = prior$bSig),
         alpha = stats::rgamma(K, prior$aA, rate = prior$bA))
  })

#' @describeIn componentDensity MSNBurr density closure (stable reference form).
setMethod("componentDensity", "MSNBurrUvSpec",
  function(spec, ...) {
    function(x, params) dmsnburr(x, params[["mu"]], params[["sigma"]],
                                 params[["alpha"]])
  })

#' @describeIn componentDensity MSNBurr-IIa density closure (stable form).
setMethod("componentDensity", "MSNBurr2aUvSpec",
  function(spec, ...) {
    function(x, params) dmsnburr2a(x, params[["mu"]], params[["sigma"]],
                                   params[["alpha"]])
  })

# --- kernels --------------------------------------------------------------------

.msnburrFixedKCode <- function(densName) {
  tmpl <- sprintf("{
    for (i in 1:n) {
      z[i] ~ dcat(weights[1:K])
      y[i] ~ %s(muTilde[z[i]], sigmaTilde[z[i]], alphaTilde[z[i]])
    }
    weights[1:K] ~ ddirch(alphaVec[1:K])
    for (j in 1:K) {
      muTilde[j]    ~ dnorm(mu0, sd = muSd)
      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
      alphaTilde[j] ~ dgamma(shape = aA, rate = bA)
    }
  }", densName)
  str2lang(tmpl)
}

.msnburrDPMCode <- function(densName) {
  tmpl <- sprintf("{
    for (i in 1:n) {
      y[i] ~ %s(muTilde[xi[i]], sigmaTilde[xi[i]], alphaTilde[xi[i]])
    }
    xi[1:n] ~ dCRP(alpha, size = n)
    for (j in 1:L) {
      muTilde[j]    ~ dnorm(mu0, sd = muSd)
      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
      alphaTilde[j] ~ dgamma(shape = aA, rate = bA)
    }
    alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
  }", densName)
  str2lang(tmpl)
}

.msnburrMc <- function(code, alloc) {
  mons <- c(alloc, "muTilde", "sigmaTilde", "alphaTilde",
            if (alloc == "z") "weights" else "alpha")
  list(code = code, monitors = mons,
       paramNodes = c(mu = "muTilde", sigma = "sigmaTilde",
                      alpha = "alphaTilde"),
       allocNode = alloc)
}

#' @describeIn buildModelCode MSNBurr finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("MSNBurrUvSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...)
    .msnburrMc(.msnburrFixedKCode("dMSNBurr_k"), "z"))

#' @describeIn buildModelCode MSNBurr-IIa finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("MSNBurr2aUvSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...)
    .msnburrMc(.msnburrFixedKCode("dMSNBurr2a_k"), "z"))

#' @describeIn buildModelCode MSNBurr DPM mixture.
#' @export
setMethod("buildModelCode", signature("MSNBurrUvSpec", "DPMEngine"),
  function(spec, engine, n, L, ...)
    .msnburrMc(.msnburrDPMCode("dMSNBurr_k"), "xi"))

#' @describeIn buildModelCode MSNBurr-IIa DPM mixture.
#' @export
setMethod("buildModelCode", signature("MSNBurr2aUvSpec", "DPMEngine"),
  function(spec, engine, n, L, ...)
    .msnburrMc(.msnburrDPMCode("dMSNBurr2a_k"), "xi"))

# --- engine-facing methods ------------------------------------------------------

.msnburrConstants <- function(prior, n)
  list(n = n, mu0 = prior$mu0, muSd = prior$muSd,
       aSig = prior$aSig, bSig = prior$bSig, aA = prior$aA, bA = prior$bA)

#' @describeIn buildConstants MSNBurr constants.
setMethod("buildConstants", "MSNBurrUvSpec",
  function(spec, prior, n, ...) .msnburrConstants(prior, n))

#' @describeIn buildConstants MSNBurr-IIa constants.
setMethod("buildConstants", "MSNBurr2aUvSpec",
  function(spec, prior, n, ...) .msnburrConstants(prior, n))

#' @describeIn buildDataList Univariate response.
setMethod("buildDataList", "MSNBurrUvSpec",
  function(spec, data, ...) list(y = as.numeric(data)))

#' @describeIn buildDataList Univariate response.
setMethod("buildDataList", "MSNBurr2aUvSpec",
  function(spec, data, ...) list(y = as.numeric(data)))

.msnburrInits <- function(prior, data, count, initMethod, initRatio) {
  y <- as.numeric(data); n <- length(y)
  nUnique <- length(unique(y))
  k0 <- max(1L, min(as.integer(floor(initRatio * count)),
                    as.integer(ceiling(sqrt(n)))))
  k0 <- min(k0, max(1L, nUnique))
  alloc <- rep(1L, n)
  centers <- mean(y); sds <- stats::sd(y)
  if (!is.finite(sds) || sds <= 0) sds <- 1
  if (identical(initMethod, "kmeans") && k0 >= 2L && nUnique >= k0) {
    km <- tryCatch(stats::kmeans(y, centers = k0, nstart = 5L),
                   error = function(e) NULL)
    if (!is.null(km)) {
      alloc <- as.integer(km$cluster)
      centers <- as.numeric(km$centers)
      sds <- vapply(seq_len(k0), function(j) {
        s <- stats::sd(y[alloc == j])
        if (!is.finite(s) || s <= 0) 1 else s
      }, numeric(1))
    }
  }
  muInit  <- rep(prior$mu0, count)
  sigInit <- rep(prior$bSig / (prior$aSig - 1), count)
  alInit  <- rep(1, count)                 # logistic-shaped start
  occ <- sort(unique(alloc))
  for (idx in seq_along(occ)) {
    j <- occ[idx]
    if (length(centers) >= idx) muInit[j] <- centers[idx]
    if (length(sds) >= idx && is.finite(sds[idx]) && sds[idx] > 0)
      sigInit[j] <- sds[idx]
  }
  list(alloc = alloc,
       params = list(muTilde = muInit, sigmaTilde = sigInit,
                     alphaTilde = alInit))
}

#' @describeIn componentInits Dispersed k-means start for MSNBurr.
setMethod("componentInits", "MSNBurrUvSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...)
    .msnburrInits(prior, data, count, initMethod, .initRatioArg(...)))

#' @describeIn componentInits Dispersed k-means start for MSNBurr-IIa.
setMethod("componentInits", "MSNBurr2aUvSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...)
    .msnburrInits(prior, data, count, initMethod, .initRatioArg(...)))

#' @describeIn extractParamTraces Parse mu / sigma / alpha traces.
setMethod("extractParamTraces", "MSNBurrUvSpec",
  function(spec, samples, L, ...) {
    list(mu    = .nodeToArray(samples, "muTilde",    L),
         sigma = .nodeToArray(samples, "sigmaTilde", L),
         alpha = .nodeToArray(samples, "alphaTilde", L))
  })

#' @describeIn extractParamTraces Parse mu / sigma / alpha traces.
setMethod("extractParamTraces", "MSNBurr2aUvSpec",
  function(spec, samples, L, ...) {
    list(mu    = .nodeToArray(samples, "muTilde",    L),
         sigma = .nodeToArray(samples, "sigmaTilde", L),
         alpha = .nodeToArray(samples, "alphaTilde", L))
  })

.msnburrRelabel <- function(paramTrace, idx, occList, perms, modalK, weights) {
  muTr <- paramTrace$mu[idx, , drop = FALSE]
  sgTr <- paramTrace$sigma[idx, , drop = FALSE]
  alTr <- paramTrace$alpha[idx, , drop = FALSE]
  m <- length(idx)
  muRe <- matrix(NA_real_, m, modalK)
  sgRe <- matrix(NA_real_, m, modalK)
  alRe <- matrix(NA_real_, m, modalK)
  for (r in seq_len(m)) {
    occ <- occList[[r]]; pr <- perms[r, ]
    muRe[r, ] <- muTr[r, occ][pr]
    sgRe[r, ] <- sgTr[r, occ][pr]
    alRe[r, ] <- alTr[r, occ][pr]
  }
  q <- function(M, p) apply(M, 2L, stats::quantile, probs = p, names = FALSE)
  summ <- data.frame(
    component = seq_len(modalK),
    weight    = colMeans(weights),
    mu_mean   = colMeans(muRe),
    mu_lwr    = q(muRe, 0.025), mu_upr = q(muRe, 0.975),
    sigma_mean = colMeans(sgRe),
    sigma_lwr  = q(sgRe, 0.025), sigma_upr = q(sgRe, 0.975),
    alpha_mean = colMeans(alRe),
    alpha_lwr  = q(alRe, 0.025), alpha_upr = q(alRe, 0.975)
  )
  list(mu = muRe, sigma = sgRe, alpha = alRe, summary = summ)
}

#' @describeIn relabelComponents Permute mu / sigma / alpha and summarise.
setMethod("relabelComponents", "MSNBurrUvSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...)
    .msnburrRelabel(paramTrace, idx, occList, perms, modalK, weights))

#' @describeIn relabelComponents Permute mu / sigma / alpha and summarise.
setMethod("relabelComponents", "MSNBurr2aUvSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...)
    .msnburrRelabel(paramTrace, idx, occList, perms, modalK, weights))
