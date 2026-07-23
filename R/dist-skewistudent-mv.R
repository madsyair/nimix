## dist-skewistudent-mv.R ---------------------------------------------------------
## Ferreira-Steel skew multivariate independent-Student (FS 2007, Sec 5.2):
## eta = A' eps + mu with independent FS-skew Student-t margins, per-dimension
## degrees of freedom nu_j. Closed-form density (no lambda augmentation), the
## model most supported by the data in FS's own application. Same design as
## dist-skewnormal-mv.R: A = chol(Sigma) (O = I; see that header), harmonised
## gamma convention, Sigma ~ inverse-Wishart. nu_j is a stochastic node
## truncated below at 2 so component variances exist (dCRP rejects
## deterministic per-component nodes); FS additionally recommend
## nu > max(3, m-1) in their improper-prior regression setting -- with the
## proper priors used here the truncation at 2 suffices, but raising it via
## control$nuLower is supported.

#' @include class-DistributionSpec.R
#' @include dist-skewnormal-mv.R
NULL

#' Ferreira-Steel skew multivariate independent-Student
#'
#' Density and RNG for the FS skew multivariate independent-Student: FS-skew
#' Student-t margins with per-dimension \code{nu}, transformed by the
#' upper-triangular Cholesky factor of \code{Sigma}. \code{gamma = 1} gives the
#' symmetric independent-Student; \code{nu -> Inf} recovers
#' \code{\link{dskewmvn}}.
#'
#' @inheritParams skewnormal-mv-distribution
#' @param nu Positive numeric vector of per-dimension degrees of freedom.
#' @return \code{dskewmvit} numeric vector; \code{rskewmvit} an \code{n x d}
#'   matrix.
#' @references Ferreira & Steel (2007), Statistica Sinica 17, 505--529.
#' @name skewistudent-mv-distribution
NULL

#' @rdname skewistudent-mv-distribution
#' @export
dskewmvit <- function(x, mu, Sigma, gamma, nu, log = FALSE) {
  X <- if (is.matrix(x)) x else matrix(x, nrow = 1L)
  d <- length(mu)
  stopifnot(ncol(X) == d, length(gamma) == d, length(nu) == d,
            all(gamma > 0), all(nu > 0))
  U  <- chol(Sigma)
  Ui <- backsolve(U, diag(d))
  ldet <- sum(log(diag(U)))
  E <- sweep(X, 2L, mu) %*% Ui
  G <- matrix(gamma, nrow(E), d, byrow = TRUE)
  S <- ifelse(E < 0, E * G, E / G)
  lp <- numeric(nrow(E))
  for (j in seq_len(d))
    lp <- lp + log(2) - log(gamma[j] + 1 / gamma[j]) +
      stats::dt(S[, j], df = nu[j], log = TRUE)
  lp <- lp - ldet
  if (log) lp else exp(lp)
}

#' @rdname skewistudent-mv-distribution
#' @export
rskewmvit <- function(n, mu, Sigma, gamma, nu) {
  d <- length(mu)
  stopifnot(length(gamma) == d, length(nu) == d, all(gamma > 0), all(nu > 0))
  U <- chol(Sigma)
  G <- matrix(gamma, n, d, byrow = TRUE)
  W <- abs(matrix(stats::rt(n * d, df = rep(nu, each = n)), n, d))
  pos <- matrix(stats::runif(n * d), n, d) < G^2 / (1 + G^2)
  Eps <- ifelse(pos, W * G, -W / G)
  sweep(Eps %*% U, 2L, mu, "+")
}

#' Skew multivariate independent-Student mixture components (Ferreira-Steel)
#'
#' Heavy-tailed multivariate component: FS-skew Student-t margins with
#' per-dimension degrees of freedom \eqn{\nu_j} (truncated below at 2),
#' per-dimension skewness \code{gamma}, and \eqn{A = \mathrm{chol}(\Sigma)}
#' upper triangular (orthogonal factor fixed at the identity; see
#' \code{\link{SkewNormalMvSpec-class}}). \code{gamma = 1} is the symmetric
#' independent-Student; \eqn{\nu \to \infty} recovers the skew multivariate
#' Normal. Non-conjugate.
#'
#' @references Ferreira & Steel (2007), Statistica Sinica 17, 505--529.
#' @keywords internal
#' @rdname SkewIStudentMvSpec-class
#' @export
setClass("SkewIStudentMvSpec", contains = "DistributionSpec",
  prototype = prototype(name = "skewistudent-mv",
                        paramNames = c("mu", "Sigma", "gamma", "nu"),
                        dataDim = NA_integer_))

#' @rdname SkewIStudentMvSpec-class
#' @export
SkewIStudentMvSpec <- function() methods::new("SkewIStudentMvSpec")

#' @describeIn defaultPrior Data-scaled skew-mv-IStudent prior.
setMethod("defaultPrior", "SkewIStudentMvSpec",
  function(spec, data, control = list(), ...) {
    p <- getMethod("defaultPrior", "SkewNormalMvSpec")(spec, data, control, ...)
    p$aNu <- if (!is.null(control$aNuShape)) control$aNuShape else 2
    p$bNu <- if (!is.null(control$bNuRate))  control$bNuRate  else 0.15
    p$nuLower <- if (!is.null(control$nuLower)) control$nuLower else 2
    p
  })

#' @describeIn validateParams Skew-mv-IStudent hyperparameter checks.
setMethod("validateParams", "SkewIStudentMvSpec",
  function(spec, params, ...) {
    getMethod("validateParams", "SkewNormalMvSpec")(spec, params, ...)
    stopifnot(params$aNu > 0, params$bNu > 0, params$nuLower >= 2)
    invisible(TRUE)
  })

#' @describeIn simulateParams Draw skew-mv-IStudent component parameters.
setMethod("simulateParams", "SkewIStudentMvSpec",
  function(spec, prior, nClust, ...) {
    p <- getMethod("simulateParams", "SkewNormalMvSpec")(spec, prior, nClust, ...)
    p$nu <- lapply(seq_len(nClust), function(k)
      prior$nuLower + stats::rgamma(prior$d, prior$aNu, rate = prior$bNu))
    p
  })

#' @describeIn componentDensity Skew-mv-IStudent density closure.
setMethod("componentDensity", "SkewIStudentMvSpec",
  function(spec, ...) {
    function(x, params) dskewmvit(x, params[["mu"]], params[["Sigma"]],
                                  params[["gamma"]], params[["nu"]])
  })

.skewMvITFixedKCode <- function() nimble::nimbleCode({
  for (i in 1:n) {
    z[i] ~ dcat(weights[1:K])
    muObs[i, 1:d]       <- muTilde[z[i], 1:d]
    SigObs[i, 1:d, 1:d] <- SigTilde[z[i], 1:d, 1:d]
    gamObs[i, 1:d]      <- gamTilde[z[i], 1:d]
    nuObs[i, 1:d]       <- nuTilde[z[i], 1:d]
    y[i, 1:d] ~ dSkewMvIT_k(muObs[i, 1:d], SigObs[i, 1:d, 1:d],
                            gamObs[i, 1:d], nuObs[i, 1:d])
  }
  weights[1:K] ~ ddirch(alphaVec[1:K])
  for (j in 1:K) {
    SigTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
    muTilde[j, 1:d] ~ dmnorm(mu0[1:d], cov = covMu[1:d, 1:d])
    for (r in 1:d) {
      gamTilde[j, r] ~ dlnorm(0, sdlog = gamLogSd)
      nuTilde[j, r]  ~ T(dgamma(shape = aNu, rate = bNu), nuLower, )
    }
  }
})

.skewMvITDPMCode <- function() nimble::nimbleCode({
  for (i in 1:n) {
    muObs[i, 1:d]       <- muTilde[xi[i], 1:d]
    SigObs[i, 1:d, 1:d] <- SigTilde[xi[i], 1:d, 1:d]
    gamObs[i, 1:d]      <- gamTilde[xi[i], 1:d]
    nuObs[i, 1:d]       <- nuTilde[xi[i], 1:d]
    y[i, 1:d] ~ dSkewMvIT_k(muObs[i, 1:d], SigObs[i, 1:d, 1:d],
                            gamObs[i, 1:d], nuObs[i, 1:d])
  }
  xi[1:n] ~ dCRP(alpha, size = n)
  for (j in 1:L) {
    SigTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
    muTilde[j, 1:d] ~ dmnorm(mu0[1:d], cov = covMu[1:d, 1:d])
    for (r in 1:d) {
      gamTilde[j, r] ~ dlnorm(0, sdlog = gamLogSd)
      nuTilde[j, r]  ~ T(dgamma(shape = aNu, rate = bNu), nuLower, )
    }
  }
  alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
})

.skewMvITMc <- function(code, alloc) list(
  code = code,
  monitors = c(alloc, "muTilde", "SigTilde", "gamTilde", "nuTilde",
               if (alloc == "z") "weights" else "alpha"),
  paramNodes = c(mu = "muTilde", Sigma = "SigTilde", gamma = "gamTilde",
                 nu = "nuTilde"),
  allocNode = alloc)

#' @describeIn buildModelCode Skew-mv-IStudent finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("SkewIStudentMvSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) .skewMvITMc(.skewMvITFixedKCode(), "z"))

#' @describeIn buildModelCode Skew-mv-IStudent DPM mixture.
#' @export
setMethod("buildModelCode", signature("SkewIStudentMvSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) .skewMvITMc(.skewMvITDPMCode(), "xi"))

#' @describeIn buildConstants Skew-mv-IStudent constants.
setMethod("buildConstants", "SkewIStudentMvSpec",
  function(spec, prior, n, ...)
    list(n = n, d = prior$d, mu0 = prior$mu0, covMu = prior$covMu,
         S0 = prior$S0, df0 = prior$df0, gamLogSd = prior$gamLogSd,
         aNu = prior$aNu, bNu = prior$bNu, nuLower = prior$nuLower))

#' @describeIn buildDataList Multivariate response matrix.
setMethod("buildDataList", "SkewIStudentMvSpec",
  function(spec, data, ...) list(y = as.matrix(data)))

#' @describeIn componentInits k-means start (gamma at 1, nu at 8).
setMethod("componentInits", "SkewIStudentMvSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    init <- getMethod("componentInits", "SkewNormalMvSpec")(
      spec, prior, data, count, initMethod, ...)
    init$params$nuTilde <- matrix(8, count, prior$d)
    init
  })

#' @describeIn extractParamTraces Parse mu / Sigma / gamma / nu traces.
setMethod("extractParamTraces", "SkewIStudentMvSpec",
  function(spec, samples, L, d = NULL, ...) {
    if (is.null(d)) stop("extractParamTraces(SkewIStudentMvSpec) needs 'd'.",
                         call. = FALSE)
    list(mu    = .nodeToArray(samples, "muTilde",  c(L, d)),
         Sigma = .nodeToArray(samples, "SigTilde", c(L, d, d)),
         gamma = .nodeToArray(samples, "gamTilde", c(L, d)),
         nu    = .nodeToArray(samples, "nuTilde",  c(L, d)),
         d = d)
  })

#' @describeIn relabelComponents Permute mv (mu, Sigma, gamma, nu), summarise.
setMethod("relabelComponents", "SkewIStudentMvSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...) {
    out <- getMethod("relabelComponents", "SkewNormalMvSpec")(
      spec, paramTrace, idx, occList, perms, modalK, weights, ...)
    d <- paramTrace$d; nuTr <- paramTrace$nu; m <- length(idx)
    nuRe <- array(NA_real_, c(m, modalK, d))
    for (t in seq_len(m)) {
      r <- idx[t]; occ <- occList[[t]]; ord <- perms[t, ]
      for (k in seq_len(modalK)) nuRe[t, k, ] <- nuTr[r, occ[ord[k]], ]
    }
    nuMean <- apply(nuRe, c(2L, 3L), mean)
    for (j in seq_len(d))
      out$summary[[paste0("nu_", j, "_mean")]] <- nuMean[, j]
    out$nu <- nuRe
    out
  })
