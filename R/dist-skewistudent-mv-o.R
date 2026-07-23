## dist-skewistudent-mv-o.R ------------------------------------------------------
## Ferreira-Steel skew multivariate independent-Student with the orthogonal
## factor O of A = OU ESTIMATED via the Householder angle theta (FS 2007,
## Lemma 1 and Appendix A). This is to "skewistudent-mv" what
## "skewnormal-mv-o" is to "skewnormal-mv": the skewness *directions* are
## learned rather than tied to the coordinate axes (FS Sec 3.3).
##
## The angle machinery is shared with the skew-Normal-O family: for m = 2,
## O = I - 2 v v' with v = (sin theta, cos theta), so O11 = cos(2 theta),
## O21 = -sin(2 theta), |O| = -1, and FS's identifiability restriction (8)
## reduces to theta in (-pi/8, pi/8). Because |O| = -1 always, O = I is not in
## the restricted set: theta = 0 gives O = diag(1, -1), which equals
## "skewistudent-mv" with gamma_2 replaced by 1/gamma_2.
##
## Identifiability: BETTER here than in the Gaussian case, and for a reason
## worth stating. In skewnormal-mv-o the density is theta-invariant at
## gamma = 1, because eps is then spherical Normal and Lemma 1 says only
## Sigma = A'A is identified; theta is recoverable only through the skewness.
## Independent Student margins are NOT spherical, so the skew-IStudent density
## depends on theta even at gamma = 1 -- verified numerically, and the profile
## likelihood recovers theta from symmetric data. Letting nu -> Inf restores
## sphericity and with it the invariance, exactly as the theory predicts.
##
## The near-mirror secondary mode of the theta likelihood still exists, so
## theta remains a large-sample quantity: read its trace, and do not
## over-interpret a sign estimated from a small component. thetaTilde is
## slice-sampled and grid-initialised.
##
## Scope: m = 2, matching skewnormal-mv-o.

#' @include class-DistributionSpec.R
#' @include dist-skewistudent-mv.R
#' @include dist-skewnormal-mv-o.R
NULL

#' Ferreira-Steel skew multivariate independent-Student with estimated O
#'
#' Density and random generation for the FS skew multivariate
#' independent-Student in which the orthogonal factor of \eqn{A = OU} is
#' estimated through the Householder angle \code{theta}. Bivariate only:
#' \code{theta} lies in \eqn{(-\pi/8, \pi/8)}, exactly FS's restriction (8).
#'
#' @inheritParams skewnormal-mv-o-distribution
#' @param nu Positive numeric vector of per-dimension degrees of freedom.
#' @return \code{dskewmvito} a numeric vector; \code{rskewmvito} an
#'   \code{n x 2} matrix.
#' @references Ferreira & Steel (2007), Statistica Sinica 17, 505--529.
#' @name skewistudent-mv-o-distribution
NULL

#' @rdname skewistudent-mv-o-distribution
#' @export
dskewmvito <- function(x, mu, Sigma, gamma, nu, theta, log = FALSE) {
  X <- if (is.matrix(x)) x else matrix(x, nrow = 1L)
  stopifnot(length(mu) == 2L, ncol(X) == 2L, length(gamma) == 2L,
            length(nu) == 2L, all(gamma > 0), all(nu > 0))
  U  <- chol(Sigma)
  Ui <- backsolve(U, diag(2))
  W  <- sweep(X, 2L, mu) %*% Ui
  E  <- W %*% .householderO(theta)
  G  <- matrix(gamma, nrow(E), 2L, byrow = TRUE)
  S  <- ifelse(E < 0, E * G, E / G)
  lp <- numeric(nrow(E))
  for (j in 1:2)
    lp <- lp + log(2) - log(gamma[j] + 1 / gamma[j]) +
      stats::dt(S[, j], df = nu[j], log = TRUE)
  lp <- lp - sum(log(diag(U)))
  if (log) lp else exp(lp)
}

#' @rdname skewistudent-mv-o-distribution
#' @export
rskewmvito <- function(n, mu, Sigma, gamma, nu, theta) {
  stopifnot(length(mu) == 2L, length(gamma) == 2L, length(nu) == 2L,
            all(gamma > 0), all(nu > 0))
  U <- chol(Sigma)
  G <- matrix(gamma, n, 2L, byrow = TRUE)
  W <- abs(matrix(stats::rt(n * 2L, df = rep(nu, each = n)), n, 2L))
  pos <- matrix(stats::runif(n * 2L), n, 2L) < G^2 / (1 + G^2)
  Eps <- ifelse(pos, W * G, -W / G)
  sweep(Eps %*% .householderO(theta) %*% U, 2L, mu, "+")
}

#' Skew mv independent-Student components with estimated orthogonal factor
#'
#' As \code{\link{SkewIStudentMvSpec-class}}, but the orthogonal factor of
#' \eqn{A = OU} is estimated through the Householder angle \code{theta} with a
#' uniform prior on \eqn{(-\pi/8, \pi/8)} (FS restriction (8)). Bivariate data
#' only. \code{theta} is identified only through the skewness.
#'
#' @references Ferreira & Steel (2007), Statistica Sinica 17, 505--529.
#' @keywords internal
#' @rdname SkewIStudentMvOSpec-class
#' @export
setClass("SkewIStudentMvOSpec", contains = "SkewIStudentMvSpec",
  prototype = prototype(name = "skewistudent-mv-o",
                        paramNames = c("mu", "Sigma", "gamma", "nu", "theta"),
                        dataDim = NA_integer_))

#' @rdname SkewIStudentMvOSpec-class
#' @export
SkewIStudentMvOSpec <- function() methods::new("SkewIStudentMvOSpec")

#' @describeIn defaultPrior Adds the Householder angle bound.
setMethod("defaultPrior", "SkewIStudentMvOSpec",
  function(spec, data, control = list(), ...) {
    p <- getMethod("defaultPrior", "SkewIStudentMvSpec")(spec, data, control, ...)
    if (p$d != 2L)
      stop("distribution = 'skewistudent-mv-o' currently supports m = 2 only; ",
           "got d = ", p$d, ". Use 'skewistudent-mv' for higher dimensions.",
           call. = FALSE)
    p$thetaBound <- .thetaBound
    p
  })

#' @describeIn validateParams Checks the Householder bound.
setMethod("validateParams", "SkewIStudentMvOSpec",
  function(spec, params, ...) {
    getMethod("validateParams", "SkewIStudentMvSpec")(spec, params, ...)
    stopifnot(params$d == 2L, params$thetaBound > 0,
              params$thetaBound <= pi / 8 + 1e-12)
    invisible(TRUE)
  })

#' @describeIn simulateParams Draw components, theta uniform on Theta^2.
setMethod("simulateParams", "SkewIStudentMvOSpec",
  function(spec, prior, nClust, ...) {
    p <- getMethod("simulateParams", "SkewIStudentMvSpec")(spec, prior, nClust, ...)
    p$theta <- as.list(stats::runif(nClust, -prior$thetaBound, prior$thetaBound))
    p
  })

#' @describeIn componentDensity Density closure with nu and the angle.
setMethod("componentDensity", "SkewIStudentMvOSpec",
  function(spec, ...) {
    function(x, params) dskewmvito(x, params[["mu"]], params[["Sigma"]],
                                   params[["gamma"]], params[["nu"]],
                                   params[["theta"]])
  })

.skewMvITOFixedKCode <- function() nimble::nimbleCode({
  for (i in 1:n) {
    z[i] ~ dcat(weights[1:K])
    muObs[i, 1:d]       <- muTilde[z[i], 1:d]
    SigObs[i, 1:d, 1:d] <- SigTilde[z[i], 1:d, 1:d]
    gamObs[i, 1:d]      <- gamTilde[z[i], 1:d]
    nuObs[i, 1:d]       <- nuTilde[z[i], 1:d]
    thetaObs[i]         <- thetaTilde[z[i]]
    y[i, 1:d] ~ dSkewMvITO_k(muObs[i, 1:d], SigObs[i, 1:d, 1:d],
                             gamObs[i, 1:d], nuObs[i, 1:d], thetaObs[i])
  }
  weights[1:K] ~ ddirch(alphaVec[1:K])
  for (j in 1:K) {
    SigTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
    muTilde[j, 1:d] ~ dmnorm(mu0[1:d], cov = covMu[1:d, 1:d])
    for (r in 1:d) {
      gamTilde[j, r] ~ dlnorm(0, sdlog = gamLogSd)
      nuTilde[j, r]  ~ T(dgamma(shape = aNu, rate = bNu), nuLower, )
    }
    thetaTilde[j] ~ dunif(-thetaBound, thetaBound)
  }
})

.skewMvITODPMCode <- function() nimble::nimbleCode({
  for (i in 1:n) {
    muObs[i, 1:d]       <- muTilde[xi[i], 1:d]
    SigObs[i, 1:d, 1:d] <- SigTilde[xi[i], 1:d, 1:d]
    gamObs[i, 1:d]      <- gamTilde[xi[i], 1:d]
    nuObs[i, 1:d]       <- nuTilde[xi[i], 1:d]
    thetaObs[i]         <- thetaTilde[xi[i]]
    y[i, 1:d] ~ dSkewMvITO_k(muObs[i, 1:d], SigObs[i, 1:d, 1:d],
                             gamObs[i, 1:d], nuObs[i, 1:d], thetaObs[i])
  }
  xi[1:n] ~ dCRP(alpha, size = n)
  for (j in 1:L) {
    SigTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
    muTilde[j, 1:d] ~ dmnorm(mu0[1:d], cov = covMu[1:d, 1:d])
    for (r in 1:d) {
      gamTilde[j, r] ~ dlnorm(0, sdlog = gamLogSd)
      nuTilde[j, r]  ~ T(dgamma(shape = aNu, rate = bNu), nuLower, )
    }
    thetaTilde[j] ~ dunif(-thetaBound, thetaBound)
  }
  alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
})

.skewMvITOMc <- function(code, alloc) list(
  code = code,
  monitors = c(alloc, "muTilde", "SigTilde", "gamTilde", "nuTilde",
               "thetaTilde", if (alloc == "z") "weights" else "alpha"),
  paramNodes = c(mu = "muTilde", Sigma = "SigTilde", gamma = "gamTilde",
                 nu = "nuTilde", theta = "thetaTilde"),
  allocNode = alloc)

#' @describeIn buildModelCode Skew-mv-IStudent-O finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("SkewIStudentMvOSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) .skewMvITOMc(.skewMvITOFixedKCode(), "z"))

#' @describeIn buildModelCode Skew-mv-IStudent-O DPM mixture.
#' @export
setMethod("buildModelCode", signature("SkewIStudentMvOSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) .skewMvITOMc(.skewMvITODPMCode(), "xi"))

#' @describeIn buildConstants Skew-mv-IStudent-O constants.
setMethod("buildConstants", "SkewIStudentMvOSpec",
  function(spec, prior, n, ...) {
    cn <- getMethod("buildConstants", "SkewIStudentMvSpec")(spec, prior, n, ...)
    cn$thetaBound <- prior$thetaBound
    cn
  })

#' @describeIn componentInits k-means start; theta grid-initialised per cluster.
setMethod("componentInits", "SkewIStudentMvOSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    init <- getMethod("componentInits", "SkewIStudentMvSpec")(
      spec, prior, data, count, initMethod, ...)
    Y <- as.matrix(data)
    thInit <- rep(0, count)
    alloc <- init$alloc
    for (j in sort(unique(alloc))) {
      Yj <- Y[alloc == j, , drop = FALSE]
      if (nrow(Yj) >= 10L && j <= count) {
        nuJ <- init$params$nuTilde[j, ]
        thInit[j] <- .initThetaGrid(
          Yj, init$params$muTilde[j, ],
          matrix(init$params$SigTilde[j, , ], 2, 2), prior$thetaBound,
          llFun = function(Y, mu, Sigma, g, th)
            sum(dskewmvito(Y, mu, Sigma, g, nuJ, th, log = TRUE)))
      }
    }
    init$params$thetaTilde <- thInit
    init
  })

#' @describeIn customizeSamplers Slice-sample the Householder angles.
#' @export
setMethod("customizeSamplers", "SkewIStudentMvOSpec",
  function(spec, conf, model, ...) {
    nodes <- model$expandNodeNames("thetaTilde")
    if (length(nodes) == 0L) return(invisible(conf))
    conf$removeSamplers("thetaTilde")
    for (nd in nodes) conf$addSampler(target = nd, type = "slice")
    invisible(conf)
  })

#' @describeIn extractParamTraces Parse mu / Sigma / gamma / nu / theta traces.
setMethod("extractParamTraces", "SkewIStudentMvOSpec",
  function(spec, samples, L, d = NULL, ...) {
    tr <- getMethod("extractParamTraces", "SkewIStudentMvSpec")(
      spec, samples, L, d = d, ...)
    tr$theta <- .nodeToArray(samples, "thetaTilde", L)
    tr
  })

#' @describeIn relabelComponents Permute mv params plus nu and theta.
setMethod("relabelComponents", "SkewIStudentMvOSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...) {
    out <- getMethod("relabelComponents", "SkewIStudentMvSpec")(
      spec, paramTrace, idx, occList, perms, modalK, weights, ...)
    thTr <- paramTrace$theta; m <- length(idx)
    thRe <- matrix(NA_real_, m, modalK)
    for (t in seq_len(m)) {
      r <- idx[t]; occ <- occList[[t]]; ord <- perms[t, ]
      thRe[t, ] <- thTr[r, occ][ord]
    }
    q <- function(M, p) apply(M, 2L, stats::quantile, probs = p, names = FALSE)
    out$summary$theta_mean <- colMeans(thRe)
    out$summary$theta_lwr  <- q(thRe, 0.025)
    out$summary$theta_upr  <- q(thRe, 0.975)
    out$theta <- thRe
    out
  })
