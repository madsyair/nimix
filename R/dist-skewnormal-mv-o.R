## dist-skewnormal-mv-o.R ---------------------------------------------------------
## Ferreira-Steel skew multivariate Normal with the orthogonal factor O of
## A = OU ESTIMATED rather than fixed (FS 2007, Lemma 1 and Appendix A).
##
## Why O matters. Under symmetry only Sigma = A'A is identified, so the
## orthogonal factor is irrelevant. Under skewness it is not: O fixes the
## *directions* along which the FS mechanism skews (FS Sec 3.3). The sibling
## family "skewnormal-mv" holds O fixed; this one learns it.
##
## Parameterisation (m = 2 only). FS reparameterise O by Householder angles:
## O = I - 2 v v' with v = (sin theta, cos theta). In closed form
##   O11 = cos(2 theta), O21 = -sin(2 theta), and |O| = -1 always.
## The identifiability restriction (8), O11 > |O21| > 0, is therefore exactly
##   theta in Theta^2 = (-pi/8, pi/8),
## which we impose as the (proper, uniform) prior support. Note that |O| = -1
## means O = I is NOT a member of FS's restricted set: theta = 0 gives
## O = diag(1, -1), which coincides with "skewnormal-mv" after replacing
## gamma_2 by 1/gamma_2. The O = I family is thus nested here at theta = 0, up
## to that reflection of gamma.
##
## Identifiability caveats, stated plainly.
##  (1) When gamma = 1 the density is invariant in theta (it collapses to
##      dmnorm for every theta), so theta is identified ONLY through the
##      skewness. Expect a diffuse posterior on theta for near-symmetric
##      components.
##  (2) Even with clear skewness, the likelihood in theta has a secondary,
##      near-mirror mode. With mu, Sigma and gamma free it can sit within a
##      couple of log-likelihood units of the true mode at moderate sample
##      sizes: in our checks, at n = 150 per component the mirror mode was only
##      1.65 log-lik units worse (and both beat the log-lik at the true
##      parameters), and chains reliably settled on the wrong sign; at n = 500
##      per component theta was recovered with 95% intervals covering the truth.
##      Treat theta as a large-sample quantity, inspect its trace, and do not
##      over-interpret a sign estimated from small components. thetaTilde is
##      slice-sampled and grid-initialised for this reason.
##
## Scope: m = 2. General m needs O = O_theta^m x ... x O_theta^2 with
## m(m-1)/2 angles and the full chain restriction (8); planned separately.

#' @include class-DistributionSpec.R
#' @include dist-skewnormal-mv.R
NULL

.thetaBound <- pi / 8   # = 0.3926991; Theta^2 for m = 2

.householderO <- function(theta) {
  v <- c(sin(theta), cos(theta))
  diag(2) - 2 * tcrossprod(v)
}

#' Ferreira-Steel skew multivariate Normal with estimated orthogonal factor
#'
#' Density and random generation for the FS skew multivariate Normal in which
#' the orthogonal factor of \eqn{A = OU} is estimated via the Householder angle
#' \code{theta}. Only \code{m = 2} is supported: \code{theta} must lie in
#' \eqn{(-\pi/8, \pi/8)}, which is exactly FS's identifiability restriction (8).
#'
#' @inheritParams skewnormal-mv-distribution
#' @param theta Householder angle in \eqn{(-\pi/8, \pi/8)}.
#' @return \code{dskewmvno} a numeric vector; \code{rskewmvno} an
#'   \code{n x 2} matrix.
#' @references Ferreira & Steel (2007), Statistica Sinica 17, 505--529.
#' @name skewnormal-mv-o-distribution
NULL

#' @rdname skewnormal-mv-o-distribution
#' @export
dskewmvno <- function(x, mu, Sigma, gamma, theta, log = FALSE) {
  X <- if (is.matrix(x)) x else matrix(x, nrow = 1L)
  d <- length(mu)
  stopifnot(d == 2L, ncol(X) == 2L, length(gamma) == 2L, all(gamma > 0))
  U  <- chol(Sigma)
  Ui <- backsolve(U, diag(2))
  W  <- sweep(X, 2L, mu) %*% Ui            # rows: w' = (x - mu)' U^{-1}
  E  <- W %*% .householderO(theta)         # eps' = w' O' ; O symmetric
  G  <- matrix(gamma, nrow(E), 2L, byrow = TRUE)
  S  <- ifelse(E < 0, E * G, E / G)
  lp <- rowSums(matrix(log(2) - log(gamma + 1 / gamma), nrow(E), 2L,
                       byrow = TRUE) + stats::dnorm(S, log = TRUE)) -
    sum(log(diag(U)))
  if (log) lp else exp(lp)
}

#' @rdname skewnormal-mv-o-distribution
#' @export
rskewmvno <- function(n, mu, Sigma, gamma, theta) {
  stopifnot(length(mu) == 2L, length(gamma) == 2L, all(gamma > 0))
  U <- chol(Sigma)
  G <- matrix(gamma, n, 2L, byrow = TRUE)
  W <- abs(matrix(stats::rnorm(n * 2L), n, 2L))
  pos <- matrix(stats::runif(n * 2L), n, 2L) < G^2 / (1 + G^2)
  Eps <- ifelse(pos, W * G, -W / G)
  # eta' = eps' O U + mu'   (eta = A' eps + mu, A = OU, O symmetric)
  sweep(Eps %*% .householderO(theta) %*% U, 2L, mu, "+")
}

#' Skew multivariate Normal components with estimated orthogonal factor
#'
#' As \code{\link{SkewNormalMvSpec-class}}, but the orthogonal factor of
#' \eqn{A = OU} is estimated through the Householder angle \code{theta} with a
#' uniform prior on \eqn{(-\pi/8, \pi/8)} (FS restriction (8)). Bivariate data
#' only. \code{theta} is identified only through the skewness: at
#' \code{gamma = 1} the density is invariant in \code{theta}.
#'
#' @references Ferreira & Steel (2007), Statistica Sinica 17, 505--529.
#' @keywords internal
#' @rdname SkewNormalMvOSpec-class
#' @export
setClass("SkewNormalMvOSpec", contains = "SkewNormalMvSpec",
  prototype = prototype(name = "skewnormal-mv-o",
                        paramNames = c("mu", "Sigma", "gamma", "theta"),
                        dataDim = NA_integer_))

#' @rdname SkewNormalMvOSpec-class
#' @export
SkewNormalMvOSpec <- function() methods::new("SkewNormalMvOSpec")

#' @describeIn defaultPrior Adds the Householder angle bound to the mv prior.
setMethod("defaultPrior", "SkewNormalMvOSpec",
  function(spec, data, control = list(), ...) {
    p <- getMethod("defaultPrior", "SkewNormalMvSpec")(spec, data, control, ...)
    if (p$d != 2L)
      stop("distribution = 'skewnormal-mv-o' currently supports m = 2 only; ",
           "got d = ", p$d, ". Use 'skewnormal-mv' for higher dimensions.",
           call. = FALSE)
    p$thetaBound <- .thetaBound
    p
  })

#' @describeIn validateParams Checks the Householder bound.
setMethod("validateParams", "SkewNormalMvOSpec",
  function(spec, params, ...) {
    getMethod("validateParams", "SkewNormalMvSpec")(spec, params, ...)
    stopifnot(params$d == 2L, params$thetaBound > 0,
              params$thetaBound <= pi / 8 + 1e-12)
    invisible(TRUE)
  })

#' @describeIn simulateParams Draw components, theta uniform on Theta^2.
setMethod("simulateParams", "SkewNormalMvOSpec",
  function(spec, prior, K, ...) {
    p <- getMethod("simulateParams", "SkewNormalMvSpec")(spec, prior, K, ...)
    p$theta <- as.list(stats::runif(K, -prior$thetaBound, prior$thetaBound))
    p
  })

#' @describeIn componentDensity Density closure with the Householder angle.
setMethod("componentDensity", "SkewNormalMvOSpec",
  function(spec, ...) {
    function(x, params) dskewmvno(x, params[["mu"]], params[["Sigma"]],
                                  params[["gamma"]], params[["theta"]])
  })

.skewMvNOFixedKCode <- function() nimble::nimbleCode({
  for (i in 1:n) {
    z[i] ~ dcat(weights[1:K])
    muObs[i, 1:d]       <- muTilde[z[i], 1:d]
    SigObs[i, 1:d, 1:d] <- SigTilde[z[i], 1:d, 1:d]
    gamObs[i, 1:d]      <- gamTilde[z[i], 1:d]
    thetaObs[i]         <- thetaTilde[z[i]]
    y[i, 1:d] ~ dSkewMvNO_k(muObs[i, 1:d], SigObs[i, 1:d, 1:d],
                            gamObs[i, 1:d], thetaObs[i])
  }
  weights[1:K] ~ ddirch(alphaVec[1:K])
  for (j in 1:K) {
    SigTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
    muTilde[j, 1:d] ~ dmnorm(mu0[1:d], cov = covMu[1:d, 1:d])
    for (r in 1:d) gamTilde[j, r] ~ dlnorm(0, sdlog = gamLogSd)
    thetaTilde[j] ~ dunif(-thetaBound, thetaBound)
  }
})

.skewMvNODPMCode <- function() nimble::nimbleCode({
  for (i in 1:n) {
    muObs[i, 1:d]       <- muTilde[xi[i], 1:d]
    SigObs[i, 1:d, 1:d] <- SigTilde[xi[i], 1:d, 1:d]
    gamObs[i, 1:d]      <- gamTilde[xi[i], 1:d]
    thetaObs[i]         <- thetaTilde[xi[i]]
    y[i, 1:d] ~ dSkewMvNO_k(muObs[i, 1:d], SigObs[i, 1:d, 1:d],
                            gamObs[i, 1:d], thetaObs[i])
  }
  xi[1:n] ~ dCRP(alpha, size = n)
  for (j in 1:L) {
    SigTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
    muTilde[j, 1:d] ~ dmnorm(mu0[1:d], cov = covMu[1:d, 1:d])
    for (r in 1:d) gamTilde[j, r] ~ dlnorm(0, sdlog = gamLogSd)
    thetaTilde[j] ~ dunif(-thetaBound, thetaBound)
  }
  alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
})

.skewMvNOMc <- function(code, alloc) list(
  code = code,
  monitors = c(alloc, "muTilde", "SigTilde", "gamTilde", "thetaTilde",
               if (alloc == "z") "weights" else "alpha"),
  paramNodes = c(mu = "muTilde", Sigma = "SigTilde", gamma = "gamTilde",
                 theta = "thetaTilde"),
  allocNode = alloc)

#' @describeIn buildModelCode Skew-mv-Normal-O finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("SkewNormalMvOSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) .skewMvNOMc(.skewMvNOFixedKCode(), "z"))

#' @describeIn buildModelCode Skew-mv-Normal-O DPM mixture.
#' @export
setMethod("buildModelCode", signature("SkewNormalMvOSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) .skewMvNOMc(.skewMvNODPMCode(), "xi"))

#' @describeIn buildConstants Skew-mv-Normal-O constants.
setMethod("buildConstants", "SkewNormalMvOSpec",
  function(spec, prior, n, ...) {
    cn <- getMethod("buildConstants", "SkewNormalMvSpec")(spec, prior, n, ...)
    cn$thetaBound <- prior$thetaBound
    cn
  })

# Coarse per-cluster starting value for theta. Starting at theta = 0 leaves the
# chain to find its own way across Theta^2, and the FS likelihood has a
# secondary local mode at the boundary of Theta^2 that an adaptive random walk
# can get stuck in. Instead we grid over theta and, for each candidate, obtain
# gamma in closed form from the FS identity P(eps_j > 0) = gamma_j^2/(1+gamma_j^2)
# -- i.e. gamma_j = sqrt(p_j / (1 - p_j)) with p_j the empirical right-tail
# fraction of eps_j -- then keep the theta with the highest log-likelihood.
# `llFun(Y, mu, Sigma, gamma, theta)` must return the family's log-likelihood.
# It is a parameter rather than hard-wired to the Gaussian density: grid-initialising
# a heavy-tailed family with a Gaussian profile picks the wrong angle, because
# outliers dominate the Gaussian fit.
.initThetaGrid <- function(Y, mu, Sigma, bound, ngrid = 31L,
                           llFun = function(Y, mu, Sigma, g, th)
                             sum(dskewmvno(Y, mu, Sigma, g, th, log = TRUE))) {
  ths <- seq(-bound * 0.95, bound * 0.95, length.out = ngrid)
  U <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(U) || nrow(Y) < 5L) return(0)
  Ui <- backsolve(U, diag(2))
  W <- sweep(Y, 2L, mu) %*% Ui
  best <- 0; bestLL <- -Inf
  for (th in ths) {
    E <- W %*% .householderO(th)
    p <- pmin(pmax(colMeans(E > 0), 0.02), 0.98)
    g <- sqrt(p / (1 - p))
    ll <- tryCatch(llFun(Y, mu, Sigma, g, th), error = function(e) -Inf)
    if (is.finite(ll) && ll > bestLL) { bestLL <- ll; best <- th }
  }
  best
}

#' @describeIn componentInits k-means start; theta grid-initialised per cluster.
setMethod("componentInits", "SkewNormalMvOSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    init <- getMethod("componentInits", "SkewNormalMvSpec")(
      spec, prior, data, count, initMethod, ...)
    Y <- as.matrix(data)
    thInit <- rep(0, count)
    alloc <- init$alloc
    for (j in sort(unique(alloc))) {
      Yj <- Y[alloc == j, , drop = FALSE]
      if (nrow(Yj) >= 10L && j <= count)
        thInit[j] <- .initThetaGrid(Yj, init$params$muTilde[j, ],
                                    matrix(init$params$SigTilde[j, , ], 2, 2),
                                    prior$thetaBound)
    }
    init$params$thetaTilde <- thInit
    init
  })

#' @describeIn customizeSamplers Slice-sample the Householder angles: the FS
#'   likelihood in theta is bounded and can be multimodal near the edge of
#'   Theta^2, where an adaptive random walk mixes poorly.
#' @export
setMethod("customizeSamplers", "SkewNormalMvOSpec",
  function(spec, conf, model, ...) {
    nodes <- model$expandNodeNames("thetaTilde")
    if (length(nodes) == 0L) return(invisible(conf))
    conf$removeSamplers("thetaTilde")
    for (nd in nodes) conf$addSampler(target = nd, type = "slice")
    invisible(conf)
  })

#' @describeIn extractParamTraces Parse mu / Sigma / gamma / theta traces.
setMethod("extractParamTraces", "SkewNormalMvOSpec",
  function(spec, samples, L, d = NULL, ...) {
    tr <- getMethod("extractParamTraces", "SkewNormalMvSpec")(
      spec, samples, L, d = d, ...)
    tr$theta <- .nodeToArray(samples, "thetaTilde", L)
    tr
  })

#' @describeIn relabelComponents Permute mv params plus theta, summarise.
setMethod("relabelComponents", "SkewNormalMvOSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...) {
    out <- getMethod("relabelComponents", "SkewNormalMvSpec")(
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
