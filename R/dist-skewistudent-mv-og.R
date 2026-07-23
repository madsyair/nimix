## dist-skewistudent-mv-og.R -------------------------------------------------------
## FS skew multivariate independent-Student with the orthogonal factor O
## estimated for GENERAL m (>= 3). Reached through
## distribution = "skewistudent-mv-o" when the data have more than two columns;
## m = 2 keeps its own validated implementation (SkewIStudentMvOSpec).
##
## Everything about the angles and about restriction (8) is shared with the
## Gaussian general-m family: m(m-1)/2 Householder angles sampled on the FS
## box, and (8) applied as a post-hoc canonicalisation of each posterior draw,
## because among the signed row permutations of A with |P| = +1 exactly one PO
## satisfies (8). The canonicalisation carries nu along with gamma: a row
## permutation permutes the margins, so nu*_i = nu_{perm(i)}, while a sign flip
## leaves nu alone (the Student kernel is symmetric) and inverts gamma.
##
## One genuine difference from the Gaussian case: theta remains identified even
## at gamma = 1, because independent Student margins are not spherical. See
## dist-skewistudent-mv-o.R.

#' @include class-DistributionSpec.R
#' @include dist-skewistudent-mv.R
#' @include dist-skewnormal-mv-og.R
NULL

#' Skew mv independent-Student components with estimated O, general dimension
#'
#' As \code{\link{SkewIStudentMvOSpec-class}} but for any \eqn{m \ge 3}, using
#' the general Householder parameterisation and the canonicalisation of
#' \code{\link{canonicaliseO}}. Reached via
#' \code{distribution = "skewistudent-mv-o"} for data with more than two
#' columns. Experimental; see \code{\link{SkewNormalMvOGenSpec-class}} for how
#' to read \code{gamma} and \code{O} after canonicalisation.
#'
#' @references Ferreira & Steel (2007), Statistica Sinica 17, 505--529.
#' @keywords internal
#' @rdname SkewIStudentMvOGenSpec-class
#' @export
setClass("SkewIStudentMvOGenSpec", contains = "SkewIStudentMvSpec",
  prototype = prototype(name = "skewistudent-mv-og",
                        paramNames = c("mu", "Sigma", "gamma", "nu", "theta"),
                        dataDim = NA_integer_))

#' @rdname SkewIStudentMvOGenSpec-class
#' @export
SkewIStudentMvOGenSpec <- function() methods::new("SkewIStudentMvOGenSpec")

#' @describeIn defaultPrior Adds the Householder angle box.
setMethod("defaultPrior", "SkewIStudentMvOGenSpec",
  function(spec, data, control = list(), ...) {
    p <- getMethod("defaultPrior", "SkewIStudentMvSpec")(spec, data, control, ...)
    if (p$d < 2L) stop("skewistudent-mv-o needs at least 2 dimensions.",
                       call. = FALSE)
    box <- .angleBox(p$d)
    p$nAng <- .nAngles(p$d)
    p$thLower <- box$lower
    p$thUpper <- box$upper
    p
  })

#' @describeIn validateParams Checks the angle box.
setMethod("validateParams", "SkewIStudentMvOGenSpec",
  function(spec, params, ...) {
    getMethod("validateParams", "SkewIStudentMvSpec")(spec, params, ...)
    stopifnot(params$nAng == .nAngles(params$d),
              length(params$thLower) == params$nAng,
              all(params$thUpper > params$thLower))
    invisible(TRUE)
  })

#' @describeIn simulateParams Draw components, angles uniform on the box.
setMethod("simulateParams", "SkewIStudentMvOGenSpec",
  function(spec, prior, nClust, ...) {
    p <- getMethod("simulateParams", "SkewIStudentMvSpec")(spec, prior, nClust, ...)
    p$theta <- lapply(seq_len(nClust), function(k)
      stats::runif(prior$nAng, prior$thLower, prior$thUpper))
    p
  })

#' @describeIn componentDensity Density closure, general O and per-margin nu.
setMethod("componentDensity", "SkewIStudentMvOGenSpec",
  function(spec, ...) {
    function(x, params) {
      m <- length(params[["mu"]])
      O <- orthogonalFactor(params[["theta"]], m)
      X <- if (is.matrix(x)) x else matrix(x, nrow = 1L)
      U <- chol(params[["Sigma"]]); Ui <- backsolve(U, diag(m))
      E <- (sweep(X, 2L, params[["mu"]]) %*% Ui) %*% t(O)
      g <- params[["gamma"]]; nu <- params[["nu"]]
      G <- matrix(g, nrow(E), m, byrow = TRUE)
      S <- ifelse(E < 0, E * G, E / G)
      lp <- numeric(nrow(E))
      for (j in seq_len(m))
        lp <- lp + log(2) - log(g[j] + 1 / g[j]) +
          stats::dt(S[, j], df = nu[j], log = TRUE)
      exp(lp - sum(log(diag(U))))
    }
  })

.skewMvITOGFixedKCode <- function() nimble::nimbleCode({
  for (i in 1:n) {
    z[i] ~ dcat(weights[1:K])
    muObs[i, 1:d]       <- muTilde[z[i], 1:d]
    SigObs[i, 1:d, 1:d] <- SigTilde[z[i], 1:d, 1:d]
    gamObs[i, 1:d]      <- gamTilde[z[i], 1:d]
    nuObs[i, 1:d]       <- nuTilde[z[i], 1:d]
    thObs[i, 1:nAng]    <- thetaTilde[z[i], 1:nAng]
    y[i, 1:d] ~ dSkewMvITOG_k(muObs[i, 1:d], SigObs[i, 1:d, 1:d],
                              gamObs[i, 1:d], nuObs[i, 1:d], thObs[i, 1:nAng])
  }
  weights[1:K] ~ ddirch(alphaVec[1:K])
  for (j in 1:K) {
    SigTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
    muTilde[j, 1:d] ~ dmnorm(mu0[1:d], cov = covMu[1:d, 1:d])
    for (r in 1:d) {
      gamTilde[j, r] ~ dlnorm(0, sdlog = gamLogSd)
      nuTilde[j, r]  ~ T(dgamma(shape = aNu, rate = bNu), nuLower, )
    }
    for (a in 1:nAng) thetaTilde[j, a] ~ dunif(thLower[a], thUpper[a])
  }
})

.skewMvITOGDPMCode <- function() nimble::nimbleCode({
  for (i in 1:n) {
    muObs[i, 1:d]       <- muTilde[xi[i], 1:d]
    SigObs[i, 1:d, 1:d] <- SigTilde[xi[i], 1:d, 1:d]
    gamObs[i, 1:d]      <- gamTilde[xi[i], 1:d]
    nuObs[i, 1:d]       <- nuTilde[xi[i], 1:d]
    thObs[i, 1:nAng]    <- thetaTilde[xi[i], 1:nAng]
    y[i, 1:d] ~ dSkewMvITOG_k(muObs[i, 1:d], SigObs[i, 1:d, 1:d],
                              gamObs[i, 1:d], nuObs[i, 1:d], thObs[i, 1:nAng])
  }
  xi[1:n] ~ dCRP(alpha, size = n)
  for (j in 1:L) {
    SigTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
    muTilde[j, 1:d] ~ dmnorm(mu0[1:d], cov = covMu[1:d, 1:d])
    for (r in 1:d) {
      gamTilde[j, r] ~ dlnorm(0, sdlog = gamLogSd)
      nuTilde[j, r]  ~ T(dgamma(shape = aNu, rate = bNu), nuLower, )
    }
    for (a in 1:nAng) thetaTilde[j, a] ~ dunif(thLower[a], thUpper[a])
  }
  alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
})

.skewMvITOGMc <- function(code, alloc) list(
  code = code,
  monitors = c(alloc, "muTilde", "SigTilde", "gamTilde", "nuTilde",
               "thetaTilde", if (alloc == "z") "weights" else "alpha"),
  paramNodes = c(mu = "muTilde", Sigma = "SigTilde", gamma = "gamTilde",
                 nu = "nuTilde", theta = "thetaTilde"),
  allocNode = alloc)

#' @describeIn buildModelCode General-m skew-mv-IStudent-O finite mixture.
#' @export
setMethod("buildModelCode", signature("SkewIStudentMvOGenSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) .skewMvITOGMc(.skewMvITOGFixedKCode(), "z"))

#' @describeIn buildModelCode General-m skew-mv-IStudent-O DPM mixture.
#' @export
setMethod("buildModelCode", signature("SkewIStudentMvOGenSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) .skewMvITOGMc(.skewMvITOGDPMCode(), "xi"))

#' @describeIn buildConstants General-m constants (angle box included).
setMethod("buildConstants", "SkewIStudentMvOGenSpec",
  function(spec, prior, n, ...) {
    cn <- getMethod("buildConstants", "SkewIStudentMvSpec")(spec, prior, n, ...)
    cn$nAng <- prior$nAng
    cn$thLower <- prior$thLower
    cn$thUpper <- prior$thUpper
    cn
  })

#' @describeIn componentInits k-means start; angles at the box midpoint.
setMethod("componentInits", "SkewIStudentMvOGenSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    init <- getMethod("componentInits", "SkewIStudentMvSpec")(
      spec, prior, data, count, initMethod, ...)
    mid <- (prior$thLower + prior$thUpper) / 2
    mid[abs(mid) < 1e-8] <- 0.05        # sin(0) = 0 degenerates O's first column
    init$params$thetaTilde <- matrix(rep(mid, each = count), count, prior$nAng)
    init
  })

#' @describeIn customizeSamplers Slice-sample the Householder angles.
#' @export
setMethod("customizeSamplers", "SkewIStudentMvOGenSpec",
  function(spec, conf, model, ...) {
    nodes <- model$expandNodeNames("thetaTilde")
    if (length(nodes) == 0L) return(invisible(conf))
    conf$removeSamplers("thetaTilde")
    for (nd in nodes) conf$addSampler(target = nd, type = "slice")
    invisible(conf)
  })

#' @describeIn extractParamTraces Parse mu / Sigma / gamma / nu / angle traces.
setMethod("extractParamTraces", "SkewIStudentMvOGenSpec",
  function(spec, samples, L, d = NULL, ...) {
    tr <- getMethod("extractParamTraces", "SkewIStudentMvSpec")(
      spec, samples, L, d = d, ...)
    tr$theta <- .nodeToArray(samples, "thetaTilde", c(L, .nAngles(d)))
    tr
  })

#' @describeIn relabelComponents Permute components, then canonicalise each
#'   draw's orthogonal factor, carrying gamma and nu with the permutation.
setMethod("relabelComponents", "SkewIStudentMvOGenSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...) {
    out <- getMethod("relabelComponents", "SkewIStudentMvSpec")(
      spec, paramTrace, idx, occList, perms, modalK, weights, ...)
    d <- paramTrace$d; nAng <- .nAngles(d)
    thTr <- paramTrace$theta; m <- length(idx)
    thRe <- array(NA_real_, c(m, modalK, nAng))
    for (t in seq_len(m)) {
      r <- idx[t]; occ <- occList[[t]]; ord <- perms[t, ]
      for (k in seq_len(modalK)) thRe[t, k, ] <- thTr[r, occ[ord[k]], ]
    }
    ORe   <- array(NA_real_, c(m, modalK, d, d))
    gamRe <- out$gamma
    nuRe  <- out$nu
    nCanon <- 0L
    for (t in seq_len(m)) for (k in seq_len(modalK)) {
      O <- orthogonalFactor(thRe[t, k, ], d)
      cn <- canonicaliseO(O, gamRe[t, k, ], nuRe[t, k, ])
      ORe[t, k, , ] <- cn$O
      gamRe[t, k, ] <- cn$gamma
      nuRe[t, k, ]  <- cn$nu
      nCanon <- nCanon + as.integer(cn$canonical)
    }
    gamMean <- apply(gamRe, c(2L, 3L), mean)
    gamLwr  <- apply(gamRe, c(2L, 3L), stats::quantile, 0.025, names = FALSE)
    gamUpr  <- apply(gamRe, c(2L, 3L), stats::quantile, 0.975, names = FALSE)
    nuMean  <- apply(nuRe, c(2L, 3L), mean)
    for (j in seq_len(d)) {
      out$summary[[paste0("gamma_", j, "_mean")]] <- gamMean[, j]
      out$summary[[paste0("gamma_", j, "_lwr")]]  <- gamLwr[, j]
      out$summary[[paste0("gamma_", j, "_upr")]]  <- gamUpr[, j]
      out$summary[[paste0("nu_", j, "_mean")]]    <- nuMean[, j]
    }
    out$gamma <- gamRe
    out$nu <- nuRe
    out$theta <- thRe
    out$O <- ORe
    out$O_mean <- apply(ORe, c(2L, 3L, 4L), mean)
    out$canonicalFraction <- nCanon / (m * modalK)
    out
  })
