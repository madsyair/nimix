## dist-skewnormal-mv-og.R --------------------------------------------------------
## FS skew multivariate Normal with the orthogonal factor O estimated for
## GENERAL m (>= 3). The m = 2 case keeps its own, already-validated
## implementation (SkewNormalMvOSpec); nimixClust() routes on the data
## dimension, so users write distribution = "skewnormal-mv-o" either way.
##
## Sampling and identifiability. O is parameterised by m(m-1)/2 Householder
## angles (Lemma 2). The angles are sampled on the box Theta^j of FS Appendix
## A.1 -- WITHOUT imposing restriction (8) during sampling. FS state that the
## box already confines O to O_m; it does not. Sampling angles uniformly from
## the box, the fraction of draws obeying (8) is 0.245 (m = 2), 0.069 (m = 3)
## and 0.007 (m = 4). Constraining the chain to a 0.7% slice of its own prior
## would mix badly and would need an interior starting point.
##
## Instead, (8) is applied as a CANONICALISATION of the posterior draws, which
## is what it really is: among the signed row permutations P of A with
## |P| = +1, exactly one PO satisfies (8) (verified exhaustively for
## m = 2, 3, 4). The m! 2^m row ambiguity of A is label switching in the
## dimension index, and this package already prefers post-hoc relabelling to
## ordering constraints. The map leaves Sigma untouched, sends eps -> P eps,
## and adjusts gamma by gamma_i -> gamma_{perm(i)} or 1/gamma_{perm(i)}
## according to the row sign, because p(-e | gamma) = p(e | 1/gamma). The
## density is invariant under it (checked to 1e-14).
##
## Cost: canonicalisation searches m! 2^(m-1) signed permutations per draw and
## component, so it is cheap for m <= 4 and grows quickly beyond that.
##
## Reading the output, honestly.
##  * gamma and O are reported AFTER canonicalisation, so they are not directly
##    comparable to the values used to simulate data unless the simulating
##    (O, gamma) is canonicalised too. This is expected, not a defect.
##  * `O_mean` is an elementwise posterior mean. It is a convenient summary but
##    is NOT itself an orthogonal matrix (the set is not convex); use the O
##    draws for anything that needs orthogonality.
##  * The mirror modes of the angle likelihood documented for m = 2 proliferate
##    with m. Partition and location recovery are robust; individual angles and
##    the resulting O should be treated as large-sample quantities, and read
##    with wider intervals than one would allow mu or Sigma. This family is
##    marked experimental for that reason.
##  * Budget more iterations than for the fixed-O families. Each component
##    carries m(m-1)/2 slice-sampled angles, and under the MRF engine those
##    slow the Potts sweep noticeably: on an 8x8 rook grid with m = 3 and
##    heavy tails, 1200 iterations recovered the regions poorly while 3000
##    recovered them exactly. That is mixing, not misspecification.

#' @include class-DistributionSpec.R
#' @include dist-skewnormal-mv.R
#' @include skew-mv-o-general.R
NULL

# Box Theta^j (FS Appendix A.1), flattened over j = 2..m.
.angleBox <- function(m) {
  lo <- numeric(0); hi <- numeric(0)
  for (j in 2:m) {
    if (j == 2L) { lo <- c(lo, -pi / 2); hi <- c(hi, pi / 2) }
    else {
      lo <- c(lo, 0);                       hi <- c(hi, pi / 2)          # first
      if (j > 3L) { lo <- c(lo, rep(-pi / 2, j - 3L))
                    hi <- c(hi, rep(pi / 2, j - 3L)) }                   # middle
      lo <- c(lo, -pi);                     hi <- c(hi, pi)              # last
    }
  }
  list(lower = lo, upper = hi)
}

#' Skew multivariate Normal components with estimated O, general dimension
#'
#' As \code{\link{SkewNormalMvOSpec-class}} but for any \eqn{m \ge 3}: the
#' orthogonal factor of \eqn{A = OU} is parameterised by \eqn{m(m-1)/2}
#' Householder angles, sampled on the FS angle box, and FS's identifiability
#' restriction (8) is applied as a post-hoc canonicalisation of the posterior
#' draws (see \code{\link{canonicaliseO}}). Reached via
#' \code{distribution = "skewnormal-mv-o"} when the data have more than two
#' columns.
#'
#' @references Ferreira & Steel (2007), Statistica Sinica 17, 505--529.
#' @keywords internal
#' @rdname SkewNormalMvOGenSpec-class
#' @export
setClass("SkewNormalMvOGenSpec", contains = "SkewNormalMvSpec",
  prototype = prototype(name = "skewnormal-mv-og",
                        paramNames = c("mu", "Sigma", "gamma", "theta"),
                        dataDim = NA_integer_))

#' @rdname SkewNormalMvOGenSpec-class
#' @export
SkewNormalMvOGenSpec <- function() methods::new("SkewNormalMvOGenSpec")

#' @describeIn defaultPrior Adds the Householder angle box.
setMethod("defaultPrior", "SkewNormalMvOGenSpec",
  function(spec, data, control = list(), ...) {
    p <- getMethod("defaultPrior", "SkewNormalMvSpec")(spec, data, control, ...)
    if (p$d < 2L) stop("skewnormal-mv-o needs at least 2 dimensions.",
                       call. = FALSE)
    box <- .angleBox(p$d)
    p$nAng <- .nAngles(p$d)
    p$thLower <- box$lower
    p$thUpper <- box$upper
    p
  })

#' @describeIn validateParams Checks the angle box.
setMethod("validateParams", "SkewNormalMvOGenSpec",
  function(spec, params, ...) {
    getMethod("validateParams", "SkewNormalMvSpec")(spec, params, ...)
    stopifnot(params$nAng == .nAngles(params$d),
              length(params$thLower) == params$nAng,
              all(params$thUpper > params$thLower))
    invisible(TRUE)
  })

#' @describeIn simulateParams Draw components, angles uniform on the box.
setMethod("simulateParams", "SkewNormalMvOGenSpec",
  function(spec, prior, K, ...) {
    p <- getMethod("simulateParams", "SkewNormalMvSpec")(spec, prior, K, ...)
    p$theta <- lapply(seq_len(K), function(k)
      stats::runif(prior$nAng, prior$thLower, prior$thUpper))
    p
  })

#' @describeIn componentDensity Density closure with a general orthogonal factor.
setMethod("componentDensity", "SkewNormalMvOGenSpec",
  function(spec, ...) {
    function(x, params) {
      m <- length(params[["mu"]])
      O <- orthogonalFactor(params[["theta"]], m)
      X <- if (is.matrix(x)) x else matrix(x, nrow = 1L)
      U <- chol(params[["Sigma"]]); Ui <- backsolve(U, diag(m))
      E <- (sweep(X, 2L, params[["mu"]]) %*% Ui) %*% t(O)
      g <- params[["gamma"]]
      G <- matrix(g, nrow(E), m, byrow = TRUE)
      S <- ifelse(E < 0, E * G, E / G)
      exp(rowSums(matrix(log(2) - log(g + 1 / g), nrow(E), m, byrow = TRUE) +
                    stats::dnorm(S, log = TRUE)) - sum(log(diag(U))))
    }
  })

.skewMvNOGFixedKCode <- function() nimble::nimbleCode({
  for (i in 1:n) {
    z[i] ~ dcat(weights[1:K])
    muObs[i, 1:d]       <- muTilde[z[i], 1:d]
    SigObs[i, 1:d, 1:d] <- SigTilde[z[i], 1:d, 1:d]
    gamObs[i, 1:d]      <- gamTilde[z[i], 1:d]
    thObs[i, 1:nAng]    <- thetaTilde[z[i], 1:nAng]
    y[i, 1:d] ~ dSkewMvNOG_k(muObs[i, 1:d], SigObs[i, 1:d, 1:d],
                             gamObs[i, 1:d], thObs[i, 1:nAng])
  }
  weights[1:K] ~ ddirch(alphaVec[1:K])
  for (j in 1:K) {
    SigTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
    muTilde[j, 1:d] ~ dmnorm(mu0[1:d], cov = covMu[1:d, 1:d])
    for (r in 1:d) gamTilde[j, r] ~ dlnorm(0, sdlog = gamLogSd)
    for (a in 1:nAng) thetaTilde[j, a] ~ dunif(thLower[a], thUpper[a])
  }
})

.skewMvNOGDPMCode <- function() nimble::nimbleCode({
  for (i in 1:n) {
    muObs[i, 1:d]       <- muTilde[xi[i], 1:d]
    SigObs[i, 1:d, 1:d] <- SigTilde[xi[i], 1:d, 1:d]
    gamObs[i, 1:d]      <- gamTilde[xi[i], 1:d]
    thObs[i, 1:nAng]    <- thetaTilde[xi[i], 1:nAng]
    y[i, 1:d] ~ dSkewMvNOG_k(muObs[i, 1:d], SigObs[i, 1:d, 1:d],
                             gamObs[i, 1:d], thObs[i, 1:nAng])
  }
  xi[1:n] ~ dCRP(alpha, size = n)
  for (j in 1:L) {
    SigTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
    muTilde[j, 1:d] ~ dmnorm(mu0[1:d], cov = covMu[1:d, 1:d])
    for (r in 1:d) gamTilde[j, r] ~ dlnorm(0, sdlog = gamLogSd)
    for (a in 1:nAng) thetaTilde[j, a] ~ dunif(thLower[a], thUpper[a])
  }
  alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
})

.skewMvNOGMc <- function(code, alloc) list(
  code = code,
  monitors = c(alloc, "muTilde", "SigTilde", "gamTilde", "thetaTilde",
               if (alloc == "z") "weights" else "alpha"),
  paramNodes = c(mu = "muTilde", Sigma = "SigTilde", gamma = "gamTilde",
                 theta = "thetaTilde"),
  allocNode = alloc)

#' @describeIn buildModelCode General-m skew-mv-Normal-O finite mixture.
#' @export
setMethod("buildModelCode", signature("SkewNormalMvOGenSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) .skewMvNOGMc(.skewMvNOGFixedKCode(), "z"))

#' @describeIn buildModelCode General-m skew-mv-Normal-O DPM mixture.
#' @export
setMethod("buildModelCode", signature("SkewNormalMvOGenSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) .skewMvNOGMc(.skewMvNOGDPMCode(), "xi"))

#' @describeIn buildConstants General-m constants (angle box included).
setMethod("buildConstants", "SkewNormalMvOGenSpec",
  function(spec, prior, n, ...) {
    cn <- getMethod("buildConstants", "SkewNormalMvSpec")(spec, prior, n, ...)
    cn$nAng <- prior$nAng
    cn$thLower <- prior$thLower
    cn$thUpper <- prior$thUpper
    cn
  })

#' @describeIn componentInits k-means start; angles at the box midpoint.
setMethod("componentInits", "SkewNormalMvOGenSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    init <- getMethod("componentInits", "SkewNormalMvSpec")(
      spec, prior, data, count, initMethod, ...)
    mid <- (prior$thLower + prior$thUpper) / 2
    # nudge off exact zeros: sin(0) = 0 makes the first column of O degenerate
    mid[abs(mid) < 1e-8] <- 0.05
    init$params$thetaTilde <- matrix(rep(mid, each = count), count, prior$nAng)
    init
  })

#' @describeIn customizeSamplers Slice-sample the Householder angles.
#' @export
setMethod("customizeSamplers", "SkewNormalMvOGenSpec",
  function(spec, conf, model, ...) {
    nodes <- model$expandNodeNames("thetaTilde")
    if (length(nodes) == 0L) return(invisible(conf))
    conf$removeSamplers("thetaTilde")
    for (nd in nodes) conf$addSampler(target = nd, type = "slice")
    invisible(conf)
  })

#' @describeIn extractParamTraces Parse mu / Sigma / gamma / angle traces.
setMethod("extractParamTraces", "SkewNormalMvOGenSpec",
  function(spec, samples, L, d = NULL, ...) {
    tr <- getMethod("extractParamTraces", "SkewNormalMvSpec")(
      spec, samples, L, d = d, ...)
    tr$theta <- .nodeToArray(samples, "thetaTilde", c(L, .nAngles(d)))
    tr
  })

#' @describeIn relabelComponents Permute components, then canonicalise each
#'   draw's orthogonal factor via FS restriction (8), adjusting gamma with it.
setMethod("relabelComponents", "SkewNormalMvOGenSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...) {
    out <- getMethod("relabelComponents", "SkewNormalMvSpec")(
      spec, paramTrace, idx, occList, perms, modalK, weights, ...)
    d <- paramTrace$d; nAng <- .nAngles(d)
    thTr <- paramTrace$theta; m <- length(idx)
    thRe <- array(NA_real_, c(m, modalK, nAng))
    for (t in seq_len(m)) {
      r <- idx[t]; occ <- occList[[t]]; ord <- perms[t, ]
      for (k in seq_len(modalK)) thRe[t, k, ] <- thTr[r, occ[ord[k]], ]
    }
    # canonicalise: (O, gamma) -> unique representative obeying restriction (8)
    ORe   <- array(NA_real_, c(m, modalK, d, d))
    gamRe <- out$gamma
    nCanon <- 0L
    for (t in seq_len(m)) for (k in seq_len(modalK)) {
      O <- orthogonalFactor(thRe[t, k, ], d)
      cn <- canonicaliseO(O, gamRe[t, k, ])
      ORe[t, k, , ] <- cn$O
      gamRe[t, k, ] <- cn$gamma
      nCanon <- nCanon + as.integer(cn$canonical)
    }
    gamMean <- apply(gamRe, c(2L, 3L), mean)
    gamLwr  <- apply(gamRe, c(2L, 3L), stats::quantile, 0.025, names = FALSE)
    gamUpr  <- apply(gamRe, c(2L, 3L), stats::quantile, 0.975, names = FALSE)
    for (j in seq_len(d)) {
      out$summary[[paste0("gamma_", j, "_mean")]] <- gamMean[, j]
      out$summary[[paste0("gamma_", j, "_lwr")]]  <- gamLwr[, j]
      out$summary[[paste0("gamma_", j, "_upr")]]  <- gamUpr[, j]
    }
    out$gamma <- gamRe
    out$theta <- thRe
    out$O <- ORe
    out$O_mean <- apply(ORe, c(2L, 3L, 4L), mean)
    out$canonicalFraction <- nCanon / (m * modalK)
    out
  })
