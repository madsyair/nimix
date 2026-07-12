## dist-skewnormal-mv.R ----------------------------------------------------------
## Ferreira-Steel skew multivariate Normal (Ferreira & Steel 2007, Statistica
## Sinica 17:505-529). Construction: eta = A' eps + mu with independent FS-skew
## Normal components eps_j (skewness gamma_j; gamma = 1 recovers dmnorm).
##
## Design decisions (documented honestly):
##  * A = chol(Sigma), upper triangular -- i.e. the orthogonal factor O of
##    A = OU (FS Lemma 1) is fixed at the identity. FS Sec 3.3 warns that O
##    matters under skewness: fixing O = I ties the skewness directions to the
##    coordinate axes after the triangular transform. This keeps the model
##    closed-form, identifiable without restriction (8), and dCRP-compatible;
##    a Householder-parameterised O is planned as a separate increment.
##  * Parameterisation is (mu, Sigma, gamma) with Sigma ~ inverse-Wishart --
##    the same stochastic matrix node the Normal-mv DPM uses, which dCRP is
##    proven to cluster (a deterministic per-component node would break it).
##  * gamma follows the harmonised FS convention shared by all nimix skew
##    families: gamma_j > 1 skews dimension j right along its basic axis, and
##    the prior log(gamma_j) ~ N(0, s^2) treats both directions symmetrically
##    (FS use s = 1).

#' @include class-DistributionSpec.R
#' @include dist-normal-mv.R
NULL

# --- Reference R implementation (used by componentDensity, ppCheck, tests) -----

#' Ferreira-Steel skew multivariate Normal
#'
#' Density and random generation for the Ferreira-Steel skew multivariate
#' Normal with location \code{mu}, scatter \code{Sigma} (its upper-triangular
#' Cholesky factor is the FS transformation \code{A}), and per-dimension
#' skewness \code{gamma} in the shared Fernandez-Steel convention
#' (\code{gamma = 1} symmetric; \code{gamma > 1} skews right).
#'
#' @param x Numeric vector of length \code{d}, or an \code{n x d} matrix.
#' @param n Integer number of draws.
#' @param mu Numeric location vector (length \code{d}).
#' @param Sigma Positive-definite \code{d x d} scatter matrix.
#' @param gamma Positive numeric vector of per-dimension FS skewness
#'   parameters (length \code{d}).
#' @param log Logical; return the log-density?
#' @return \code{dskewmvn} a numeric vector of (log-)densities;
#'   \code{rskewmvn} an \code{n x d} matrix of draws.
#' @references
#' Ferreira, J. T. A. S. & Steel, M. F. J. (2007). A new class of skewed
#' multivariate distributions with applications to regression analysis.
#' Statistica Sinica 17, 505--529.
#' @name skewnormal-mv-distribution
NULL

#' @rdname skewnormal-mv-distribution
#' @export
dskewmvn <- function(x, mu, Sigma, gamma, log = FALSE) {
  X <- if (is.matrix(x)) x else matrix(x, nrow = 1L)
  d <- length(mu)
  stopifnot(ncol(X) == d, length(gamma) == d, all(gamma > 0))
  U  <- chol(Sigma)                       # Sigma = t(U) %*% U, A = U
  Ui <- backsolve(U, diag(d))             # U^{-1}
  ldet <- sum(log(diag(U)))
  E <- sweep(X, 2L, mu) %*% Ui            # rows: eps' = (x - mu)' A^{-1}
  G <- matrix(gamma, nrow(E), d, byrow = TRUE)
  S <- ifelse(E < 0, E * G, E / G)
  lp <- rowSums(matrix(log(2) - log(gamma + 1 / gamma), nrow(E), d,
                       byrow = TRUE) +
                stats::dnorm(S, log = TRUE)) - ldet
  if (log) lp else exp(lp)
}

#' @rdname skewnormal-mv-distribution
#' @export
rskewmvn <- function(n, mu, Sigma, gamma) {
  d <- length(mu)
  stopifnot(length(gamma) == d, all(gamma > 0))
  U <- chol(Sigma)
  G <- matrix(gamma, n, d, byrow = TRUE)
  W <- abs(matrix(stats::rnorm(n * d), n, d))
  pos <- matrix(stats::runif(n * d), n, d) < G^2 / (1 + G^2)
  Eps <- ifelse(pos, W * G, -W / G)
  sweep(Eps %*% U, 2L, mu, "+")           # eta' = eps' U + mu'  (A' eps + mu)
}

# --- S4 spec --------------------------------------------------------------------

#' Skew multivariate Normal mixture components (Ferreira-Steel)
#'
#' Multivariate component family \eqn{\eta = A^T \epsilon + \mu} with independent
#' Fernandez-Steel skew-Normal margins for \code{eps} and \code{A} the
#' upper-triangular Cholesky factor of \code{Sigma} (the orthogonal factor of
#' FS Lemma 1 is fixed at the identity; see the file header for what that
#' implies). \code{gamma = 1} in every dimension recovers the multivariate
#' Normal exactly. Non-conjugate.
#'
#' @references
#' Ferreira, J. T. A. S. & Steel, M. F. J. (2007). Statistica Sinica 17,
#' 505--529.
#' @keywords internal
#' @rdname SkewNormalMvSpec-class
#' @export
setClass("SkewNormalMvSpec", contains = "DistributionSpec",
  prototype = prototype(name = "skewnormal-mv",
                        paramNames = c("mu", "Sigma", "gamma"),
                        dataDim = NA_integer_))

#' @rdname SkewNormalMvSpec-class
#' @export
SkewNormalMvSpec <- function() methods::new("SkewNormalMvSpec")

#' @describeIn defaultPrior Data-scaled skew-mv-Normal prior.
setMethod("defaultPrior", "SkewNormalMvSpec",
  function(spec, data, control = list(), ...) {
    Y <- as.matrix(data); d <- ncol(Y)
    v <- apply(Y, 2L, stats::var); v[!is.finite(v) | v <= 0] <- 1
    cLoc <- if (!is.null(control$cLoc)) control$cLoc else 2
    df0  <- if (!is.null(control$df0)) control$df0 else d + 2
    S0   <- if (!is.null(control$S0)) control$S0 else diag(v, d) * (df0 - d - 1)
    gamLogSd <- if (!is.null(control$gamLogSd)) control$gamLogSd else 1  # FS s=1
    list(mu0 = colMeans(Y), covMu = diag(cLoc^2 * v, d),
         S0 = S0, df0 = df0, gamLogSd = gamLogSd, d = d)
  })

#' @describeIn validateParams Skew-mv-Normal hyperparameter checks.
setMethod("validateParams", "SkewNormalMvSpec",
  function(spec, params, ...) {
    d <- params$d
    stopifnot(length(params$mu0) == d, all(is.finite(params$mu0)),
              is.matrix(params$covMu), nrow(params$covMu) == d,
              is.matrix(params$S0), nrow(params$S0) == d,
              params$df0 > d + 1, params$gamLogSd > 0)
    invisible(TRUE)
  })

#' @describeIn simulateParams Draw skew-mv-Normal component parameters.
setMethod("simulateParams", "SkewNormalMvSpec",
  function(spec, prior, K, ...) {
    d <- prior$d
    list(mu = lapply(seq_len(K), function(k)
           as.numeric(prior$mu0 + chol(prior$covMu) %*% stats::rnorm(d))),
         Sigma = lapply(seq_len(K), function(k)
           prior$S0 / (prior$df0 - d - 1)),
         gamma = lapply(seq_len(K), function(k)
           stats::rlnorm(d, 0, prior$gamLogSd)))
  })

#' @describeIn componentDensity Skew-mv-Normal density closure.
setMethod("componentDensity", "SkewNormalMvSpec",
  function(spec, ...) {
    function(x, params) dskewmvn(x, params[["mu"]], params[["Sigma"]],
                                 params[["gamma"]])
  })

.skewMvNFixedKCode <- function() nimble::nimbleCode({
  for (i in 1:n) {
    z[i] ~ dcat(weights[1:K])
    muObs[i, 1:d]       <- muTilde[z[i], 1:d]
    SigObs[i, 1:d, 1:d] <- SigTilde[z[i], 1:d, 1:d]
    gamObs[i, 1:d]      <- gamTilde[z[i], 1:d]
    y[i, 1:d] ~ dSkewMvN_k(muObs[i, 1:d], SigObs[i, 1:d, 1:d], gamObs[i, 1:d])
  }
  weights[1:K] ~ ddirch(alphaVec[1:K])
  for (j in 1:K) {
    SigTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
    muTilde[j, 1:d] ~ dmnorm(mu0[1:d], cov = covMu[1:d, 1:d])
    for (r in 1:d) gamTilde[j, r] ~ dlnorm(0, sdlog = gamLogSd)
  }
})

.skewMvNDPMCode <- function() nimble::nimbleCode({
  for (i in 1:n) {
    muObs[i, 1:d]       <- muTilde[xi[i], 1:d]
    SigObs[i, 1:d, 1:d] <- SigTilde[xi[i], 1:d, 1:d]
    gamObs[i, 1:d]      <- gamTilde[xi[i], 1:d]
    y[i, 1:d] ~ dSkewMvN_k(muObs[i, 1:d], SigObs[i, 1:d, 1:d], gamObs[i, 1:d])
  }
  xi[1:n] ~ dCRP(alpha, size = n)
  for (j in 1:L) {
    SigTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
    muTilde[j, 1:d] ~ dmnorm(mu0[1:d], cov = covMu[1:d, 1:d])
    for (r in 1:d) gamTilde[j, r] ~ dlnorm(0, sdlog = gamLogSd)
  }
  alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
})

.skewMvNMc <- function(code, alloc) list(
  code = code,
  monitors = c(alloc, "muTilde", "SigTilde", "gamTilde",
               if (alloc == "z") "weights" else "alpha"),
  paramNodes = c(mu = "muTilde", Sigma = "SigTilde", gamma = "gamTilde"),
  allocNode = alloc)

#' @describeIn buildModelCode Skew-mv-Normal finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("SkewNormalMvSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) .skewMvNMc(.skewMvNFixedKCode(), "z"))

#' @describeIn buildModelCode Skew-mv-Normal DPM mixture.
#' @export
setMethod("buildModelCode", signature("SkewNormalMvSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) .skewMvNMc(.skewMvNDPMCode(), "xi"))

#' @describeIn buildConstants Skew-mv-Normal constants.
setMethod("buildConstants", "SkewNormalMvSpec",
  function(spec, prior, n, ...)
    list(n = n, d = prior$d, mu0 = prior$mu0, covMu = prior$covMu,
         S0 = prior$S0, df0 = prior$df0, gamLogSd = prior$gamLogSd))

#' @describeIn buildDataList Multivariate response matrix.
setMethod("buildDataList", "SkewNormalMvSpec",
  function(spec, data, ...) list(y = as.matrix(data)))

#' @describeIn componentInits Dispersed k-means start (gamma at 1: symmetric).
setMethod("componentInits", "SkewNormalMvSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    Y <- as.matrix(data); n <- nrow(Y); d <- ncol(Y)
    nUnique <- nrow(unique(Y))
    initRatio <- .initRatioArg(...)
    k0 <- max(1L, min(as.integer(floor(initRatio * count)),
                      as.integer(ceiling(sqrt(n)))))
    k0 <- min(k0, max(1L, nUnique))
    alloc <- rep(1L, n); centers <- matrix(colMeans(Y), 1, d)
    if (identical(initMethod, "kmeans") && k0 >= 2L && nUnique >= k0) {
      km <- tryCatch(stats::kmeans(Y, centers = k0, nstart = 5L),
                     error = function(e) NULL)
      if (!is.null(km)) { alloc <- as.integer(km$cluster); centers <- km$centers }
    }
    priorMeanCov <- prior$S0 / (prior$df0 - d - 1)
    muInit  <- matrix(rep(prior$mu0, count), count, d, byrow = TRUE)
    SigInit <- array(NA_real_, c(count, d, d))
    for (j in seq_len(count)) SigInit[j, , ] <- priorMeanCov
    occ <- sort(unique(alloc))
    for (idx in seq_along(occ)) {
      j <- occ[idx]
      if (nrow(centers) >= idx) muInit[j, ] <- centers[idx, ]
      Yj <- Y[alloc == j, , drop = FALSE]
      if (nrow(Yj) > d + 1) {
        Cj <- stats::cov(Yj)
        if (all(is.finite(Cj)) && all(eigen(Cj, TRUE, TRUE)$values > 1e-8))
          SigInit[j, , ] <- Cj
      }
    }
    list(alloc = alloc,
         params = list(muTilde = muInit, SigTilde = SigInit,
                       gamTilde = matrix(1, count, d)))
  })

#' @describeIn extractParamTraces Parse mu / Sigma / gamma traces.
setMethod("extractParamTraces", "SkewNormalMvSpec",
  function(spec, samples, L, d = NULL, ...) {
    if (is.null(d)) stop("extractParamTraces(SkewNormalMvSpec) needs 'd'.",
                         call. = FALSE)
    list(mu    = .nodeToArray(samples, "muTilde",  c(L, d)),
         Sigma = .nodeToArray(samples, "SigTilde", c(L, d, d)),
         gamma = .nodeToArray(samples, "gamTilde", c(L, d)),
         d = d)
  })

#' @describeIn relabelComponents Permute mv (mu, Sigma, gamma) and summarise.
setMethod("relabelComponents", "SkewNormalMvSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...) {
    d <- paramTrace$d
    muTr <- paramTrace$mu; SigTr <- paramTrace$Sigma; gamTr <- paramTrace$gamma
    m <- length(idx)
    muRe  <- array(NA_real_, c(m, modalK, d))
    SigRe <- array(NA_real_, c(m, modalK, d, d))
    gamRe <- array(NA_real_, c(m, modalK, d))
    for (t in seq_len(m)) {
      r <- idx[t]; occ <- occList[[t]]; ord <- perms[t, ]
      for (k in seq_len(modalK)) {
        j <- occ[ord[k]]
        muRe[t, k, ]    <- muTr[r, j, ]
        SigRe[t, k, , ] <- matrix(SigTr[r, j, , ], d, d)
        gamRe[t, k, ]   <- gamTr[r, j, ]
      }
    }
    muMean <- apply(muRe, c(2L, 3L), mean)
    muMed  <- apply(muRe, c(2L, 3L), stats::median)
    muLwr  <- apply(muRe, c(2L, 3L), stats::quantile, 0.025, names = FALSE)
    muUpr  <- apply(muRe, c(2L, 3L), stats::quantile, 0.975, names = FALSE)
    gamMean <- apply(gamRe, c(2L, 3L), mean)
    gamLwr  <- apply(gamRe, c(2L, 3L), stats::quantile, 0.025, names = FALSE)
    gamUpr  <- apply(gamRe, c(2L, 3L), stats::quantile, 0.975, names = FALSE)
    SigMean <- apply(SigRe, c(2L, 3L, 4L), mean)
    summ <- data.frame(component = seq_len(modalK), weight = colMeans(weights))
    for (j in seq_len(d)) {
      summ[[paste0("mu_", j, "_mean")]] <- muMean[, j]
      summ[[paste0("mu_", j, "_med")]]  <- muMed[, j]
      summ[[paste0("mu_", j, "_lwr")]]  <- muLwr[, j]
      summ[[paste0("mu_", j, "_upr")]]  <- muUpr[, j]
    }
    for (j in seq_len(d)) {
      summ[[paste0("gamma_", j, "_mean")]] <- gamMean[, j]
      summ[[paste0("gamma_", j, "_lwr")]]  <- gamLwr[, j]
      summ[[paste0("gamma_", j, "_upr")]]  <- gamUpr[, j]
    }
    for (j in seq_len(d))
      summ[[paste0("var_", j, "_mean")]] <-
        vapply(seq_len(modalK), function(k) SigMean[k, j, j], numeric(1))
    list(mu = muRe, Sigma = SigRe, gamma = gamRe,
         mu_mean = muMean, Sigma_mean = SigMean, summary = summ)
  })
