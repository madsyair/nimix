## ---------------------------------------------------------------------------
## dist-glm-reg.R
##
## GLM mixture-of-regression components: each mixture component is a generalized
## linear model with its own coefficient vector beta_k, a linear predictor
## eta = X beta_k, and a link function g(mu) = eta.
##   Poisson : log link,   y_i ~ Poisson(exp(X beta_k))
##   Binomial: logit link, y_i ~ Binomial(size, plogis(X beta_k))
## Unlike the Normal-linear component there is no separate dispersion s2 (the
## variance is determined by the mean). The GLM likelihood is NOT conjugate to
## the normal coefficient prior, so the cluster coefficients are updated
## non-conjugately (NIMBLE's CRP cluster wrapper handles this); we deliberately
## avoid Polya-Gamma augmentation to keep the dependency set unchanged and to
## avoid the partition-mixing penalty that latent-variable augmentation brings.
##
## The coefficient prior is a Zellner g-prior on the linear-predictor scale,
## beta_k ~ N(0, g (X'X)^{-1}), consistent with NormalRegSpec.
##
## The per-observation deterministic node betaObs[i, ] <- betaTilde[xi[i], ]
## resolves the dynamic cluster lookup before it enters inprod(), the same
## indirection NormalRegSpec uses.
## ---------------------------------------------------------------------------

# ----- shared helpers for beta-only regression specs -----------------------

.glmRegPrior <- function(y, control, family) {
  X <- control$X
  if (is.null(X) || !is.matrix(X))
    stop("GLM regression defaultPrior needs the design matrix in control$X.",
         call. = FALSE)
  n <- length(y); p <- ncol(X)
  g <- if (!is.null(control$g)) control$g else n
  ridge <- 1e-6 * mean(diag(crossprod(X)))
  B0 <- g * solve(crossprod(X) + diag(ridge, p))
  list(b0 = rep(0, p), B0 = B0, p = p, X = X, g = g,
       coefNames = if (!is.null(colnames(X))) colnames(X) else paste0("b", seq_len(p)),
       terms = control$terms)
}

# GLM mixture-of-regression start. A k-means split on the raw response cannot
# separate components that differ only by the SIGN of their slope (their
# regression lines cross and overlap), and copying one global fit to every
# cluster leaves them starting identical. Instead we (1) spread the per-cluster
# slopes around the pooled fit so the components start with distinct, opposite
# tilts, then (2) take one hard E-step: each observation is allocated to the
# component whose current coefficients give it the highest likelihood. This
# "dispersed starts + one classification step" is the standard robust
# initialisation for mixtures of regressions (Grun & Leisch 2008, FlexMix).
.glmRegInits <- function(y, X, count, family, size = NULL, initRatio = .DEFAULT_INIT_RATIO) {
  n <- length(y); p <- ncol(X)
  yfit <- if (is.null(size)) y else cbind(y, size - y)
  bGlobal <- tryCatch(stats::glm.fit(X, yfit, family = family)$coefficients,
                      error = function(e) rep(0, p))
  bGlobal[!is.finite(bGlobal)] <- 0

  # Treat a leading all-ones column as the intercept; spread the slopes only.
  hasInt   <- ncol(X) >= 1L && all(abs(X[, 1] - 1) < 1e-12)
  slopeIdx <- if (hasInt && p >= 2L) 2:p else seq_len(p)
  sScale   <- max(abs(bGlobal[slopeIdx]), 1)   # at least order 1 on link scale

  # Seed up to kUse clusters only. For the DPM, count = L = K_max is a hard
  # truncation, so capping at initRatio * count (default 0.8) keeps headroom for
  # cluster moves before the chain settles. The floor of 2 preserves the hard
  # E-step's ability to separate components for a small fixed K (e.g. K = 2,
  # where floor(0.8 * 2) = 1 would otherwise collapse the start; kept regardless of initRatio).
  kUse <- min(count, max(2L, as.integer(floor(initRatio * count))))

  betaMat  <- matrix(rep(bGlobal, each = count), nrow = count)   # count x p
  spread   <- if (kUse >= 2L) seq(-1, 1, length.out = kUse) else 0
  for (k in seq_len(kUse))
    betaMat[k, slopeIdx] <- bGlobal[slopeIdx] + spread[k] * sScale

  xiInit <- rep(1L, n)
  if (kUse >= 2L) {
    eta <- X %*% t(betaMat[seq_len(kUse), , drop = FALSE])  # n x kUse
    ll  <- matrix(0, n, kUse)
    for (k in seq_len(kUse)) {
      ll[, k] <- if (is.null(size))
        stats::dpois(y, exp(eta[, k]), log = TRUE)
      else
        stats::dbinom(y, size, stats::plogis(eta[, k]), log = TRUE)
    }
    ll[!is.finite(ll)] <- -1e6
    xiInit <- max.col(ll, ties.method = "first")
  }

  # small jitter so no two clusters are exactly equal after the E-step
  betaMat <- betaMat + matrix(stats::rnorm(count * p, 0, 0.05), count, p)
  list(alloc = as.integer(xiInit), params = list(betaTilde = betaMat))
}

.betaExtract <- function(samples, L, p, coefNames)
  list(beta = .nodeToArray(samples, "betaTilde", c(L, p)), p = p,
       coefNames = coefNames)

.betaRelabel <- function(paramTrace, idx, occList, perms, modalK, weights) {
  p <- paramTrace$p; betaTr <- paramTrace$beta; coefNm <- paramTrace$coefNames
  m <- length(idx)
  betaRe <- array(NA_real_, dim = c(m, modalK, p))
  for (t in seq_len(m)) {
    r <- idx[t]; occ <- occList[[t]]; ord <- perms[t, ]
    for (k in seq_len(modalK)) betaRe[t, k, ] <- betaTr[r, occ[ord[k]], ]
  }
  betaMean <- apply(betaRe, c(2L, 3L), mean)
  summ <- data.frame(component = seq_len(modalK), weight = colMeans(weights))
  for (j in seq_len(p)) summ[[coefNm[j]]] <- betaMean[, j]
  list(beta = betaRe, beta_mean = betaMean, coefNames = coefNm, summary = summ)
}

# ===========================================================================
# Poisson regression (log link)
# ===========================================================================

#' Poisson GLM regression component (log link)
#' @slot name Fixed to \code{"poisson-reg"}.
#' @slot paramNames \code{c("beta")}.
#' @export
setClass("PoissonRegSpec", contains = "DistributionSpec",
  prototype = prototype(name = "poisson-reg", paramNames = "beta",
                        dataDim = 1L))

#' Construct a Poisson regression component spec
#' @return A \code{\linkS4class{PoissonRegSpec}}.
#' @examples
#' spec <- PoissonRegSpec()
#' @export
PoissonRegSpec <- function() new("PoissonRegSpec")

#' @describeIn isRegressionSpec Poisson regression is a regression spec.
#' @export
setMethod("isRegressionSpec", "PoissonRegSpec", function(spec, ...) TRUE)

#' @describeIn linkInv Log link inverse (\code{exp}).
#' @export
setMethod("linkInv", "PoissonRegSpec",
  function(spec, eta, prior = NULL, ...) exp(eta))

#' @describeIn defaultPrior g-prior on the coefficients.
#' @export
setMethod("defaultPrior", "PoissonRegSpec",
  function(spec, data, control = list(), ...)
    .glmRegPrior(as.numeric(data), control, stats::poisson()))

#' @describeIn validateParams Validate the coefficient prior.
#' @export
setMethod("validateParams", "PoissonRegSpec",
  function(spec, params, ...) {
    if (is.null(params$B0) || is.null(params$p))
      stop("Poisson regression prior needs B0 and p.", call. = FALSE)
    invisible(TRUE)
  })

#' @describeIn componentDensity Poisson pmf at the fitted rate.
setMethod("componentDensity", "PoissonRegSpec",
  function(spec, ...) function(x, params)
    stats::dpois(round(params[["y"]]), lambda = exp(params[["eta"]])))

#' @describeIn buildModelCode Poisson GLM regression DPM code (log link).
#' @export
setMethod("buildModelCode", signature("PoissonRegSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        betaObs[i, 1:p] <- betaTilde[xi[i], 1:p]
        log(mu[i]) <- inprod(X[i, 1:p], betaObs[i, 1:p])
        y[i] ~ dpois(mu[i])
      }
      xi[1:n] ~ dCRP(alpha, size = n)
      for (j in 1:L) betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = B0[1:p, 1:p])
      alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
    })
    list(code = code, monitors = c("xi", "betaTilde", "alpha"),
         paramNodes = c(beta = "betaTilde"), allocNode = "xi")
  })

#' @describeIn buildModelCode Poisson GLM regression fixed-K code (log link).
#' @export
setMethod("buildModelCode", signature("PoissonRegSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        z[i] ~ dcat(weights[1:K])
        betaObs[i, 1:p] <- betaTilde[z[i], 1:p]
        log(mu[i]) <- inprod(X[i, 1:p], betaObs[i, 1:p])
        y[i] ~ dpois(mu[i])
      }
      weights[1:K] ~ ddirch(alphaVec[1:K])
      for (j in 1:K) betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = B0[1:p, 1:p])
    })
    list(code = code, monitors = c("z", "betaTilde", "weights"),
         paramNodes = c(beta = "betaTilde"), allocNode = "z")
  })

#' @describeIn buildConstants Poisson regression constants.
setMethod("buildConstants", "PoissonRegSpec",
  function(spec, prior, n, ...)
    list(n = n, p = prior$p, X = prior$X, b0 = prior$b0, B0 = prior$B0))

#' @describeIn buildDataList Response and design matrix.
setMethod("buildDataList", "PoissonRegSpec",
  function(spec, data, ...) list(y = as.numeric(data)))

#' @describeIn componentInits Global-GLM start with k-means allocation.
setMethod("componentInits", "PoissonRegSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...)
    .glmRegInits(as.numeric(data), prior$X, count, stats::poisson(),
                 initRatio = .initRatioArg(...)))

#' @describeIn extractParamTraces Parse coefficient traces.
setMethod("extractParamTraces", "PoissonRegSpec",
  function(spec, samples, L, d = NULL, prior = NULL, ...)
    .betaExtract(samples, L, prior$p, prior$coefNames))

#' @describeIn relabelComponents Permute coefficients and summarise.
setMethod("relabelComponents", "PoissonRegSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...)
    .betaRelabel(paramTrace, idx, occList, perms, modalK, weights))

# ===========================================================================
# Binomial regression (logit link)
# ===========================================================================

#' Binomial GLM regression component (logit link)
#' @slot name Fixed to \code{"binomial-reg"}.
#' @slot paramNames \code{c("beta")}.
#' @export
setClass("BinomialRegSpec", contains = "DistributionSpec",
  prototype = prototype(name = "binomial-reg", paramNames = "beta",
                        dataDim = 1L))

#' Construct a Binomial regression component spec
#' @return A \code{\linkS4class{BinomialRegSpec}}.
#' @examples
#' spec <- BinomialRegSpec()
#' @export
BinomialRegSpec <- function() new("BinomialRegSpec")

#' @describeIn isRegressionSpec Binomial regression is a regression spec.
#' @export
setMethod("isRegressionSpec", "BinomialRegSpec", function(spec, ...) TRUE)

#' @describeIn linkInv Logit link inverse (\code{size * plogis}).
#' @export
setMethod("linkInv", "BinomialRegSpec",
  function(spec, eta, prior = NULL, ...) {
    sz <- if (!is.null(prior$size)) prior$size else 1
    sz * stats::plogis(eta)
  })

#' @describeIn defaultPrior g-prior on the coefficients; needs \code{size}.
#' @export
setMethod("defaultPrior", "BinomialRegSpec",
  function(spec, data, control = list(), ...) {
    if (is.null(control$size))
      stop("Binomial regression needs the trials in prior = list(size = ...).",
           call. = FALSE)
    pr <- .glmRegPrior(as.numeric(data), control, stats::binomial())
    pr$size <- as.integer(control$size)
    pr
  })

#' @describeIn validateParams Validate the coefficient prior and \code{size}.
#' @export
setMethod("validateParams", "BinomialRegSpec",
  function(spec, params, ...) {
    if (is.null(params$B0) || is.null(params$p))
      stop("Binomial regression prior needs B0 and p.", call. = FALSE)
    if (is.null(params$size) || params$size < 1)
      stop("Binomial regression needs size >= 1.", call. = FALSE)
    invisible(TRUE)
  })

#' @describeIn componentDensity Binomial pmf at the fitted probability.
setMethod("componentDensity", "BinomialRegSpec",
  function(spec, size = NULL, ...) function(x, params) {
    sz <- if (!is.null(params[["size"]])) params[["size"]] else size
    stats::dbinom(round(params[["y"]]), size = sz,
                  prob = stats::plogis(params[["eta"]]))
  })

#' @describeIn buildModelCode Binomial GLM regression DPM code (logit link).
#' @export
setMethod("buildModelCode", signature("BinomialRegSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        betaObs[i, 1:p] <- betaTilde[xi[i], 1:p]
        logit(pp[i]) <- inprod(X[i, 1:p], betaObs[i, 1:p])
        y[i] ~ dbin(pp[i], size)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
      for (j in 1:L) betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = B0[1:p, 1:p])
      alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
    })
    list(code = code, monitors = c("xi", "betaTilde", "alpha"),
         paramNodes = c(beta = "betaTilde"), allocNode = "xi")
  })

#' @describeIn buildModelCode Binomial GLM regression fixed-K code (logit link).
#' @export
setMethod("buildModelCode", signature("BinomialRegSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        z[i] ~ dcat(weights[1:K])
        betaObs[i, 1:p] <- betaTilde[z[i], 1:p]
        logit(pp[i]) <- inprod(X[i, 1:p], betaObs[i, 1:p])
        y[i] ~ dbin(pp[i], size)
      }
      weights[1:K] ~ ddirch(alphaVec[1:K])
      for (j in 1:K) betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = B0[1:p, 1:p])
    })
    list(code = code, monitors = c("z", "betaTilde", "weights"),
         paramNodes = c(beta = "betaTilde"), allocNode = "z")
  })

#' @describeIn buildConstants Binomial regression constants plus \code{size}.
setMethod("buildConstants", "BinomialRegSpec",
  function(spec, prior, n, ...)
    list(n = n, p = prior$p, X = prior$X, b0 = prior$b0, B0 = prior$B0,
         size = prior$size))

#' @describeIn buildDataList Response and design matrix.
setMethod("buildDataList", "BinomialRegSpec",
  function(spec, data, ...) list(y = as.numeric(data)))

#' @describeIn componentInits Global-GLM start with k-means allocation.
setMethod("componentInits", "BinomialRegSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...)
    .glmRegInits(as.numeric(data), prior$X, count, stats::binomial(),
                 size = prior$size, initRatio = .initRatioArg(...)))

#' @describeIn extractParamTraces Parse coefficient traces.
setMethod("extractParamTraces", "BinomialRegSpec",
  function(spec, samples, L, d = NULL, prior = NULL, ...)
    .betaExtract(samples, L, prior$p, prior$coefNames))

#' @describeIn relabelComponents Permute coefficients and summarise.
setMethod("relabelComponents", "BinomialRegSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...)
    .betaRelabel(paramTrace, idx, occList, perms, modalK, weights))
