#' @include class-DistributionSpec.R
NULL

## ---------------------------------------------------------------------------
## dist-mv-reg.R
##
## Multivariate-RESPONSE mixture regression: the response y_i is a d-vector and
## each component is a multivariate linear regression with its own coefficient
## matrix B_k (p x d), mean X_i B_k, and a d x d error covariance Sigma_k.
##   NormalMvRegSpec      : y_i ~ N_d(X_i B_k, Sigma_k)
##   StudentTMvRegSpec    : y_i ~ multivariate-t(X_i B_k, Sigma_k, df)   (direct)
##   NormalGammaMvRegSpec : omega_i ~ Gamma(df/2, df/2),
##                          y_i ~ N_d(X_i B_k, Sigma_k / omega_i)        (mixture)
## The coefficient matrix is indexed dynamically by the cluster label and bound
## to a per-observation deterministic node before the matrix product, the same
## indirection NormalMvSpec / NormalRegSpec use.
##
## The cluster prior is the conjugate matrix-Normal-Inverse-Wishart (Backlund &
## Hobert 2020): Sigma_k ~ inverse-Wishart, and each coefficient row has prior
## covariance v0[l] * Sigma_k with a g-prior among-row scale
## v0 = g * diag((X'X)^{-1}). Conditional on Sigma the coefficients are Gaussian,
## so NIMBLE keeps the collapsed CRP cluster sampler. The Student-t direct kernel
## reuses the user-defined dmvt_nimix density (non-conjugate); the Normal-Gamma
## route keeps conjugate cluster updates and slice-samples omega.
##
## References: Zellner (1976); Fernandez & Steel (1999) (small-df pitfalls);
## Backlund & Hobert (2020) (Gibbs for multivariate regression, scale-mixture
## errors, conjugate NIW).
## ---------------------------------------------------------------------------

# ----- shared helpers ------------------------------------------------------

.mvRegPrior <- function(Y, control) {
  X <- control$X
  if (is.null(X) || !is.matrix(X))
    stop("Multivariate regression defaultPrior needs control$X.", call. = FALSE)
  Y <- as.matrix(Y); n <- nrow(Y); p <- ncol(X); d <- ncol(Y)
  ridge <- 1e-6 * mean(diag(crossprod(X)))
  Bols <- solve(crossprod(X) + diag(ridge, p), crossprod(X, Y))   # p x d OLS
  resid <- Y - X %*% Bols
  Sg <- stats::cov(resid); if (any(!is.finite(Sg))) Sg <- diag(d)
  df0 <- d + 2
  # matrix-Normal-Inverse-Wishart coefficient prior (Backlund & Hobert 2020):
  # each coefficient row l has prior covariance v0[l] * Sigma, with the among-row
  # scale a g-prior, v0 = g * diag((X'X)^{-1}); this keeps the cluster updates
  # conjugate (B | Sigma Gaussian, Sigma inverse-Wishart).
  ridge2 <- 1e-6 * mean(diag(crossprod(X)))
  g <- if (!is.null(control$g)) control$g else n
  v0 <- g * diag(solve(crossprod(X) + diag(ridge2, p)))
  list(p = p, d = d, X = X,
       mb0 = rep(0, d),
       v0 = v0,
       df0 = df0, S0 = Sg * (df0 - d - 1),
       coefNames = if (!is.null(colnames(X))) colnames(X) else paste0("b", seq_len(p)),
       respNames = if (!is.null(colnames(Y))) colnames(Y) else paste0("y", seq_len(d)),
       terms = control$terms, Bols = Bols)
}

.mvRegInits <- function(Y, prior, count, initRatio = .DEFAULT_INIT_RATIO) {
  Y <- as.matrix(Y); n <- nrow(Y); p <- prior$p; d <- prior$d
  # Dispersed k-means start, capped at initRatio * count (default 0.8) to leave
  # the cap: for the DPM, count = L = K_max is a hard truncation, and early CRP
  # sweeps can briefly occupy more clusters than the modal K before merging
  # down. Seeding right at the ceiling left no room for that transient.
  k0 <- max(1L, min(as.integer(floor(initRatio * count)), as.integer(ceiling(sqrt(n)))))
  xiInit <- rep(1L, n)
  if (k0 >= 2L) {
    km <- tryCatch(stats::kmeans(Y, centers = k0, nstart = 5L),
                   error = function(e) NULL)
    if (!is.null(km)) xiInit <- as.integer(km$cluster)
  }
  betaArr <- array(0, dim = c(count, p, d))
  for (j in seq_len(count)) betaArr[j, , ] <- prior$Bols
  covArr <- array(rep(diag(d), each = count), dim = c(count, d, d))
  list(alloc = xiInit, params = list(betaTilde = betaArr, covTilde = covArr))
}

.mvRegConstants <- function(prior, n)
  list(n = n, p = prior$p, d = prior$d, X = prior$X, mb0 = prior$mb0,
       v0 = prior$v0, df0 = prior$df0, S0 = prior$S0)

.mvRegExtract <- function(samples, L, p, d, coefNames, respNames)
  list(beta = .nodeToArray(samples, "betaTilde", c(L, p, d)),
       cov  = .nodeToArray(samples, "covTilde",  c(L, d, d)),
       p = p, d = d, coefNames = coefNames, respNames = respNames)

.mvRegRelabel <- function(paramTrace, idx, occList, perms, modalK, weights) {
  p <- paramTrace$p; d <- paramTrace$d
  betaTr <- paramTrace$beta; m <- length(idx)
  betaRe <- array(NA_real_, dim = c(m, modalK, p, d))
  for (t in seq_len(m)) {
    r <- idx[t]; occ <- occList[[t]]; ord <- perms[t, ]
    for (k in seq_len(modalK)) betaRe[t, k, , ] <- betaTr[r, occ[ord[k]], , ]
  }
  betaMean <- apply(betaRe, c(2L, 3L, 4L), mean)          # modalK x p x d
  summ <- data.frame(component = seq_len(modalK), weight = colMeans(weights))
  for (l in seq_len(p)) for (jj in seq_len(d))
    summ[[paste0(paramTrace$coefNames[l], ":", paramTrace$respNames[jj])]] <-
      betaMean[, l, jj]
  list(beta = betaRe, beta_mean = betaMean, p = p, d = d,
       coefNames = paramTrace$coefNames, respNames = paramTrace$respNames,
       summary = summ)
}

# ===========================================================================
# Normal multivariate-response regression
# ===========================================================================

#' Multivariate-response Normal mixture regression
#' @slot name Fixed to \code{"normal-mv-reg"}.
#' @export
setClass("NormalMvRegSpec", contains = "DistributionSpec",
  prototype = prototype(name = "normal-mv-reg",
                        paramNames = c("beta", "Sigma"), dataDim = NA_integer_))

#' Construct a multivariate-response Normal regression spec
#' @return A \code{\linkS4class{NormalMvRegSpec}}.
#' @examples
#' spec <- NormalMvRegSpec()
#' @export
NormalMvRegSpec <- function() new("NormalMvRegSpec")

#' @describeIn isRegressionSpec Multivariate-response regression spec.
#' @export
setMethod("isRegressionSpec", "NormalMvRegSpec", function(spec, ...) TRUE)

#' @describeIn defaultPrior Inverse-Wishart on Sigma + matrix-normal coefficients.
#' @export
setMethod("defaultPrior", "NormalMvRegSpec",
  function(spec, data, control = list(), ...) .mvRegPrior(data, control))

#' @describeIn validateParams Validate the multivariate-regression prior.
#' @export
setMethod("validateParams", "NormalMvRegSpec",
  function(spec, params, ...) {
    if (is.null(params$S0) || is.null(params$p) || is.null(params$d))
      stop("Multivariate regression prior needs S0, p, d.", call. = FALSE)
    if (params$df0 <= params$d + 1)
      stop("df0 must exceed d + 1.", call. = FALSE)
    invisible(TRUE)
  })

#' @describeIn componentDensity Multivariate normal density of the residual.
setMethod("componentDensity", "NormalMvRegSpec",
  function(spec, ...) function(x, params)
    .dmvnorm(as.numeric(params[["resid"]]), rep(0, length(params[["resid"]])),
             params[["Sigma"]]))

#' @describeIn buildModelCode Multivariate-response Normal regression DPM code.
#' @export
setMethod("buildModelCode", signature("NormalMvRegSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        betaObs[i, 1:p, 1:d] <- betaTilde[xi[i], 1:p, 1:d]
        mu[i, 1:d] <- (X[i, 1:p] %*% betaObs[i, 1:p, 1:d])[1, 1:d]
        covObs[i, 1:d, 1:d] <- covTilde[xi[i], 1:d, 1:d]
        y[i, 1:d] ~ dmnorm(mu[i, 1:d], cov = covObs[i, 1:d, 1:d])
      }
      xi[1:n] ~ dCRP(alpha, size = n)
      for (j in 1:L) {
        covTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
        for (l in 1:p) {
          covBeta[j, l, 1:d, 1:d] <- v0[l] * covTilde[j, 1:d, 1:d]
          betaTilde[j, l, 1:d] ~ dmnorm(mb0[1:d], cov = covBeta[j, l, 1:d, 1:d])
        }
      }
      alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
    })
    list(code = code, monitors = c("xi", "betaTilde", "covTilde", "alpha"),
         paramNodes = c(beta = "betaTilde", cov = "covTilde"), allocNode = "xi")
  })

#' @describeIn buildModelCode Multivariate-response Normal regression fixed-K.
#' @export
setMethod("buildModelCode", signature("NormalMvRegSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        z[i] ~ dcat(weights[1:K])
        betaObs[i, 1:p, 1:d] <- betaTilde[z[i], 1:p, 1:d]
        mu[i, 1:d] <- (X[i, 1:p] %*% betaObs[i, 1:p, 1:d])[1, 1:d]
        covObs[i, 1:d, 1:d] <- covTilde[z[i], 1:d, 1:d]
        y[i, 1:d] ~ dmnorm(mu[i, 1:d], cov = covObs[i, 1:d, 1:d])
      }
      weights[1:K] ~ ddirch(alphaVec[1:K])
      for (j in 1:K) {
        covTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
        for (l in 1:p) {
          covBeta[j, l, 1:d, 1:d] <- v0[l] * covTilde[j, 1:d, 1:d]
          betaTilde[j, l, 1:d] ~ dmnorm(mb0[1:d], cov = covBeta[j, l, 1:d, 1:d])
        }
      }
    })
    list(code = code, monitors = c("z", "betaTilde", "covTilde", "weights"),
         paramNodes = c(beta = "betaTilde", cov = "covTilde"), allocNode = "z")
  })

#' @describeIn buildConstants Multivariate-regression constants.
setMethod("buildConstants", "NormalMvRegSpec",
  function(spec, prior, n, ...) .mvRegConstants(prior, n))

#' @describeIn buildDataList Matrix response.
setMethod("buildDataList", "NormalMvRegSpec",
  function(spec, data, ...) list(y = as.matrix(data)))

#' @describeIn componentInits Global multivariate-OLS start, k-means allocation.
setMethod("componentInits", "NormalMvRegSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...)
    .mvRegInits(data, prior, count, initRatio = .initRatioArg(...)))

#' @describeIn extractParamTraces Parse coefficient and covariance traces.
setMethod("extractParamTraces", "NormalMvRegSpec",
  function(spec, samples, L, d = NULL, prior = NULL, ...)
    .mvRegExtract(samples, L, prior$p, prior$d, prior$coefNames, prior$respNames))

#' @describeIn relabelComponents Permute coefficient matrices and summarise.
setMethod("relabelComponents", "NormalMvRegSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...)
    .mvRegRelabel(paramTrace, idx, occList, perms, modalK, weights))

# ===========================================================================
# Student-t multivariate-response regression (direct multivariate-t)
# ===========================================================================

#' Multivariate-response Student-t mixture regression (direct density)
#' @slot name Fixed to \code{"student-t-mv-reg"}.
#' @export
setClass("StudentTMvRegSpec", contains = "NormalMvRegSpec",
  prototype = prototype(name = "student-t-mv-reg",
                        paramNames = c("beta", "Sigma"), dataDim = NA_integer_))

#' Construct a multivariate-response Student-t regression spec
#' @return A \code{\linkS4class{StudentTMvRegSpec}}.
#' @examples
#' spec <- StudentTMvRegSpec()
#' @export
StudentTMvRegSpec <- function() new("StudentTMvRegSpec")

#' @describeIn defaultPrior Inverse-Wishart + matrix-normal plus a fixed \code{df}.
#' @export
setMethod("defaultPrior", "StudentTMvRegSpec",
  function(spec, data, control = list(), ...) {
    base <- callNextMethod()
    df <- if (!is.null(control$df)) control$df else 5
    if (df <= 2) stop("df must exceed 2.", call. = FALSE)
    base$df <- df; base
  })

#' @describeIn buildConstants Multivariate-regression constants plus \code{df}.
setMethod("buildConstants", "StudentTMvRegSpec",
  function(spec, prior, n, ...) { cs <- callNextMethod(); cs$df <- prior$df; cs })

#' @describeIn buildModelCode Multivariate-response Student-t regression DPM code
#'   (direct multivariate-t kernel).
#' @export
setMethod("buildModelCode", signature("StudentTMvRegSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        betaObs[i, 1:p, 1:d] <- betaTilde[xi[i], 1:p, 1:d]
        mu[i, 1:d] <- (X[i, 1:p] %*% betaObs[i, 1:p, 1:d])[1, 1:d]
        covObs[i, 1:d, 1:d] <- covTilde[xi[i], 1:d, 1:d]
        y[i, 1:d] ~ dmvt_nimix(mu[i, 1:d], covObs[i, 1:d, 1:d], df)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
      for (j in 1:L) {
        covTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
        for (l in 1:p) {
          covBeta[j, l, 1:d, 1:d] <- v0[l] * covTilde[j, 1:d, 1:d]
          betaTilde[j, l, 1:d] ~ dmnorm(mb0[1:d], cov = covBeta[j, l, 1:d, 1:d])
        }
      }
      alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
    })
    list(code = code, monitors = c("xi", "betaTilde", "covTilde", "alpha"),
         paramNodes = c(beta = "betaTilde", cov = "covTilde"), allocNode = "xi")
  })

#' @describeIn buildModelCode Multivariate-response Student-t regression fixed-K.
#' @export
setMethod("buildModelCode", signature("StudentTMvRegSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        z[i] ~ dcat(weights[1:K])
        betaObs[i, 1:p, 1:d] <- betaTilde[z[i], 1:p, 1:d]
        mu[i, 1:d] <- (X[i, 1:p] %*% betaObs[i, 1:p, 1:d])[1, 1:d]
        covObs[i, 1:d, 1:d] <- covTilde[z[i], 1:d, 1:d]
        y[i, 1:d] ~ dmvt_nimix(mu[i, 1:d], covObs[i, 1:d, 1:d], df)
      }
      weights[1:K] ~ ddirch(alphaVec[1:K])
      for (j in 1:K) {
        covTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
        for (l in 1:p) {
          covBeta[j, l, 1:d, 1:d] <- v0[l] * covTilde[j, 1:d, 1:d]
          betaTilde[j, l, 1:d] ~ dmnorm(mb0[1:d], cov = covBeta[j, l, 1:d, 1:d])
        }
      }
    })
    list(code = code, monitors = c("z", "betaTilde", "covTilde", "weights"),
         paramNodes = c(beta = "betaTilde", cov = "covTilde"), allocNode = "z")
  })

# ===========================================================================
# Normal-Gamma multivariate-response regression (conjugate scale mixture)
# ===========================================================================

#' Multivariate-response Normal-Gamma mixture regression (scale mixture)
#' @slot name Fixed to \code{"normal-gamma-mv-reg"}.
#' @export
setClass("NormalGammaMvRegSpec", contains = "NormalMvRegSpec",
  prototype = prototype(name = "normal-gamma-mv-reg",
                        paramNames = c("beta", "Sigma"), dataDim = NA_integer_))

#' Construct a multivariate-response Normal-Gamma regression spec
#' @return A \code{\linkS4class{NormalGammaMvRegSpec}}.
#' @examples
#' spec <- NormalGammaMvRegSpec()
#' @export
NormalGammaMvRegSpec <- function() new("NormalGammaMvRegSpec")

#' @describeIn defaultPrior Inverse-Wishart + matrix-normal plus a fixed \code{df}.
#' @export
setMethod("defaultPrior", "NormalGammaMvRegSpec",
  function(spec, data, control = list(), ...) {
    base <- callNextMethod()
    df <- if (!is.null(control$df)) control$df else 5
    if (df <= 2) stop("df must exceed 2.", call. = FALSE)
    base$df <- df; base
  })

#' @describeIn buildConstants Multivariate-regression constants plus \code{df}.
setMethod("buildConstants", "NormalGammaMvRegSpec",
  function(spec, prior, n, ...) { cs <- callNextMethod(); cs$df <- prior$df; cs })

#' @describeIn componentInits Global multivariate-OLS start plus unit \code{omega}.
setMethod("componentInits", "NormalGammaMvRegSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    ci <- callNextMethod(); ci$params$omega <- rep(1, nrow(as.matrix(data)))
    ci
  })

#' @describeIn customizeSamplers Slice-sample the latent precision multipliers.
#' @export
setMethod("customizeSamplers", "NormalGammaMvRegSpec",
  function(spec, conf, model, ...) .omegaToSlice(conf, model))

#' @describeIn buildModelCode Multivariate-response Normal-Gamma regression DPM
#'   code (scale mixture; conjugate cluster updates).
#' @export
setMethod("buildModelCode", signature("NormalGammaMvRegSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        omega[i] ~ dgamma(shape = df / 2, rate = df / 2)
        betaObs[i, 1:p, 1:d] <- betaTilde[xi[i], 1:p, 1:d]
        mu[i, 1:d] <- (X[i, 1:p] %*% betaObs[i, 1:p, 1:d])[1, 1:d]
        covObs[i, 1:d, 1:d] <- covTilde[xi[i], 1:d, 1:d] / omega[i]
        y[i, 1:d] ~ dmnorm(mu[i, 1:d], cov = covObs[i, 1:d, 1:d])
      }
      xi[1:n] ~ dCRP(alpha, size = n)
      for (j in 1:L) {
        covTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
        for (l in 1:p) {
          covBeta[j, l, 1:d, 1:d] <- v0[l] * covTilde[j, 1:d, 1:d]
          betaTilde[j, l, 1:d] ~ dmnorm(mb0[1:d], cov = covBeta[j, l, 1:d, 1:d])
        }
      }
      alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
    })
    list(code = code,
         monitors = c("xi", "betaTilde", "covTilde", "alpha"),
         paramNodes = c(beta = "betaTilde", cov = "covTilde"), allocNode = "xi")
  })

#' @describeIn buildModelCode Multivariate-response Normal-Gamma regression
#'   fixed-K code.
#' @export
setMethod("buildModelCode", signature("NormalGammaMvRegSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        z[i] ~ dcat(weights[1:K])
        omega[i] ~ dgamma(shape = df / 2, rate = df / 2)
        betaObs[i, 1:p, 1:d] <- betaTilde[z[i], 1:p, 1:d]
        mu[i, 1:d] <- (X[i, 1:p] %*% betaObs[i, 1:p, 1:d])[1, 1:d]
        covObs[i, 1:d, 1:d] <- covTilde[z[i], 1:d, 1:d] / omega[i]
        y[i, 1:d] ~ dmnorm(mu[i, 1:d], cov = covObs[i, 1:d, 1:d])
      }
      weights[1:K] ~ ddirch(alphaVec[1:K])
      for (j in 1:K) {
        covTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
        for (l in 1:p) {
          covBeta[j, l, 1:d, 1:d] <- v0[l] * covTilde[j, 1:d, 1:d]
          betaTilde[j, l, 1:d] ~ dmnorm(mb0[1:d], cov = covBeta[j, l, 1:d, 1:d])
        }
      }
    })
    list(code = code,
         monitors = c("z", "betaTilde", "covTilde", "weights"),
         paramNodes = c(beta = "betaTilde", cov = "covTilde"), allocNode = "z")
  })
