#' @include dist-normal-reg.R
NULL

## ---------------------------------------------------------------------------
## dist-heavytail-reg.R
##
## Univariate heavy-tailed mixture regression: each component is a linear
## regression with Student-t errors, in the two equivalent parameterisations
## already used for clustering.
##   StudentTRegSpec   : y_i ~ t(x_i' beta_k, tau_k, df)               (direct)
##   NormalGammaRegSpec: omega_i ~ Gamma(df/2, df/2),
##                       y_i ~ N(x_i' beta_k, s2_k / omega_i)          (mixture)
## Both marginalise to the same Student-t regression. The Normal-Gamma route
## keeps the conjugate Normal-Inverse-Gamma g-prior on (beta, s2) conditional on
## omega -- a pure Gibbs sampler (Geweke 1993) -- and exposes per-observation
## robustness weights; the direct route is non-conjugate but stores no omega.
## df is a fixed hyperparameter (> 2).
##
## Both inherit the NIG g-prior, the coefficient/scale trace handling and the
## k-means start of NormalRegSpec; only the kernel (and, for Normal-Gamma, the
## omega augmentation + slice sampler) differ.
##
## References: Geweke (1993) for the Student-t linear model via Gibbs
## data-augmentation; West (1987), Lange, Little & Taylor (1989) for the
## scale-mixture / robustness-weight interpretation.
## ---------------------------------------------------------------------------

# ===========================================================================
# Student-t regression (direct t density)
# ===========================================================================

#' Student-t mixture regression component (direct t density)
#' @slot name Fixed to \code{"student-t-reg"}.
#' @export
setClass("StudentTRegSpec", contains = "NormalRegSpec",
  prototype = prototype(name = "student-t-reg",
                        paramNames = c("beta", "s2"), dataDim = 1L))

#' Construct a Student-t regression component spec
#' @return A \code{\linkS4class{StudentTRegSpec}}.
#' @examples
#' spec <- StudentTRegSpec()
#' @export
StudentTRegSpec <- function() new("StudentTRegSpec")

#' @describeIn defaultPrior NIG g-prior plus a fixed \code{df} (default 4, > 2).
#' @export
setMethod("defaultPrior", "StudentTRegSpec",
  function(spec, data, control = list(), ...) {
    base <- callNextMethod()
    df <- if (!is.null(control$df)) control$df else 4
    if (df <= 2) stop("df must exceed 2 for a finite error variance.",
                      call. = FALSE)
    base$df <- df
    base
  })

#' @describeIn buildConstants NIG g-prior constants plus \code{df}.
setMethod("buildConstants", "StudentTRegSpec",
  function(spec, prior, n, ...) { cs <- callNextMethod(); cs$df <- prior$df; cs })

#' @describeIn buildModelCode Student-t regression DPM code (direct t density;
#'   the scale enters as precision \code{tau = 1 / s2}).
#' @export
setMethod("buildModelCode", signature("StudentTRegSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        betaObs[i, 1:p] <- betaTilde[xi[i], 1:p]
        mu[i] <- inprod(X[i, 1:p], betaObs[i, 1:p])
        tauObs[i] <- 1 / s2Tilde[xi[i]]
        y[i] ~ dt(mu[i], tauObs[i], df)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
      for (j in 1:L) {
        s2Tilde[j] ~ dinvgamma(shape = nu0, scale = s0)
        covBeta[j, 1:p, 1:p] <- s2Tilde[j] * B0[1:p, 1:p]
        betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = covBeta[j, 1:p, 1:p])
      }
      alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
    })
    list(code = code, monitors = c("xi", "betaTilde", "s2Tilde", "alpha"),
         paramNodes = c(beta = "betaTilde", s2 = "s2Tilde"), allocNode = "xi")
  })

#' @describeIn buildModelCode Student-t regression fixed-K code.
#' @export
setMethod("buildModelCode", signature("StudentTRegSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        z[i] ~ dcat(weights[1:K])
        betaObs[i, 1:p] <- betaTilde[z[i], 1:p]
        mu[i] <- inprod(X[i, 1:p], betaObs[i, 1:p])
        tauObs[i] <- 1 / s2Tilde[z[i]]
        y[i] ~ dt(mu[i], tauObs[i], df)
      }
      weights[1:K] ~ ddirch(alphaVec[1:K])
      for (j in 1:K) {
        s2Tilde[j] ~ dinvgamma(shape = nu0, scale = s0)
        covBeta[j, 1:p, 1:p] <- s2Tilde[j] * B0[1:p, 1:p]
        betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = covBeta[j, 1:p, 1:p])
      }
    })
    list(code = code, monitors = c("z", "betaTilde", "s2Tilde", "weights"),
         paramNodes = c(beta = "betaTilde", s2 = "s2Tilde"), allocNode = "z")
  })

# ===========================================================================
# Normal-Gamma regression (conjugate scale-mixture)
# ===========================================================================

#' Normal-Gamma mixture regression component (conjugate scale-mixture)
#' @slot name Fixed to \code{"normal-gamma-reg"}.
#' @export
setClass("NormalGammaRegSpec", contains = "NormalRegSpec",
  prototype = prototype(name = "normal-gamma-reg",
                        paramNames = c("beta", "s2"), dataDim = 1L))

#' Construct a Normal-Gamma regression component spec
#' @return A \code{\linkS4class{NormalGammaRegSpec}}.
#' @examples
#' spec <- NormalGammaRegSpec()
#' @export
NormalGammaRegSpec <- function() new("NormalGammaRegSpec")

#' @describeIn defaultPrior NIG g-prior plus a fixed \code{df} (default 4, > 2).
#' @export
setMethod("defaultPrior", "NormalGammaRegSpec",
  function(spec, data, control = list(), ...) {
    base <- callNextMethod()
    df <- if (!is.null(control$df)) control$df else 4
    if (df <= 2) stop("df must exceed 2 for a finite error variance.",
                      call. = FALSE)
    base$df <- df
    base
  })

#' @describeIn buildConstants NIG g-prior constants plus \code{df}.
setMethod("buildConstants", "NormalGammaRegSpec",
  function(spec, prior, n, ...) { cs <- callNextMethod(); cs$df <- prior$df; cs })

#' @describeIn componentInits k-means start (inherited) plus unit \code{omega}.
setMethod("componentInits", "NormalGammaRegSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    ci <- callNextMethod(); ci$params$omega <- rep(1, length(as.numeric(data)))
    ci
  })

#' @describeIn customizeSamplers Slice-sample the latent precision multipliers.
#' @export
setMethod("customizeSamplers", "NormalGammaRegSpec",
  function(spec, conf, model, ...) .omegaToSlice(conf, model))

#' @describeIn buildModelCode Normal-Gamma regression DPM code (scale mixture;
#'   conjugate coefficients conditional on \code{omega}).
#' @export
setMethod("buildModelCode", signature("NormalGammaRegSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        omega[i] ~ dgamma(shape = df / 2, rate = df / 2)
        betaObs[i, 1:p] <- betaTilde[xi[i], 1:p]
        mu[i] <- inprod(X[i, 1:p], betaObs[i, 1:p])
        s2Obs[i] <- s2Tilde[xi[i]] / omega[i]
        y[i] ~ dnorm(mu[i], var = s2Obs[i])
      }
      xi[1:n] ~ dCRP(alpha, size = n)
      for (j in 1:L) {
        s2Tilde[j] ~ dinvgamma(shape = nu0, scale = s0)
        covBeta[j, 1:p, 1:p] <- s2Tilde[j] * B0[1:p, 1:p]
        betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = covBeta[j, 1:p, 1:p])
      }
      alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
    })
    list(code = code, monitors = c("xi", "betaTilde", "s2Tilde", "alpha"),
         paramNodes = c(beta = "betaTilde", s2 = "s2Tilde"), allocNode = "xi")
  })

#' @describeIn buildModelCode Normal-Gamma regression fixed-K code.
#' @export
setMethod("buildModelCode", signature("NormalGammaRegSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        z[i] ~ dcat(weights[1:K])
        omega[i] ~ dgamma(shape = df / 2, rate = df / 2)
        betaObs[i, 1:p] <- betaTilde[z[i], 1:p]
        mu[i] <- inprod(X[i, 1:p], betaObs[i, 1:p])
        s2Obs[i] <- s2Tilde[z[i]] / omega[i]
        y[i] ~ dnorm(mu[i], var = s2Obs[i])
      }
      weights[1:K] ~ ddirch(alphaVec[1:K])
      for (j in 1:K) {
        s2Tilde[j] ~ dinvgamma(shape = nu0, scale = s0)
        covBeta[j, 1:p, 1:p] <- s2Tilde[j] * B0[1:p, 1:p]
        betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = covBeta[j, 1:p, 1:p])
      }
    })
    list(code = code, monitors = c("z", "betaTilde", "s2Tilde", "weights"),
         paramNodes = c(beta = "betaTilde", s2 = "s2Tilde"), allocNode = "z")
  })
