#' @include dist-normal-mv.R
NULL

## ---------------------------------------------------------------------------
## dist-normal-gamma-mv.R
##
## Multivariate Normal-Gamma component: the conjugate scale-mixture route to a
## multivariate Student-t component. A per-observation latent precision
## multiplier
##     omega_i ~ Gamma(df/2, df/2),   y_i ~ N_d(mu_k, Sigma_k / omega_i)
## marginalises (over omega) to y_i ~ multivariate-t(mu_k, Sigma_k, df), with df
## a fixed hyperparameter. Conditional on omega the kernel is Gaussian, so the
## Normal-Inverse-Wishart cluster prior stays conjugate and NIMBLE keeps its
## collapsed CRP cluster sampler -- the same advantage as the univariate
## Normal-Gamma over the direct multivariate-t density. The latent omega_i
## double as robustness weights.
##
## This inherits the NIW cluster prior, the cov+dinvwish parameterisation, and
## the mu/cov trace handling of NormalMvSpec; it only adds df + omega and swaps
## the predictive density to the multivariate-t marginal. As for NormalMvSpec,
## a multivariate distribution parameter may not be an expression, so the
## omega-scaled covariance is bound to a deterministic node.
##
## SAME marginal as StudentTMvSpec, reached by a cheaper (conjugate) route.
##
## References: Andrews & Mallows (1974); Backlund & Hobert (2020) for the Gibbs
## sampler under a conjugate Normal-Inverse-Wishart prior.
## ---------------------------------------------------------------------------

#' Multivariate Normal-Gamma (scale-mixture multivariate-t) component
#'
#' Conjugate scale-mixture representation of a multivariate Student-t component:
#' identical marginal to \code{\linkS4class{StudentTMvSpec}}, with conjugate
#' Normal-Inverse-Wishart cluster updates. \code{df} is a fixed hyperparameter.
#'
#' @slot name Fixed to \code{"normal-gamma-mv"}.
#' @slot paramNames \code{c("mu", "Sigma")}.
#'
#' @references
#' Andrews, D.F., & Mallows, C.L. (1974). Scale mixtures of normal
#' distributions. \emph{JRSS-B}, 36(1), 99--102.
#' \doi{10.1111/j.2517-6161.1974.tb00989.x}
#'
#' Backlund, E., & Hobert, J.P. (2020). [Gibbs sampling for multivariate linear
#' regression with errors that are scale mixtures of normals under a conjugate
#' Normal-Inverse-Wishart prior.]
#'
#' @seealso \code{\link{StudentTMvSpec}} for the direct (non-conjugate) path.
#' @export
setClass(
  "NormalGammaMvSpec",
  contains = "NormalMvSpec",
  prototype = prototype(
    name = "normal-gamma-mv",
    paramNames = c("mu", "Sigma"),
    dataDim = NA_integer_
  )
)

#' Construct a multivariate Normal-Gamma component spec
#' @return A \code{\linkS4class{NormalGammaMvSpec}}.
#' @examples
#' spec <- NormalGammaMvSpec()
#' @export
NormalGammaMvSpec <- function() new("NormalGammaMvSpec")

#' @describeIn defaultPrior Normal-Inverse-Wishart prior plus a fixed \code{df}
#'   (default 5, must exceed 2 for a finite component covariance).
#' @export
setMethod("defaultPrior", "NormalGammaMvSpec",
  function(spec, data, control = list(), ...) {
    base <- callNextMethod()                       # NIW: d, mu0, kappa0, df0, S0
    df <- if (!is.null(control$df)) control$df else 5
    if (df <= 2)
      stop("df must exceed 2 so the component has a finite covariance.",
           call. = FALSE)
    base$df <- df
    base
  }
)

#' @describeIn validateParams Validate the NIW prior and the fixed \code{df}.
#' @export
setMethod("validateParams", "NormalGammaMvSpec",
  function(spec, params, ...) {
    callNextMethod()
    if (is.null(params$df) || params$df <= 2)
      stop("df must be > 2 (finite covariance).", call. = FALSE)
    invisible(TRUE)
  }
)

#' @describeIn componentDensity Multivariate Student-t marginal density.
setMethod("componentDensity", "NormalGammaMvSpec",
  function(spec, df = 5, ...) {
    function(x, params) {
      dfv <- if (!is.null(params[["df"]])) params[["df"]] else df
      .dmvt(as.numeric(x), params[["mu"]], params[["Sigma"]], dfv)
    }
  }
)

#' @describeIn buildModelCode Multivariate Normal-Gamma scale-mixture DPM code
#'   with a per-observation latent precision multiplier \code{omega}.
#' @export
setMethod("buildModelCode", signature("NormalGammaMvSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        omega[i] ~ dgamma(shape = df / 2, rate = df / 2)
        muObs[i, 1:d]       <- muTilde[xi[i], 1:d]
        covObs[i, 1:d, 1:d] <- covTilde[xi[i], 1:d, 1:d] / omega[i]
        y[i, 1:d] ~ dmnorm(muObs[i, 1:d], cov = covObs[i, 1:d, 1:d])
      }
      xi[1:n] ~ dCRP(alpha, size = n)
      for (j in 1:L) {
        covTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
        covMu[j, 1:d, 1:d] <- covTilde[j, 1:d, 1:d] / kappa0
        muTilde[j, 1:d] ~ dmnorm(mu0[1:d], cov = covMu[j, 1:d, 1:d])
      }
      alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
    })
    list(code = code,
         monitors  = c("xi", "muTilde", "covTilde", "alpha"),
         paramNodes = c(mu = "muTilde", cov = "covTilde"),
         allocNode  = "xi")
  }
)

#' @describeIn buildConstants Normal-Inverse-Wishart constants plus \code{df}.
setMethod("buildConstants", "NormalGammaMvSpec",
  function(spec, prior, n, ...) {
    cs <- callNextMethod(); cs$df <- prior$df; cs
  }
)

#' @describeIn componentInits k-means start (inherited) plus unit \code{omega}.
setMethod("componentInits", "NormalGammaMvSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    ci <- callNextMethod()
    ci$params$omega <- rep(1, .nObs(data))
    ci
  }
)

#' @describeIn customizeSamplers Slice-sample the latent precision multipliers.
#' @export
setMethod("customizeSamplers", "NormalGammaMvSpec",
  function(spec, conf, model, ...) .omegaToSlice(conf, model))
