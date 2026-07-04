#' @include dist-normal-uv.R
NULL

## ---------------------------------------------------------------------------
## dist-normal-gamma.R
##
## Univariate Normal-Gamma component: a conjugate scale-mixture representation of
## a Student-t component. A per-observation latent precision multiplier
##     omega_i ~ Gamma(df/2, df/2),   y_i ~ N(mu_k, s2_k / omega_i)
## gives, after integrating omega out, exactly y_i ~ t(mu_k, s2_k, df). df is a
## fixed hyperparameter. The advantage over the direct Student-t density is
## cost: conditional on omega the kernel is Gaussian, so the cluster parameters
## (mu_k, s2_k) keep the conjugate Normal-Inverse-Gamma updates that NIMBLE
## samples with its collapsed CRP sampler, instead of the non-conjugate update a
## raw t density needs. The latent omega_i double as robustness weights -- a
## small posterior omega_i flags an observation the model treats as an outlier.
##
## This shares the Normal-Inverse-Gamma cluster prior and the mu/s2 trace
## handling of the plain Gaussian component, so it inherits from NormalUvSpec
## and only overrides the model code (adding omega), the constants (adding df),
## the initial values (adding omega), and the predictive density (which is the
## Student-t marginal, not Gaussian).
##
## IMPORTANT: this is the SAME marginal distribution as StudentTUvSpec, reached
## by a different (cheaper) sampling route. It is also unrelated to the
## "Normal-Gamma" shrinkage prior on mixture means; that is a different idea.
##
## References: Andrews & Mallows (1974); West (1987); Lange, Little & Taylor
## (1989); Geweke (1993).
## ---------------------------------------------------------------------------

#' Univariate Normal-Gamma (scale-mixture Student-t) component specification
#'
#' A conjugate scale-mixture representation of a univariate Student-t component:
#' identical marginal to \code{\linkS4class{StudentTUvSpec}}, but with conjugate
#' cluster updates because the kernel is Gaussian conditional on a latent
#' per-observation precision multiplier. The degrees of freedom \code{df} are a
#' fixed hyperparameter.
#'
#' @slot name Fixed to \code{"normal-gamma"}.
#' @slot paramNames \code{c("mu", "s2")}.
#'
#' @references
#' Andrews, D.F., & Mallows, C.L. (1974). Scale mixtures of normal
#' distributions. \emph{JRSS-B}, 36(1), 99--102.
#' \doi{10.1111/j.2517-6161.1974.tb00989.x}
#'
#' West, M. (1987). On scale mixtures of normal distributions.
#' \emph{Biometrika}, 74(3), 646--648. \doi{10.1093/biomet/74.3.646}
#'
#' Lange, K.L., Little, R.J.A., & Taylor, J.M.G. (1989). Robust statistical
#' modeling using the t distribution. \emph{JASA}, 84(408), 881--896.
#' \doi{10.1080/01621459.1989.10478852}
#'
#' @seealso \code{\link{StudentTUvSpec}} for the direct (non-conjugate) path to
#'   the same marginal.
#' @export
setClass(
  "NormalGammaUvSpec",
  contains = "NormalUvSpec",
  prototype = prototype(
    name = "normal-gamma",
    paramNames = c("mu", "s2"),
    dataDim = 1L
  )
)

#' Construct a univariate Normal-Gamma component spec
#' @return A \code{\linkS4class{NormalGammaUvSpec}}.
#' @examples
#' spec <- NormalGammaUvSpec()
#' @export
NormalGammaUvSpec <- function() new("NormalGammaUvSpec")

# --- defaultPrior: inherit the NIG prior, add fixed df ---------------------

#' @describeIn defaultPrior Normal-Inverse-Gamma prior plus a fixed \code{df}
#'   (degrees of freedom, default 4, must exceed 2 for a finite component
#'   variance).
#' @export
setMethod("defaultPrior", "NormalGammaUvSpec",
  function(spec, data, control = list(), ...) {
    base <- callNextMethod()                       # NIG: mu0, kappa0, nu0, s0
    df <- if (!is.null(control$df)) control$df else 4
    if (df <= 2)
      stop("df must exceed 2 so the component has a finite variance.",
           call. = FALSE)
    base$df <- df
    base
  }
)

#' @describeIn validateParams Validate the NIG prior and the fixed \code{df}.
#' @export
setMethod("validateParams", "NormalGammaUvSpec",
  function(spec, params, ...) {
    callNextMethod()
    if (is.null(params$df) || params$df <= 2)
      stop("df must be > 2 (finite variance).", call. = FALSE)
    invisible(TRUE)
  }
)

# --- componentDensity: the Student-t MARGINAL ------------------------------

#' @describeIn componentDensity Student-t marginal density (location \code{mu},
#'   scale \eqn{\sqrt{s2}}, \code{df}).
setMethod("componentDensity", "NormalGammaUvSpec",
  function(spec, df = 4, ...) {
    function(x, params) {
      sigma <- sqrt(params[["s2"]])
      dfv <- if (!is.null(params[["df"]])) params[["df"]] else df
      stats::dt((x - params[["mu"]]) / sigma, df = dfv) / sigma
    }
  }
)

# --- buildModelCode: NormalGammaUvSpec x DPMEngine -------------------------

#' @describeIn buildModelCode Normal-Gamma scale-mixture DPM model code (dCRP)
#'   with a per-observation latent precision multiplier \code{omega}.
#' @export
setMethod("buildModelCode", signature("NormalGammaUvSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        omega[i] ~ dgamma(shape = df / 2, rate = df / 2)
        y[i] ~ dnorm(muTilde[xi[i]], var = s2Tilde[xi[i]] / omega[i])
      }
      xi[1:n] ~ dCRP(alpha, size = n)
      for (j in 1:L) {
        s2Tilde[j] ~ dinvgamma(shape = nu0, scale = s0)
        muTilde[j] ~ dnorm(mu0, var = s2Tilde[j] / kappa0)
      }
      alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
    })
    list(code = code,
         monitors  = c("xi", "muTilde", "s2Tilde", "alpha"),
         paramNodes = c(mu = "muTilde", s2 = "s2Tilde"),
         allocNode  = "xi")
  }
)

# --- Engine-facing overrides: add df to constants, omega to inits ----------

#' @describeIn buildConstants Normal-Inverse-Gamma constants plus \code{df}.
setMethod("buildConstants", "NormalGammaUvSpec",
  function(spec, prior, n, ...) {
    cs <- callNextMethod()
    cs$df <- prior$df
    cs
  }
)

#' @describeIn componentInits k-means start (inherited) plus unit \code{omega}.
setMethod("componentInits", "NormalGammaUvSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    ci <- callNextMethod()
    ci$params$omega <- rep(1, .nObs(data))
    ci
  }
)

#' @describeIn customizeSamplers Slice-sample the latent precision multipliers,
#'   which mixes the partition better than the default random walk.
#' @export
setMethod("customizeSamplers", "NormalGammaUvSpec",
  function(spec, conf, model, ...) .omegaToSlice(conf, model))
