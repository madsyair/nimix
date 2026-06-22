## ---------------------------------------------------------------------------
## dist-student-t-mv.R
##
## Multivariate Student-t component evaluated DIRECTLY as a multivariate-t
## density. NIMBLE has no built-in multivariate-t, so nimix supplies one as a
## user-defined distribution (registered with NIMBLE at load time) -- the same
## mechanism registerDistribution() exposes to package users. df is a fixed
## hyperparameter. Because the t density is not conjugate to the cluster prior,
## the cluster parameters are updated non-conjugately; the conjugate route to
## the SAME marginal is NormalGammaMvSpec.
##
## The cluster prior is the Normal-Inverse-Wishart base measure of NormalMvSpec
## (mu_k, Sigma_k), reused by inheritance; only the kernel and the predictive
## density differ.
## ---------------------------------------------------------------------------

#' Multivariate Student-t log density (nimbleFunction)
#'
#' User-defined NIMBLE distribution for a location-scale multivariate-t with
#' location \code{mu}, scale matrix \code{cov}, and degrees of freedom
#' \code{df}. Registered with NIMBLE when the package loads.
#'
#' @param x,mu Numeric vectors (observation, location).
#' @param cov Scale matrix.
#' @param df Degrees of freedom.
#' @param log Return the log density?
#' @param n Number of draws (always 1 for \code{rmvt_nimix}).
#' @return A density (or log density), or a single draw.
#' @keywords internal
#' @export
dmvt_nimix <- nimble::nimbleFunction(
  run = function(x = double(1), mu = double(1), cov = double(2),
                 df = double(0), log = integer(0, default = 0)) {
    returnType(double(0))
    d <- length(x)
    ch <- chol(cov)                                  # cov = ch^T ch (upper)
    logdet <- 2 * sum(log(diag(ch)))
    z <- forwardsolve(t(ch), x - mu)                 # solve ch^T z = x - mu
    quad <- inprod(z, z)                             # (x-mu)^T cov^{-1} (x-mu)
    ll <- lgamma((df + d) / 2) - lgamma(df / 2) -
      (d / 2) * log(df * 3.141592653589793) -
      0.5 * logdet - ((df + d) / 2) * log(1 + quad / df)
    if (log) return(ll)
    return(exp(ll))
  }
)

#' @rdname dmvt_nimix
#' @keywords internal
#' @export
rmvt_nimix <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(1), cov = double(2),
                 df = double(0)) {
    returnType(double(1))
    g <- rgamma(1, shape = df / 2, rate = df / 2)
    out <- rmnorm_chol(1, mu, chol(cov / g), prec_param = 0)
    return(out)
  }
)

#' Multivariate Student-t component specification (direct density)
#'
#' Direct multivariate-t kernel via a user-defined NIMBLE distribution. Same
#' marginal as \code{\linkS4class{NormalGammaMvSpec}} but non-conjugate. \code{df}
#' is a fixed hyperparameter.
#'
#' @slot name Fixed to \code{"student-t-mv"}.
#' @slot paramNames \code{c("mu", "Sigma")}.
#'
#' @references
#' Lange, K.L., Little, R.J.A., & Taylor, J.M.G. (1989). Robust statistical
#' modeling using the t distribution. \emph{JASA}, 84(408), 881--896.
#' \doi{10.1080/01621459.1989.10478852}
#'
#' @seealso \code{\link{NormalGammaMvSpec}} for the conjugate path.
#' @export
setClass(
  "StudentTMvSpec",
  contains = "NormalMvSpec",
  prototype = prototype(
    name = "student-t-mv",
    paramNames = c("mu", "Sigma"),
    dataDim = NA_integer_
  )
)

#' Construct a multivariate Student-t component spec
#' @return A \code{\linkS4class{StudentTMvSpec}}.
#' @examples
#' spec <- StudentTMvSpec()
#' @export
StudentTMvSpec <- function() new("StudentTMvSpec")

#' @describeIn defaultPrior Normal-Inverse-Wishart prior plus a fixed \code{df}
#'   (default 5, must exceed 2 for a finite component covariance).
#' @export
setMethod("defaultPrior", "StudentTMvSpec",
  function(spec, data, control = list(), ...) {
    base <- callNextMethod()
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
setMethod("validateParams", "StudentTMvSpec",
  function(spec, params, ...) {
    callNextMethod()
    if (is.null(params$df) || params$df <= 2)
      stop("df must be > 2 (finite covariance).", call. = FALSE)
    invisible(TRUE)
  }
)

#' @describeIn componentDensity Multivariate Student-t density.
setMethod("componentDensity", "StudentTMvSpec",
  function(spec, df = 5, ...) {
    function(x, params) {
      dfv <- if (!is.null(params[["df"]])) params[["df"]] else df
      .dmvt(as.numeric(x), params[["mu"]], params[["Sigma"]], dfv)
    }
  }
)

#' @describeIn buildModelCode Multivariate Student-t DPM code using the
#'   user-defined \code{dmvt_nimix} kernel with a Normal-Inverse-Wishart cluster
#'   prior.
#' @export
setMethod("buildModelCode", signature("StudentTMvSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        muObs[i, 1:d]       <- muTilde[xi[i], 1:d]
        covObs[i, 1:d, 1:d] <- covTilde[xi[i], 1:d, 1:d]
        y[i, 1:d] ~ dmvt_nimix(muObs[i, 1:d], covObs[i, 1:d, 1:d], df)
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
setMethod("buildConstants", "StudentTMvSpec",
  function(spec, prior, n, ...) {
    cs <- callNextMethod(); cs$df <- prior$df; cs
  }
)
