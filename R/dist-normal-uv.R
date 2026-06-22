## ---------------------------------------------------------------------------
## dist-normal-uv.R
##
## Univariate Gaussian component (the v0.1.0 starting distribution).
##
## Cluster prior is a conjugate Normal-Inverse-Gamma:
##     s2_j      ~ InvGamma(shape = nu0, scale = s0)
##     mu_j | s2 ~ Normal(mu0, var = s2_j / kappa0)
## This conjugacy lets NIMBLE assign efficient conjugate CRP samplers
## (exploit conjugacy for good mixing) and is
## one of the relationships NIMBLE recognises for dCRP with both mean and
## variance unknown.
##
## Defaults are DATA-SCALED, not vague: mu0 = mean(y),
## prior sd of mu ~= c * sd(y), and prior mean of s2 ~= var(y).
## ---------------------------------------------------------------------------

#' Univariate Gaussian component specification
#'
#' @slot name Fixed to \code{"normal-uv"}.
#' @slot paramNames \code{c("mu", "s2")}.
#' @slot priorSpec Filled by \code{\link{defaultPrior}}.
#'
#' @references
#' Frühwirth-Schnatter, S. (2006). \emph{Finite Mixture and Markov Switching
#' Models}. Springer. \doi{10.1007/978-0-387-35768-3}
#'
#' Neal, R.M. (2000). Markov chain sampling methods for Dirichlet process
#' mixture models. \emph{JCGS}, 9(2), 249--265.
#' \doi{10.1080/10618600.2000.10474879}
#'
#' @seealso \code{\link{nimixClust}}
#' @export
setClass(
  "NormalUvSpec",
  contains = "DistributionSpec",
  prototype = prototype(
    name = "normal-uv",
    paramNames = c("mu", "s2"),
    dataDim = 1L
  )
)

#' Construct a univariate Gaussian component spec
#' @return A \code{\linkS4class{NormalUvSpec}}.
#' @examples
#' spec <- NormalUvSpec()
#' @export
NormalUvSpec <- function() new("NormalUvSpec")

# --- defaultPrior ----------------------------------------------------------

#' @describeIn defaultPrior Data-scaled Normal-Inverse-Gamma prior.
#'
#' Control overrides: \code{cLoc} (location spread multiplier; prior sd of
#' \code{mu} ~ \code{cLoc * sd(data)}, default 2), \code{nu0} (InvGamma shape,
#' default 3, must exceed 2 for a finite prior variance), and \code{concPrior}.
#' @export
setMethod("defaultPrior", "NormalUvSpec",
  function(spec, data, control = list(), ...) {
    data <- as.numeric(data)
    vy <- stats::var(data)
    if (!is.finite(vy) || vy <= 0) vy <- 1
    cLoc <- if (!is.null(control$cLoc)) control$cLoc else 2
    nu0  <- if (!is.null(control$nu0))  control$nu0  else 3
    if (nu0 <= 2)
      stop("nu0 must exceed 2 so the prior variance of s2 is finite ",
 "(avoid singular/degenerate priors).",
           call. = FALSE)
    # var(mu | s2) = s2 / kappa0; target prior sd of mu ~= cLoc * sd(y) at
    # s2 ~= var(y)  ==>  kappa0 = 1 / cLoc^2.
    kappa0 <- 1 / (cLoc^2)
    # InvGamma(nu0, s0) has mean s0 / (nu0 - 1); target mean = var(y).
    s0 <- vy * (nu0 - 1)
    list(
      mu0    = mean(data),
      kappa0 = kappa0,
      nu0    = nu0,
      s0     = s0,
      cLoc   = cLoc
    )
  }
)

# --- validateParams --------------------------------------------------------

#' @describeIn validateParams Validate a Normal-Inverse-Gamma prior list.
#' @export
setMethod("validateParams", "NormalUvSpec",
  function(spec, params, ...) {
    req <- c("mu0", "kappa0", "nu0", "s0")
    miss <- setdiff(req, names(params))
    if (length(miss))
      stop("NormalUvSpec prior is missing: ", paste(miss, collapse = ", "),
           call. = FALSE)
    if (params$kappa0 <= 0) stop("kappa0 must be > 0.", call. = FALSE)
    if (params$nu0 <= 2)
 stop("nu0 must be > 2 (finite prior variance).",
           call. = FALSE)
    if (params$s0 <= 0) stop("s0 (InvGamma scale) must be > 0.", call. = FALSE)
    invisible(TRUE)
  }
)

# --- simulateParams (for inits / recovery) ---------------------------------

#' @describeIn simulateParams Draw (mu, s2) from the NIG prior.
setMethod("simulateParams", "NormalUvSpec",
  function(spec, prior, nClust, ...) {
    s2 <- 1 / stats::rgamma(nClust, shape = prior$nu0, rate = prior$s0)
    mu <- stats::rnorm(nClust, mean = prior$mu0,
                       sd = sqrt(s2 / prior$kappa0))
    list(mu = mu, s2 = s2)
  }
)

# --- componentDensity ------------------------------------------------------

#' @describeIn componentDensity Gaussian density for posterior predictive use.
setMethod("componentDensity", "NormalUvSpec",
  function(spec, ...) {
    function(x, params) stats::dnorm(x, mean = params[["mu"]],
                                     sd = sqrt(params[["s2"]]))
  }
)

# --- buildModelCode: NormalUvSpec x DPMEngine ------------------------------
# This is the extensibility seam: a new
# (distribution, engine) pairing is one method, not a rewrite.

#' @describeIn buildModelCode Univariate Gaussian DPM model code (dCRP).
#'
#' Builds the NIMBLE code for
#' \deqn{y_i \sim N(\mu_{\xi_i}, s^2_{\xi_i}), \quad
#'       \xi_{1:n} \sim CRP(\alpha, n),}
#' with a conjugate Normal-Inverse-Gamma cluster prior and a Gamma hyperprior
#' on \eqn{\alpha}. Cluster-parameter vectors have length \code{L = K_max}
#' (NIMBLE's exact truncation; the sampler stays proper as long as the number
#' of occupied clusters is strictly below \code{L}).
#'
#' @references
#' Neal, R.M. (2000) \doi{10.1080/10618600.2000.10474879};
#' de Valpine et al. (2017) \doi{10.1080/10618600.2016.1172487}.
#' @export
setMethod("buildModelCode", signature("NormalUvSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        y[i] ~ dnorm(muTilde[xi[i]], var = s2Tilde[xi[i]])
      }
      xi[1:n] ~ dCRP(alpha, size = n)
      for (j in 1:L) {
        s2Tilde[j] ~ dinvgamma(shape = nu0, scale = s0)
        muTilde[j] ~ dnorm(mu0, var = s2Tilde[j] / kappa0)
      }
      alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
    })
    list(
      code      = code,
      monitors  = c("xi", "muTilde", "s2Tilde", "alpha"),
      paramNodes = c(mu = "muTilde", s2 = "s2Tilde"),
      allocNode  = "xi"
    )
  }
)

# --- Engine-facing methods (dimension-agnostic DPM orchestration) ----------
# These move the univariate-specific bits of constants/inits/trace-parsing and
# relabelling OUT of the engine and relabel() core and onto the spec, so adding
# NormalMvSpec (v0.2.0) needs no engine edits.

#' @describeIn buildConstants Univariate Normal-Inverse-Gamma constants.
setMethod("buildConstants", "NormalUvSpec",
  function(spec, prior, n, ...) {
    list(n = n,
         mu0 = prior$mu0, kappa0 = prior$kappa0,
         nu0 = prior$nu0, s0 = prior$s0)
  }
)

#' @describeIn buildDataList Univariate data vector.
setMethod("buildDataList", "NormalUvSpec",
  function(spec, data, ...) list(y = as.numeric(data))
)

#' @describeIn componentInits k-means dispersed start for univariate DPM.
#'
#' A k-means allocation gives a dispersed initial partition that shortens
#' burn-in.
setMethod("componentInits", "NormalUvSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    y <- as.numeric(data)
    n <- length(y)
    nUnique <- length(unique(y))
    k0 <- max(1L, min(count - 1L, as.integer(ceiling(sqrt(n)))))
    k0 <- min(k0, max(1L, nUnique))

    xiInit <- rep(1L, n)
    centers <- mean(y)
    vars <- stats::var(y)
    if (!is.finite(vars) || vars <= 0) vars <- 1

    if (identical(initMethod, "kmeans") && k0 >= 2L && nUnique >= k0) {
      km <- tryCatch(stats::kmeans(y, centers = k0, nstart = 5L),
                     error = function(e) NULL)
      if (!is.null(km)) {
        xiInit  <- as.integer(km$cluster)
        centers <- as.numeric(km$centers)
        vars <- vapply(seq_len(k0), function(j) {
          v <- stats::var(y[xiInit == j])
          if (!is.finite(v) || v <= 0) prior$s0 / (prior$nu0 - 1) else v
        }, numeric(1))
      }
    }

    muInit <- rep(prior$mu0, count)
    s2Init <- rep(prior$s0 / (prior$nu0 - 1), count)
    occ <- sort(unique(xiInit))
    for (idx in seq_along(occ)) {
      j <- occ[idx]
      muInit[j] <- if (length(centers) >= idx) centers[idx] else prior$mu0
      s2Init[j] <- if (length(vars) >= idx && is.finite(vars[idx]))
        vars[idx] else s2Init[j]
    }
    list(alloc = xiInit, params = list(muTilde = muInit, s2Tilde = s2Init))
  }
)

#' @describeIn extractParamTraces Parse muTilde / s2Tilde traces.
setMethod("extractParamTraces", "NormalUvSpec",
  function(spec, samples, L, ...) {
    list(
      mu = .nodeToArray(samples, "muTilde", L),
      s2 = .nodeToArray(samples, "s2Tilde", L)
    )
  }
)

#' @describeIn relabelComponents Permute (mu, s2) and summarise components.
setMethod("relabelComponents", "NormalUvSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...) {
    muTr <- paramTrace$mu[idx, , drop = FALSE]
    s2Tr <- paramTrace$s2[idx, , drop = FALSE]
    m <- length(idx)

    # Align each retained iteration's occupied clusters to columns 1..modalK,
    # then apply the label permutation.
    muRe <- matrix(NA_real_, m, modalK)
    s2Re <- matrix(NA_real_, m, modalK)
    for (r in seq_len(m)) {
      occ <- occList[[r]]
      muRe[r, ] <- muTr[r, occ][perms[r, ]]
      s2Re[r, ] <- s2Tr[r, occ][perms[r, ]]
    }

    q <- function(M, p) apply(M, 2L, stats::quantile, probs = p, names = FALSE)
    summ <- data.frame(
      component = seq_len(modalK),
      weight    = colMeans(weights),
      mu_mean   = colMeans(muRe),
      mu_lwr    = q(muRe, 0.025),
      mu_upr    = q(muRe, 0.975),
      s2_mean   = colMeans(s2Re),
      s2_lwr    = q(s2Re, 0.025),
      s2_upr    = q(s2Re, 0.975)
    )
    list(mu = muRe, s2 = s2Re, summary = summ)
  }
)
