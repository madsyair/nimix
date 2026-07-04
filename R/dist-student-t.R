#' @include class-DistributionSpec.R
NULL

## ---------------------------------------------------------------------------
## dist-student-t.R
##
## Univariate Student-t component. The kernel is evaluated DIRECTLY as a t
## density, y_i ~ t(mu_k, tau_k, df), with df a fixed hyperparameter. Because
## the t density is not conjugate to the cluster prior, the cluster parameters
## are updated non-conjugately. This is the "direct" heavy-tail path; the
## NormalGammaUvSpec gives the SAME marginal distribution through a conjugate
## scale-mixture augmentation, which is cheaper per iteration. The two are
## offered as distinct choices on purpose -- identical marginal, different
## sampling cost -- and must not be described as different distributions.
##
## Cluster prior: mu_k ~ Normal(mu0, muSd^2) and tau_k ~ Gamma(aTau, bTau),
## where tau is the precision. Defaults are data-scaled: mu0 = mean(y),
## muSd ~= cLoc * sd(y), and E[tau] ~= 1 / var(y).
##
## Reference: Lange, Little & Taylor (1989) for the t as a robust model.
## ---------------------------------------------------------------------------

#' Univariate Student-t component specification
#'
#' @slot name Fixed to \code{"student-t"}.
#' @slot paramNames \code{c("mu", "tau")} (location and precision; the reported
#'   scale is \eqn{\sigma = \tau^{-1/2}}).
#' @slot dataDim \code{1L}.
#'
#' @references
#' Lange, K.L., Little, R.J.A., & Taylor, J.M.G. (1989). Robust statistical
#' modeling using the t distribution. \emph{JASA}, 84(408), 881--896.
#' \doi{10.1080/01621459.1989.10478852}
#'
#' @seealso \code{\link{NormalGammaUvSpec}} for the conjugate scale-mixture path
#'   to the same marginal distribution.
#' @export
setClass(
  "StudentTUvSpec",
  contains = "DistributionSpec",
  prototype = prototype(
    name = "student-t",
    paramNames = c("mu", "tau"),
    dataDim = 1L
  )
)

#' Construct a univariate Student-t component spec
#' @return A \code{\linkS4class{StudentTUvSpec}}.
#' @examples
#' spec <- StudentTUvSpec()
#' @export
StudentTUvSpec <- function() new("StudentTUvSpec")

# --- defaultPrior ----------------------------------------------------------

#' @describeIn defaultPrior Data-scaled Normal location / Gamma precision prior
#'   for the Student-t component.
#'
#' Control overrides: \code{cLoc} (location spread multiplier, default 2),
#' \code{df} (degrees of freedom, a fixed hyperparameter, default 4, must
#' exceed 2 for a finite component variance).
#' @export
setMethod("defaultPrior", "StudentTUvSpec",
  function(spec, data, control = list(), ...) {
    y <- as.numeric(data)
    vy <- stats::var(y)
    if (!is.finite(vy) || vy <= 0) vy <- 1
    cLoc <- if (!is.null(control$cLoc)) control$cLoc else 2
    df   <- if (!is.null(control$df))   control$df   else 4
    if (df <= 2)
      stop("df must exceed 2 so the component has a finite variance.",
           call. = FALSE)
    # E[tau] = aTau / bTau; target 1 / var(y) with a weak shape.
    aTau <- 2
    bTau <- aTau * vy
    list(mu0 = mean(y), muSd = cLoc * sqrt(vy),
         aTau = aTau, bTau = bTau, df = df, cLoc = cLoc)
  }
)

# --- validateParams --------------------------------------------------------

#' @describeIn validateParams Validate the Student-t prior list.
#' @export
setMethod("validateParams", "StudentTUvSpec",
  function(spec, params, ...) {
    req <- c("mu0", "muSd", "aTau", "bTau", "df")
    miss <- setdiff(req, names(params))
    if (length(miss))
      stop("StudentTUvSpec prior is missing: ", paste(miss, collapse = ", "),
           call. = FALSE)
    if (params$muSd <= 0) stop("muSd must be > 0.", call. = FALSE)
    if (params$aTau <= 0 || params$bTau <= 0)
      stop("aTau and bTau must be > 0.", call. = FALSE)
    if (params$df <= 2) stop("df must be > 2 (finite variance).", call. = FALSE)
    invisible(TRUE)
  }
)

# --- simulateParams --------------------------------------------------------

#' @describeIn simulateParams Draw (mu, tau) from the location/precision prior.
setMethod("simulateParams", "StudentTUvSpec",
  function(spec, prior, nClust, ...) {
    list(mu  = stats::rnorm(nClust, prior$mu0, prior$muSd),
         tau = stats::rgamma(nClust, shape = prior$aTau, rate = prior$bTau))
  }
)

# --- componentDensity ------------------------------------------------------

#' @describeIn componentDensity Location-scale Student-t density.
setMethod("componentDensity", "StudentTUvSpec",
  function(spec, df = 4, ...) {
    function(x, params) {
      sigma <- 1 / sqrt(params[["tau"]])
      dfv <- if (!is.null(params[["df"]])) params[["df"]] else df
      stats::dt((x - params[["mu"]]) / sigma, df = dfv) / sigma
    }
  }
)

# --- buildModelCode: StudentTUvSpec x DPMEngine ----------------------------

#' @describeIn buildModelCode Univariate Student-t DPM model code (dCRP). The
#'   t density is evaluated directly; df is a constant.
#' @export
setMethod("buildModelCode", signature("StudentTUvSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        y[i] ~ dt(mu = muTilde[xi[i]], tau = tauTilde[xi[i]], df = df)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
      for (j in 1:L) {
        tauTilde[j] ~ dgamma(shape = aTau, rate = bTau)
        muTilde[j]  ~ dnorm(mu0, sd = muSd)
      }
      alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
    })
    list(code = code,
         monitors  = c("xi", "muTilde", "tauTilde", "alpha"),
         paramNodes = c(mu = "muTilde", tau = "tauTilde"),
         allocNode  = "xi")
  }
)

# --- Engine-facing methods -------------------------------------------------

#' @describeIn buildConstants Student-t constants (location/precision prior + df).
setMethod("buildConstants", "StudentTUvSpec",
  function(spec, prior, n, ...) {
    list(n = n, mu0 = prior$mu0, muSd = prior$muSd,
         aTau = prior$aTau, bTau = prior$bTau, df = prior$df)
  }
)

#' @describeIn buildDataList Univariate data vector.
setMethod("buildDataList", "StudentTUvSpec",
  function(spec, data, ...) list(y = as.numeric(data))
)

#' @describeIn componentInits k-means dispersed start (location + precision).
setMethod("componentInits", "StudentTUvSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    y <- as.numeric(data); n <- length(y)
    nUnique <- length(unique(y))
    # Dispersed k-means start, capped at initRatio * count (default 0.8, tunable
    # via mcmcControl$initRatio). For the DPM, count = L = K_max is a hard
    # truncation, and early CRP sweeps can briefly occupy more clusters than the
    # modal K before merging; seeding right at the ceiling leaves no headroom.
    initRatio <- .initRatioArg(...)
    k0 <- max(1L, min(as.integer(floor(initRatio * count)), as.integer(ceiling(sqrt(n)))))
    k0 <- min(k0, max(1L, nUnique))

    xiInit <- rep(1L, n)
    centers <- mean(y)
    vars <- stats::var(y); if (!is.finite(vars) || vars <= 0) vars <- 1
    if (identical(initMethod, "kmeans") && k0 >= 2L && nUnique >= k0) {
      km <- tryCatch(stats::kmeans(y, centers = k0, nstart = 5L),
                     error = function(e) NULL)
      if (!is.null(km)) {
        xiInit  <- as.integer(km$cluster)
        centers <- as.numeric(km$centers)
        vars <- vapply(seq_len(k0), function(j) {
          v <- stats::var(y[xiInit == j]); if (!is.finite(v) || v <= 0) 1 else v
        }, numeric(1))
      }
    }
    muInit  <- rep(prior$mu0, count)
    tauInit <- rep(prior$aTau / prior$bTau, count)
    occ <- sort(unique(xiInit))
    for (idx in seq_along(occ)) {
      j <- occ[idx]
      if (length(centers) >= idx) muInit[j] <- centers[idx]
      if (length(vars) >= idx && is.finite(vars[idx]) && vars[idx] > 0)
        tauInit[j] <- 1 / vars[idx]
    }
    list(alloc = xiInit, params = list(muTilde = muInit, tauTilde = tauInit))
  }
)

#' @describeIn extractParamTraces Parse muTilde / tauTilde traces.
setMethod("extractParamTraces", "StudentTUvSpec",
  function(spec, samples, L, ...) {
    list(mu  = .nodeToArray(samples, "muTilde",  L),
         tau = .nodeToArray(samples, "tauTilde", L))
  }
)

#' @describeIn relabelComponents Permute (mu, tau) and summarise; the reported
#'   scale is \eqn{\sigma = \tau^{-1/2}}.
setMethod("relabelComponents", "StudentTUvSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...) {
    muTr  <- paramTrace$mu[idx, , drop = FALSE]
    tauTr <- paramTrace$tau[idx, , drop = FALSE]
    m <- length(idx)
    muRe  <- matrix(NA_real_, m, modalK)
    sigRe <- matrix(NA_real_, m, modalK)
    for (r in seq_len(m)) {
      occ <- occList[[r]]; pr <- perms[r, ]
      muRe[r, ]  <- muTr[r, occ][pr]
      sigRe[r, ] <- (1 / sqrt(tauTr[r, occ]))[pr]
    }
    q <- function(M, p) apply(M, 2L, stats::quantile, probs = p, names = FALSE)
    summ <- data.frame(
      component = seq_len(modalK),
      weight    = colMeans(weights),
      mu_mean   = colMeans(muRe),
      mu_lwr    = q(muRe, 0.025), mu_upr = q(muRe, 0.975),
      sigma_mean = colMeans(sigRe),
      sigma_lwr  = q(sigRe, 0.025), sigma_upr = q(sigRe, 0.975)
    )
    list(mu = muRe, sigma = sigRe, summary = summ)
  }
)
