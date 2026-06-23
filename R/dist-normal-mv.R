## ---------------------------------------------------------------------------
## dist-normal-mv.R
##
## Multivariate Gaussian component (the v0.2.0 distribution).
##
## Conjugate cluster base measure is Normal-Inverse-Wishart on (mean, COV):
##     covTilde_j         ~ InverseWishart(S0, df0)      # covariance matrix
##     muTilde_j | cov    ~ N(mu0, covTilde_j / kappa0)
## which is the standard conjugate prior for a multivariate normal with unknown
## mean and covariance (Frühwirth-Schnatter 2006, Section
## 10.6). It is mathematically equivalent to a Normal-Wishart prior on the
## PRECISION; we use the covariance + inverse-Wishart parameterisation because
## (a) it is the pairing NIMBLE recognises as conjugate for dmnorm with cov =
## (so the dCRP sampler assigns the collapsed conjugate CRP_cluster_wrapper
## updates -- mixing empty components are not sampled),
## and (b) the precision + dwish path triggers a Cholesky-lifting failure under
## the dynamic CRP indexing in NIMBLE 1.4.x (verified empirically). The
## conjugacy fact itself is native NIMBLE (0.5): nimix only
## supplies the S4 wiring, the data-scaled prior, and the relabelling.
##
## NOTE on native vs nimix (0.5): dCRP, the conjugate CRP
## samplers, dmnorm and dinvwish are all native NIMBLE. nimix contributes the
## DistributionSpec architecture, the DATA-SCALED Normal-Inverse-Wishart default
## prior, and the multivariate relabelling/summaries.
##
## Reference for the multivariate-normal mixture base measure and recovery:
## Zhang, Chan, Wu & Chen (2004);
## Dellaportas & Papageorgiou (2006);
## Frühwirth-Schnatter (2006).
##
## Defaults are DATA-SCALED: mu0 = colMeans(data), the prior mean
## of the cluster covariance equals cov(data), and the mean's prior dispersion
## is cLoc-scaled. df0 > d + 1 is ENFORCED so prior draws on empty components
## have a finite, non-singular covariance.
## ---------------------------------------------------------------------------

#' Multivariate Gaussian component specification
#'
#' @slot name Fixed to \code{"normal-mv"}.
#' @slot paramNames \code{c("mu", "Sigma")}.
#' @slot dataDim \code{NA_integer_}: the actual dimension \eqn{d} is taken from
#'   the data at fit time and carried in the prior list.
#'
#' @references
#' Zhang, Z., Chan, K.L., Wu, Y., & Chen, C. (2004). Learning a multivariate
#' Gaussian mixture model with the reversible jump MCMC algorithm.
#' \emph{Statistics and Computing}, 14, 343--355.
#' \doi{10.1023/B:STCO.0000039481.32735.0c}
#'
#' Dellaportas, P., & Papageorgiou, I. (2006). Multivariate mixtures of normals
#' with unknown number of components. \emph{Statistics and Computing}, 16,
#' 57--68. \doi{10.1007/s11222-006-5338-6}
#'
#' @seealso \code{\link{nimixClust}}, \code{\linkS4class{NormalUvSpec}}
#' @export
setClass(
  "NormalMvSpec",
  contains = "DistributionSpec",
  prototype = prototype(
    name = "normal-mv",
    paramNames = c("mu", "Sigma"),
    dataDim = NA_integer_
  )
)

#' Construct a multivariate Gaussian component spec
#' @return A \code{\linkS4class{NormalMvSpec}}.
#' @examples
#' spec <- NormalMvSpec()
#' @export
NormalMvSpec <- function() new("NormalMvSpec")

# --- defaultPrior ----------------------------------------------------------

#' @describeIn defaultPrior Data-scaled Normal-Inverse-Wishart prior (multivariate).
#'
#' Control overrides: \code{cLoc} (mean-dispersion multiplier; prior covariance
#' of \code{mu} is \code{Sigma / kappa0} with \code{kappa0 = 1 / cLoc^2},
#' default \code{cLoc = 2}) and \code{df0} (inverse-Wishart degrees of freedom,
#' default \code{d + 2}, must exceed \code{d + 1} for a finite, non-singular
#' prior covariance on empty components).
#' @export
setMethod("defaultPrior", "NormalMvSpec",
  function(spec, data, control = list(), ...) {
    if (!is.matrix(data)) data <- as.matrix(data)
    d  <- ncol(data)
    Sy <- stats::cov(data)
    if (any(!is.finite(Sy))) Sy <- diag(d)

    cLoc <- if (!is.null(control$cLoc)) control$cLoc else 2
    df0  <- if (!is.null(control$df0))  control$df0  else (d + 2)
    if (df0 <= d + 1)
      stop("df0 must exceed d + 1 (here d = ", d, ") so the prior covariance ",
           "on empty components is finite and non-singular ",
 ".", call. = FALSE)

    kappa0 <- 1 / (cLoc^2)
    # Covariance ~ InverseWishart(S0, df0) has E[Sigma] = S0 / (df0 - d - 1).
    # Target E[Sigma] = cov(data)  ==>  S0 = cov(data) * (df0 - d - 1).
    S0 <- Sy * (df0 - d - 1)
    ev <- tryCatch(min(eigen((S0 + t(S0)) / 2, symmetric = TRUE,
                             only.values = TRUE)$values),
                   error = function(e) NA_real_)
    if (is.na(ev) || ev <= 0)
      S0 <- S0 + diag(1e-6 * mean(diag(S0)) + 1e-8, d)

    list(mu0 = colMeans(data), kappa0 = kappa0, df0 = df0, S0 = S0,
         d = d, cLoc = cLoc)
  }
)

# --- validateParams --------------------------------------------------------

#' @describeIn validateParams Validate a Normal-Inverse-Wishart prior list and
#'   enforce the dimension invariant \eqn{\dim(\mu_0) = \dim(S_0) = d} and
#' \eqn{df_0 > d + 1} (.b encapsulation,
#' ).
#' @export
setMethod("validateParams", "NormalMvSpec",
  function(spec, params, ...) {
    req <- c("mu0", "kappa0", "df0", "S0", "d")
    miss <- setdiff(req, names(params))
    if (length(miss))
      stop("NormalMvSpec prior is missing: ", paste(miss, collapse = ", "),
           call. = FALSE)
    d <- params$d
    if (params$kappa0 <= 0) stop("kappa0 must be > 0.", call. = FALSE)
    if (params$df0 <= d + 1)
 stop("df0 must exceed d + 1 (avoid singular prior draws ",
           "on empty components).", call. = FALSE)
    if (length(params$mu0) != d)
      stop("mu0 must have length d = ", d, ".", call. = FALSE)
    if (!is.matrix(params$S0) || !all(dim(params$S0) == c(d, d)))
      stop("S0 must be a ", d, " x ", d, " matrix.", call. = FALSE)
    ev <- min(eigen((params$S0 + t(params$S0)) / 2, symmetric = TRUE,
                    only.values = TRUE)$values)
    if (ev <= 0) stop("S0 must be positive definite.", call. = FALSE)
    invisible(TRUE)
  }
)

# --- simulateParams (for inits / recovery) ---------------------------------

#' @describeIn simulateParams Draw (mu, Sigma) per cluster from the
#'   Normal-Inverse-Wishart prior. Returns \code{mu} (nClust x d) and
#'   \code{Sigma} (d x d x nClust).
setMethod("simulateParams", "NormalMvSpec",
  function(spec, prior, nClust, ...) {
    d  <- prior$d
    # Sigma ~ InverseWishart(S0, df0): if W ~ Wishart(scale = S0^{-1}, df0)
    # then Sigma = W^{-1} ~ InverseWishart(S0, df0). stats::rWishart uses the
    # SCALE parameterisation (E[W] = df * scale).
    scaleMat <- solve(prior$S0)
    mu  <- matrix(NA_real_, nClust, d)
    Sig <- array(NA_real_, dim = c(d, d, nClust))
    for (j in seq_len(nClust)) {
      W <- stats::rWishart(1, df = prior$df0, Sigma = scaleMat)[, , 1]
      covj <- solve(W)
      Sig[, , j] <- covj
      mu[j, ] <- .rmvnorm1(prior$mu0, covj / prior$kappa0)
    }
    list(mu = mu, Sigma = Sig)
  }
)

# --- componentDensity ------------------------------------------------------

#' @describeIn componentDensity Multivariate normal density for predictive use.
setMethod("componentDensity", "NormalMvSpec",
  function(spec, ...) {
    function(x, params) .dmvnorm(as.numeric(x), params[["mu"]], params[["Sigma"]])
  }
)

# --- buildModelCode: NormalMvSpec x DPMEngine ------------------------------

#' @describeIn buildModelCode Multivariate Gaussian DPM model code (dCRP) with a
#'   conjugate Normal-Inverse-Wishart cluster base measure.
#'
#' Builds NIMBLE code for
#' \deqn{y_i \sim N_d(\mu_{\xi_i}, \Sigma_{\xi_i}), \quad
#'       \xi_{1:n} \sim CRP(\alpha, n),}
#' with \eqn{\Sigma_j \sim InvWishart(S_0, df_0)} and
#' \eqn{\mu_j \sim N_d(\mu_0, \Sigma_j / \kappa_0)}. The dimension \eqn{d} is a
#' constant so the index ranges \code{1:d} expand at model-build time.
#'
#' @references
#' Zhang et al. (2004) \doi{10.1023/B:STCO.0000039481.32735.0c};
#' de Valpine et al. (2017) \doi{10.1080/10618600.2016.1172487}.
#' @export
setMethod("buildModelCode", signature("NormalMvSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) {
    # Conjugate Normal-Inverse-Wishart DPM. dmnorm(cov=) + dinvwish is the
    # pairing NIMBLE recognises as conjugate for dCRP, enabling the collapsed
 # CRP_cluster_wrapper updates (Neal 2000, Algorithm 2).
    #
    # Two NIMBLE requirements shape this code:
    #  (1) a multivariate distribution parameter may not be an *expression*, so
    #      the scaled mean covariance is bound to a deterministic node covMu;
    #  (2) a dynamically indexed multivariate parameter (covTilde[xi[i], , ])
    #      cannot be fed straight to dmnorm, so per-observation deterministic
    #      nodes muObs[i, ] / covObs[i, , ] resolve the cluster lookup first --
    #      the same indirection the NIMBLE BNP examples use for the univariate
    #      kernel (muTilde[xi[i]] -> mu[i]).
    #
    # Cluster-parameter arrays have length L = K_max (NIMBLE's exact truncation;
    # the sampler is proper as long as the number of occupied clusters stays
    # strictly below L -- NIMBLE warns if a draw would exceed it).
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        muObs[i, 1:d]        <- muTilde[xi[i], 1:d]
        covObs[i, 1:d, 1:d]  <- covTilde[xi[i], 1:d, 1:d]
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
    list(
      code       = code,
      monitors   = c("xi", "muTilde", "covTilde", "alpha"),
      paramNodes = c(mu = "muTilde", cov = "covTilde"),
      allocNode  = "xi"
    )
  }
)

# --- Engine-facing methods -------------------------------------------------

#' @describeIn buildConstants Multivariate Normal-Inverse-Wishart constants
#'   (includes the dimension \code{d}, mean vector \code{mu0} and scale matrix
#'   \code{S0}).
setMethod("buildConstants", "NormalMvSpec",
  function(spec, prior, n, ...) {
    list(n = n, d = prior$d,
         mu0 = prior$mu0, kappa0 = prior$kappa0,
         df0 = prior$df0, S0 = prior$S0)
  }
)

#' @describeIn buildDataList Multivariate data matrix (one row per observation).
setMethod("buildDataList", "NormalMvSpec",
  function(spec, data, ...) list(y = as.matrix(data))
)

#' @describeIn componentInits k-means dispersed start for the multivariate DPM.
setMethod("componentInits", "NormalMvSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    Y <- as.matrix(data)
    n <- nrow(Y); d <- ncol(Y)
    nUnique <- nrow(unique(Y))
    # Dispersed k-means start, but capped at 0.8 * count to leave headroom below
    # the cap: for the DPM, count = L = K_max is a hard truncation, and early CRP
    # sweeps can briefly occupy more clusters than the modal K before merging
    # down. Seeding right at the ceiling left no room for that transient.
    k0 <- max(1L, min(as.integer(floor(0.8 * count)), as.integer(ceiling(sqrt(n)))))
    k0 <- min(k0, max(1L, nUnique))

    priorMeanCov <- prior$S0 / (prior$df0 - d - 1)

    xiInit  <- rep(1L, n)
    centers <- matrix(prior$mu0, nrow = 1L)
    covList <- list(priorMeanCov)

    if (identical(initMethod, "kmeans") && k0 >= 2L && nUnique >= k0) {
      km <- tryCatch(stats::kmeans(Y, centers = k0, nstart = 5L),
                     error = function(e) NULL)
      if (!is.null(km)) {
        xiInit  <- as.integer(km$cluster)
        centers <- km$centers
        covList <- lapply(seq_len(k0), function(j) {
          rows <- which(xiInit == j)
          if (length(rows) > d) {
            cv <- stats::cov(Y[rows, , drop = FALSE])
            cv + diag(1e-6 * mean(diag(cv)) + 1e-8, d)
          } else priorMeanCov
        })
      }
    }

    muInit  <- matrix(rep(prior$mu0, each = count), nrow = count)  # count x d
    covInit <- array(0, dim = c(count, d, d))
    for (j in seq_len(count)) covInit[j, , ] <- priorMeanCov
    occ <- sort(unique(xiInit))
    for (idx in seq_along(occ)) {
      j <- occ[idx]
      if (idx <= nrow(centers)) muInit[j, ] <- centers[idx, ]
      if (idx <= length(covList)) covInit[j, , ] <- covList[[idx]]
    }
    list(alloc = xiInit, params = list(muTilde = muInit, covTilde = covInit))
  }
)

#' @describeIn extractParamTraces Parse muTilde (L x d) and covTilde
#'   (L x d x d) traces into arrays.
setMethod("extractParamTraces", "NormalMvSpec",
  function(spec, samples, L, d = NULL, ...) {
    if (is.null(d)) stop("extractParamTraces(NormalMvSpec) needs 'd'.",
                         call. = FALSE)
    list(
      mu  = .nodeToArray(samples, "muTilde",  c(L, d)),
      cov = .nodeToArray(samples, "covTilde", c(L, d, d)),
      d   = d
    )
  }
)

#' @describeIn relabelComponents Permute multivariate (mu, Sigma) and summarise.
setMethod("relabelComponents", "NormalMvSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...) {
    d     <- paramTrace$d
    muTr  <- paramTrace$mu                      # (iter x L x d)
    covTr <- paramTrace$cov                     # (iter x L x d x d)
    m <- length(idx)

    muRe  <- array(NA_real_, dim = c(m, modalK, d))
    covRe <- array(NA_real_, dim = c(m, modalK, d, d))
    for (t in seq_len(m)) {
      r   <- idx[t]
      occ <- occList[[t]]
      ord <- perms[t, ]                         # permutation of occupied order
      for (k in seq_len(modalK)) {
        j <- occ[ord[k]]
        muRe[t, k, ] <- muTr[r, j, ]
        covRe[t, k, , ] <- matrix(covTr[r, j, , ], d, d)
      }
    }

    muMean  <- apply(muRe,  c(2L, 3L),     mean)            # modalK x d
    SigMean <- apply(covRe, c(2L, 3L, 4L), mean)            # modalK x d x d
    varMean <- t(vapply(seq_len(modalK),
                        function(k) diag(matrix(SigMean[k, , ], d, d)),
                        numeric(d)))                        # modalK x d

    summ <- data.frame(component = seq_len(modalK),
                       weight = colMeans(weights))
    for (j in seq_len(d)) summ[[paste0("mu_", j)]]  <- muMean[, j]
    for (j in seq_len(d)) summ[[paste0("var_", j)]] <- varMean[, j]

    list(mu = muRe, Sigma = covRe, mu_mean = muMean,
         Sigma_mean = SigMean, summary = summ)
  }
)
