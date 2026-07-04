## ---------------------------------------------------------------------------
## dist-poisson-binomial.R
##
## Discrete count components, proving the DistributionSpec contract is not tied
## to continuous data. Both use conjugate cluster priors that NIMBLE recognises
## for the collapsed CRP sampler:
##   Poisson : y_i ~ Poisson(lambda_k),   lambda_k ~ Gamma(a0, b0)
##   Binomial: y_i ~ Binomial(size, p_k), p_k      ~ Beta(a0, b0)
## The Binomial number of trials `size` is a known constant (supplied in the
## prior). Defaults are data-scaled: E[lambda] ~= mean(y) for Poisson and
## E[p] ~= mean(y)/size for Binomial.
## ---------------------------------------------------------------------------

# ===========================================================================
# Poisson
# ===========================================================================

#' Poisson count component specification
#'
#' @slot name Fixed to \code{"poisson"}.
#' @slot paramNames \code{c("lambda")}.
#' @export
setClass(
  "PoissonSpec",
  contains = "DistributionSpec",
  prototype = prototype(name = "poisson", paramNames = "lambda", dataDim = 1L)
)

#' Construct a Poisson component spec
#' @return A \code{\linkS4class{PoissonSpec}}.
#' @examples
#' spec <- PoissonSpec()
#' @export
PoissonSpec <- function() new("PoissonSpec")

#' @describeIn defaultPrior Data-scaled Gamma prior on the Poisson rate
#'   (\code{E[lambda] ~= mean(y)}).
#' @export
setMethod("defaultPrior", "PoissonSpec",
  function(spec, data, control = list(), ...) {
    y <- as.numeric(data)
    my <- mean(y); if (!is.finite(my) || my <= 0) my <- 1
    a0 <- if (!is.null(control$a0)) control$a0 else 2
    list(a0 = a0, b0 = a0 / my)
  }
)

#' @describeIn validateParams Validate the Gamma rate prior.
#' @export
setMethod("validateParams", "PoissonSpec",
  function(spec, params, ...) {
    if (is.null(params$a0) || is.null(params$b0) ||
        params$a0 <= 0 || params$b0 <= 0)
      stop("Poisson prior needs a0 > 0 and b0 > 0.", call. = FALSE)
    invisible(TRUE)
  }
)

#' @describeIn simulateParams Draw lambda from the Gamma prior.
setMethod("simulateParams", "PoissonSpec",
  function(spec, prior, nClust, ...)
    list(lambda = stats::rgamma(nClust, shape = prior$a0, rate = prior$b0))
)

#' @describeIn componentDensity Poisson pmf.
setMethod("componentDensity", "PoissonSpec",
  function(spec, ...) function(x, params)
    stats::dpois(round(x), lambda = params[["lambda"]])
)

#' @describeIn buildModelCode Poisson DPM model code (Gamma-Poisson conjugate).
#' @export
setMethod("buildModelCode", signature("PoissonSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) y[i] ~ dpois(lambda[xi[i]])
      xi[1:n] ~ dCRP(alpha, size = n)
      for (j in 1:L) lambda[j] ~ dgamma(shape = a0, rate = b0)
      alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
    })
    list(code = code, monitors = c("xi", "lambda", "alpha"),
         paramNodes = c(lambda = "lambda"), allocNode = "xi")
  }
)

#' @describeIn buildConstants Poisson Gamma-prior constants.
setMethod("buildConstants", "PoissonSpec",
  function(spec, prior, n, ...) list(n = n, a0 = prior$a0, b0 = prior$b0)
)

#' @describeIn buildDataList Count data vector.
setMethod("buildDataList", "PoissonSpec",
  function(spec, data, ...) list(y = as.numeric(data))
)

#' @describeIn componentInits k-means start on counts.
setMethod("componentInits", "PoissonSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...)
    .countInits(as.numeric(data), count, initMethod, prior$a0 / prior$b0,
                "lambda", transform = function(m) pmax(m, 1e-2),
                initRatio = .initRatioArg(...))
)

#' @describeIn extractParamTraces Parse lambda traces.
setMethod("extractParamTraces", "PoissonSpec",
  function(spec, samples, L, ...)
    list(lambda = .nodeToArray(samples, "lambda", L))
)

#' @describeIn relabelComponents Permute lambda and summarise.
setMethod("relabelComponents", "PoissonSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...)
    .scalarRelabel(paramTrace$lambda, idx, occList, perms, modalK, weights,
                   "lambda")
)

# ===========================================================================
# Binomial
# ===========================================================================

#' Binomial count component specification
#'
#' @slot name Fixed to \code{"binomial"}.
#' @slot paramNames \code{c("prob")}.
#' @export
setClass(
  "BinomialSpec",
  contains = "DistributionSpec",
  prototype = prototype(name = "binomial", paramNames = "prob", dataDim = 1L)
)

#' Construct a Binomial component spec
#' @return A \code{\linkS4class{BinomialSpec}}.
#' @examples
#' spec <- BinomialSpec()
#' @export
BinomialSpec <- function() new("BinomialSpec")

#' @describeIn defaultPrior Data-scaled Beta prior on the success probability.
#'   Requires the number of trials in \code{control$size}.
#' @export
setMethod("defaultPrior", "BinomialSpec",
  function(spec, data, control = list(), ...) {
    y <- as.numeric(data)
    size <- control$size
    if (is.null(size))
      stop("Binomial needs the number of trials in prior = list(size = ...).",
           call. = FALSE)
    if (any(y > size)) stop("Some counts exceed size.", call. = FALSE)
    phat <- mean(y) / size
    phat <- min(max(phat, 1e-3), 1 - 1e-3)
    list(a0 = max(0.5, 2 * phat), b0 = max(0.5, 2 * (1 - phat)),
         size = as.integer(size))
  }
)

#' @describeIn validateParams Validate the Beta prior and \code{size}.
#' @export
setMethod("validateParams", "BinomialSpec",
  function(spec, params, ...) {
    if (is.null(params$a0) || is.null(params$b0) ||
        params$a0 <= 0 || params$b0 <= 0)
      stop("Binomial prior needs a0 > 0 and b0 > 0.", call. = FALSE)
    if (is.null(params$size) || params$size < 1)
      stop("Binomial needs size >= 1.", call. = FALSE)
    invisible(TRUE)
  }
)

#' @describeIn simulateParams Draw prob from the Beta prior.
setMethod("simulateParams", "BinomialSpec",
  function(spec, prior, nClust, ...)
    list(prob = stats::rbeta(nClust, prior$a0, prior$b0))
)

#' @describeIn componentDensity Binomial pmf.
setMethod("componentDensity", "BinomialSpec",
  function(spec, size = NULL, ...) function(x, params) {
    sz <- if (!is.null(params[["size"]])) params[["size"]] else size
    stats::dbinom(round(x), size = sz, prob = params[["prob"]])
  }
)

#' @describeIn buildModelCode Binomial DPM model code (Beta-Binomial conjugate).
#' @export
setMethod("buildModelCode", signature("BinomialSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) y[i] ~ dbin(prob[xi[i]], size)
      xi[1:n] ~ dCRP(alpha, size = n)
      for (j in 1:L) prob[j] ~ dbeta(a0, b0)
      alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
    })
    list(code = code, monitors = c("xi", "prob", "alpha"),
         paramNodes = c(prob = "prob"), allocNode = "xi")
  }
)

#' @describeIn buildConstants Binomial Beta-prior constants plus \code{size}.
setMethod("buildConstants", "BinomialSpec",
  function(spec, prior, n, ...)
    list(n = n, a0 = prior$a0, b0 = prior$b0, size = prior$size)
)

#' @describeIn buildDataList Count data vector.
setMethod("buildDataList", "BinomialSpec",
  function(spec, data, ...) list(y = as.numeric(data))
)

#' @describeIn componentInits k-means start on proportions.
setMethod("componentInits", "BinomialSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...)
    .countInits(as.numeric(data), count, initMethod,
                prior$a0 / (prior$a0 + prior$b0), "prob",
                transform = function(m) min(max(m / prior$size, 1e-3),
                                            1 - 1e-3),
                initRatio = .initRatioArg(...))
)

#' @describeIn extractParamTraces Parse prob traces.
setMethod("extractParamTraces", "BinomialSpec",
  function(spec, samples, L, ...)
    list(prob = .nodeToArray(samples, "prob", L))
)

#' @describeIn relabelComponents Permute prob and summarise.
setMethod("relabelComponents", "BinomialSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...)
    .scalarRelabel(paramTrace$prob, idx, occList, perms, modalK, weights,
                   "prob")
)

# ===========================================================================
# Shared helpers for single-parameter discrete components
# ===========================================================================

# k-means start producing an allocation plus one scalar parameter per slot.
.countInits <- function(y, count, initMethod, priorMean, nodeName,
                        transform = identity, initRatio = .DEFAULT_INIT_RATIO) {
  n <- length(y); nUnique <- length(unique(y))
  # Dispersed k-means start, capped at initRatio * count (default 0.8) to leave
  # the cap: for the DPM, count = L = K_max is a hard truncation, and early CRP
  # sweeps can briefly occupy more clusters than the modal K before merging
  # down. Seeding right at the ceiling left no room for that transient.
  k0 <- max(1L, min(as.integer(floor(initRatio * count)), as.integer(ceiling(sqrt(n)))))
  k0 <- min(k0, max(1L, nUnique))
  xiInit <- rep(1L, n); centers <- mean(y)
  if (identical(initMethod, "kmeans") && k0 >= 2L && nUnique >= k0) {
    km <- tryCatch(stats::kmeans(y, centers = k0, nstart = 5L),
                   error = function(e) NULL)
    if (!is.null(km)) { xiInit <- as.integer(km$cluster)
                        centers <- as.numeric(km$centers) }
  }
  par <- rep(transform(priorMean), count)
  occ <- sort(unique(xiInit))
  for (idx in seq_along(occ))
    if (length(centers) >= idx) par[occ[idx]] <- transform(centers[idx])
  out <- list(alloc = xiInit, params = stats::setNames(list(par), nodeName))
  out
}

# Permute a single scalar parameter across retained draws and summarise.
.scalarRelabel <- function(trMat, idx, occList, perms, modalK, weights, nm) {
  tr <- trMat[idx, , drop = FALSE]; m <- length(idx)
  re <- matrix(NA_real_, m, modalK)
  for (r in seq_len(m)) { occ <- occList[[r]]; re[r, ] <- tr[r, occ][perms[r, ]] }
  q <- function(M, p) apply(M, 2L, stats::quantile, probs = p, names = FALSE)
  summ <- data.frame(component = seq_len(modalK), weight = colMeans(weights))
  summ[[paste0(nm, "_mean")]] <- colMeans(re)
  summ[[paste0(nm, "_lwr")]]  <- q(re, 0.025)
  summ[[paste0(nm, "_upr")]]  <- q(re, 0.975)
  stats::setNames(list(re, summ), c(nm, "summary"))
}
