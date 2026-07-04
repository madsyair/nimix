#' @include class-FitResult.R
#' @include class-DistributionSpec.R
NULL

## workflow-ppc.R ---------------------------------------------------------------
## Posterior predictive checking (Bayesian workflow; Gelman et al. 2020,
## "Bayesian workflow", Section 6; Gelman, Meng & Stern 1996 for posterior
## predictive p-values). Replicated data are drawn CONDITIONALLY on each
## posterior draw's fitted allocation ("mixed" posterior predictive): y*_i is
## simulated from the component that observation i occupies in that draw. This
## makes the check label-invariant by construction and valid for every engine,
## including the MRF (the replicate inherits the spatially smoothed partition).

# Per-family replicate simulators: parse THAT draw's component parameters from
# the monitored samples and simulate one replicate data set.

#' @keywords internal
setGeneric(".ppcSimulate", function(spec, samples, draw, alloc, prior)
  standardGeneric(".ppcSimulate"))

setMethod(".ppcSimulate", "DistributionSpec",
  function(spec, samples, draw, alloc, prior)
    stop("ppCheck() is not yet available for '", spec@name, "' components.",
         call. = FALSE))

.ppcCols <- function(samples, node) {
  cols <- grep(paste0("^", node, "\\["), colnames(samples))
  if (!length(cols)) stop("ppCheck(): node '", node, "' is not monitored.",
                          call. = FALSE)
  cols
}

setMethod(".ppcSimulate", "NormalUvSpec",
  function(spec, samples, draw, alloc, prior) {
    mu <- samples[draw, .ppcCols(samples, "muTilde")]
    s2 <- samples[draw, .ppcCols(samples, "s2Tilde")]
    stats::rnorm(length(alloc), mu[alloc], sqrt(s2[alloc]))
  })

## Normal-Gamma (uv): marginally Student-t(df) around muTilde with scale
## sqrt(s2Tilde) -- the scale-mixture identity (West 1987; Geweke 1993).
setMethod(".ppcSimulate", "NormalGammaUvSpec",
  function(spec, samples, draw, alloc, prior) {
    mu <- samples[draw, .ppcCols(samples, "muTilde")]
    s2 <- samples[draw, .ppcCols(samples, "s2Tilde")]
    df <- if (!is.null(prior$df)) prior$df else 4
    mu[alloc] + sqrt(s2[alloc]) * stats::rt(length(alloc), df)
  })

setMethod(".ppcSimulate", "StudentTUvSpec",
  function(spec, samples, draw, alloc, prior) {
    mu  <- samples[draw, .ppcCols(samples, "muTilde")]
    tau <- samples[draw, .ppcCols(samples, "tauTilde")]
    df  <- if (!is.null(prior$df)) prior$df else 4
    mu[alloc] + (1 / sqrt(tau[alloc])) * stats::rt(length(alloc), df)
  })

setMethod(".ppcSimulate", "PoissonSpec",
  function(spec, samples, draw, alloc, prior) {
    lam <- samples[draw, .ppcCols(samples, "lambda")]
    stats::rpois(length(alloc), lam[alloc])
  })

setMethod(".ppcSimulate", "BinomialSpec",
  function(spec, samples, draw, alloc, prior) {
    pr <- samples[draw, .ppcCols(samples, "prob")]
    size <- if (!is.null(prior$size)) prior$size else 1
    stats::rbinom(length(alloc), size, pr[alloc])
  })

setMethod(".ppcSimulate", "NormalMvSpec",
  function(spec, samples, draw, alloc, prior) {
    d <- prior$d
    K <- max(alloc)
    muC  <- samples[draw, .ppcCols(samples, "muTilde")]
    covC <- samples[draw, .ppcCols(samples, "covTilde")]
    n <- length(alloc)
    out <- matrix(0, n, d)
    muM <- matrix(muC, ncol = d)                       # K x d (column-major)
    covA <- array(covC, dim = c(K, d, d))
    for (k in seq_len(K)) {
      idx <- which(alloc == k)
      if (!length(idx)) next
      ch <- chol(covA[k, , ])
      z <- matrix(stats::rnorm(length(idx) * d), ncol = d)
      out[idx, ] <- rep(muM[k, ], each = length(idx)) + z %*% ch
    }
    out
  })

#' Posterior predictive check
#'
#' Simulates replicated data sets from the fitted mixture -- conditionally on
#' each posterior draw's allocation, which makes the check label-invariant and
#' valid for every engine including the spatial MRF -- and compares test
#' statistics of the replicates against the observed data via posterior
#' predictive p-values (Gelman, Meng & Stern 1996; Gelman et al. 2020,
#' Section 6). Extreme p-values (near 0 or 1) flag aspects of the data the
#' model fails to reproduce; values near 0.5 indicate no evidence of misfit
#' for that statistic.
#'
#' @param fit A \code{\linkS4class{FitResult}} from \code{\link{nimixClust}}.
#' @param nrep Number of replicated data sets (posterior draws used), thinned
#'   evenly from the retained draws. Default 200.
#' @param statistics Named list of test-statistic functions, or a character
#'   vector naming built-ins from \code{mean}, \code{sd}, \code{min},
#'   \code{max}, \code{skew}. Multivariate data applies each statistic
#'   column-wise and reports per-column results.
#' @param seed RNG seed for the replicate simulation.
#' @return An object of class \code{nimixPPC}: a data frame of observed value,
#'   replicate mean, and posterior predictive p-value per statistic, printed
#'   with guidance.
#' @references
#' Gelman, A., Meng, X.-L., & Stern, H. (1996). Posterior predictive
#' assessment of model fitness via realized discrepancies.
#' \emph{Statistica Sinica}, 6, 733--807.
#'
#' Gelman, A., et al. (2020). Bayesian workflow. \emph{arXiv:2011.01808},
#' Section 6.
#' @export
ppCheck <- function(fit, nrep = 200, statistics = c("mean", "sd", "min", "max"),
                    seed = 1L) {
  if (!methods::is(fit, "FitResult"))
    stop("ppCheck() expects a FitResult.", call. = FALSE)
  if (methods::is(fit@distSpec, "RegressionMixModel") ||
      length(fit@prior$X))
    stop("ppCheck() currently supports clustering fits; regression ",
         "predictive checks are planned.", call. = FALSE)
  y <- fit@data
  spec <- fit@distSpec
  samples <- fit@mcmcSamples
  alloc <- fit@clusterAllocation
  m <- nrow(alloc)
  nrep <- min(nrep, m)
  draws <- unique(round(seq(1L, m, length.out = nrep)))

  skew <- function(v) {
    v <- v[is.finite(v)]
    mean((v - mean(v))^3) / stats::sd(v)^3
  }
  builtins <- list(mean = mean, sd = stats::sd, min = min, max = max,
                   skew = skew)
  statFns <- if (is.character(statistics)) builtins[statistics] else statistics
  if (any(vapply(statFns, is.null, logical(1))))
    stop("Unknown statistic; built-ins are: ",
         paste(names(builtins), collapse = ", "), call. = FALSE)

  applyStats <- function(dat) {
    if (is.matrix(dat)) {
      unlist(lapply(names(statFns), function(nm) {
        v <- apply(dat, 2L, statFns[[nm]])
        stats::setNames(v, paste0(nm, "_", seq_along(v)))
      }))
    } else {
      vapply(names(statFns), function(nm) statFns[[nm]](dat), numeric(1))
    }
  }

  obs <- applyStats(y)
  repMat <- matrix(NA_real_, length(draws), length(obs),
                   dimnames = list(NULL, names(obs)))
  set.seed(seed)
  for (r in seq_along(draws)) {
    dr <- draws[r]
    yrep <- .ppcSimulate(spec, samples, dr, as.integer(alloc[dr, ]), fit@prior)
    repMat[r, ] <- applyStats(yrep)
  }
  pval <- colMeans(sweep(repMat, 2L, obs, ">="))
  out <- data.frame(statistic = names(obs), observed = as.numeric(obs),
                    repMean = colMeans(repMat), ppp = pval, row.names = NULL)
  class(out) <- c("nimixPPC", class(out))
  out
}

#' @export
print.nimixPPC <- function(x, ...) {
  cat("Posterior predictive check (", nrow(x), " statistics):\n", sep = "")
  df <- x; class(df) <- "data.frame"
  df$observed <- signif(df$observed, 4)
  df$repMean  <- signif(df$repMean, 4)
  df$ppp      <- round(df$ppp, 3)
  print(df, row.names = FALSE)
  flag <- df$statistic[df$ppp < 0.05 | df$ppp > 0.95]
  if (length(flag))
    cat("\nStatistics with extreme posterior predictive p-values (",
        paste(flag, collapse = ", "),
        "): the model does not reproduce these aspects of the data well.\n",
        sep = "")
  else
    cat("\nNo statistic shows extreme posterior predictive p-values.\n")
  invisible(x)
}

# --- posterior predictive simulators (workflow) ---------------------------------

setMethod(".ppcSimulate", "MSNBurrUvSpec",
  function(spec, samples, draw, alloc, prior) {
    mu <- samples[draw, .ppcCols(samples, "muTilde")]
    sg <- samples[draw, .ppcCols(samples, "sigmaTilde")]
    al <- samples[draw, .ppcCols(samples, "alphaTilde")]
    rmsnburr(length(alloc), mu[alloc], sg[alloc], al[alloc])
  })

setMethod(".ppcSimulate", "MSNBurr2aUvSpec",
  function(spec, samples, draw, alloc, prior) {
    mu <- samples[draw, .ppcCols(samples, "muTilde")]
    sg <- samples[draw, .ppcCols(samples, "sigmaTilde")]
    al <- samples[draw, .ppcCols(samples, "alphaTilde")]
    rmsnburr2a(length(alloc), mu[alloc], sg[alloc], al[alloc])
  })

setMethod(".ppcSimulate", "GMSNBurrUvSpec",
  function(spec, samples, draw, alloc, prior) {
    mu <- samples[draw, .ppcCols(samples, "muTilde")]
    sg <- samples[draw, .ppcCols(samples, "sigmaTilde")]
    al <- samples[draw, .ppcCols(samples, "alphaTilde")]
    th <- samples[draw, .ppcCols(samples, "thetaTilde")]
    rgmsnburr(length(alloc), mu[alloc], sg[alloc], al[alloc], th[alloc])
  })
