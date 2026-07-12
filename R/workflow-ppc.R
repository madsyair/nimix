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

# Column lookup for a monitored node. The grep costs ~3e-05 s per call
# (measured), so this memoisation is about clarity and about very large
# monitor sets, not about speed: ppCheck()/posteriorPredict() attach an
# environment to `samples` once, and every .ppcSimulate call across the draw
# loop then resolves each node's columns exactly once instead of nrep times.
# The attribute is optional, so direct .ppcCols(samples, node) calls (as in
# tests) keep working unchanged.
.ppcCols <- function(samples, node) {
  cache <- attr(samples, "ppcColsCache", exact = TRUE)
  if (!is.null(cache) && !is.null(cache[[node]])) return(cache[[node]])
  cols <- grep(paste0("^", node, "\\["), colnames(samples))
  if (!length(cols)) stop("ppCheck(): node '", node, "' is not monitored.",
                          call. = FALSE)
  if (!is.null(cache)) cache[[node]] <- cols
  cols
}

.withPpcColsCache <- function(samples) {
  attr(samples, "ppcColsCache") <- new.env(parent = emptyenv())
  samples
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
    muC  <- samples[draw, .ppcCols(samples, "muTilde")]
    covC <- samples[draw, .ppcCols(samples, "covTilde")]
    # The number of monitored components is fixed by the number of columns, NOT
    # by max(alloc): when a component is empty on this draw, max(alloc) < Kmon,
    # and reshaping covC to c(max(alloc), d, d) would truncate and shift the
    # covariance entries -- yielding a non-symmetric (or, worse, still-PD but
    # wrong) matrix silently. Derive the dimension from the monitored count.
    Kmon <- length(muC) / d
    if (length(covC) != Kmon * d * d)
      stop("ppc: covariance trace length ", length(covC), " is not Kmon*d*d = ",
           Kmon * d * d, "; monitored components/dimension are inconsistent.",
           call. = FALSE)
    n <- length(alloc)
    out <- matrix(0, n, d)
    muM  <- matrix(muC, ncol = d)                      # Kmon x d (column-major)
    covA <- array(covC, dim = c(Kmon, d, d))
    for (k in seq_len(max(alloc))) {
      idx <- which(alloc == k)
      if (!length(idx)) next
      ch <- chol(matrix(covA[k, , ], d, d))
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
#' @param store_yrep If \code{TRUE}, attach the simulated replicates as
#'   attribute \code{"yrep"} (with \code{"y"} and \code{"draws"}) on the
#'   result -- the inputs graphical PPC functions consume.
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
                    seed = 1L, store_yrep = FALSE) {
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
  # Replicates are expensive to produce and are exactly what graphical PPC
  # (e.g. bayesplot::ppc_dens_overlay) consumes; discarding them forces users
  # to re-simulate. Storage is opt-in because yrep is nrep x n (x d).
  yrepStore <- if (store_yrep) vector("list", length(draws)) else NULL
  samples <- .withPpcColsCache(samples)   # node columns resolved once, not per draw
  set.seed(seed)
  for (r in seq_along(draws)) {
    dr <- draws[r]
    yrep <- .ppcSimulate(spec, samples, dr, as.integer(alloc[dr, ]), fit@prior)
    if (store_yrep) yrepStore[[r]] <- yrep
    repMat[r, ] <- applyStats(yrep)
  }
  pval <- colMeans(sweep(repMat, 2L, obs, ">="))
  out <- data.frame(statistic = names(obs), observed = as.numeric(obs),
                    repMean = colMeans(repMat), ppp = pval, row.names = NULL)
  class(out) <- c("nimixPPC", class(out))
  if (store_yrep) {
    yrepOut <- if (is.matrix(yrepStore[[1L]]))
      aperm(simplify2array(yrepStore), c(3L, 1L, 2L))   # nrep x n x d
    else do.call(rbind, yrepStore)                       # nrep x n
    attr(out, "yrep")  <- yrepOut
    attr(out, "y")     <- y
    attr(out, "draws") <- draws
  }
  out
}

#' Posterior predictive replicates
#'
#' Simulates replicated data sets from the fitted mixture, conditionally on
#' each retained draw's parameters and allocations, and returns them --
#' the input that graphical posterior predictive checking (for example
#' \code{bayesplot::ppc_dens_overlay(y, yrep)}) consumes. \code{\link{ppCheck}}
#' computes summary-statistic tail probabilities from the same replicates;
#' this function exposes the replicates themselves.
#'
#' @param fit A \code{FitResult} from a clustering fit.
#' @param ndraws Number of posterior draws to use (thinned evenly).
#' @param seed RNG seed for the replicate noise.
#' @return For univariate data, an \code{ndraws x n} matrix; for multivariate
#'   data, an \code{ndraws x n x d} array. The attribute \code{"draws"} records
#'   which posterior iterations were used.
#' @seealso \code{\link{ppCheck}} for tail-probability summaries.
#' @export
posteriorPredict <- function(fit, ndraws = 100, seed = 1L) {
  if (!methods::is(fit, "FitResult"))
    stop("posteriorPredict() expects a FitResult.", call. = FALSE)
  if (length(fit@prior$X))
    stop("posteriorPredict() currently supports clustering fits.",
         call. = FALSE)
  spec <- fit@distSpec
  samples <- fit@mcmcSamples
  alloc <- fit@clusterAllocation
  m <- nrow(alloc)
  ndraws <- min(ndraws, m)
  draws <- unique(round(seq(1L, m, length.out = ndraws)))
  samples <- .withPpcColsCache(samples)   # node columns resolved once, not per draw
  set.seed(seed)
  reps <- lapply(draws, function(dr)
    .ppcSimulate(spec, samples, dr, as.integer(alloc[dr, ]), fit@prior))
  out <- if (is.matrix(reps[[1L]]))
    aperm(simplify2array(reps), c(3L, 1L, 2L))
  else do.call(rbind, reps)
  attr(out, "draws") <- draws
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

setMethod(".ppcSimulate", "SEPUvSpec",
  function(spec, samples, draw, alloc, prior) {
    mu <- samples[draw, .ppcCols(samples, "muTilde")]
    sg <- samples[draw, .ppcCols(samples, "sigmaTilde")]
    nu <- samples[draw, .ppcCols(samples, "nuTilde")]
    rsep(length(alloc), mu[alloc], sg[alloc], nu[alloc])
  })

setMethod(".ppcSimulate", "LEPUvSpec",
  function(spec, samples, draw, alloc, prior) {
    mu <- samples[draw, .ppcCols(samples, "muTilde")]
    sg <- samples[draw, .ppcCols(samples, "sigmaTilde")]
    nu <- samples[draw, .ppcCols(samples, "nuTilde")]
    rlep(length(alloc), mu[alloc], sg[alloc], nu[alloc])
  })

setMethod(".ppcSimulate", "FSSNUvSpec",
  function(spec, samples, draw, alloc, prior) {
    mu <- samples[draw, .ppcCols(samples, "muTilde")]
    sg <- samples[draw, .ppcCols(samples, "sigmaTilde")]
    al <- samples[draw, .ppcCols(samples, "alphaTilde")]
    rfssn(length(alloc), mu[alloc], sg[alloc], al[alloc])
  })

setMethod(".ppcSimulate", "FOSSEPUvSpec",
  function(spec, samples, draw, alloc, prior) {
    mu <- samples[draw, .ppcCols(samples, "muTilde")]
    sg <- samples[draw, .ppcCols(samples, "sigmaTilde")]
    al <- samples[draw, .ppcCols(samples, "alphaTilde")]
    th <- samples[draw, .ppcCols(samples, "thetaTilde")]
    rfossep(length(alloc), mu[alloc], sg[alloc], al[alloc], th[alloc])
  })

setMethod(".ppcSimulate", "FSSTUvSpec",
  function(spec, samples, draw, alloc, prior) {
    mu <- samples[draw, .ppcCols(samples, "muTilde")]
    sg <- samples[draw, .ppcCols(samples, "sigmaTilde")]
    al <- samples[draw, .ppcCols(samples, "alphaTilde")]
    nu <- samples[draw, .ppcCols(samples, "nuTilde")]
    rfsst(length(alloc), mu[alloc], sg[alloc], al[alloc], nu[alloc])
  })

setMethod(".ppcSimulate", "JFSTUvSpec",
  function(spec, samples, draw, alloc, prior) {
    mu <- samples[draw, .ppcCols(samples, "muTilde")]
    sg <- samples[draw, .ppcCols(samples, "sigmaTilde")]
    al <- samples[draw, .ppcCols(samples, "alphaTilde")]
    th <- samples[draw, .ppcCols(samples, "thetaTilde")]
    rjfst(length(alloc), mu[alloc], sg[alloc], al[alloc], th[alloc])
  })

setMethod(".ppcSimulate", "SkewNormalMvSpec",
  function(spec, samples, draw, alloc, prior) {
    d <- prior$d
    K <- max(alloc)
    muC  <- samples[draw, .ppcCols(samples, "muTilde")]
    sigC <- samples[draw, .ppcCols(samples, "SigTilde")]
    gamC <- samples[draw, .ppcCols(samples, "gamTilde")]
    muM  <- matrix(muC, ncol = d)
    sigA <- array(sigC, dim = c(length(muC) / d, d, d))
    gamM <- matrix(gamC, ncol = d)
    n <- length(alloc)
    out <- matrix(0, n, d)
    for (k in seq_len(K)) {
      idx <- which(alloc == k)
      if (length(idx))
        out[idx, ] <- rskewmvn(length(idx), muM[k, ],
                               matrix(sigA[k, , ], d, d), gamM[k, ])
    }
    out
  })

setMethod(".ppcSimulate", "SkewIStudentMvSpec",
  function(spec, samples, draw, alloc, prior) {
    d <- prior$d; K <- max(alloc)
    muC  <- samples[draw, .ppcCols(samples, "muTilde")]
    sigC <- samples[draw, .ppcCols(samples, "SigTilde")]
    gamC <- samples[draw, .ppcCols(samples, "gamTilde")]
    nuC  <- samples[draw, .ppcCols(samples, "nuTilde")]
    muM  <- matrix(muC, ncol = d)
    sigA <- array(sigC, dim = c(length(muC) / d, d, d))
    gamM <- matrix(gamC, ncol = d); nuM <- matrix(nuC, ncol = d)
    out <- matrix(0, length(alloc), d)
    for (k in seq_len(K)) {
      idx <- which(alloc == k)
      if (length(idx))
        out[idx, ] <- rskewmvit(length(idx), muM[k, ],
                                matrix(sigA[k, , ], d, d), gamM[k, ], nuM[k, ])
    }
    out
  })

setMethod(".ppcSimulate", "SkewIStudentMvSpec",
  function(spec, samples, draw, alloc, prior) {
    d <- prior$d; K <- max(alloc)
    muC  <- samples[draw, .ppcCols(samples, "muTilde")]
    sigC <- samples[draw, .ppcCols(samples, "SigTilde")]
    gamC <- samples[draw, .ppcCols(samples, "gamTilde")]
    nuC  <- samples[draw, .ppcCols(samples, "nuTilde")]
    muM  <- matrix(muC, ncol = d)
    sigA <- array(sigC, dim = c(length(muC) / d, d, d))
    gamM <- matrix(gamC, ncol = d)
    nuM  <- matrix(nuC, ncol = d)
    out <- matrix(0, length(alloc), d)
    for (k in seq_len(K)) {
      idx <- which(alloc == k)
      if (length(idx))
        out[idx, ] <- rskewmvit(length(idx), muM[k, ],
                                matrix(sigA[k, , ], d, d), gamM[k, ], nuM[k, ])
    }
    out
  })

setMethod(".ppcSimulate", "SkewIStudentMvSpec",
  function(spec, samples, draw, alloc, prior) {
    d <- prior$d
    K <- max(alloc)
    muM  <- matrix(samples[draw, .ppcCols(samples, "muTilde")],  ncol = d)
    gamM <- matrix(samples[draw, .ppcCols(samples, "gamTilde")], ncol = d)
    nuM  <- matrix(samples[draw, .ppcCols(samples, "nuTilde")],  ncol = d)
    sigA <- array(samples[draw, .ppcCols(samples, "SigTilde")],
                  dim = c(nrow(muM), d, d))
    out <- matrix(0, length(alloc), d)
    for (k in seq_len(K)) {
      idx <- which(alloc == k)
      if (length(idx))
        out[idx, ] <- rskewmvit(length(idx), muM[k, ],
                                matrix(sigA[k, , ], d, d), gamM[k, ], nuM[k, ])
    }
    out
  })

setMethod(".ppcSimulate", "SkewNormalMvOSpec",
  function(spec, samples, draw, alloc, prior) {
    d <- prior$d
    K <- max(alloc)
    muM  <- matrix(samples[draw, .ppcCols(samples, "muTilde")],  ncol = d)
    gamM <- matrix(samples[draw, .ppcCols(samples, "gamTilde")], ncol = d)
    thV  <- samples[draw, .ppcCols(samples, "thetaTilde")]
    sigA <- array(samples[draw, .ppcCols(samples, "SigTilde")],
                  dim = c(nrow(muM), d, d))
    out <- matrix(0, length(alloc), d)
    for (k in seq_len(K)) {
      idx <- which(alloc == k)
      if (length(idx))
        out[idx, ] <- rskewmvno(length(idx), muM[k, ],
                                matrix(sigA[k, , ], d, d), gamM[k, ], thV[k])
    }
    out
  })

setMethod(".ppcSimulate", "SkewIStudentMvOSpec",
  function(spec, samples, draw, alloc, prior) {
    d <- prior$d
    K <- max(alloc)
    muM  <- matrix(samples[draw, .ppcCols(samples, "muTilde")],  ncol = d)
    gamM <- matrix(samples[draw, .ppcCols(samples, "gamTilde")], ncol = d)
    nuM  <- matrix(samples[draw, .ppcCols(samples, "nuTilde")],  ncol = d)
    thV  <- samples[draw, .ppcCols(samples, "thetaTilde")]
    sigA <- array(samples[draw, .ppcCols(samples, "SigTilde")],
                  dim = c(nrow(muM), d, d))
    out <- matrix(0, length(alloc), d)
    for (k in seq_len(K)) {
      idx <- which(alloc == k)
      if (length(idx))
        out[idx, ] <- rskewmvito(length(idx), muM[k, ],
                                 matrix(sigA[k, , ], d, d), gamM[k, ],
                                 nuM[k, ], thV[k])
    }
    out
  })

setMethod(".ppcSimulate", "SkewNormalMvOGenSpec",
  function(spec, samples, draw, alloc, prior) {
    d <- prior$d; K <- max(alloc); nAng <- .nAngles(d)
    muM  <- matrix(samples[draw, .ppcCols(samples, "muTilde")],  ncol = d)
    gamM <- matrix(samples[draw, .ppcCols(samples, "gamTilde")], ncol = d)
    thM  <- matrix(samples[draw, .ppcCols(samples, "thetaTilde")], ncol = nAng)
    sigA <- array(samples[draw, .ppcCols(samples, "SigTilde")],
                  dim = c(nrow(muM), d, d))
    out <- matrix(0, length(alloc), d)
    for (k in seq_len(K)) {
      idx <- which(alloc == k)
      if (!length(idx)) next
      O <- orthogonalFactor(thM[k, ], d)
      U <- chol(matrix(sigA[k, , ], d, d))
      g <- gamM[k, ]; nk <- length(idx)
      G <- matrix(g, nk, d, byrow = TRUE)
      W <- abs(matrix(stats::rnorm(nk * d), nk, d))
      pos <- matrix(stats::runif(nk * d), nk, d) < G^2 / (1 + G^2)
      Eps <- ifelse(pos, W * G, -W / G)
      out[idx, ] <- sweep(Eps %*% O %*% U, 2L, muM[k, ], "+")
    }
    out
  })

setMethod(".ppcSimulate", "SkewIStudentMvOGenSpec",
  function(spec, samples, draw, alloc, prior) {
    d <- prior$d; K <- max(alloc); nAng <- .nAngles(d)
    muM  <- matrix(samples[draw, .ppcCols(samples, "muTilde")],  ncol = d)
    gamM <- matrix(samples[draw, .ppcCols(samples, "gamTilde")], ncol = d)
    nuM  <- matrix(samples[draw, .ppcCols(samples, "nuTilde")],  ncol = d)
    thM  <- matrix(samples[draw, .ppcCols(samples, "thetaTilde")], ncol = nAng)
    sigA <- array(samples[draw, .ppcCols(samples, "SigTilde")],
                  dim = c(nrow(muM), d, d))
    out <- matrix(0, length(alloc), d)
    for (k in seq_len(K)) {
      idx <- which(alloc == k)
      if (!length(idx)) next
      O <- orthogonalFactor(thM[k, ], d)
      U <- chol(matrix(sigA[k, , ], d, d))
      g <- gamM[k, ]; nu <- nuM[k, ]; nk <- length(idx)
      G <- matrix(g, nk, d, byrow = TRUE)
      W <- abs(matrix(stats::rt(nk * d, df = rep(nu, each = nk)), nk, d))
      pos <- matrix(stats::runif(nk * d), nk, d) < G^2 / (1 + G^2)
      Eps <- ifelse(pos, W * G, -W / G)
      out[idx, ] <- sweep(Eps %*% O %*% U, 2L, muM[k, ], "+")
    }
    out
  })
