#' @include class-FitResult.R
NULL

## model-selection.R -------------------------------------------------------------
## Predictive model selection and ensembling for fitted mixtures (Bayesian
## workflow, comparison layer). Everything is built on the pointwise
## log-likelihood matrix (draws x observations):
##
##   loglik[s, i] = log sum_k pi_k^s f(y_i | phi_k^s)
##
## with pi_k^s the occupied-cluster proportions in draw s and phi_k^s that
## draw's component parameters. This is the mixture density marginal over the
## latent label, so it is label-invariant and valid for every engine
## (FixedK, DPM, MRF). WAIC is computed natively (Watanabe 2010; Gelman,
## Hwang & Vehtari 2014); PSIS-LOO and stacking / Pseudo-BMA+ weights use the
## 'loo' package when available (Vehtari, Gelman & Gabry 2017; Yao, Vehtari,
## Simpson & Gelman 2018).

# log(mean(exp(v))) computed stably.
.logMeanExp <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return(-Inf)
  mx <- max(v)
  mx + log(mean(exp(v - mx)))
}

# Pointwise log-likelihood matrix (nDrawsUsed x n). Clustering fits only.
.pointwiseLogLik <- function(fit, maxDraws = 1000L) {
  spec <- fit@distSpec
  if (isRegressionSpec(spec))
    stop("Predictive model selection currently supports clustering fits; ",
         "regression comparison is planned.", call. = FALSE)
  d <- .dataDimOf(fit@data)
  alloc <- fit@clusterAllocation
  S <- nrow(alloc); n <- ncol(alloc)
  use <- if (S > maxDraws)
    as.integer(round(seq(1, S, length.out = maxDraws))) else seq_len(S)
  dfun <- componentDensity(spec, df = fit@prior$df)
  ll <- matrix(-Inf, length(use), n)

  if (d == 1L) {
    xs <- as.numeric(fit@data)
    keys <- names(fit@paramTrace)
    keys <- keys[vapply(keys, function(k) is.matrix(fit@paramTrace[[k]]),
                        logical(1))]
    for (ri in seq_along(use)) {
      r <- use[ri]
      occ <- sort(unique(alloc[r, ]))
      w <- tabulate(match(alloc[r, ], occ), nbins = length(occ)) / n
      mix <- numeric(n)
      for (kk in seq_along(occ)) {
        j <- occ[kk]
        params <- lapply(fit@paramTrace[keys], function(M) M[r, j])
        if (!is.null(fit@prior$df))   params$df   <- fit@prior$df
        if (!is.null(fit@prior$size)) params$size <- fit@prior$size
        mix <- mix + w[kk] * dfun(xs, params)
      }
      ll[ri, ] <- log(mix)
    }
  } else {
    X <- as.matrix(fit@data)
    muTr <- fit@paramTrace$mu; covTr <- fit@paramTrace$cov
    q <- nrow(X)
    for (ri in seq_along(use)) {
      r <- use[ri]
      occ <- sort(unique(alloc[r, ]))
      w <- tabulate(match(alloc[r, ], occ), nbins = length(occ)) / n
      mix <- numeric(q)
      for (kk in seq_along(occ)) {
        j <- occ[kk]
        Sig <- matrix(covTr[r, j, , ], d, d); muj <- muTr[r, j, ]
        mix <- mix + w[kk] * vapply(seq_len(q),
          function(t) dfun(X[t, ], list(mu = muj, Sigma = Sig,
                                        df = fit@prior$df)), numeric(1))
      }
      ll[ri, ] <- log(mix)
    }
  }
  ll
}

# Native WAIC from a pointwise log-likelihood matrix.
.waicFromLL <- function(ll) {
  lppd_i  <- apply(ll, 2L, .logMeanExp)
  pwaic_i <- apply(ll, 2L, stats::var)
  elpd_i  <- lppd_i - pwaic_i
  list(elpd_waic = sum(elpd_i), p_waic = sum(pwaic_i),
       waic = -2 * sum(elpd_i), pointwise = elpd_i,
       se = sqrt(length(elpd_i)) * stats::sd(elpd_i))
}

#' WAIC for a fitted mixture
#'
#' Widely Applicable Information Criterion (Watanabe 2010; Gelman, Hwang &
#' Vehtari 2014) from the label-invariant pointwise mixture log-likelihood.
#' Lower WAIC (equivalently higher \code{elpd_waic}) indicates better expected
#' out-of-sample predictive fit. Useful for choosing K or the component family.
#'
#' @param fit A clustering \code{\linkS4class{FitResult}}.
#' @param maxDraws Cap on posterior draws used (thinned evenly). Default 1000.
#' @return A list with \code{waic}, \code{elpd_waic}, \code{p_waic}, and
#'   \code{se} (standard error of elpd).
#' @references Watanabe, S. (2010). JMLR 11, 3571--3594.
#' Gelman, A., Hwang, J., & Vehtari, A. (2014). Stat. Comput. 24, 997--1016.
#' @export
nimixWAIC <- function(fit, maxDraws = 1000L) {
  ll <- .pointwiseLogLik(fit, maxDraws)
  w <- .waicFromLL(ll)
  structure(list(waic = w$waic, elpd_waic = w$elpd_waic,
                 p_waic = w$p_waic, se = w$se, n = ncol(ll)),
            class = "nimixIC")
}

#' PSIS-LOO for a fitted mixture
#'
#' Pareto-smoothed importance-sampling leave-one-out cross-validation
#' (Vehtari, Gelman & Gabry 2017) via the \pkg{loo} package, on the
#' label-invariant pointwise mixture log-likelihood. Requires \pkg{loo};
#' if it is not installed, use \code{\link{nimixWAIC}} instead.
#'
#' @param fit A clustering \code{\linkS4class{FitResult}}.
#' @param maxDraws Cap on posterior draws used. Default 1000.
#' @return A \code{loo} object (see \code{loo::loo}); its \code{estimates}
#'   carry \code{elpd_loo}, \code{p_loo}, and \code{looic}, and high Pareto-k
#'   values flag observations where the approximation is unreliable.
#' @references Vehtari, A., Gelman, A., & Gabry, J. (2017).
#' Stat. Comput. 27, 1413--1432.
#' @export
nimixLOO <- function(fit, maxDraws = 1000L) {
  if (!requireNamespace("loo", quietly = TRUE))
    stop("nimixLOO() needs the 'loo' package. Install it, or use nimixWAIC().",
         call. = FALSE)
  ll <- .pointwiseLogLik(fit, maxDraws)
  loo::loo(ll)
}

#' Compare mixture models by predictive fit
#'
#' Ranks several fitted mixtures by WAIC (native) and, when \pkg{loo} is
#' available, by PSIS-LOO. Models must be fitted to the same data. Use to
#' choose K, or to compare component families (e.g. Normal vs Student-t vs
#' MSNBurr) on the same data.
#'
#' @param ... Two or more clustering \code{\linkS4class{FitResult}} objects,
#'   or a single named list of them.
#' @param maxDraws Cap on posterior draws used per fit. Default 1000.
#' @return A data.frame ordered best-first, with WAIC, elpd, and (if available)
#'   LOO columns, plus \code{dWAIC} relative to the best model.
#' @export
modelSelect <- function(..., maxDraws = 1000L) {
  fits <- list(...)
  if (length(fits) == 1L && is.list(fits[[1]]) &&
      !methods::is(fits[[1]], "FitResult")) fits <- fits[[1]]
  if (length(fits) < 2L)
    stop("modelSelect() needs at least two fitted models.", call. = FALSE)
  nm <- names(fits)
  if (is.null(nm)) nm <- paste0("model", seq_along(fits))
  ns <- vapply(fits, function(f) ncol(f@clusterAllocation), integer(1))
  if (length(unique(ns)) != 1L)
    stop("All models must be fitted to the same data (n differs).",
         call. = FALSE)
  hasLoo <- requireNamespace("loo", quietly = TRUE)
  rows <- lapply(seq_along(fits), function(i) {
    ll <- .pointwiseLogLik(fits[[i]], maxDraws)
    w <- .waicFromLL(ll)
    r <- data.frame(model = nm[i], modalK = length(unique(
      apply(fits[[i]]@clusterAllocation, 2L, function(c)
        as.integer(names(which.max(table(c))))))),
      elpd_waic = w$elpd_waic, p_waic = w$p_waic, waic = w$waic,
      se_waic = w$se, stringsAsFactors = FALSE)
    if (hasLoo) {
      lo <- suppressWarnings(loo::loo(ll))
      r$elpd_loo <- lo$estimates["elpd_loo", "Estimate"]
      r$looic    <- lo$estimates["looic", "Estimate"]
    }
    r
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$waic), , drop = FALSE]
  out$dWAIC <- out$waic - min(out$waic)
  rownames(out) <- NULL
  class(out) <- c("nimixModelSelect", "data.frame")
  out
}

#' @export
print.nimixModelSelect <- function(x, ...) {
  cat("Mixture model comparison (best first; lower WAIC/LOOIC = better):\n")
  df <- x; class(df) <- "data.frame"
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], function(v) round(v, 1))
  print(df, row.names = FALSE)
  invisible(x)
}

#' Ensemble several fitted mixtures
#'
#' Combines several fitted mixtures into a single weighted predictive model,
#' rather than selecting one winner. Weights come from Bayesian stacking or
#' Pseudo-BMA+ (Yao et al. 2018) via \pkg{loo} when available, or from WAIC
#' (Akaike-style weights) natively. Stacking is the default and is the most
#' robust to model misspecification. Models must share the same data.
#'
#' @param ... Two or more clustering \code{\linkS4class{FitResult}} objects,
#'   or a single named list of them.
#' @param method Weighting scheme: \code{"stacking"} or \code{"pseudobma"}
#'   (both need \pkg{loo}), or \code{"waic"} (native, no dependency).
#' @param maxDraws Cap on posterior draws used per fit. Default 1000.
#' @return A \code{nimixEnsemble} object carrying the fits and their weights,
#'   with a \code{predict} method for the weighted predictive density.
#' @references Yao, Y., Vehtari, A., Simpson, D., & Gelman, A. (2018).
#' Bayesian Analysis 13(3), 917--1007.
#' @export
ensembleFit <- function(..., method = c("stacking", "pseudobma", "waic"),
                        maxDraws = 1000L) {
  method <- match.arg(method)
  fits <- list(...)
  if (length(fits) == 1L && is.list(fits[[1]]) &&
      !methods::is(fits[[1]], "FitResult")) fits <- fits[[1]]
  if (length(fits) < 2L)
    stop("ensembleFit() needs at least two fitted models.", call. = FALSE)
  nm <- names(fits); if (is.null(nm)) nm <- paste0("model", seq_along(fits))
  ns <- vapply(fits, function(f) ncol(f@clusterAllocation), integer(1))
  if (length(unique(ns)) != 1L)
    stop("All models must be fitted to the same data (n differs).",
         call. = FALSE)
  lls <- lapply(fits, .pointwiseLogLik, maxDraws = maxDraws)

  if (method %in% c("stacking", "pseudobma")) {
    if (!requireNamespace("loo", quietly = TRUE))
      stop("method = '", method, "' needs the 'loo' package; ",
           "use method = 'waic' for a native alternative.", call. = FALSE)
    loos <- lapply(lls, function(ll) suppressWarnings(loo::loo(ll)))
    wm <- if (method == "stacking") "stacking" else "pseudobma"
    w <- as.numeric(loo::loo_model_weights(loos, method = wm))
  } else {
    # Akaike-style WAIC weights: proportional to exp(-0.5 * dWAIC)
    waics <- vapply(lls, function(ll) .waicFromLL(ll)$waic, numeric(1))
    dw <- waics - min(waics)
    w <- exp(-0.5 * dw); w <- w / sum(w)
  }
  names(w) <- nm
  structure(list(fits = fits, weights = w, method = method, names = nm),
            class = "nimixEnsemble")
}

#' @export
print.nimixEnsemble <- function(x, ...) {
  cat("Mixture ensemble (", length(x$fits), " models, weighting = ",
      x$method, "):\n", sep = "")
  wdf <- data.frame(model = x$names, weight = round(x$weights, 3),
                    row.names = NULL)
  print(wdf[order(-wdf$weight), ], row.names = FALSE)
  invisible(x)
}

#' Weighted predictive density from a mixture ensemble
#'
#' @param object A \code{nimixEnsemble} from \code{\link{ensembleFit}}.
#' @param newdata Points at which to evaluate the density (univariate) or a
#'   matrix of rows (multivariate). Defaults to each model's own data grid.
#' @param ... Unused.
#' @return A data.frame of evaluation points and the ensemble-weighted density.
#' @name predict-nimixEnsemble
#' @exportMethod predict
setOldClass("nimixEnsemble")

#' @rdname predict-nimixEnsemble
setMethod("predict", "nimixEnsemble", function(object, newdata = NULL, ...) {
  preds <- lapply(object$fits, function(f)
    if (is.null(newdata)) predict(f) else predict(f, newdata))
  base <- preds[[1]]
  dens <- Reduce(`+`, Map(function(p, w) w * p$density, preds, object$weights))
  base$density <- dens
  base
})
