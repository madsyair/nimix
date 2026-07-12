#' @include class-FitResult.R
NULL

## ---------------------------------------------------------------------------
## methods-summary-plot.R
##
## summary() / plot() / predict() for FitResult.
##
## summary() ALWAYS reports relabelled component estimates: if relabel() has not
## been run it is run automatically (presenting raw per-component posterior
## summaries before relabelling is not meaningful). It also reports a simple
## ESS-based mixing diagnostic and warns on poor mixing. Cross-chain Rhat
## requires multiple chains and is deferred to v0.5.0 (multi-chain support).
##
## plot()/predict() are dimension-aware: they dispatch their per-component
## density through the spec's componentDensity() so univariate and multivariate
## share one code path.
## ---------------------------------------------------------------------------

#' @describeIn FitResult Posterior summary (relabelled component estimates,
#'   posterior of the number of clusters, and a mixing diagnostic).
#' @param object A \code{FitResult}.
#' @param ... Passed to \code{\link{relabel}} when relabelling is triggered.
#' @export
setMethod("summary", "FitResult", function(object, ...) {
  if (length(object@relabeled) == 0L) {
    message("Relabelling MCMC output before summarising ",
 "(label switching)...")
    object <- relabel(object, ...)
  }
  rl <- object@relabeled

  Ktab <- prop.table(table(object@Kposterior))

  diag   <- object@diagnostics
  nch    <- if (!is.null(diag$nchains)) diag$nchains else 1L
  nDraws <- nrow(object@mcmcSamples)
  # Fall back to a single-chain computation for results produced before the
  # diagnostics slot was populated.
  essK <- if (!is.null(diag$essK)) diag$essK else tryCatch(
    coda::effectiveSize(coda::as.mcmc(as.numeric(object@Kposterior))),
    error = function(e) NA_real_)
  essAlpha <- if (!is.null(diag$essAlpha)) diag$essAlpha else tryCatch(
    coda::effectiveSize(coda::as.mcmc(object@mcmcSamples[, "alpha"])),
    error = function(e) NA_real_)
  RhatK <- diag$RhatK

  if (nch > 1L && is.finite(RhatK) && RhatK > 1.1)
    warning("Cross-chain split-Rhat for the cluster count is ", round(RhatK, 3),
            " (> 1.1): the ", nch, " chains disagree on the number of clusters ",
            "and have likely not converged. Run longer or add chains.",
            call. = FALSE)
  if (is.finite(essK) && essK < 0.05 * nDraws)
    warning("Low effective sample size for the cluster count (ESS = ",
            round(essK), " of ", nDraws, " draws): the chain may be mixing ",
            "poorly across partitions. ",
            "Consider a longer run or k-means initialisation.", call. = FALSE)

  cat("nimix mixture summary (engine: ", object@engineUsed,
      ", distribution: ", object@distSpec@name, ")\n", sep = "")
  cat("Observations: ", .nObs(object@data), " (dimension d = ",
      .dataDimOf(object@data), ")\n", sep = "")
  cat("Relabelling: ", rl$method, " conditioned on modal K = ", rl$modalK,
      " (", rl$nDraws, " draws)\n\n", sep = "")
  cat("Posterior of number of occupied clusters:\n")
  print(round(Ktab, 3))
  cat("\nRelabelled component estimates (posterior mean; ",
      "CIs for univariate):\n", sep = "")
  print(format(rl$summary, digits = 3), row.names = FALSE)
  if (nch > 1L) {
    cat("\nConvergence (", nch, " chains, ", nDraws, " pooled draws):\n", sep = "")
    cat("  split-Rhat(#clusters) = ", round(RhatK, 3),
        "  ESS(#clusters) = ", round(essK), "\n", sep = "")
    if (!is.null(diag$RhatAlpha) && is.finite(diag$RhatAlpha))
      cat("  split-Rhat(alpha)     = ", round(diag$RhatAlpha, 3),
          "  ESS(alpha) = ", round(essAlpha), "\n", sep = "")
    if (!is.null(diag$RhatBeta) && is.finite(diag$RhatBeta))
      cat("  split-Rhat(beta MRF)  = ", round(diag$RhatBeta, 3),
          "  ESS(beta) = ", round(diag$essBeta),
          "  E[beta|y] = ", round(diag$betaMean, 3), "\n", sep = "")
    if (!is.null(diag$functionals)) {
      cat("  Convergence over label-invariant functionals",
          " (Vehtari et al. 2021; aim Rhat < 1.01):\n", sep = "")
      ft <- diag$functionals
      for (r in seq_len(nrow(ft)))
        cat(sprintf("    %-8s Rhat %-6s folded %-6s bulkESS %-6s tailESS %s\n",
            ft$functional[r], round(ft$Rhat[r], 3), round(ft$foldedRhat[r], 3),
            round(ft$bulkESS[r]), round(ft$tailESS[r])))
    }
    cat("  (Rhat computed on label-invariant quantities; ",
        "target Rhat < 1.1)\n", sep = "")
  } else {
    cat("\nMixing diagnostic (single chain): ESS(#clusters) = ",
        round(essK), sep = "")
    if (is.finite(essAlpha)) cat(", ESS(alpha) = ", round(essAlpha), sep = "")
    cat("\n  Set mcmcControl$nchains > 1 for cross-chain split-Rhat.\n", sep = "")
  }

  invisible(list(Kposterior = Ktab, components = rl$summary,
                 ess = c(alpha = essAlpha, K = essK),
                 diagnostics = diag, object = object))
})

#' @describeIn FitResult Diagnostic and result plots.
#' @param x A \code{FitResult}.
#' @param y Ignored.
#' @param type One of \code{"K"} (posterior of #clusters), \code{"trace_raw"}
#'   (raw cluster-parameter traces; zig-zags reveal label switching),
#'   \code{"trace_relabeled"} (traces after relabelling), \code{"density"}
#'   (univariate clustering only: data histogram with posterior predictive
#'   overlay), \code{"cluster"} (multivariate clustering, \eqn{d \ge 2}: scatter
#'   coloured by MAP cluster), or \code{"fitted"} (regression only: observed
#'   response vs posterior predictive mean).
#' @export
#' @return Invisibly, a tidy data frame of exactly what was drawn (e.g.
#'   \code{iteration}/\code{component}/\code{value} for traces,
#'   \code{x}/\code{density} for the predictive density), so the plot can be
#'   reproduced with ggplot2 or any other system without nimix depending on
#'   them.
setMethod("plot", signature(x = "FitResult", y = "missing"),
  function(x, y, type = c("K", "trace_raw", "trace_relabeled", "density",
                          "cluster", "fitted"), ...) {
    type <- match.arg(type)
    d <- .dataDimOf(x@data)
    isReg <- isRegressionSpec(x@distSpec)

    # Primary location-trace as an (iter x K) matrix, whatever the component:
    # cluster means for clustering, first coefficient for regression.
    .locTrace <- function(pt, relabeled = FALSE) {
      src <- if (relabeled) x@relabeled else pt
      if (!is.null(src$mu))   return(if (d == 1L) src$mu else src$mu[, , 1L])
      if (!is.null(src$beta)) return(src$beta[, , 1L])
      stop("no location trace available.", call. = FALSE)
    }

    # Every branch draws with base graphics AND returns (invisibly) the tidy
    # data frame it plotted. That keeps `graphics` as the only plotting
    # dependency while letting users replot with ggplot2/lattice/plotly from
    # the returned data -- no need for nimix to carry those packages.
    out <- NULL
    if (type == "K") {
      tab <- prop.table(table(x@Kposterior))
      graphics::barplot(tab,
                        xlab = "number of occupied clusters",
                        ylab = "posterior probability",
                        main = "Posterior of K")
      out <- data.frame(K = as.integer(names(tab)),
                        probability = as.numeric(tab))

    } else if (type == "trace_raw") {
      lab <- if (isReg) "first coefficient (raw)" else
        if (d == 1L) "muTilde (raw)" else "muTilde[, , dim 1] (raw)"
      tr <- .locTrace(x@paramTrace)
      graphics::matplot(tr, type = "l", lty = 1,
                        xlab = "iteration", ylab = lab,
                        main = "Raw cluster-parameter traces (look for switching)")
      out <- data.frame(iteration = rep(seq_len(nrow(tr)), ncol(tr)),
                        component = rep(seq_len(ncol(tr)), each = nrow(tr)),
                        value = as.numeric(tr))

    } else if (type == "trace_relabeled") {
      if (length(x@relabeled) == 0L) x <- relabel(x)
      lab <- if (isReg) "first coefficient (relabelled)" else
        if (d == 1L) "mu (relabelled)" else "mu[, , dim 1] (relabelled)"
      tr <- .locTrace(x@paramTrace, relabeled = TRUE)
      graphics::matplot(tr, type = "l",
                        lty = 1, xlab = "iteration (modal-K draws)", ylab = lab,
                        main = "Relabelled component traces")
      out <- data.frame(iteration = rep(seq_len(nrow(tr)), ncol(tr)),
                        component = rep(seq_len(ncol(tr)), each = nrow(tr)),
                        value = as.numeric(tr))

    } else if (type == "density") {
      if (isReg)
        stop("type = 'density' is for clustering; for a regression fit use ",
             "type = 'fitted'.", call. = FALSE)
      if (d != 1L)
        stop("type = 'density' is univariate only; for d >= 2 use ",
             "type = 'cluster'.", call. = FALSE)
      xs <- seq(min(x@data), max(x@data), length.out = 256)
      pp <- predict(x, newdata = xs)
      graphics::hist(x@data, breaks = "FD", freq = FALSE, col = "grey90",
                     border = "grey60", main = "Posterior predictive density",
                     xlab = "y")
      graphics::lines(pp$x, pp$density, lwd = 2)
      out <- data.frame(x = pp$x, density = pp$density)

    } else if (type == "cluster") {
      if (isReg || d < 2L)
        stop("type = 'cluster' needs multivariate clustering (d >= 2).",
             call. = FALSE)
      X <- as.matrix(x@data)
      # MAP allocation per observation = most frequent cluster across draws,
      # recoded to consecutive integers for colouring (labels themselves are
      # not identified; this is purely a visual partition).
      mapAlloc <- apply(x@clusterAllocation, 2L, function(col) {
        tb <- table(col); as.integer(names(tb)[which.max(tb)])
      })
      grp <- as.integer(factor(mapAlloc))
      graphics::plot.default(X[, 1L], X[, 2L], col = grp, pch = 19,
                     xlab = "dimension 1", ylab = "dimension 2",
                     main = "Data coloured by MAP cluster (dims 1-2)")
      out <- data.frame(X, cluster = grp)
      names(out) <- c(paste0("dim", seq_len(ncol(X))), "cluster")

    } else if (type == "fitted") {
      if (!isReg)
        stop("type = 'fitted' is for regression fits (nimixReg).",
             call. = FALSE)
      pr <- predict(x)
      graphics::plot.default(pr$.fitted, x@data, pch = 19, col = "grey40",
                     xlab = "posterior predictive mean", ylab = "observed y",
                     main = "Observed vs fitted (mixture of regressions)")
      graphics::abline(0, 1, lty = 2)
      out <- data.frame(fitted = pr$.fitted, observed = as.numeric(x@data))
    }
    invisible(out)
  }
)

#' @describeIn FitResult Posterior predictive density at new points.
#'
#' Returns the posterior predictive density averaged over MCMC draws, using the
#' occupied-cluster mixture in each draw (weights = cluster sizes / n). This
#' quantity is label-invariant, so no relabelling is required. For multivariate
#' fits the density is evaluated at the supplied rows of \code{newdata}. For a
#' regression fit (\code{nimixReg}) it instead returns the posterior predictive
#' \emph{mean} \eqn{E[y \mid x]} per row of \code{newdata} (column
#' \code{.fitted}). To keep evaluation tractable the draws are subsampled to at
#' most \code{maxDraws}.
#' @param newdata Points at which to evaluate: a numeric vector (univariate
#'   clustering), a matrix with \code{d} columns (multivariate clustering), or a
#'   data frame of predictors (regression). Defaults to the training data.
#' @param maxDraws Integer cap on the number of posterior draws used (default
#'   500); draws are thinned uniformly if exceeded.
#' @export
setMethod("predict", "FitResult", function(object, newdata, maxDraws = 500L,
                                           ...) {
  if (isRegressionSpec(object@distSpec))
    return(.predictReg(object, if (missing(newdata)) NULL else newdata,
                       maxDraws))

  d <- .dataDimOf(object@data)
  alloc <- object@clusterAllocation
  n <- ncol(alloc)
  m <- nrow(alloc)
  dfun <- componentDensity(object@distSpec, df = object@prior$df)

  use <- if (m > maxDraws)
    as.integer(round(seq(1, m, length.out = maxDraws))) else seq_len(m)
  mUse <- length(use)

  if (d == 1L) {
    if (missing(newdata)) newdata <- object@data
    xs <- as.numeric(newdata)
    # Per-cluster scalar parameter traces (mu, s2 or mu, tau, ...). Build the
    # params list generically from the trace keys and evaluate the component
    # density, so univariate Normal, Student-t, and Normal-Gamma all work.
    keys <- names(object@paramTrace)
    keys <- keys[vapply(keys, function(k) is.matrix(object@paramTrace[[k]]),
                        logical(1))]
    dfval <- object@prior$df
    dens <- numeric(length(xs))
    for (r in use) {
      occ <- sort(unique(alloc[r, ]))
      w <- tabulate(match(alloc[r, ], occ), nbins = length(occ)) / n
      for (k in seq_along(occ)) {
        j <- occ[k]
        params <- lapply(object@paramTrace[keys], function(M) M[r, j])
        if (!is.null(dfval)) params$df <- dfval
        if (!is.null(object@prior$size)) params$size <- object@prior$size
        dens <- dens + w[k] * dfun(xs, params)
      }
    }
    return(data.frame(x = xs, density = dens / mUse))
  }

  # multivariate
  if (missing(newdata)) newdata <- object@data
  X <- as.matrix(newdata)
  q <- nrow(X)
  muTr  <- object@paramTrace$mu      # iter x L x d
  covTr <- object@paramTrace$cov     # iter x L x d x d
  dens <- numeric(q)
  for (r in use) {
    occ <- sort(unique(alloc[r, ]))
    w <- tabulate(match(alloc[r, ], occ), nbins = length(occ)) / n
    for (k in seq_along(occ)) {
      j <- occ[k]
      Sig <- matrix(covTr[r, j, , ], d, d)
      muj <- muTr[r, j, ]
      dens <- dens + w[k] *
        vapply(seq_len(q), function(t) dfun(X[t, ],
                 list(mu = muj, Sigma = Sig, df = object@prior$df)),
               numeric(1))
    }
  }
  out <- as.data.frame(X)
  out$density <- dens / mUse
  out
})

# Posterior predictive MEAN E[y | x] for a mixture-of-regressions fit. For each
# draw the prediction is the CRP-weighted mixture of cluster linear predictors;
# we average over (subsampled) draws. Label-invariant, so no relabelling needed.
.predictReg <- function(object, newdata, maxDraws = 500L) {
  prior <- object@prior
  if (is.null(newdata)) {
    Xnew <- prior$X
    nd <- NULL
  } else {
    tt <- stats::delete.response(prior$terms)
    mf <- stats::model.frame(tt, newdata, na.action = stats::na.pass)
    Xnew <- stats::model.matrix(tt, mf)
    nd <- newdata
  }
  alloc <- object@clusterAllocation
  n <- ncol(alloc); m <- nrow(alloc)
  betaTr <- object@paramTrace$beta            # iter x L x p   (x d if mv-response)
  q <- nrow(Xnew)
  isMv <- length(dim(betaTr)) == 4L

  use <- if (m > maxDraws)
    as.integer(round(seq(1, m, length.out = maxDraws))) else seq_len(m)

  if (isMv) {
    d <- dim(betaTr)[4]
    fitted <- matrix(0, q, d)
    for (r in use) {
      occ <- sort(unique(alloc[r, ]))
      w <- tabulate(match(alloc[r, ], occ), nbins = length(occ)) / n
      for (k in seq_along(occ))
        fitted <- fitted + w[k] * (Xnew %*% betaTr[r, occ[k], , ])
    }
    fitted <- fitted / length(use)
    rn <- object@prior$respNames
    if (is.null(rn)) rn <- paste0("y", seq_len(d))
    out <- if (!is.null(nd)) as.data.frame(nd) else as.data.frame(Xnew)
    for (jj in seq_len(d)) out[[paste0(".fitted.", rn[jj])]] <- fitted[, jj]
    return(out)
  }

  fitted <- numeric(q)
  for (r in use) {
    occ <- sort(unique(alloc[r, ]))
    w <- tabulate(match(alloc[r, ], occ), nbins = length(occ)) / n
    for (k in seq_along(occ)) {
      j <- occ[k]
      fitted <- fitted + w[k] *
        linkInv(object@distSpec, as.numeric(Xnew %*% betaTr[r, j, ]),
                prior = prior)
    }
  }
  fitted <- fitted / length(use)
  out <- if (!is.null(nd)) as.data.frame(nd) else as.data.frame(Xnew)
  out$.fitted <- fitted
  out
}
