# workflow-prior-predictive.R
#
# Prior predictive check (Gelman et al. 2020, "Bayesian workflow", section 2.4):
# simulate whole datasets from the prior -- parameters from simulateParams(),
# mixture weights from the Dirichlet prior, data through the same per-family
# simulators the posterior predictive check uses -- and put the observed data
# beside them. If the observed data sit far outside what the prior can produce,
# the prior (usually its scale) deserves a second look before any MCMC is run.
#
# Reuse is deliberate: simulateParams() draws component parameters from the
# prior for every family that supports it, and .ppcSimulate() knows how to turn
# parameters into data. The bridge between them is the naming convention
# `<name>` in simulateParams -> monitor column `<name>Tilde[j]`, which holds
# across the univariate families.

#' Prior predictive check for a mixture model
#'
#' Simulates \code{nsim} datasets from the prior implied by
#' \code{distribution}, \code{K}, and the data-scaled default prior (optionally
#' overridden), then compares summary statistics of the observed data with
#' their prior predictive distributions. Run this \emph{before} fitting: if the
#' observed statistics fall far outside the prior predictive range, the prior
#' scale is off for these data.
#'
#' @param data Numeric vector of observations (univariate; the multivariate
#'   families are not yet covered).
#' @param K Number of mixture components to simulate under.
#' @param distribution Component family name, as in \code{\link{nimixClust}}.
#' @param nsim Number of prior predictive datasets. Default 200.
#' @param prior Optional named list overriding entries of the data-scaled
#'   default prior (as in \code{nimixClust}).
#' @param conc Dirichlet concentration for the mixture weights. Default 1
#'   (uniform on the simplex), matching the fixed-K engine default.
#' @param seed Optional RNG seed.
#' @return An object of class \code{nimixPriorPred}: a list with the observed
#'   statistics (\code{obs}), the matrix of simulated statistics (\code{sim},
#'   \code{nsim} rows), the tail probability of each observed statistic under
#'   the prior predictive (\code{pTail}, two-sided), and the simulated datasets'
#'   summary. Its \code{print} method flags statistics with
#'   \code{pTail < 0.05}, and its \code{plot} method overlays prior predictive
#'   densities on the observed data.
#' @references Gelman, A., et al. (2020). Bayesian workflow.
#'   \emph{arXiv:2011.01808}.
#' @examples
#' \dontrun{
#' y <- c(rnorm(80, -2), rnorm(120, 3))
#' pp <- priorPredictive(y, K = 2, distribution = "normal")
#' pp          # flags any statistic the prior cannot reach
#' plot(pp)    # observed density over prior predictive draws
#' }
#' @export
priorPredictive <- function(data, K, distribution = "normal", nsim = 200L,
                            prior = NULL, conc = 1, seed = NULL) {
  if (!is.numeric(data) || is.matrix(data))
    stop("priorPredictive() currently covers univariate data only.",
         call. = FALSE)
  y <- as.numeric(data); n <- length(y)
  if (n < 2L) stop("Need at least two observations.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)

  spec <- .selectClusterSpec(distribution, isMv = FALSE, d = 1L)
  pr <- defaultPrior(spec, y)
  if (!is.null(prior)) pr[names(prior)] <- prior

  statFn <- function(v) c(mean = mean(v), sd = stats::sd(v),
                          min = min(v), max = max(v),
                          skew = mean(((v - mean(v)) / stats::sd(v))^3))
  obs <- statFn(y)

  sims <- matrix(NA_real_, nsim, length(obs),
                 dimnames = list(NULL, names(obs)))
  yreps <- vector("list", min(nsim, 20L))   # keep a few full draws for plot()
  for (s in seq_len(nsim)) {
    params <- simulateParams(spec, pr, K)
    # Build a one-row samples matrix in monitor-column format so the same
    # per-family simulator as ppCheck() can be reused verbatim.
    cols <- unlist(lapply(names(params), function(nm)
      stats::setNames(as.numeric(params[[nm]]),
                      sprintf("%sTilde[%d]", nm, seq_len(K)))))
    samples <- matrix(cols, nrow = 1L, dimnames = list(NULL, names(cols)))
    w <- stats::rgamma(K, conc); w <- w / sum(w)
    alloc <- sample.int(K, n, replace = TRUE, prob = w)
    yr <- .ppcSimulate(spec, samples, 1L, alloc, pr)
    sims[s, ] <- statFn(yr)
    if (s <= length(yreps)) yreps[[s]] <- yr
  }

  # Two-sided tail probability of each observed statistic under the prior
  # predictive: min(P(sim <= obs), P(sim >= obs)) * 2, capped at 1.
  pTail <- vapply(seq_along(obs), function(j) {
    lo <- mean(sims[, j] <= obs[j]); hi <- mean(sims[, j] >= obs[j])
    min(1, 2 * min(lo, hi))
  }, numeric(1))
  names(pTail) <- names(obs)

  structure(list(obs = obs, sim = sims, pTail = pTail, yreps = yreps, y = y,
                 distribution = distribution, K = K, n = n, nsim = nsim),
            class = "nimixPriorPred")
}

#' @method print nimixPriorPred
#' @export
print.nimixPriorPred <- function(x, ...) {
  cat("Prior predictive check (", x$distribution, ", K = ", x$K,
      ", nsim = ", x$nsim, ")\n\n", sep = "")
  qs <- apply(x$sim, 2L, stats::quantile, probs = c(.025, .975), na.rm = TRUE)
  df <- data.frame(observed = round(x$obs, 3),
                   prior_2.5 = round(qs[1, ], 3),
                   prior_97.5 = round(qs[2, ], 3),
                   p_tail = round(x$pTail, 3),
                   flag = ifelse(x$pTail < 0.05, "<-- outside", ""))
  print(df)
  if (any(x$pTail < 0.05))
    cat("\nFlagged statistics fall in the far tail of what this prior can",
        "produce;\nrevisit the prior scale (or the distribution choice)",
        "before fitting.\n")
  else
    cat("\nAll observed statistics are comfortably within the prior",
        "predictive range.\n")
  invisible(x)
}

#' @method plot nimixPriorPred
#' @export
plot.nimixPriorPred <- function(x, ...) {
  dens <- lapply(x$yreps, stats::density)
  dObs <- stats::density(x$y)
  xr <- range(c(vapply(dens, function(d) range(d$x), numeric(2)), dObs$x))
  yr <- c(0, 1.1 * max(dObs$y, vapply(dens, function(d) max(d$y), numeric(1))))
  plot(NULL, xlim = xr, ylim = yr, xlab = "y", ylab = "density",
       main = sprintf("Prior predictive check (%s, K = %d)",
                      x$distribution, x$K), ...)
  for (d in dens)
    graphics::lines(d, col = grDevices::adjustcolor("steelblue", 0.30))
  graphics::lines(dObs, col = "black", lwd = 2)
  graphics::legend("topright",
                   legend = c("observed data", "prior predictive draws"),
                   col = c("black", "steelblue"), lwd = c(2, 1), bty = "n")
  invisible(x)
}
