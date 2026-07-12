## bayesplot-adaptors.R -------------------------------------------------------
## Interoperability with the bayesplot ecosystem, WITHOUT taking it on as a
## dependency. bayesplot's mcmc_* functions natively accept a plain 3-D array
## of dimension iterations x chains x parameters, and its ppc_* functions take
## (y, yrep) directly -- so the adaptors below return base R objects and nimix
## keeps bayesplot in Suggests only (used by examples and vignettes, guarded by
## requireNamespace()).
##
## The one piece of real substance here is the safety guard. Component
## parameters (muTilde, gamTilde, ...) are NOT valid inputs for convergence
## diagnostics before relabelling: the mixture likelihood is invariant to label
## permutations, so "muTilde[1]" refers to different components in different
## chains (and different iterations). R-hat computed on such a trace looks
## legitimate and means nothing. drawsArray() therefore refuses to hand out
## component draws unless relabel() has been run, and clearly says why --
## label-INVARIANT functionals (number of occupied clusters, allocation
## entropy, the DP concentration alpha) are the quantities that are safe on
## raw draws, and they are what the default returns.

#' @include class-FitResult.R
NULL

#' Posterior draws as an iterations x chains x parameters array
#'
#' Returns a plain 3-D array in the layout \code{bayesplot}'s \code{mcmc_*}
#' functions accept natively (\code{iterations x chains x parameters}), so no
#' extra packages are required to hand a nimix fit to
#' \code{bayesplot::mcmc_trace()}, \code{mcmc_rhat_hist()} and friends.
#'
#' Two views are available, and the distinction is statistical, not cosmetic:
#' \describe{
#'   \item{\code{"invariant"} (default)}{Label-invariant functionals -- the
#'     number of occupied clusters, the allocation entropy, and (for DPM fits)
#'     the concentration parameter \code{alpha}. These are meaningful on raw
#'     draws and retain the per-chain structure, so cross-chain diagnostics
#'     like R-hat apply.}
#'   \item{\code{"components"}}{Per-component parameters \emph{after}
#'     \code{\link{relabel}}. Refused if \code{relabel()} has not been run:
#'     under label switching, \code{muTilde[1]} names different components in
#'     different chains, and an R-hat computed on it looks valid while meaning
#'     nothing. Because relabelling conditions on the modal cluster count,
#'     chains lose equal lengths; the draws are therefore pooled into a single
#'     chain, suitable for posterior density/interval plots but not for
#'     cross-chain R-hat.}
#' }
#'
#' @param fit A \code{FitResult}.
#' @param params \code{"invariant"} or \code{"components"}.
#' @return A numeric array \code{iterations x chains x parameters} with
#'   dimnames on the parameter margin.
#' @seealso \code{\link{ppcData}}, \code{\link{relabel}}, \code{\link{psm}}.
#' @export
drawsArray <- function(fit, params = c("invariant", "components")) {
  if (!methods::is(fit, "FitResult"))
    stop("drawsArray() expects a FitResult.", call. = FALSE)
  params <- match.arg(params)

  if (params == "invariant") {
    alloc <- fit@clusterAllocation
    m <- nrow(alloc)
    K <- as.numeric(fit@Kposterior)
    ent <- apply(alloc, 1L, function(z) {
      p <- tabulate(z) / length(z); p <- p[p > 0]; -sum(p * log(p))
    })
    cols <- list(K = K, entropy = ent)
    if ("alpha" %in% colnames(fit@mcmcSamples))
      cols$alpha <- as.numeric(fit@mcmcSamples[, "alpha"])
    M <- do.call(cbind, cols)
    cid <- fit@diagnostics$chainId
    if (is.null(cid)) cid <- rep(1L, m)          # fits predating chainId
    nch <- length(unique(cid))
    nit <- m / nch
    if (nit != floor(nit))
      stop("Chains have unequal lengths; cannot form a draws array.",
           call. = FALSE)
    out <- array(NA_real_, c(nit, nch, ncol(M)),
                 dimnames = list(NULL, NULL, colnames(M)))
    for (ch in seq_len(nch)) out[, ch, ] <- M[cid == unique(cid)[ch], ]
    return(out)
  }

  # params == "components"
  rl <- fit@relabeled
  if (is.null(rl) || !length(rl))
    stop("drawsArray(params = \"components\") needs relabel() first.\n",
         "  Component draws are not identified before relabelling: the ",
         "mixture likelihood is invariant to label permutations, so ",
         "'muTilde[1]' refers to different components in different chains ",
         "and iterations. Convergence diagnostics computed on raw component ",
         "traces look valid and mean nothing. Run fit <- relabel(fit), or ",
         "use params = \"invariant\" for label-invariant functionals.",
         call. = FALSE)
  skip <- c("method", "modalK", "nDraws", "idx", "permutations", "summary",
            "O", "O_mean", "canonicalFraction")
  keep <- setdiff(names(rl), skip)
  colsL <- list()
  for (nm in keep) {
    v <- rl[[nm]]
    if (is.matrix(v)) {                            # m x K
      for (k in seq_len(ncol(v)))
        colsL[[sprintf("%s[%d]", nm, k)]] <- v[, k]
    } else if (is.array(v) && length(dim(v)) == 3L) {   # m x K x d
      for (k in seq_len(dim(v)[2L])) for (j in seq_len(dim(v)[3L]))
        colsL[[sprintf("%s[%d,%d]", nm, k, j)]] <- v[, k, j]
    }
  }
  if (!length(colsL))
    stop("No array-valued relabelled parameters found.", call. = FALSE)
  M <- do.call(cbind, colsL)
  # Conditioning on the modal K leaves unequal draws per chain, so the chain
  # dimension is honestly collapsed to 1 (see roxygen above).
  array(M, c(nrow(M), 1L, ncol(M)),
        dimnames = list(NULL, NULL, colnames(M)))
}

#' Observed data and posterior predictive replicates for graphical PPC
#'
#' Packages \code{y} and \code{yrep} in the shapes \code{bayesplot}'s
#' \code{ppc_*} functions consume, e.g.
#' \code{bayesplot::ppc_dens_overlay(d$y, d$yrep)}.
#'
#' @param fit A \code{FitResult} from a clustering fit.
#' @param ndraws Number of replicate draws (rows of \code{yrep}).
#' @param margin For multivariate fits, which data dimension to extract
#'   (bayesplot's PPC graphics are univariate). Ignored for univariate fits.
#' @param seed RNG seed passed to \code{\link{posteriorPredict}}.
#' @return A list with \code{y} (numeric vector) and \code{yrep}
#'   (\code{ndraws x n} matrix).
#' @seealso \code{\link{posteriorPredict}}, \code{\link{ppCheck}}.
#' @export
ppcData <- function(fit, ndraws = 100, margin = 1L, seed = 1L) {
  if (!methods::is(fit, "FitResult"))
    stop("ppcData() expects a FitResult.", call. = FALSE)
  yrep <- posteriorPredict(fit, ndraws = ndraws, seed = seed)
  y <- fit@data
  if (length(dim(yrep)) == 3L) {
    d <- dim(yrep)[3L]
    if (margin < 1L || margin > d)
      stop("margin must be in 1..", d, " for this fit.", call. = FALSE)
    list(y = as.numeric(y[, margin]), yrep = yrep[, , margin])
  } else {
    list(y = as.numeric(y), yrep = yrep)
  }
}
