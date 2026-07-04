## ---------------------------------------------------------------------------
## utils-labelswitch.R
##
## Label switching correction. Raw per-component
## posterior summaries are meaningless because the mixture likelihood is
## invariant to relabelling. We:
##   1. condition on the modal number of OCCUPIED clusters (relabelling
##      algorithms assume a fixed K; under a DPM K varies across iterations);
##   2. recode allocations to 1..K (the permutation derivation depends ONLY on
##      the allocations, hence is identical for univariate and multivariate);
##   3. call label.switching (ECR-ITERATIVE-1 by default, the pivot-free
##      iterative ECR of Papastamoulis & Iliopoulos 2010) to obtain the
##      per-iteration label permutations;
##   4. delegate parameter permutation + the component summary to the spec via
##      relabelComponents() -- so multivariate covariance handling lives in
## NormalMvSpec, not here (extensibility).
##
## We WRAP the label.switching package's ECR algorithm rather than
## re-implementing it (reference Papastamoulis 2016).
## ---------------------------------------------------------------------------

#' Correct label switching in a fitted mixture
#'
#' Post-hoc relabelling of MCMC output so that per-component posterior
#' summaries are meaningful. Conditions on the modal number of occupied
#' clusters, then applies an algorithm from the \pkg{label.switching} package.
#'
#' @param fit A \code{\linkS4class{FitResult}}.
#' @param method One of \code{"ECR-ITERATIVE-1"} (default) or \code{"ECR"}.
#'   The default is the pivot-free iterative ECR; \code{"ECR"} uses the
#'   highest-posterior allocation as the pivot.
#' @param ... Reserved.
#' @return The \code{fit} with its \code{relabeled} slot populated.
#'
#' @references
#' Papastamoulis, P., & Iliopoulos, G. (2010). An artificial allocations based
#' solution to the label switching problem. \emph{JCGS}, 19(2), 313--331.
#' \doi{10.1198/jcgs.2010.09008}
#'
#' Papastamoulis, P. (2016). label.switching: An R package for dealing with the
#' label switching problem in MCMC outputs. \emph{JSS, Code Snippets}, 69(1).
#' \doi{10.18637/jss.v069.c01}
#' @export
setGeneric("relabel", function(fit, method = "ECR-ITERATIVE-1", ...) {
  standardGeneric("relabel")
})

#' @describeIn relabel Relabelling for a fitted result.
#' @export
setMethod("relabel", "FitResult",
  function(fit, method = "ECR-ITERATIVE-1", ...) {
    method <- match.arg(method, c("ECR-ITERATIVE-1", "ECR"))

    Kpost <- fit@Kposterior
    kt <- sort(table(Kpost), decreasing = TRUE)
    modalK <- as.integer(names(kt)[1])
    idx <- which(Kpost == modalK)
    m <- length(idx)

    if (m < 20L)
      warning("Only ", m, " posterior draws have the modal #clusters (K = ",
              modalK, "); relabelled summaries may be unstable. Consider a ",
              "longer chain or a different K_max.", call. = FALSE)

    alloc <- fit@clusterAllocation[idx, , drop = FALSE]
    n <- ncol(alloc)

    # Recode each retained iteration's occupied labels to 1..modalK and record
    # the (sorted) occupied labels so the spec can align its parameters. The
    # per-row sort(unique)/match is vectorised: a presence matrix (one tabulate
    # over combined row-label ids), a row-wise cumulative sum giving each
    # label's rank among the row's occupied labels, then one matrix-index
    # lookup. Reproduces match(alloc[t, ], sort(unique(alloc[t, ]))) exactly.
    Lmax <- max(alloc)
    pres <- .rowPresence(alloc, Lmax)
    rk  <- matrix(0L, nrow = m, ncol = Lmax)
    acc <- integer(m)
    for (l in seq_len(Lmax)) {                # Lmax (small) vectorised passes
      acc <- acc + pres[, l]
      rk[, l] <- acc
    }
    z <- matrix(rk[cbind(rep(seq_len(m), times = n), as.vector(alloc))],
                nrow = m, ncol = n)
    storage.mode(z) <- "integer"
    occList <- lapply(seq_len(m), function(t) which(pres[t, ]))

    # Derive label permutations (the wrapped label.switching ECR algorithm).
    if (modalK == 1L) {
      perms <- matrix(1L, nrow = m, ncol = 1L)
    } else {
      perms <- switch(method,
        "ECR-ITERATIVE-1" = label.switching::ecr.iterative.1(
            z = z, K = modalK)$permutations,
        "ECR" = label.switching::ecr(
            zpivot = z[1L, ], z = z, K = modalK)$permutations
      )
    }

    # Mixing weights from relabelled occupancy (distribution-independent),
    # vectorised: one tabulate for all rows, one matrix-index permute.
    wMat <- matrix(tabulate(seq_len(m) + m * (as.vector(z) - 1L),
                            nbins = m * modalK), nrow = m) / n
    wPerm <- matrix(wMat[cbind(rep(seq_len(m), times = modalK),
                               as.vector(perms))],
                    nrow = m, ncol = modalK)

    # Delegate parameter permutation + component summary to the spec.
    comp <- relabelComponents(fit@distSpec, fit@paramTrace, idx, occList,
                              perms, modalK, wPerm)

    # Add a posterior median next to every scalar-per-component mean, so the
    # summary reports mean, median, and a 95% credible interval. Works for any
    # family whose relabelled draws are returned as m x modalK matrices;
    # array-valued parameters are summarised by their own spec method.
    comp <- .augmentSummaryMedian(comp, m, modalK)

    fit@relabeled <- c(
      list(method = method, modalK = modalK, nDraws = m, idx = idx,
           permutations = perms, weight = wPerm),
      comp
    )
    fit
  }
)

# Insert a posterior median column after each scalar-per-component mean.
.augmentSummaryMedian <- function(comp, nDraws, modalK) {
  summ <- comp$summary
  if (is.null(summ) || !is.data.frame(summ)) return(comp)
  for (nm in names(comp)) {
    M <- comp[[nm]]
    if (is.matrix(M) && nrow(M) == nDraws && ncol(M) == modalK &&
        paste0(nm, "_mean") %in% names(summ) &&
        !paste0(nm, "_med") %in% names(summ)) {
      summ[[paste0(nm, "_med")]] <- apply(M, 2L, stats::median)
    }
  }
  ord <- character(0)
  for (cn in names(summ)) {
    if (grepl("_med$", cn)) next
    ord <- c(ord, cn)
    medcn <- sub("_mean$", "_med", cn)
    if (grepl("_mean$", cn) && medcn %in% names(summ)) ord <- c(ord, medcn)
  }
  comp$summary <- summ[, ord, drop = FALSE]
  comp
}
