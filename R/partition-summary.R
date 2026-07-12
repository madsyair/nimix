## partition-summary.R --------------------------------------------------------
## Label-free partition summaries for mixture fits.
##
## relabel() summarises COMPONENT PARAMETERS, and to do so it must condition on
## the iterations whose number of occupied clusters equals the posterior mode
## (a varying K has no fixed set of components to align). On DPM fits that can
## discard a substantial share of the draws -- measured at 34% on a two-cluster
## example. The functions here summarise the PARTITION instead, and need
## neither a fixed K nor any labelling, so they use every draw:
##
##  * psm(fit)              posterior similarity matrix  S_ij = P(z_i = z_j | y)
##  * binderPartition(fit)  point-estimate partition minimising the expected
##                          Binder loss over the sampled partitions
##
## The Binder (1978) loss between a candidate partition c and the truth z
## counts pairwise disagreements. Its posterior expectation depends on the
## draws only through S, and minimising it over candidates is equivalent to
## Dahl's (2006) least-squares criterion:
##     c* = argmin_c  sum_{i<j} ( I(c_i = c_j) - S_ij )^2 .
## Searching over the sampled partitions themselves (Dahl 2006) keeps the
## optimum inside the support actually visited by the chain.
##
## Both are complements to relabel(), not replacements: relabel() answers
## "what are the component parameters", these answer "which observations
## belong together".

#' @include class-FitResult.R
NULL

#' Posterior similarity matrix
#'
#' Returns the \eqn{n \times n} matrix \eqn{S_{ij} = \Pr(z_i = z_j \mid y)},
#' estimated as the fraction of posterior draws in which observations \eqn{i}
#' and \eqn{j} share a cluster. The quantity is invariant to label
#' permutations and to the number of occupied clusters, so \emph{every} draw
#' contributes -- unlike \code{\link{relabel}}, which must condition on the
#' modal number of clusters before component parameters can be aligned.
#'
#' @param fit A \code{FitResult}.
#' @return A symmetric matrix with unit diagonal.
#' @references Binder, D. A. (1978), Biometrika 65, 31--38.
#' @seealso \code{\link{binderPartition}} for a point-estimate partition,
#'   \code{\link{relabel}} for component-parameter summaries.
#' @export
psm <- function(fit) {
  if (!methods::is(fit, "FitResult"))
    stop("psm() expects a FitResult.", call. = FALSE)
  A <- fit@clusterAllocation
  m <- nrow(A); n <- ncol(A)
  S <- matrix(0, n, n)
  for (t in seq_len(m)) {
    z <- A[t, ]
    S <- S + outer(z, z, "==")
  }
  S <- S / m
  dimnames(S) <- NULL
  S
}

#' Binder-loss point partition (Dahl's least-squares criterion)
#'
#' Selects, among the partitions actually visited by the chain, the one
#' minimising the posterior expected Binder loss -- equivalently, the draw
#' whose pairwise co-clustering matrix is closest (in squared error) to the
#' posterior similarity matrix (Dahl 2006). All draws inform the similarity
#' matrix; none are discarded.
#'
#' @param fit A \code{FitResult}.
#' @param S Optional precomputed \code{\link{psm}} matrix.
#' @return A list with \code{partition} (integer vector, labels recoded to
#'   \code{1..K}), \code{K}, \code{draw} (index of the selected iteration),
#'   \code{score} (the least-squares criterion value), and \code{psm}.
#' @references Binder, D. A. (1978), Biometrika 65, 31--38. Dahl, D. B.
#'   (2006), in \emph{Bayesian Inference for Gene Expression and Proteomics},
#'   Cambridge University Press, 201--218.
#' @seealso \code{\link{psm}}, \code{\link{relabel}}.
#' @export
binderPartition <- function(fit, S = NULL) {
  if (!methods::is(fit, "FitResult"))
    stop("binderPartition() expects a FitResult.", call. = FALSE)
  if (is.null(S)) S <- psm(fit)
  A <- fit@clusterAllocation
  m <- nrow(A)
  best <- Inf; bestT <- 1L
  for (t in seq_len(m)) {
    z <- A[t, ]
    D <- outer(z, z, "==") - S
    sc <- sum(D * D)
    if (sc < best) { best <- sc; bestT <- t }
  }
  zBest <- A[bestT, ]
  zBest <- as.integer(factor(zBest, levels = unique(zBest)))
  list(partition = zBest, K = length(unique(zBest)),
       draw = bestT, score = best, psm = S)
}
