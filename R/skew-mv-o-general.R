## skew-mv-o-general.R -----------------------------------------------------------
## Machinery for the Ferreira-Steel orthogonal factor O when m > 2.
##
## Decomposition (FS 2007, Appendix A.1-A.2, after Golub & Van Loan 1989).
## Write O = O_{theta^m} x ... x O_{theta^2}, where theta^j has j-1 angles,
##   v^j_1 = sin(theta^j_1),
##   v^j_i = prod_{l<i} cos(theta^j_l) * sin(theta^j_i)   for i < j,
##   v^j_j = prod_{l<j} cos(theta^j_l),
## H_{theta^j} = I_j - 2 v^j (v^j)' / (v^j)'v^j, and
## O_{theta^j} = blockdiag(I_{m-j}, H_{theta^j}). Each H has determinant -1, so
## |O| = (-1)^(m-1) = (-1)^(m+1), as FS require. In total m(m-1)/2 angles.
##
## Restriction (8) is a CANONICALISATION, not a sampling constraint. FS say
## that confining theta^j to the angle box Theta^j puts O in O_m; it does not.
## Sampling angles uniformly from the box and testing (8) directly, the fraction
## of draws that satisfy it is 0.245 (m = 2), 0.069 (m = 3), 0.007 (m = 4).
## Constraining the sampler to a 0.7% slice of its own prior would mix badly.
##
## What is true, and what we use: among the signed row permutations P of A with
## |P| = +1, exactly one PO satisfies (8) -- verified exhaustively for
## m = 2, 3, 4. The m! 2^m ambiguity of A's rows is therefore label switching in
## the dimension index, and nimix already prefers post-hoc relabelling over
## ordering constraints. So we sample the angles unconstrained on the box and
## map each posterior draw to its unique representative.
##
## Under A -> PA: Sigma = A'A and U = chol(Sigma) are unchanged, eps -> P eps,
## and since the FS mechanism obeys p(-e | gamma) = p(e | 1/gamma),
##   gamma*_i = gamma_{perm(i)}      if row i keeps its sign,
##   gamma*_i = 1 / gamma_{perm(i)}  if row i is negated,
##   nu*_i    = nu_{perm(i)}.
## For m = 2 the (8)-satisfying set is exactly theta in (-pi/8, pi/8), which is
## the prior support already used by SkewNormalMvOSpec: the general treatment
## reduces to the bivariate one rather than replacing it.

#' @include class-DistributionSpec.R
NULL

# Number of Householder angles for dimension m.
.nAngles <- function(m) as.integer(m * (m - 1L) / 2L)

# Split a flat angle vector into the per-j blocks theta^j, j = 2..m.
.angleBlocks <- function(theta, m) {
  out <- vector("list", m)
  pos <- 1L
  for (j in 2:m) {
    len <- j - 1L
    out[[j]] <- theta[pos:(pos + len - 1L)]
    pos <- pos + len
  }
  out
}

# Unit vector in polar coordinates (FS Appendix A.1).
.vPolar <- function(th) {
  j <- length(th) + 1L
  v <- numeric(j)
  v[1L] <- sin(th[1L])
  if (j > 2L)
    for (i in 2:(j - 1L)) v[i] <- prod(cos(th[1:(i - 1L)])) * sin(th[i])
  v[j] <- prod(cos(th))
  v
}

.householder <- function(th) {
  v <- .vPolar(th)
  diag(length(v)) - 2 * outer(v, v) / sum(v * v)
}

#' Orthogonal factor from Householder angles
#'
#' Builds \eqn{O = O_{\theta^m} \cdots O_{\theta^2}} from a flat vector of
#' \eqn{m(m-1)/2} Householder angles (Ferreira & Steel 2007, Appendix A).
#'
#' @param theta Numeric vector of \code{m * (m - 1) / 2} angles.
#' @param m Dimension.
#' @return An \code{m x m} orthogonal matrix with determinant \eqn{(-1)^{m+1}}.
#' @references Ferreira & Steel (2007), Statistica Sinica 17, 505--529.
#' @export
orthogonalFactor <- function(theta, m) {
  stopifnot(m >= 2L, length(theta) == .nAngles(m))
  blocks <- .angleBlocks(theta, m)
  O <- diag(m)
  for (j in m:2) {
    B <- diag(m)
    B[(m - j + 1L):m, (m - j + 1L):m] <- .householder(blocks[[j]])
    O <- O %*% B
  }
  O
}

# FS restriction (8) on the first column of O:
#   O11 > -O_m1 > -O_(m-1)1 > ... > |O21| > 0,  |O| = (-1)^(m+1).
.restriction8 <- function(O) {
  m <- nrow(O)
  c1 <- O[, 1L]
  if (!(c1[1L] > 0)) return(FALSE)
  chain <- c(c1[1L], if (m > 2L) -c1[m:3L], abs(c1[2L]))
  all(diff(chain) < 0) && abs(c1[2L]) > 1e-12
}

# All signed row permutations P with |P| = +1 (these preserve |O|).
.signedPerms <- function(m) {
  perms <- .permList(m)
  signs <- as.matrix(expand.grid(rep(list(c(1, -1)), m)))
  out <- list()
  for (p in perms) for (r in seq_len(nrow(signs))) {
    P <- matrix(0, m, m)
    for (i in seq_len(m)) P[i, p[i]] <- signs[r, i]
    if (abs(det(P) - 1) < 1e-9) out[[length(out) + 1L]] <- P
  }
  out
}

.permList <- function(m) {
  if (m == 1L) return(list(1L))
  out <- list()
  for (i in seq_len(m))
    for (rest in .permList(m - 1L)) {
      others <- setdiff(seq_len(m), i)
      out[[length(out) + 1L]] <- c(i, others[rest])
    }
  out
}

#' Canonical representative of an FS orthogonal factor
#'
#' Maps \code{(O, gamma, nu)} to the unique signed row permutation satisfying
#' Ferreira & Steel's identifiability restriction (8). Exactly one such
#' representative exists. \code{Sigma} is invariant under the map.
#'
#' @param O Orthogonal matrix with determinant \eqn{(-1)^{m+1}}.
#' @param gamma Positive skewness vector, length \code{m}.
#' @param nu Optional positive degrees-of-freedom vector, length \code{m}.
#' @return A list with the canonical \code{O}, \code{gamma}, and (if supplied)
#'   \code{nu}. If no representative is found (a measure-zero event, e.g. a zero
#'   in the first column of \code{O}), the inputs are returned unchanged with
#'   \code{canonical = FALSE}.
#' @references Ferreira & Steel (2007), Statistica Sinica 17, 505--529.
#' @export
canonicaliseO <- function(O, gamma, nu = NULL) {
  m <- nrow(O)
  stopifnot(length(gamma) == m, is.null(nu) || length(nu) == m)
  for (P in .signedPerms(m)) {
    Oc <- P %*% O
    if (.restriction8(Oc)) {
      # row i of P picks column perm(i) of the old basis with sign s_i
      perm <- apply(P, 1L, function(r) which(r != 0))
      sgn  <- vapply(seq_len(m), function(i) sign(P[i, perm[i]]), numeric(1))
      g <- gamma[perm]
      g[sgn < 0] <- 1 / g[sgn < 0]      # p(-e | gamma) = p(e | 1/gamma)
      out <- list(O = Oc, gamma = g, canonical = TRUE)
      if (!is.null(nu)) out$nu <- nu[perm]
      return(out)
    }
  }
  out <- list(O = O, gamma = gamma, canonical = FALSE)
  if (!is.null(nu)) out$nu <- nu
  out
}
