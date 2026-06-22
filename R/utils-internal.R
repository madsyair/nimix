## ---------------------------------------------------------------------------
## utils-internal.R
##
## Small internal helpers shared across specs/engine. None are exported.
##
## Two groups:
##   (1) data-shape helpers so callers never branch on vector-vs-matrix inline;
##   (2) multivariate-normal math (density + draws) implemented locally so we do
## NOT add a dependency on mvtnorm/MASS (no new
##       dependencies without explicit justification). The Wishart draw reuses
##       stats::rWishart, which ships with base R.
## ---------------------------------------------------------------------------

# Number of observations regardless of univariate (vector) / multivariate
# (matrix, one row per observation) storage.
.nObs <- function(data) if (is.matrix(data)) nrow(data) else length(data)

# Data dimension d (1 for a univariate vector).
.dataDimOf <- function(data) if (is.matrix(data)) ncol(data) else 1L

# Parse the MCMC columns for a NIMBLE node (e.g. "muTilde", "covTilde") into a
# dense array of shape c(nIter, dimSizes), using the bracket indices in the
# column names. Robust to NIMBLE's column ordering because each column is placed
# by its own parsed indices rather than by position.
.nodeColInfo <- function(samples, node) {
  cn  <- colnames(samples)
  pat <- paste0("^", node, "\\[")
  hit <- grep(pat, cn)
  if (!length(hit)) return(NULL)
  inside <- sub(paste0("^", node, "\\[(.*)\\]$"), "\\1", cn[hit])
  idx <- lapply(strsplit(inside, ","), function(s) as.integer(trimws(s)))
  list(cols = hit, idx = idx)
}

.nodeToArray <- function(samples, node, dimSizes) {
  info <- .nodeColInfo(samples, node)
  m <- nrow(samples)
  if (is.null(info))
    stop("node '", node, "' not found in MCMC samples.", call. = FALSE)
  arr <- array(NA_real_, dim = c(m, dimSizes))
  for (t in seq_along(info$cols)) {
    ix <- info$idx[[t]]
    mIdx <- cbind(seq_len(m),
                  matrix(ix, nrow = m, ncol = length(ix), byrow = TRUE))
    arr[mIdx] <- samples[, info$cols[t]]
  }
  arr
}

# Cholesky with a tiny ridge fallback so near-singular sample covariances (e.g.
# a cluster with barely more members than dimensions) do
# not abort a draw or a density evaluation.
.cholSafe <- function(S) {
  d <- nrow(S)
  out <- tryCatch(chol(S), error = function(e) NULL)
  if (!is.null(out)) return(out)
  ridge <- 1e-6 * mean(diag(S)) + 1e-8
  chol(S + diag(ridge, d))
}

# One draw from N(mu, Sigma) via its Cholesky factor.
.rmvnorm1 <- function(mu, Sigma) {
  d <- length(mu)
  R <- .cholSafe(Sigma)
  as.numeric(mu + t(R) %*% stats::rnorm(d))
}

# Multivariate normal density (not log) at a single point x given mean mu and
# covariance Sigma. Uses the Cholesky factor for a stable determinant/quadratic
# form. Standard MVN density (see Frühwirth-Schnatter 2006,
# ).
.dmvnorm <- function(x, mu, Sigma) {
  d <- length(mu)
  R <- .cholSafe(Sigma)                 # Sigma = R'R, R upper-triangular
  z <- backsolve(R, x - mu, transpose = TRUE)
  logdet <- 2 * sum(log(diag(R)))
  exp(-0.5 * (d * log(2 * pi) + logdet + sum(z^2)))
}

# Multivariate Student-t density (location mu, scale matrix Sigma, df). This is
# the marginal of the multivariate Normal-Gamma scale mixture and the kernel of
# the direct multivariate-t component.
.dmvt <- function(x, mu, Sigma, df) {
  d <- length(mu)
  R <- .cholSafe(Sigma)
  z <- backsolve(R, x - mu, transpose = TRUE)
  logdet <- 2 * sum(log(diag(R)))
  quad <- sum(z^2)
  ll <- lgamma((df + d) / 2) - lgamma(df / 2) - (d / 2) * log(df * pi) -
    0.5 * logdet - ((df + d) / 2) * log1p(quad / df)
  exp(ll)
}

# Replace the default (random-walk) sampler on the latent precision multipliers
# omega with per-node slice samplers. Marginalising omega is what makes the
# direct Student-t mix well; when omega is kept explicit (the scale-mixture
# route) a slice sampler mixes the partition markedly better than the default
# random walk (van Dyk & Meng 2001 on the cost of augmentation).
.omegaToSlice <- function(conf, model) {
  nodes <- model$expandNodeNames("omega")
  if (length(nodes) == 0L) return(invisible(conf))
  conf$removeSamplers("omega")
  for (nd in nodes) conf$addSampler(target = nd, type = "slice")
  invisible(conf)
}

# Turn a finite-mixture model code into a single-component model for K = 1. A
# one-category Dirichlet / categorical is degenerate and NIMBLE will not build
# it, yet K = 1 is a legitimate baseline (e.g. the null model in a WAIC
# comparison across K). We drop the allocation layer -- the stochastic
# declarations z[i] ~ dcat(...) and weights[1:K] ~ ddirch(...) -- and leave the
# per-observation lookups betaTilde[z[i], ...] in place; z is supplied as the
# constant 1 vector so those lookups fold to the single component. All other
# statements (component priors, latent omega for the scale mixtures, the
# deterministic mean/covariance nodes) are preserved untouched.
.stripMixtureLayer <- function(code) {
  isAllocDecl <- function(stmt) {
    # TRUE for `z[...] ~ ...` or `weights[...] ~ ...`
    if (!is.call(stmt) || length(stmt) < 3L) return(FALSE)
    if (!identical(stmt[[1]], as.name("~"))) return(FALSE)
    lhs <- stmt[[2]]
    is.call(lhs) && identical(lhs[[1]], as.name("[")) &&
      as.character(lhs[[2]]) %in% c("z", "weights")
  }
  filterBlock <- function(block) {
    kept <- list(block[[1]])                       # the `{`
    for (k in seq.int(2L, length(block))) {
      stmt <- block[[k]]
      if (is.call(stmt) && identical(stmt[[1]], as.name("for"))) {
        stmt[[4]] <- filterBlock(stmt[[4]])         # recurse into loop body
        kept[[length(kept) + 1L]] <- stmt
      } else if (!isAllocDecl(stmt)) {
        kept[[length(kept) + 1L]] <- stmt
      }
    }
    as.call(kept)
  }
  filterBlock(code)
}
