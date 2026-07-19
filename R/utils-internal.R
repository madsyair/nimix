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
  # Fast path for vector nodes (e.g. the allocation xi[1..n]): a single block
  # column assignment instead of one indexed write per column.
  if (length(dimSizes) == 1L && all(lengths(info$idx) == 1L)) {
    arr[, unlist(info$idx, use.names = FALSE)] <- samples[, info$cols]
    return(arr)
  }
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
# DPM truncation default. K_max is the dCRP truncation level L: the sampler
# errors if the occupied-cluster count ever needs to exceed it, so the default
# must sit comfortably above the expected number of clusters (which grows with
# n). This gives generous, bounded headroom; an explicit K_max overrides it.
.defaultTruncation <- function(n) {
  k <- min(40L, max(20L, as.integer(ceiling(n / 10))))
  as.integer(min(k, n))
}

# Default headroom ratio for the dispersed cluster initialisation: the k-means /
# E-step start seeds at most initRatio * count clusters (count = the truncation
# level L = K_max for the DPM, or K for the fixed-K engine), leaving headroom for
# the early transient before the chain settles. 0.8 guarantees >= 0.2 * K_max
# free slots; users can override via mcmcControl$initRatio.
.DEFAULT_INIT_RATIO <- 0.8

# Above this ratio the dispersed start crowds the truncation: still allowed (an
# advanced user may want it), but warned, because for the DPM it can breach the
# truncation by leaving almost no room for the early CRP transient.
.INIT_RATIO_WARN <- 0.95

# Resolve and validate mcmcControl$initRatio once per run (in runEngine), before
# it is threaded to componentInits. It must stay in the open interval (0, 1): a
# ratio of 0 seeds no clusters, and 1 (or more) seeds the entire cap, leaving no
# headroom at all. Ratios in [.INIT_RATIO_WARN, 1) are accepted but warned, since
# they leave little headroom below the truncation.
.resolveInitRatio <- function(mcmcControl) {
  r <- mcmcControl$initRatio
  if (is.null(r)) return(.DEFAULT_INIT_RATIO)
  if (!is.numeric(r) || length(r) != 1L || is.na(r) || r <= 0 || r >= 1)
    stop("mcmcControl$initRatio must be a single number in (0, 1); the default ",
         .DEFAULT_INIT_RATIO, " leaves headroom below the truncation. Got: ",
         paste(utils::head(r, 3L), collapse = ", "), call. = FALSE)
  if (r >= .INIT_RATIO_WARN)
    warning("mcmcControl$initRatio = ", r, " leaves little headroom below the ",
            "truncation K_max; the dispersed start may breach it for the DPM. ",
            "Consider a ratio <= 0.9, or a larger K_max.", call. = FALSE)
  as.numeric(r)
}

# Extract the (already validated) initRatio passed through `...` to a
# componentInits method, falling back to the default for direct calls (e.g. unit
# tests) that do not supply it.
.initRatioArg <- function(...) {
  r <- list(...)$initRatio
  if (is.null(r)) .DEFAULT_INIT_RATIO else r
}

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

# Run `fun()` under a fixed RNG seed, restoring the global RNG state afterwards.
# Used to make the (k-means) dispersed initialisation reproducible given the
# fit's `seed`, so that repeated and reused fits coincide bit-for-bit, without
# disturbing the caller's global random stream.
.withSeed <- function(seed, fun) {
  if (is.null(seed)) return(fun())
  hadSeed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (hadSeed) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE)
  }
  set.seed(seed)
  fun()
}

# Rank-normalized split-Rhat for a set of chains (each a numeric vector).
# Follows Vehtari, Gelman, Simpson, Carpenter & Buerkner (2021): split every
# chain in half, rank-normalize the pooled draws, then form the potential scale
# reduction factor. Returns NA when there are too few chains/draws or no
# between-chain variation (e.g. a fixed quantity). Applied only to
# label-invariant scalars (occupied-cluster count, concentration) so that
# label switching does not corrupt the diagnostic.
.splitRhat <- function(chains) {
  chains <- lapply(chains, function(x) x[is.finite(x)])
  chains <- chains[lengths(chains) >= 4L]
  if (length(chains) < 2L) return(NA_real_)
  half <- floor(min(lengths(chains)) / 2L)
  if (half < 2L) return(NA_real_)
  split <- unlist(lapply(chains, function(x)
    list(x[seq_len(half)], x[half + seq_len(half)])), recursive = FALSE)
  alld <- unlist(split)
  if (stats::var(alld) == 0) return(NA_real_)      # constant -> undefined/1
  z  <- stats::qnorm((rank(alld, ties.method = "average") - 0.5) / length(alld))
  zc <- split(z, rep(seq_along(split), each = half))
  n  <- half
  B  <- n * stats::var(vapply(zc, mean, numeric(1)))
  W  <- mean(vapply(zc, stats::var, numeric(1)))
  if (W == 0) return(NA_real_)
  sqrt((((n - 1) / n) * W + B / n) / W)
}

# Effective sample size summed across chains (coda per chain). A conservative,
# transparent multi-chain ESS: total independent information, reported next to
# Rhat (which flags between-chain disagreement the sum would otherwise hide).
.sumESS <- function(chains) {
  e <- vapply(chains, function(x) {
    v <- x[is.finite(x)]
    if (length(v) < 4L || stats::var(v) == 0) return(0)
    as.numeric(coda::effectiveSize(coda::as.mcmc(v)))
  }, numeric(1))
  sum(e)
}

# Assemble the multi-chain diagnostic list from per-chain label-invariant
# scalars. `perChainK` is a list of occupied-cluster-count vectors; `perChainAlpha`
# an optional list of concentration draws (DPM only). Rhat > 1.1 is the standard
# not-yet-converged flag (see the mixing-diagnostics guidance in the design
# notes).
.multiChainDiag <- function(perChainK, perChainAlpha = NULL, perChainBeta = NULL,
                            perChainEntropy = NULL) {
  d <- list(nchains = length(perChainK),
            RhatK   = .splitRhat(perChainK),
            essK    = .sumESS(perChainK))
  if (!is.null(perChainAlpha)) {
    d$RhatAlpha <- .splitRhat(perChainAlpha)
    d$essAlpha  <- .sumESS(perChainAlpha)
  }
  if (!is.null(perChainBeta)) {
    d$RhatBeta <- .splitRhat(perChainBeta)
    d$essBeta  <- .sumESS(perChainBeta)
    d$betaMean <- mean(unlist(perChainBeta))
  }
  if (!is.null(perChainEntropy)) {
    d$RhatEntropy <- .splitRhat(perChainEntropy)
    d$essEntropy  <- .sumESS(perChainEntropy)
  }
  ## Full Vehtari et al. (2021) table over the label-invariant functionals:
  ## rank-normalized split-Rhat, folded split-Rhat, bulk-ESS, tail-ESS.
  fns <- list(K = perChainK, alpha = perChainAlpha, beta = perChainBeta,
              entropy = perChainEntropy)
  fns <- fns[!vapply(fns, is.null, logical(1))]
  if (length(fns)) {
    d$functionals <- data.frame(
      functional = names(fns),
      Rhat       = vapply(fns, .splitRhat,  numeric(1)),
      foldedRhat = vapply(fns, .foldedRhat, numeric(1)),
      bulkESS    = vapply(fns, .bulkESS,    numeric(1)),
      tailESS    = vapply(fns, .tailESS,    numeric(1)),
      row.names  = NULL)
  }
  d
}

# Vectorised row-wise label utilities (performance, v0.5.0). For an m x n
# integer label matrix with values in 1..L, .rowPresence() returns the m x L
# logical presence matrix via a single tabulate() over combined (row, label)
# ids -- O(mn + mL) with no per-row R overhead -- and .rowDistinct() the number
# of distinct labels per row (the occupied-cluster count per posterior draw).
# Both reproduce the per-row length(unique(.)) result exactly.
.rowPresence <- function(a, L) {
  m <- nrow(a)
  comb <- seq_len(m) + m * (as.integer(a) - 1L)
  matrix(tabulate(comb, nbins = m * L) > 0L, nrow = m)
}

.rowDistinct <- function(a, L) {
  as.integer(rowSums(.rowPresence(a, L)))
}

# --- Vehtari et al. (2021) convergence suite ----------------------------------
# Rank-normalized bulk-ESS, tail-ESS, and folded split-Rhat (Vehtari, Gelman,
# Simpson, Carpenter & Buerkner 2021, Bayesian Analysis 16(2), Sections 3-4),
# complementing the rank-normalized split-Rhat above. Applied in nimix ONLY to
# label-invariant functionals (occupied K, concentration alpha, MRF beta,
# allocation entropy): per-component parameter traces are not identified under
# label switching, so chain diagnostics on them are meaningless.

.splitHalves <- function(chains) {
  chains <- lapply(chains, function(x) x[is.finite(x)])
  chains <- chains[lengths(chains) >= 4L]
  if (length(chains) < 1L) return(NULL)
  half <- floor(min(lengths(chains)) / 2L)
  if (half < 2L) return(NULL)
  unlist(lapply(chains, function(x)
    list(x[seq_len(half)], x[half + seq_len(half)])), recursive = FALSE)
}

# Bulk-ESS: ESS of the rank-normalized split chains (Vehtari et al. 2021,
# Section 4.1), summed across split chains.
.bulkESS <- function(chains) {
  sp <- .splitHalves(chains)
  if (is.null(sp)) return(NA_real_)
  alld <- unlist(sp)
  if (stats::var(alld) == 0) return(NA_real_)
  z <- stats::qnorm((rank(alld, ties.method = "average") - 0.5) / length(alld))
  zc <- split(z, rep(seq_along(sp), lengths(sp)))
  sum(vapply(zc, function(v)
    as.numeric(coda::effectiveSize(coda::as.mcmc(v))), numeric(1)))
}

# Tail-ESS: minimum of the ESS of the 5% and 95% quantile indicators
# (Vehtari et al. 2021, Section 4.3), quantiles taken on the pooled draws.
.tailESS <- function(chains) {
  sp <- .splitHalves(chains)
  if (is.null(sp)) return(NA_real_)
  alld <- unlist(sp)
  q <- stats::quantile(alld, c(0.05, 0.95), names = FALSE, type = 1)
  essInd <- function(thr, lower) {
    sum(vapply(sp, function(v) {
      ind <- if (lower) as.numeric(v <= thr) else as.numeric(v >= thr)
      if (stats::var(ind) == 0) return(0)
      as.numeric(coda::effectiveSize(coda::as.mcmc(ind)))
    }, numeric(1)))
  }
  min(essInd(q[1], TRUE), essInd(q[2], FALSE))
}

# Folded split-Rhat: split-Rhat of |x - median| (Vehtari et al. 2021,
# Section 3.2) -- sensitive to scale (variance) disagreement between chains
# that the location-based Rhat can miss.
.foldedRhat <- function(chains) {
  alld <- unlist(lapply(chains, function(x) x[is.finite(x)]))
  if (length(alld) < 8L) return(NA_real_)
  med <- stats::median(alld)
  .splitRhat(lapply(chains, function(x) abs(x - med)))
}

# Per-iteration allocation entropy H_t = -sum_k p_k(t) log p_k(t): a
# label-invariant summary of how evenly the current partition spreads mass
# across components; a standard mixture functional for chain diagnostics.
.allocEntropy <- function(zMat, K) {
  apply(zMat, 1L, function(z) {
    p <- tabulate(z, nbins = K) / length(z)
    p <- p[p > 0]
    -sum(p * log(p))
  })
}

# --- initial-allocation helpers ------------------------------------------------
#
# Shared across every componentInits method so that a new init strategy is
# written once, not eleven times. .initClusters returns a hard allocation
# vector (integer, one label per row) or NULL on failure -- exactly the
# contract the inline kmeans blocks already expect, so callers keep their
# center/variance extraction unchanged.
#
# initMethod:
#   "kmeans"  -> stats::kmeans (the default; unchanged behaviour)
#   "spread"  -> univariate only: split by |y - median(y)| quantiles, which
#                separates components by SCALE rather than location. This is
#                the one case measured to defeat k-means (heterogeneous
#                variance, overlapping means; see kajian_inisialisasi_*). For
#                a multivariate response it has no natural analogue, so it
#                falls back to kmeans rather than guessing.
#   "single"  -> handled by the caller (all-one-cluster); never reaches here.
.initClusters <- function(y, k0, initMethod = "kmeans") {
  if (k0 < 2L) return(NULL)
  isMv <- is.matrix(y) && ncol(y) > 1L

  if (identical(initMethod, "spread") && !isMv) {
    yv <- as.numeric(y)
    dev <- abs(yv - stats::median(yv))
    # k0 bands of increasing deviation from the centre: band 1 is the tight
    # core, band k0 the diffuse tail. Ties in dev are broken by rank so every
    # band is populated whenever there are enough distinct values.
    br <- stats::quantile(dev, probs = seq(0, 1, length.out = k0 + 1L),
                          names = FALSE)
    br[1L] <- -Inf; br[length(br)] <- Inf
    cl <- as.integer(cut(dev, breaks = br, labels = FALSE,
                         include.lowest = TRUE))
    if (length(unique(cl)) < k0) return(NULL)   # too few distinct: let caller fall back
    return(cl)
  }

  # default / multivariate spread fallback: kmeans
  km <- tryCatch(stats::kmeans(y, centers = k0, nstart = 5L),
                 error = function(e) NULL)
  if (is.null(km)) NULL else as.integer(km$cluster)
}
