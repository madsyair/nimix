#' @include class-EngineConfig.R
#' @include dist-skewistudent-mv-og.R
#' @include dist-skewnormal-mv-og.R
#' @include dist-skewistudent-mv-o.R
#' @include dist-skewnormal-mv-o.R
#' @include dist-skewistudent-mv.R
#' @include dist-skewnormal-mv.R
#' @include dist-jfst.R
#' @include dist-fsst.R
#' @include dist-fossep.R
#' @include dist-fssn.R
#' @include dist-lep.R
#' @include dist-gmsnburr.R
#' @include dist-msnburr.R
NULL

## engine-mrf.R ----------------------------------------------------------------
## Spatially constrained finite mixture (v0.6.0): the latent labels z[1:n]
## follow a Potts Markov random field on a user-supplied neighbourhood graph
## (SpatialWeightSpec), so neighbouring observations favour the same component.
##
## Statistical basis and an honest limitation, stated up front:
## the Potts prior p(z | beta) is proportional to exp(beta * #{edges i~j with
## z_i = z_j}) (Potts 1952; Besag 1974). Its normalising constant ("partition
## function") is intractable, BUT it depends only on beta -- so with beta held
## FIXED (this engine) the constant is absorbed and an unnormalised log-density
## yields exact MCMC for (z, theta). Bayesian estimation of beta itself
## requires that constant and is deliberately deferred to a later release.
## The per-site full conditional is tractable (Besag 1974):
##   p(z_i = k | z_-i, y, theta)  propto  f(y_i | theta_k) *
##                                        exp(beta * #{j ~ i : z_j = k}),
## which the custom sampler below sweeps with single-site Gibbs updates --
## the classical sampler for spatially variant finite mixtures (Blekas et al.
## 2005, Section III).

# --- Potts prior as a user-defined NIMBLE distribution -----------------------
# Unnormalised on purpose (see header). dPottsNimix / rPottsNimix are built and
# registered in the GLOBAL environment by .nimixDefinePotts() (registerDistrib-
# ution.R); defining them here, in the namespace frame, makes NIMBLE fail to
# find rPottsNimix during code generation for the latent label node z.

# --- single-site Gibbs sweep for z under Potts x Normal emissions ------------
# Reads muTilde / s2Tilde / y directly from the model (the node names of the
# NormalUvSpec finite-mixture kernel), so a sweep costs O(n * K + edges) with
# no per-site model$calculate round trips.

#' @keywords internal
.pottsGibbsNormalSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K    <- control$K
    nbrs <- control$nbrs          # n x maxDeg neighbour indices (0-padded)
    nDeg <- control$nDeg          # degree per site
    n    <- control$n
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z  <- model[["z"]]
    mu <- model[["muTilde"]]
    s2 <- model[["s2Tilde"]]
    y  <- model[["y"]]
    beta <- model[["beta"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) {
          for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        }
        lp[k] <- dnorm(y[i], mu[k], sqrt(s2[k]), log = TRUE) + beta * same
      }
      mx <- max(lp)
      p <- exp(lp - mx)
      p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {})
)

# --- model code: NormalUvSpec under the MRF engine ---------------------------
# Mirrors the fixed-K Gaussian kernel, with the Dirichlet-categorical
# allocation replaced by the joint Potts node. There are no mixing weights:
# spatial cohesion takes their place.


# --- engine runner ------------------------------------------------------------

#' @describeIn runEngine Spatially constrained finite mixture (Potts MRF on the
#'   labels, fixed interaction \code{beta}); univariate Gaussian components.
#' @export
setMethod("runEngine", "MRFEngine",
  function(engine, model, mcmcControl = list(), initMethod = "kmeans",
           seed = 1L, verbose = TRUE, ...) {
    spec  <- model@distSpec
    data  <- model@data
    prior <- model@prior
    K     <- model@Kmax                       # fixed number of components
    n     <- .nObs(data)
    d     <- .dataDimOf(data)

    samplerFn <- .mrfSamplerFor(spec)     # S4 dispatch; errors for unsupported
    if (K < 2L)
      stop("method = 'mrf' needs K >= 2 (a one-state Potts field is ",
           "degenerate); use method = 'fixedk' for a single component.",
           call. = FALSE)

    nReg <- nRegions(engine@spatial)
    if (nReg != n)
      stop("spatialWeights has ", nReg, " regions but the data has ", n,
           " observations; they must match one-to-one.", call. = FALSE)

    # Edge list + neighbour matrix (structural constants for dPotts +
    # sampler), taken sparsely: getEdges() is canonically ordered to match
    # the column-major which(upper.tri(A) & A > 0) of the old dense path, so
    # e1/e2/deg/nbrs -- and therefore fits under a fixed seed -- are
    # identical, without ever allocating an n x n matrix.
    E <- getEdges(engine@spatial)
    if (nrow(E) == 0L)
      stop("spatialWeights has no edges; use method = 'fixedk' for an ",
           "unstructured finite mixture.", call. = FALSE)
    e1 <- as.numeric(E[, 1L]); e2 <- as.numeric(E[, 2L])
    deg <- tabulate(c(E[, 1L], E[, 2L]), nbins = n)
    nbrs <- matrix(0, n, max(deg))
    cnt <- integer(n)
    for (k in seq_len(nrow(E))) {
      a <- E[k, 1L]; b <- E[k, 2L]
      cnt[a] <- cnt[a] + 1L; nbrs[a, cnt[a]] <- b
      cnt[b] <- cnt[b] + 1L; nbrs[b, cnt[b]] <- a
    }

    mc <- buildModelCode(spec, engine, n = n, L = K, d = d)
    dataList <- buildDataList(spec, data)
    baseConst <- buildConstants(spec, prior, n)
    constants <- c(baseConst,
                   list(K = K, nE = length(e1), e1 = e1, e2 = e2,
                        betaMax = engine@betaMax))

    initRatio <- .resolveInitRatio(mcmcControl)
    initsFn <- function(s) {
      ci <- .withSeed(s, function() componentInits(spec, prior, data, K,
                      initMethod = initMethod, initRatio = initRatio))
      c(list(z = ci$alloc, beta = engine@beta), ci$params)
    }

    # Replace whatever default sampler configureMCMC assigned to the joint
    # Potts node with the single-site Gibbs sweep.
    estBeta <- engine@estimateBeta
    betaMaxV <- engine@betaMax
    configureHook <- function(conf, rmodel) {
      conf$removeSamplers("z")
      # NIMBLE's default sampler on beta would target the UNNORMALISED Potts
      # density, which is wrong in beta; it is removed unconditionally. When
      # beta is fixed, no sampler replaces it (beta stays at its init).
      conf$removeSamplers("beta")
      conf$addSampler(target = paste0("z[1:", n, "]"),
                      type = samplerFn,
                      control = list(K = K, nbrs = nbrs,
                                     nDeg = deg, n = n, d = d,
                                     p = if (!is.null(prior$p)) prior$p else 1L,
                                     X = if (!is.null(prior$X)) prior$X else matrix(0, 1, 1),
                                     df = if (!is.null(prior$df)) prior$df else 4,
                                     size = if (!is.null(prior$size)) prior$size else 1))
      if (estBeta)
        conf$addSampler(target = "beta", type = .pottsBetaPLSampler,
                        control = list(K = K, nbrs = nbrs, nDeg = deg, n = n,
                                       betaMax = betaMaxV, propSd = 0.15))
    }

    paramDim <- if (!is.null(prior$p)) prior$p else d
    out <- .runNimbleMixture(spec, mc, constants, dataList, initsFn,
                      n = n, count = K, paramDim = paramDim, prior = prior,
                      mcmcControl = mcmcControl, seed = seed, verbose = verbose,
                      configureHook = configureHook,
                      cacheExtra = list(engine = "MRFEngine", estimateBeta = estBeta))
    if (estBeta && "beta" %in% colnames(out$mcmcSamples)) {
      bDraw <- as.numeric(out$mcmcSamples[, "beta"])
      if (mean(bDraw > 0.98 * betaMaxV) > 0.1)
        warning("The MRF interaction beta piles up near its upper bound ",
                "betaMax = ", betaMaxV, ": the field may be (near) saturated. ",
                "Consider a larger prior$betaMax, or interpret the smoothing ",
                "as effectively maximal.", call. = FALSE)
    }
    out
  }
)

# --- multivariate Gaussian sweep (v0.7.0) -------------------------------------
# Same single-site Gibbs full conditional, with multivariate-normal emissions
# under the Normal-Inverse-Wishart kernel. Component Cholesky factors are
# hoisted once per sweep, then each site costs one dmnorm_chol evaluation per
# component.

#' @keywords internal
.pottsGibbsMvNormalSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K    <- control$K
    nbrs <- control$nbrs
    nDeg <- control$nDeg
    n    <- control$n
    d    <- control$d
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z    <- model[["z"]]
    beta <- model[["beta"]]
    mu   <- model[["muTilde"]]        # K x d
    covT <- model[["covTilde"]]       # K x d x d
    y    <- model[["y"]]              # n x d
    ch <- nimArray(0, dim = c(K, d, d))
    for (k in 1:K) ch[k, 1:d, 1:d] <- chol(covT[k, 1:d, 1:d])
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) {
          for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        }
        lp[k] <- dmnorm_chol(y[i, 1:d], mu[k, 1:d], ch[k, 1:d, 1:d],
                             prec_param = 0, log = TRUE) + beta * same
      }
      mx <- max(lp)
      p <- exp(lp - mx)
      p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {})
)


# --- polymorphic sampler selection (no class-name if/else) --------------------
# S4 dispatch chooses the label sampler for the component family. Heavy-tail
# subclasses INHERIT from the Gaussian specs, so without their own methods they
# would silently dispatch to the Gaussian sweep with the wrong emission
# density; the explicit blocking methods below turn that into a clear message.

#' @keywords internal
setGeneric(".mrfSamplerFor", function(spec) standardGeneric(".mrfSamplerFor"))

setMethod(".mrfSamplerFor", "NormalUvSpec",
          function(spec) .pottsGibbsNormalSampler)
setMethod(".mrfSamplerFor", "NormalMvSpec",
          function(spec) .pottsGibbsMvNormalSampler)
.mrfUnsupported <- function(spec)
  stop("method = 'mrf' currently supports Gaussian components (univariate or ",
       "multivariate); '", spec@name, "' components are not yet available ",
       "under the spatial engine.", call. = FALSE)
setMethod(".mrfSamplerFor", "ANY", function(spec) .mrfUnsupported(spec))

# --- mixture-of-regressions sweep (v0.7.0) ------------------------------------
# Emission is the Gaussian regression density at the component's coefficient
# vector: f(y_i | k) = N(y_i; x_i' beta_k, s2_k). The design matrix X is a
# model CONSTANT (not a variable), so it is supplied through the sampler
# control rather than read from the model.

#' @keywords internal
.pottsGibbsNormalRegSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K    <- control$K
    nbrs <- control$nbrs
    nDeg <- control$nDeg
    n    <- control$n
    p    <- control$p
    X    <- control$X                 # n x p design matrix (constant)
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z    <- model[["z"]]
    beta <- model[["beta"]]
    bT   <- model[["betaTilde"]]      # K x p
    s2   <- model[["s2Tilde"]]
    y    <- model[["y"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) {
          for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        }
        mu_ik <- inprod(X[i, 1:p], bT[k, 1:p])
        lp[k] <- dnorm(y[i], mu_ik, sqrt(s2[k]), log = TRUE) + beta * same
      }
      mx <- max(lp)
      pr <- exp(lp - mx)
      pr <- pr / sum(pr)
      z[i] <- rcat(1, pr)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {})
)


setMethod(".mrfSamplerFor", "NormalRegSpec",
          function(spec) .pottsGibbsNormalRegSampler)

# --- pseudo-likelihood Metropolis update for beta (v0.7.0) --------------------
# The Potts partition function C(beta) is intractable, so the exact posterior
# of beta is doubly intractable. Following the classical route for hidden Potts
# fields, beta is updated by random-walk Metropolis against the Besag (1975)
# PSEUDO-likelihood
#   PL(beta; z) = prod_i p(z_i | z_{N(i)}, beta)
#               = prod_i exp(beta * m_i(z_i)) / sum_k exp(beta * m_i(k)),
# with m_i(k) the number of neighbours of i currently labelled k. This is an
# APPROXIMATION (slightly biased for strong interactions), documented as such;
# an exchange-algorithm refinement (Murray, Ghahramani & MacKay 2006) is a
# possible future upgrade. The uniform prior beta ~ U(0, betaMax) cancels in
# the acceptance ratio; proposals outside the support are rejected.
#
# CORRECTNESS NOTE: this must be the ONLY sampler ever assigned to the beta
# node. NIMBLE's default sampler would target the model's UNNORMALISED Potts
# log-density, which is wrong in beta (it omits -log C(beta)); the engine
# removes it unconditionally.

#' @keywords internal
.pottsBetaPLSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K       <- control$K
    nbrs    <- control$nbrs
    nDeg    <- control$nDeg
    n       <- control$n
    betaMax <- control$betaMax
    propSd  <- control$propSd
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z  <- model[["z"]]
    b0 <- model[["beta"]]
    b1 <- rnorm(1, b0, propSd)
    if (b1 < 0 | b1 > betaMax) {
      ## outside the uniform support: reject outright
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
      return()
    }
    ## neighbour same-label counts m_i(k), shared by both log-PL evaluations
    nm <- nimMatrix(0, nrow = n, ncol = K)
    for (i in 1:n) {
      if (nDeg[i] > 0) {
        for (m in 1:nDeg[i]) {
          lab <- z[nbrs[i, m]]
          nm[i, lab] <- nm[i, lab] + 1
        }
      }
    }
    lp0 <- 0; lp1 <- 0
    for (i in 1:n) {
      zi <- z[i]
      s0 <- 0; s1 <- 0
      for (k in 1:K) {
        s0 <- s0 + exp(b0 * nm[i, k])
        s1 <- s1 + exp(b1 * nm[i, k])
      }
      lp0 <- lp0 + b0 * nm[i, zi] - log(s0)
      lp1 <- lp1 + b1 * nm[i, zi] - log(s1)
    }
    logAlpha <- lp1 - lp0
    if (runif(1) < exp(min(0, logAlpha))) {
      model[["beta"]] <<- b1
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    } else {
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list(reset = function() {})
)

# --- .pottsify: derive any family's MRF kernel from its fixed-K kernel -------
# (v0.8.0, per the maintainer-approved feasibility study.) The Potts field
# only replaces the LABEL layer, so every fixed-K kernel converts mechanically:
# drop `z[i] ~ dcat(weights[...])` and `weights[...] ~ ddirch(...)`, prepend
# `beta ~ dunif(0, betaMax)` and the joint Potts node. AST transformation in
# the style of .stripMixtureLayer. Monitors lose "weights" and gain "beta".

#' @keywords internal
.pottsify <- function(mcFixedK, n) {
  isLhs <- function(stmt, var) {
    if (!is.call(stmt) || !identical(stmt[[1]], as.name("~"))) return(FALSE)
    lhs <- stmt[[2]]
    v <- if (is.call(lhs) && identical(lhs[[1]], as.name("["))) lhs[[2]] else lhs
    identical(v, as.name(var))
  }
  dropIn <- function(block) {
    kept <- list(block[[1]])
    for (k in seq.int(2L, length(block))) {
      stmt <- block[[k]]
      if (is.call(stmt) && identical(stmt[[1]], as.name("for"))) {
        stmt[[4]] <- dropIn(stmt[[4]])
        if (length(stmt[[4]]) > 1L) kept[[length(kept) + 1L]] <- stmt
      } else if (!isLhs(stmt, "z") && !isLhs(stmt, "weights")) {
        kept[[length(kept) + 1L]] <- stmt
      }
    }
    as.call(kept)
  }
  body <- dropIn(mcFixedK$code)
  header <- str2lang(paste0(
    "{ beta ~ dunif(0, betaMax)\n  z[1:", n,
    "] ~ dPottsNimix(beta, e1[1:nE], e2[1:nE]) }"))
  stmts <- c(as.list(header)[-1L], as.list(body)[-1L])
  mcFixedK$code <- as.call(c(list(as.name("{")), stmts))
  mcFixedK$monitors <- unique(c("beta", setdiff(mcFixedK$monitors, "weights")))
  mcFixedK$allocNode <- "z"
  mcFixedK
}

#' @describeIn buildModelCode Default MRF kernel for any family: the family's
#'   fixed-K kernel with the Dirichlet-categorical label layer replaced by the
#'   joint Potts node (derived mechanically; see the spatial design notes).
#' @export
setMethod("buildModelCode", signature("DistributionSpec", "MRFEngine"),
  function(spec, engine, n, L, ...) {
    fk <- buildModelCode(spec, FixedKEngine(dirichletConc = 1), n = n, L = L, ...)
    .pottsify(fk, n = n)
  }
)

# --- batch-1 sweeps (closed-form emissions; DSL calls verified empirically) ---

#' @keywords internal
.pottsGibbsPoissonSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; lam <- model[["lambda"]]; y <- model[["y"]]
    beta <- model[["beta"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        lp[k] <- dpois(y[i], lam[k], log = TRUE) + beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

#' @keywords internal
.pottsGibbsBinomialSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    size <- control$size
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; pr <- model[["prob"]]; y <- model[["y"]]
    beta <- model[["beta"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        lp[k] <- dbinom(y[i], size = size, prob = pr[k], log = TRUE) + beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

#' @keywords internal
.pottsGibbsStudentTSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    df <- control$df
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; mu <- model[["muTilde"]]; tau <- model[["tauTilde"]]
    y <- model[["y"]]; beta <- model[["beta"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        ## kernel is dt(mu, tau, df): sigma = 1/sqrt(tau)
        lp[k] <- dt_nonstandard(y[i], df, mu[k], 1 / sqrt(tau[k]), log = TRUE) +
                 beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

#' @keywords internal
.pottsGibbsPoissonRegSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    p <- control$p; X <- control$X
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; bT <- model[["betaTilde"]]; y <- model[["y"]]
    beta <- model[["beta"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        lp[k] <- dpois(y[i], exp(inprod(X[i, 1:p], bT[k, 1:p])), log = TRUE) +
                 beta * same
      }
      mx <- max(lp); pr <- exp(lp - mx); pr <- pr / sum(pr)
      z[i] <- rcat(1, pr)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

#' @keywords internal
.pottsGibbsBinomialRegSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    p <- control$p; X <- control$X; size <- control$size
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; bT <- model[["betaTilde"]]; y <- model[["y"]]
    beta <- model[["beta"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        lp[k] <- dbinom(y[i], size = size,
                        prob = ilogit(inprod(X[i, 1:p], bT[k, 1:p])),
                        log = TRUE) + beta * same
      }
      mx <- max(lp); pr <- exp(lp - mx); pr <- pr / sum(pr)
      z[i] <- rcat(1, pr)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

#' @keywords internal
.pottsGibbsStudentTRegSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    p <- control$p; X <- control$X; df <- control$df
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; bT <- model[["betaTilde"]]; s2 <- model[["s2Tilde"]]
    y <- model[["y"]]; beta <- model[["beta"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        ## kernel is dt(mu, tau = 1/s2, df): sigma = sqrt(s2)
        lp[k] <- dt_nonstandard(y[i], df, inprod(X[i, 1:p], bT[k, 1:p]),
                                sqrt(s2[k]), log = TRUE) + beta * same
      }
      mx <- max(lp); pr <- exp(lp - mx); pr <- pr / sum(pr)
      z[i] <- rcat(1, pr)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

# --- batch-1 dispatch ---------------------------------------------------------
setMethod(".mrfSamplerFor", "PoissonSpec",
          function(spec) .pottsGibbsPoissonSampler)
setMethod(".mrfSamplerFor", "BinomialSpec",
          function(spec) .pottsGibbsBinomialSampler)
setMethod(".mrfSamplerFor", "StudentTUvSpec",
          function(spec) .pottsGibbsStudentTSampler)
setMethod(".mrfSamplerFor", "PoissonRegSpec",
          function(spec) .pottsGibbsPoissonRegSampler)
setMethod(".mrfSamplerFor", "BinomialRegSpec",
          function(spec) .pottsGibbsBinomialRegSampler)
setMethod(".mrfSamplerFor", "StudentTRegSpec",
          function(spec) .pottsGibbsStudentTRegSampler)

# --- batch-2 sweeps (v0.9.0): multivariate and augmented families -------------
# Augmented (Normal-Gamma) sweeps CONDITION ON the current latent omega_i --
# the correct full conditional in the augmented joint; the documented
# augmentation mixing caveat applies and is benchmarked against the direct-t
# routes (see NEWS / knowledge patch).

#' @keywords internal
.pottsGibbsMvTSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    d <- control$d; df <- control$df
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; mu <- model[["muTilde"]]; covT <- model[["covTilde"]]
    y <- model[["y"]]; beta <- model[["beta"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        lp[k] <- dmvt_nimix(y[i, 1:d], mu[k, 1:d], covT[k, 1:d, 1:d], df,
                            log = 1) + beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

#' @keywords internal
.pottsGibbsNGSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; mu <- model[["muTilde"]]; s2 <- model[["s2Tilde"]]
    om <- model[["omega"]]; y <- model[["y"]]; beta <- model[["beta"]]
    for (i in 1:n) {
      lp <- numeric(K)
      sdi <- 1 / sqrt(om[i])
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        lp[k] <- dnorm(y[i], mu[k], sqrt(s2[k]) * sdi, log = TRUE) + beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

#' @keywords internal
.pottsGibbsNGMvSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    d <- control$d
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; mu <- model[["muTilde"]]; covT <- model[["covTilde"]]
    om <- model[["omega"]]; y <- model[["y"]]; beta <- model[["beta"]]
    ch <- nimArray(0, dim = c(K, d, d))
    for (k in 1:K) ch[k, 1:d, 1:d] <- chol(covT[k, 1:d, 1:d])
    for (i in 1:n) {
      sci <- 1 / sqrt(om[i])                 # chol(Sigma/om) = chol(Sigma)/sqrt(om)
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        lp[k] <- dmnorm_chol(y[i, 1:d], mu[k, 1:d], ch[k, 1:d, 1:d] * sci,
                             prec_param = 0, log = TRUE) + beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

#' @keywords internal
.pottsGibbsNGRegSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    p <- control$p; X <- control$X
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; bT <- model[["betaTilde"]]; s2 <- model[["s2Tilde"]]
    om <- model[["omega"]]; y <- model[["y"]]; beta <- model[["beta"]]
    for (i in 1:n) {
      lp <- numeric(K)
      sdi <- 1 / sqrt(om[i])
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        lp[k] <- dnorm(y[i], inprod(X[i, 1:p], bT[k, 1:p]),
                       sqrt(s2[k]) * sdi, log = TRUE) + beta * same
      }
      mx <- max(lp); pr <- exp(lp - mx); pr <- pr / sum(pr)
      z[i] <- rcat(1, pr)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

#' @keywords internal
.pottsGibbsMvRegSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    p <- control$p; d <- control$d; X <- control$X
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; bT <- model[["betaTilde"]]; covT <- model[["covTilde"]]
    y <- model[["y"]]; beta <- model[["beta"]]
    ch <- nimArray(0, dim = c(K, d, d))
    for (k in 1:K) ch[k, 1:d, 1:d] <- chol(covT[k, 1:d, 1:d])
    muv <- numeric(d)
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        for (dd in 1:d) muv[dd] <- inprod(X[i, 1:p], bT[k, 1:p, dd])
        lp[k] <- dmnorm_chol(y[i, 1:d], muv[1:d], ch[k, 1:d, 1:d],
                             prec_param = 0, log = TRUE) + beta * same
      }
      mx <- max(lp); pr <- exp(lp - mx); pr <- pr / sum(pr)
      z[i] <- rcat(1, pr)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

#' @keywords internal
.pottsGibbsMvTRegSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    p <- control$p; d <- control$d; X <- control$X; df <- control$df
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; bT <- model[["betaTilde"]]; covT <- model[["covTilde"]]
    y <- model[["y"]]; beta <- model[["beta"]]
    muv <- numeric(d)
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        for (dd in 1:d) muv[dd] <- inprod(X[i, 1:p], bT[k, 1:p, dd])
        lp[k] <- dmvt_nimix(y[i, 1:d], muv[1:d], covT[k, 1:d, 1:d], df,
                            log = 1) + beta * same
      }
      mx <- max(lp); pr <- exp(lp - mx); pr <- pr / sum(pr)
      z[i] <- rcat(1, pr)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

#' @keywords internal
.pottsGibbsNGMvRegSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    p <- control$p; d <- control$d; X <- control$X
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; bT <- model[["betaTilde"]]; covT <- model[["covTilde"]]
    om <- model[["omega"]]; y <- model[["y"]]; beta <- model[["beta"]]
    ch <- nimArray(0, dim = c(K, d, d))
    for (k in 1:K) ch[k, 1:d, 1:d] <- chol(covT[k, 1:d, 1:d])
    muv <- numeric(d)
    for (i in 1:n) {
      sci <- 1 / sqrt(om[i])
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        for (dd in 1:d) muv[dd] <- inprod(X[i, 1:p], bT[k, 1:p, dd])
        lp[k] <- dmnorm_chol(y[i, 1:d], muv[1:d], ch[k, 1:d, 1:d] * sci,
                             prec_param = 0, log = TRUE) + beta * same
      }
      mx <- max(lp); pr <- exp(lp - mx); pr <- pr / sum(pr)
      z[i] <- rcat(1, pr)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

setMethod(".mrfSamplerFor", "StudentTMvSpec",
          function(spec) .pottsGibbsMvTSampler)
setMethod(".mrfSamplerFor", "NormalGammaUvSpec",
          function(spec) .pottsGibbsNGSampler)
setMethod(".mrfSamplerFor", "NormalGammaMvSpec",
          function(spec) .pottsGibbsNGMvSampler)
setMethod(".mrfSamplerFor", "NormalGammaRegSpec",
          function(spec) .pottsGibbsNGRegSampler)
setMethod(".mrfSamplerFor", "NormalMvRegSpec",
          function(spec) .pottsGibbsMvRegSampler)
setMethod(".mrfSamplerFor", "StudentTMvRegSpec",
          function(spec) .pottsGibbsMvTRegSampler)
setMethod(".mrfSamplerFor", "NormalGammaMvRegSpec",
          function(spec) .pottsGibbsNGMvRegSampler)

# --- MSNBurr sweeps (neo-normal; densities are the registered stable forms) ----

#' @keywords internal
.pottsGibbsMSNBurrSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; mu <- model[["muTilde"]]; sg <- model[["sigmaTilde"]]
    al <- model[["alphaTilde"]]; y <- model[["y"]]; beta <- model[["beta"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        lp[k] <- dMSNBurr_k(y[i], mu[k], sg[k], al[k], log = 1) + beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

#' @keywords internal
.pottsGibbsMSNBurr2aSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; mu <- model[["muTilde"]]; sg <- model[["sigmaTilde"]]
    al <- model[["alphaTilde"]]; y <- model[["y"]]; beta <- model[["beta"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        lp[k] <- dMSNBurr2a_k(y[i], mu[k], sg[k], al[k], log = 1) + beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

setMethod(".mrfSamplerFor", "MSNBurrUvSpec",
          function(spec) .pottsGibbsMSNBurrSampler)
setMethod(".mrfSamplerFor", "MSNBurr2aUvSpec",
          function(spec) .pottsGibbsMSNBurr2aSampler)

# --- GMSNBurr MRF sweep (neo-normal, two shape params) -------------------------

#' @keywords internal
.pottsGibbsGMSNBurrSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; mu <- model[["muTilde"]]; sg <- model[["sigmaTilde"]]
    al <- model[["alphaTilde"]]; th <- model[["thetaTilde"]]
    y <- model[["y"]]; beta <- model[["beta"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        lp[k] <- dGMSNBurr_k(y[i], mu[k], sg[k], al[k], th[k], log = 1) +
          beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

setMethod(".mrfSamplerFor", "GMSNBurrUvSpec",
          function(spec) .pottsGibbsGMSNBurrSampler)

# --- SEP / LEP MRF sweeps (symmetric exponential power) ------------------------

#' @keywords internal
.pottsGibbsSEPSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; mu <- model[["muTilde"]]; sg <- model[["sigmaTilde"]]
    nu <- model[["nuTilde"]]; y <- model[["y"]]; beta <- model[["beta"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        lp[k] <- dSEP_k(y[i], mu[k], sg[k], nu[k], log = 1) + beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

#' @keywords internal
.pottsGibbsLEPSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; mu <- model[["muTilde"]]; sg <- model[["sigmaTilde"]]
    nu <- model[["nuTilde"]]; y <- model[["y"]]; beta <- model[["beta"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        lp[k] <- dLEP_k(y[i], mu[k], sg[k], nu[k], log = 1) + beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

setMethod(".mrfSamplerFor", "SEPUvSpec", function(spec) .pottsGibbsSEPSampler)
setMethod(".mrfSamplerFor", "LEPUvSpec", function(spec) .pottsGibbsLEPSampler)

#' @keywords internal
.pottsGibbsFSSNSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; mu <- model[["muTilde"]]; sg <- model[["sigmaTilde"]]
    al <- model[["alphaTilde"]]; y <- model[["y"]]; beta <- model[["beta"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        lp[k] <- dFSSN_k(y[i], mu[k], sg[k], al[k], log = 1) + beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

setMethod(".mrfSamplerFor", "FSSNUvSpec", function(spec) .pottsGibbsFSSNSampler)

# --- FOSSEP / FSST / JFST MRF sweeps (4-parameter neo-normal emissions) --------

#' @keywords internal
.pottsGibbsFOSSEPSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; mu <- model[["muTilde"]]; sg <- model[["sigmaTilde"]]
    al <- model[["alphaTilde"]]; th <- model[["thetaTilde"]]
    y <- model[["y"]]; beta <- model[["beta"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        lp[k] <- dFOSSEP_k(y[i], mu[k], sg[k], al[k], th[k], log = 1) + beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

#' @keywords internal
.pottsGibbsFSSTSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; mu <- model[["muTilde"]]; sg <- model[["sigmaTilde"]]
    al <- model[["alphaTilde"]]; nu <- model[["nuTilde"]]
    y <- model[["y"]]; beta <- model[["beta"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        lp[k] <- dFSST_k(y[i], mu[k], sg[k], al[k], nu[k], log = 1) + beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

#' @keywords internal
.pottsGibbsJFSTSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg; n <- control$n
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model[["z"]]; mu <- model[["muTilde"]]; sg <- model[["sigmaTilde"]]
    al <- model[["alphaTilde"]]; th <- model[["thetaTilde"]]
    y <- model[["y"]]; beta <- model[["beta"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        lp[k] <- dJFST_k(y[i], mu[k], sg[k], al[k], th[k], log = 1) + beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

setMethod(".mrfSamplerFor", "FOSSEPUvSpec", function(spec) .pottsGibbsFOSSEPSampler)
setMethod(".mrfSamplerFor", "FSSTUvSpec", function(spec) .pottsGibbsFSSTSampler)
setMethod(".mrfSamplerFor", "JFSTUvSpec", function(spec) .pottsGibbsJFSTSampler)

# --- Skew multivariate MRF sweeps (Ferreira-Steel families) --------------------
# Same Potts x emission structure as .pottsGibbsMvNormalSampler, but the
# emission is the compiled FS skew-mv kernel, which takes Sigma directly and
# forms its own Cholesky factor internally.

#' @keywords internal
.pottsGibbsSkewMvNSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg
    n <- control$n; d <- control$d
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z    <- model[["z"]]
    beta <- model[["beta"]]
    mu   <- model[["muTilde"]]      # K x d
    Sig  <- model[["SigTilde"]]     # K x d x d
    gam  <- model[["gamTilde"]]     # K x d
    y    <- model[["y"]]            # n x d
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) {
          for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        }
        lp[k] <- dSkewMvN_k(y[i, 1:d], mu[k, 1:d], Sig[k, 1:d, 1:d],
                            gam[k, 1:d], log = 1) + beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

#' @keywords internal
.pottsGibbsSkewMvITSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg
    n <- control$n; d <- control$d
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z    <- model[["z"]]
    beta <- model[["beta"]]
    mu   <- model[["muTilde"]]
    Sig  <- model[["SigTilde"]]
    gam  <- model[["gamTilde"]]
    nu   <- model[["nuTilde"]]      # K x d
    y    <- model[["y"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) {
          for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        }
        lp[k] <- dSkewMvIT_k(y[i, 1:d], mu[k, 1:d], Sig[k, 1:d, 1:d],
                             gam[k, 1:d], nu[k, 1:d], log = 1) + beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

setMethod(".mrfSamplerFor", "SkewNormalMvSpec",
          function(spec) .pottsGibbsSkewMvNSampler)
setMethod(".mrfSamplerFor", "SkewIStudentMvSpec",
          function(spec) .pottsGibbsSkewMvITSampler)

# --- Skew multivariate MRF sweeps, estimated-O variants -------------------------
# Identical Potts x emission structure; the emission kernels additionally take
# the per-component Householder angle theta.

#' @keywords internal
.pottsGibbsSkewMvNOSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg
    n <- control$n; d <- control$d
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z    <- model[["z"]]
    beta <- model[["beta"]]
    mu   <- model[["muTilde"]]
    Sig  <- model[["SigTilde"]]
    gam  <- model[["gamTilde"]]
    th   <- model[["thetaTilde"]]     # K
    y    <- model[["y"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) {
          for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        }
        lp[k] <- dSkewMvNO_k(y[i, 1:d], mu[k, 1:d], Sig[k, 1:d, 1:d],
                             gam[k, 1:d], th[k], log = 1) + beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

#' @keywords internal
.pottsGibbsSkewMvITOSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg
    n <- control$n; d <- control$d
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z    <- model[["z"]]
    beta <- model[["beta"]]
    mu   <- model[["muTilde"]]
    Sig  <- model[["SigTilde"]]
    gam  <- model[["gamTilde"]]
    nu   <- model[["nuTilde"]]
    th   <- model[["thetaTilde"]]
    y    <- model[["y"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) {
          for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        }
        lp[k] <- dSkewMvITO_k(y[i, 1:d], mu[k, 1:d], Sig[k, 1:d, 1:d],
                              gam[k, 1:d], nu[k, 1:d], th[k], log = 1) +
          beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

setMethod(".mrfSamplerFor", "SkewNormalMvOSpec",
          function(spec) .pottsGibbsSkewMvNOSampler)
setMethod(".mrfSamplerFor", "SkewIStudentMvOSpec",
          function(spec) .pottsGibbsSkewMvITOSampler)

# --- Skew multivariate MRF sweeps, general-m estimated-O variants ---------------
# Same Potts x emission structure. The number of Householder angles is a
# function of the dimension, nAng = d(d-1)/2, so it is derived in setup rather
# than passed through control.

#' @keywords internal
.pottsGibbsSkewMvNOGSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg
    n <- control$n; d <- control$d
    nAng <- as.integer(d * (d - 1) / 2)
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z    <- model[["z"]]
    beta <- model[["beta"]]
    mu   <- model[["muTilde"]]
    Sig  <- model[["SigTilde"]]
    gam  <- model[["gamTilde"]]
    th   <- model[["thetaTilde"]]      # K x nAng
    y    <- model[["y"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) {
          for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        }
        lp[k] <- dSkewMvNOG_k(y[i, 1:d], mu[k, 1:d], Sig[k, 1:d, 1:d],
                              gam[k, 1:d], th[k, 1:nAng], log = 1) + beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

#' @keywords internal
.pottsGibbsSkewMvITOGSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K <- control$K; nbrs <- control$nbrs; nDeg <- control$nDeg
    n <- control$n; d <- control$d
    nAng <- as.integer(d * (d - 1) / 2)
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z    <- model[["z"]]
    beta <- model[["beta"]]
    mu   <- model[["muTilde"]]
    Sig  <- model[["SigTilde"]]
    gam  <- model[["gamTilde"]]
    nu   <- model[["nuTilde"]]
    th   <- model[["thetaTilde"]]
    y    <- model[["y"]]
    for (i in 1:n) {
      lp <- numeric(K)
      for (k in 1:K) {
        same <- 0
        if (nDeg[i] > 0) {
          for (m in 1:nDeg[i]) if (z[nbrs[i, m]] == k) same <- same + 1
        }
        lp[k] <- dSkewMvITOG_k(y[i, 1:d], mu[k, 1:d], Sig[k, 1:d, 1:d],
                               gam[k, 1:d], nu[k, 1:d], th[k, 1:nAng],
                               log = 1) + beta * same
      }
      mx <- max(lp); p <- exp(lp - mx); p <- p / sum(p)
      z[i] <- rcat(1, p)
    }
    model[["z"]] <<- z
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {}))

setMethod(".mrfSamplerFor", "SkewNormalMvOGenSpec",
          function(spec) .pottsGibbsSkewMvNOGSampler)
setMethod(".mrfSamplerFor", "SkewIStudentMvOGenSpec",
          function(spec) .pottsGibbsSkewMvITOGSampler)
