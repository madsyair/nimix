#' @include class-DistributionSpec.R
NULL

## ---------------------------------------------------------------------------
## dist-normal-reg.R
##
## Normal-linear regression component (the v0.3.0 distribution), used by
## nimixReg() for a DPM mixture of linear regressions with CONSTANT gating
## (mixing weights come from the CRP, they do NOT depend on covariates). This
## is the baseline of; covariate-dependent
## (concomitant) gating is an explicit future opt-in, not the default, because
## it introduces extra identifiability risk.
##
## Conjugate cluster base measure is Normal-Inverse-Gamma on (beta, sigma^2):
##     s2Tilde_j           ~ InvGamma(nu0, s0)
##     betaTilde_j | s2    ~ N_p(b0, s2Tilde_j * B0)
## which is the conjugate prior for a Gaussian linear model. NIMBLE recognises
## this as conjugate for the dCRP sampler (verified: CRP_cluster_wrapper is
## assigned to both betaTilde and s2Tilde) -- native NIMBLE machinery (project
## knowledge 0.5); nimix supplies the S4 wiring, the data-scaled g-prior, and
## the relabelling.
##
## As in NormalMvSpec, a multivariate distribution parameter may
## not be an expression, so the s2-scaled coefficient covariance is bound to a
## deterministic node covBeta, and the dynamically indexed cluster coefficients
## are resolved through a per-observation deterministic node betaObs.
##
## Default prior is DATA-SCALED: a Zellner unit-information
## g-prior on beta (B0 = g * solve(crossprod(X)), g = n) and an InvGamma on the residual
## variance whose mean equals the global OLS residual variance, with nu0 >= 3 so
## the prior variance is finite and the residual variance cannot collapse to
## zero on small components.
##
## Reference: Hurn, Justel & Robert (2003).
## ---------------------------------------------------------------------------

#' Normal-linear regression component specification
#'
#' @slot name Fixed to \code{"normal-reg"}.
#' @slot paramNames \code{c("beta", "s2")}.
#' @slot dataDim \code{1L} (univariate response).
#'
#' @references
#' Hurn, M., Justel, A., & Robert, C.P. (2003). Estimating mixtures of
#' regressions. \emph{Journal of Computational and Graphical Statistics},
#' 12(1), 55--79. \doi{10.1198/1061860031329}
#'
#' @seealso \code{\link{nimixReg}}
#' @export
setClass(
  "NormalRegSpec",
  contains = "DistributionSpec",
  prototype = prototype(
    name = "normal-reg",
    paramNames = c("beta", "s2"),
    dataDim = 1L
  )
)

#' Construct a Normal-linear regression component spec
#' @return A \code{\linkS4class{NormalRegSpec}}.
#' @examples
#' spec <- NormalRegSpec()
#' @export
NormalRegSpec <- function() new("NormalRegSpec")

# --- defaultPrior ----------------------------------------------------------

#' @describeIn defaultPrior Data-scaled Normal-Inverse-Gamma g-prior for the
#'   regression component. Requires the design matrix in \code{control$X}.
#'
#' Control overrides: \code{g} (g-prior factor; prior covariance of \code{beta}
#' is \code{s2 * g * solve(crossprod(X))}, default \code{g = n}, the unit-information
#' prior) and \code{nu0} (InvGamma shape, default 3, must exceed 2 for a finite
#' prior variance).
#' @export
setMethod("defaultPrior", "NormalRegSpec",
  function(spec, data, control = list(), ...) {
    y <- as.numeric(data)
    X <- control$X
    if (is.null(X) || !is.matrix(X))
      stop("NormalRegSpec defaultPrior needs the design matrix in control$X.",
           call. = FALSE)
    n <- length(y); p <- ncol(X)

    g    <- if (!is.null(control$g))   control$g   else n
    nu0  <- if (!is.null(control$nu0)) control$nu0 else 3
    if (nu0 <= 2)
      stop("nu0 must exceed 2 so the prior variance of s2 is finite ",
 ".", call. = FALSE)

    XtX <- crossprod(X)
    ridge <- 1e-8 * mean(diag(XtX)) + 1e-10
    XtXinv <- solve(XtX + diag(ridge, p))
    B0 <- g * XtXinv

    # Scale reference for the InvGamma prior on s2. The default is the global
    # OLS residual variance, which for a mixture measures the BETWEEN-component
    # spread as much as the within-component one, so it is deliberately
    # conservative: it biases s2 upward (measured 2-5x for well-separated
    # components at moderate n, vanishing as n_j grows -- see ?nimixReg). That
    # conservatism is load-bearing, not incidental: it regularises against
    # over-splitting when K is over-specified, so the default stays.
    #
    # A user who knows the within-component scale (pilot study, literature,
    # domain knowledge) can now say so:
    #   control$s2Guess -- your prior mean for s2. This is the form to reach
    #                      for: it asks for a quantity practitioners actually
    #                      know, not a hyper-parameter of the InvGamma.
    #   control$s0      -- the raw InvGamma scale, for callers who think in
    #                      those terms (mean = s0 / (nu0 - 1)).
    # The override is deliberately ABSOLUTE only. A relative multiplier on the
    # automatic scale was considered and rejected: that scale measures the
    # wrong quantity for a mixture, so "a tenth of it" is not a statement a
    # user can defend, only a knob to twiddle. Someone with real knowledge of
    # the within-component variance should say what it is.
    if (!is.null(control$s2Guess)) {
      sigma2hat <- control$s2Guess
      if (length(sigma2hat) != 1L || !is.finite(sigma2hat) || sigma2hat <= 0)
        stop("prior$s2Guess must be a positive scalar: it is your prior mean ",
             "for the component error variance s2.", call. = FALSE)
    } else {
      bOls <- as.numeric(XtXinv %*% crossprod(X, y))
      resid <- y - as.numeric(X %*% bOls)
      sigma2hat <- sum(resid^2) / max(1, n - p)
      if (!is.finite(sigma2hat) || sigma2hat <= 0) sigma2hat <- stats::var(y)
      if (!is.finite(sigma2hat) || sigma2hat <= 0) sigma2hat <- 1
    }
    s0 <- if (!is.null(control$s0)) {
      if (!is.null(control$s2Guess))
        stop("Give either prior$s2Guess (your prior mean for s2) or ",
             "prior$s0 (the raw InvGamma scale), not both.", call. = FALSE)
      if (length(control$s0) != 1L || !is.finite(control$s0) ||
          control$s0 <= 0)
        stop("prior$s0 must be a positive scalar (the InvGamma scale; the ",
             "prior mean of s2 is s0 / (nu0 - 1)).", call. = FALSE)
      control$s0
    } else sigma2hat * (nu0 - 1)

    list(b0 = rep(0, p), B0 = B0, nu0 = nu0, s0 = s0,
         p = p, g = g, X = X)
  }
)

# --- validateParams --------------------------------------------------------

#' @describeIn validateParams Validate a Normal-Inverse-Gamma regression prior
#'   and enforce \eqn{\dim(b_0) = \dim(B_0) = p} and \eqn{nu_0 > 2}.
#' @export
setMethod("validateParams", "NormalRegSpec",
  function(spec, params, ...) {
    req <- c("b0", "B0", "nu0", "s0", "p")
    miss <- setdiff(req, names(params))
    if (length(miss))
      stop("NormalRegSpec prior is missing: ", paste(miss, collapse = ", "),
           call. = FALSE)
    p <- params$p
    if (length(params$b0) != p)
      stop("b0 must have length p = ", p, ".", call. = FALSE)
    if (!is.matrix(params$B0) || !all(dim(params$B0) == c(p, p)))
      stop("B0 must be a ", p, " x ", p, " matrix.", call. = FALSE)
    if (params$nu0 <= 2)
 stop("nu0 must be > 2 (finite prior variance).",
           call. = FALSE)
    if (params$s0 <= 0) stop("s0 (InvGamma scale) must be > 0.", call. = FALSE)
    ev <- min(eigen((params$B0 + t(params$B0)) / 2, symmetric = TRUE,
                    only.values = TRUE)$values)
    if (ev <= 0) stop("B0 must be positive definite.", call. = FALSE)
    invisible(TRUE)
  }
)

# --- simulateParams --------------------------------------------------------

#' @describeIn simulateParams Draw (beta, s2) per cluster from the NIG prior.
#'   Returns \code{beta} (nClust x p) and \code{s2} (length nClust).
setMethod("simulateParams", "NormalRegSpec",
  function(spec, prior, nClust, ...) {
    p <- prior$p
    s2 <- 1 / stats::rgamma(nClust, shape = prior$nu0, rate = prior$s0)
    beta <- matrix(NA_real_, nClust, p)
    for (j in seq_len(nClust))
      beta[j, ] <- .rmvnorm1(prior$b0, s2[j] * prior$B0)
    list(beta = beta, s2 = s2)
  }
)

# --- componentDensity ------------------------------------------------------

#' @describeIn componentDensity Gaussian density of a response given its linear
#'   predictor. \code{params} must carry \code{mu} (the fitted mean) and
#'   \code{s2}.
setMethod("componentDensity", "NormalRegSpec",
  function(spec, ...) {
    function(x, params) stats::dnorm(x, mean = params[["mu"]],
                                     sd = sqrt(params[["s2"]]))
  }
)

# --- buildModelCode: NormalRegSpec x DPMEngine -----------------------------

#' @describeIn buildModelCode DPM mixture-of-linear-regressions model code
#'   (dCRP) with a conjugate Normal-Inverse-Gamma cluster prior and constant
#'   (CRP) gating.
#'
#' Builds NIMBLE code for
#' \deqn{y_i \sim N(x_i^\top \beta_{\xi_i}, s^2_{\xi_i}), \quad
#'       \xi_{1:n} \sim CRP(\alpha, n),}
#' with \eqn{s^2_j \sim InvGamma(nu_0, s_0)} and
#' \eqn{\beta_j \mid s^2_j \sim N_p(b_0, s^2_j B_0)}. The number of predictors
#' \eqn{p} and the design matrix \eqn{X} are constants.
#'
#' @references
#' Hurn et al. (2003) \doi{10.1198/1061860031329};
#' de Valpine et al. (2017) \doi{10.1080/10618600.2016.1172487}.
#' @export
setMethod("buildModelCode", signature("NormalRegSpec", "DPMEngine"),
  function(spec, engine, n, L, ...) {
    # Conjugate NIG regression DPM. dnorm kernel with a linear predictor plus a
    # Normal-Inverse-Gamma cluster prior is recognised as conjugate for dCRP
    # (CRP_cluster_wrapper assigned to betaTilde and s2Tilde; Neal 2000, Alg 2).
    # The s2-scaled coefficient covariance and the dynamically indexed cluster
 # coefficients are resolved via deterministic nodes.
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        betaObs[i, 1:p] <- betaTilde[xi[i], 1:p]
        mu[i] <- inprod(X[i, 1:p], betaObs[i, 1:p])
        s2Obs[i] <- s2Tilde[xi[i]]
        y[i] ~ dnorm(mu[i], var = s2Obs[i])
      }
      xi[1:n] ~ dCRP(alpha, size = n)
      for (j in 1:L) {
        s2Tilde[j] ~ dinvgamma(shape = nu0, scale = s0)
        covBeta[j, 1:p, 1:p] <- s2Tilde[j] * B0[1:p, 1:p]
        betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = covBeta[j, 1:p, 1:p])
      }
      alpha ~ dgamma(shape = aAlpha, rate = bAlpha)
    })
    list(
      code       = code,
      monitors   = c("xi", "betaTilde", "s2Tilde", "alpha"),
      paramNodes = c(beta = "betaTilde", s2 = "s2Tilde"),
      allocNode  = "xi"
    )
  }
)

# --- Engine-facing methods -------------------------------------------------

#' @describeIn buildConstants Regression constants (design matrix \code{X},
#'   number of predictors \code{p}, and the NIG hyperparameters).
setMethod("buildConstants", "NormalRegSpec",
  function(spec, prior, n, ...) {
    out <- list(n = n, p = prior$p, X = prior$X,
                b0 = prior$b0, B0 = prior$B0, nu0 = prior$nu0, s0 = prior$s0)
    if (isTRUE(prior$hasRE)) {
      out$grp    <- prior$reGrp
      out$G      <- prior$reG
      out$tauMax <- prior$tauMax
      out$tauMin <- prior$tauMin
      if (isTRUE(prior$hasRESlope)) {
        out$xRE         <- prior$reSlopeX
        out$tauMaxSlope <- prior$tauMaxSlope
        out$tauMinSlope <- prior$tauMinSlope
      }
    }
    out
  }
)

#' @describeIn buildDataList Response vector for the regression mixture.
setMethod("buildDataList", "NormalRegSpec",
  function(spec, data, ...) list(y = as.numeric(data))
)

#' @describeIn componentInits k-means-on-(predictors, response) start with local OLS
#' coefficients per cluster (dispersed start).
setMethod("componentInits", "NormalRegSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    y <- as.numeric(data)
    X <- prior$X
    n <- length(y); p <- prior$p

    bOls <- as.numeric(solve(crossprod(X) + diag(1e-8, p), crossprod(X, y)))
    s2prior <- prior$s0 / (prior$nu0 - 1)

    # Dispersed k-means start, capped at initRatio * count (default 0.8, tunable via

    # the cap: for the DPM, count = L = K_max is a hard truncation, and early CRP

    # sweeps can briefly occupy more clusters than the modal K before merging

    # down. Seeding right at the ceiling left no room for that transient.

    initRatio <- .initRatioArg(...)
    k0 <- max(1L, min(as.integer(floor(initRatio * count)), as.integer(ceiling(sqrt(n)))))
    xiInit <- rep(1L, n)
    betaList <- list(bOls)
    s2List   <- list(s2prior)

    feats <- if (p > 1L) cbind(scale(X[, -1, drop = FALSE]), scale(y)) else
      cbind(scale(y))
    feats[!is.finite(feats)] <- 0
    if (!identical(initMethod, "single") && k0 >= 2L &&
        nrow(unique(feats)) >= k0) {
      km <- tryCatch(stats::kmeans(feats, centers = k0, nstart = 5L),
                     error = function(e) NULL)
      if (!is.null(km)) {
        xiInit <- as.integer(km$cluster)
        betaList <- vector("list", k0); s2List <- vector("list", k0)
        for (j in seq_len(k0)) {
          rows <- which(xiInit == j)
          if (length(rows) > p) {
            Xj <- X[rows, , drop = FALSE]; yj <- y[rows]
            bj <- tryCatch(as.numeric(solve(crossprod(Xj) + diag(1e-8, p),
                                            crossprod(Xj, yj))),
                           error = function(e) bOls)
            rj <- yj - as.numeric(Xj %*% bj)
            vj <- sum(rj^2) / max(1, length(rows) - p)
            betaList[[j]] <- bj
            s2List[[j]]   <- if (is.finite(vj) && vj > 0) vj else s2prior
          } else { betaList[[j]] <- bOls; s2List[[j]] <- s2prior }
        }
      }
    }

    betaInit <- matrix(rep(bOls, each = count), nrow = count)  # count x p
    s2Init   <- rep(s2prior, count)
    occ <- sort(unique(xiInit))
    for (idx in seq_along(occ)) {
      j <- occ[idx]
      if (idx <= length(betaList)) betaInit[j, ] <- betaList[[idx]]
      if (idx <= length(s2List))   s2Init[j]     <- s2List[[idx]]
    }
    list(alloc = xiInit, params = list(betaTilde = betaInit, s2Tilde = s2Init))
  }
)

#' @describeIn extractParamTraces Parse betaTilde (L x p) and s2Tilde (L)
#'   traces; \code{prior$coefNames} (if present) labels the coefficients.
setMethod("extractParamTraces", "NormalRegSpec",
  function(spec, samples, L, d = NULL, prior = NULL, ...) {
    p <- if (!is.null(prior$p)) prior$p else d
    if (is.null(p)) stop("extractParamTraces(NormalRegSpec) needs p.",
                         call. = FALSE)
    list(
      beta      = .nodeToArray(samples, "betaTilde", c(L, p)),
      s2        = .nodeToArray(samples, "s2Tilde",   L),
      p         = p,
      coefNames = if (!is.null(prior$coefNames)) prior$coefNames else
        paste0("b", seq_len(p))
    )
  }
)

#' @describeIn relabelComponents Permute (beta, s2) and summarise the regression
#'   components (one coefficient column per predictor).
setMethod("relabelComponents", "NormalRegSpec",
  function(spec, paramTrace, idx, occList, perms, modalK, weights, ...) {
    p        <- paramTrace$p
    betaTr   <- paramTrace$beta                 # iter x L x p
    s2Tr     <- paramTrace$s2                   # iter x L
    coefNm   <- paramTrace$coefNames
    m <- length(idx)

    betaRe <- array(NA_real_, dim = c(m, modalK, p))
    s2Re   <- matrix(NA_real_, m, modalK)
    for (t in seq_len(m)) {
      r <- idx[t]; occ <- occList[[t]]; ord <- perms[t, ]
      for (k in seq_len(modalK)) {
        j <- occ[ord[k]]
        betaRe[t, k, ] <- betaTr[r, j, ]
        s2Re[t, k]     <- s2Tr[r, j]
      }
    }

    betaMean <- apply(betaRe, c(2L, 3L), mean)             # modalK x p
    q <- function(M, pr) apply(M, 2L, stats::quantile, probs = pr, names = FALSE)
    summ <- data.frame(component = seq_len(modalK),
                       weight = colMeans(weights))
    for (j in seq_len(p)) summ[[coefNm[j]]] <- betaMean[, j]
    summ$s2_mean <- colMeans(s2Re)
    summ$s2_lwr  <- q(s2Re, 0.025)
    summ$s2_upr  <- q(s2Re, 0.975)

    list(beta = betaRe, s2 = s2Re, beta_mean = betaMean,
         coefNames = coefNm, summary = summ)
  }
)

#' @describeIn isRegressionSpec Normal-linear regression is a regression spec.
#' @export
setMethod("isRegressionSpec", "NormalRegSpec", function(spec, ...) TRUE)

# --- Conjugate Gibbs for the FixedK regression path ------------------------------
# Under the FixedK engine, NIMBLE's checkConjugacy() detects only the mixture
# weights: the dynamic indexing betaObs[i, ] <- betaTilde[z[i], ] hides the
# linear-Gaussian structure from the conjugacy checker (the same class of
# platform constraint as the multivariate Cholesky-lifting issue), so beta and
# s2 fall back to RW_block / RW in raw units. That is not merely slow: it makes
# the fit visibly non-equivariant to predictor rescaling, because a single
# adaptive block proposal cannot serve an intercept of O(1) and a slope of
# O(1e-3) at once. (The DPM path is unaffected: dCRP assigns its own conjugate
# CRP_cluster_wrapper to betaTilde and s2Tilde.)
#
# This sampler restores the exact closed-form conditionals of the
# Normal-Inverse-Gamma model (Bernardo & Smith 1994, Sec. 5.2; O'Hagan &
# Forster 2004, Ch. 11). Conditional on allocations z, for cluster j with
# members M_j = {i : z_i = j}:
#   Vn = (B0^{-1} + X_j' X_j)^{-1},  bn = Vn (B0^{-1} b0 + X_j' y_j)
#   beta_j | s2_j, z, y ~ N(bn, s2_j * Vn)
#   s2_j | z, y ~ InvGamma(nu0 + n_j / 2,
#                          s0 + (y_j' y_j + b0' B0^{-1} b0 - bn' Vn^{-1} bn) / 2)
# Empty clusters draw from the prior. Both updates are joint (beta, s2), so the
# sampler replaces BOTH default samplers for cluster j.

#' @keywords internal
.betaS2ConjSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    j     <- control$j                 # cluster index this instance updates
    X     <- control$X                 # n x p design (constant)
    yv    <- control$y                 # response (constant)
    b0    <- control$b0
    B0inv <- control$B0inv             # p x p prior precision of beta (unscaled)
    nu0   <- control$nu0
    s0    <- control$s0
    p     <- ncol(X)
    n     <- nrow(X)
    betaNode <- paste0("betaTilde[", j, ", 1:", p, "]")
    s2Node   <- paste0("s2Tilde[", j, "]")
    calcNodes <- model$getDependencies(c(betaNode, s2Node))
    B0invb0 <- (B0inv %*% asCol(b0))[, 1]
    b0Qb0 <- inprod(b0, B0invb0)
  },
  run = function() {
    z <- model[["z"]]
    # accumulate X_j'X_j, X_j'y_j, y_j'y_j, n_j over members of cluster j
    XtX <- matrix(0, p, p)
    Xty <- numeric(p)
    yty <- 0
    nj  <- 0
    for (i in 1:n) {
      if (z[i] == j) {
        nj <- nj + 1
        yty <- yty + yv[i] * yv[i]
        for (r in 1:p) {
          Xty[r] <- Xty[r] + X[i, r] * yv[i]
          for (c in 1:p) XtX[r, c] <- XtX[r, c] + X[i, r] * X[i, c]
        }
      }
    }
    Vninv <- B0inv + XtX
    Uv <- chol(Vninv)
    Vn <- inverse(Vninv)
    rhs <- B0invb0 + Xty
    bn <- (Vn %*% asCol(rhs))[, 1]
    # posterior InvGamma(shape = nu0 + nj/2, scale = s0 + quad/2)
    quad <- yty + b0Qb0 - inprod(bn, rhs)
    if (quad < 0) quad <- 0            # numeric guard, exact value is >= 0
    shape <- nu0 + nj / 2
    scale <- s0 + quad / 2
    s2new <- 1 / rgamma(1, shape = shape, rate = scale)
    # beta | s2 ~ N(bn, s2 * Vn): draw via z ~ N(0, I), beta = bn + sqrt(s2) L' z
    zdraw <- numeric(p)
    for (r in 1:p) zdraw[r] <- rnorm(1, 0, 1)
    UvVn <- chol(Vn)
    bnew <- bn + sqrt(s2new) * (t(UvVn) %*% asCol(zdraw))[, 1]
    values(model, s2Node) <<- c(s2new)
    values(model, betaNode) <<- bnew
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {})
)

# Random-effect variant of the NIG block sampler. A separate nimbleFunction
# (not a branch inside the plain one) because NIMBLE compiles every branch of
# run(): a reference to model[["b"]] would fail to compile against models
# that have no b node. Gate F4 measured why the offset is needed at all --
# conjugacy detection handles the additive *scalar* form but not the
# production inprod form, so the exact conditional must subtract the current
# random intercepts itself.
.betaS2ConjSamplerRE <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    j     <- control$j
    X     <- control$X
    yv    <- control$y
    grp   <- control$grp
    b0    <- control$b0
    B0inv <- control$B0inv
    nu0   <- control$nu0
    s0    <- control$s0
    p     <- ncol(X)
    n     <- nrow(X)
    betaNode <- paste0("betaTilde[", j, ", 1:", p, "]")
    s2Node   <- paste0("s2Tilde[", j, "]")
    calcNodes <- model$getDependencies(c(betaNode, s2Node))
    B0invb0 <- (B0inv %*% asCol(b0))[, 1]
    b0Qb0 <- inprod(b0, B0invb0)
  },
  run = function() {
    z <- model[["z"]]
    bcur <- model[["b"]]
    XtX <- matrix(0, p, p)
    Xty <- numeric(p)
    yty <- 0
    nj  <- 0
    for (i in 1:n) {
      if (z[i] == j) {
        nj <- nj + 1
        ri <- yv[i] - bcur[grp[i]]
        yty <- yty + ri * ri
        for (r in 1:p) {
          Xty[r] <- Xty[r] + X[i, r] * ri
          for (c in 1:p) XtX[r, c] <- XtX[r, c] + X[i, r] * X[i, c]
        }
      }
    }
    Vninv <- B0inv + XtX
    Vn <- inverse(Vninv)
    rhs <- B0invb0 + Xty
    bn <- (Vn %*% asCol(rhs))[, 1]
    quad <- yty + b0Qb0 - inprod(bn, rhs)
    if (quad < 0) quad <- 0
    shape <- nu0 + nj / 2
    scale <- s0 + quad / 2
    s2new <- 1 / rgamma(1, shape = shape, rate = scale)
    zdraw <- numeric(p)
    for (r in 1:p) zdraw[r] <- rnorm(1, 0, 1)
    UvVn <- chol(Vn)
    bnew <- bn + sqrt(s2new) * (t(UvVn) %*% asCol(zdraw))[, 1]
    values(model, s2Node) <<- c(s2new)
    values(model, betaNode) <<- bnew
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes,
         logProb = TRUE)
  },
  methods = list(reset = function() {})
)

# Random-intercept + random-slope variant. A THIRD separate nimbleFunction,
# for the same reason as the second (NIMBLE compiles every branch of run(),
# so a reference to sRE would break models without that node). Gate F5.2
# measured the design: the effective response subtracts BOTH offsets, and the
# slope offset enters multiplied by its covariate -- a different structure
# from the intercept-only case.
.betaS2ConjSamplerRES <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    j     <- control$j
    X     <- control$X
    yv    <- control$y
    grp   <- control$grp
    xREv  <- control$xRE
    b0    <- control$b0
    B0inv <- control$B0inv
    nu0   <- control$nu0
    s0    <- control$s0
    p     <- ncol(X)
    n     <- nrow(X)
    betaNode <- paste0("betaTilde[", j, ", 1:", p, "]")
    s2Node   <- paste0("s2Tilde[", j, "]")
    calcNodes <- model$getDependencies(c(betaNode, s2Node))
    B0invb0 <- (B0inv %*% asCol(b0))[, 1]
    b0Qb0 <- inprod(b0, B0invb0)
  },
  run = function() {
    z <- model[["z"]]
    bcur <- model[["b"]]
    scur <- model[["sRE"]]
    XtX <- matrix(0, p, p)
    Xty <- numeric(p)
    yty <- 0
    nj  <- 0
    for (i in 1:n) {
      if (z[i] == j) {
        nj <- nj + 1
        ri <- yv[i] - bcur[grp[i]] - scur[grp[i]] * xREv[i]
        yty <- yty + ri * ri
        for (r in 1:p) {
          Xty[r] <- Xty[r] + X[i, r] * ri
          for (c in 1:p) XtX[r, c] <- XtX[r, c] + X[i, r] * X[i, c]
        }
      }
    }
    Vninv <- B0inv + XtX
    Vn <- inverse(Vninv)
    rhs <- B0invb0 + Xty
    bn <- (Vn %*% asCol(rhs))[, 1]
    quad <- yty + b0Qb0 - inprod(bn, rhs)
    if (quad < 0) quad <- 0
    shape <- nu0 + nj / 2
    scale <- s0 + quad / 2
    s2new <- 1 / rgamma(1, shape = shape, rate = scale)
    zdraw <- numeric(p)
    for (r in 1:p) zdraw[r] <- rnorm(1, 0, 1)
    UvVn <- chol(Vn)
    bnew <- bn + sqrt(s2new) * (t(UvVn) %*% asCol(zdraw))[, 1]
    values(model, s2Node) <<- c(s2new)
    values(model, betaNode) <<- bnew
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes,
         logProb = TRUE)
  },
  methods = list(reset = function() {})
)

#' @describeIn customizeSamplers Replace RW samplers on \code{betaTilde} and
#'   \code{s2Tilde} with the exact Normal-Inverse-Gamma conditional (FixedK
#'   path only; the DPM path already receives the conjugate CRP wrapper).
#' @export
setMethod("customizeSamplers", "NormalRegSpec",
  function(spec, conf, model, ...) {
    # only relevant when the allocation node is z (FixedK); DPM uses xi
    if (!("z" %in% model$getVarNames())) return(invisible(conf))
    if (!all(c("betaTilde", "s2Tilde") %in% model$getVarNames()))
      return(invisible(conf))
    K <- dim(model[["betaTilde"]])[1]
    X <- model$getConstants()$X
    yv <- model[["y"]]
    cn <- model$getConstants()
    B0inv <- solve(cn$B0)
    vn <- model$getVarNames()
    hasRE  <- "b" %in% vn
    hasRES <- "sRE" %in% vn
    conf$removeSamplers("betaTilde")
    conf$removeSamplers("s2Tilde")
    for (j in seq_len(K))
      conf$addSampler(target = paste0("betaTilde[", j, ", 1:", ncol(X), "]"),
                      type = if (hasRES) .betaS2ConjSamplerRES
                             else if (hasRE) .betaS2ConjSamplerRE
                             else .betaS2ConjSampler,
                      control = c(list(j = j, X = X, y = yv, b0 = cn$b0,
                                       B0inv = B0inv, nu0 = cn$nu0,
                                       s0 = cn$s0),
                                  if (hasRE) list(grp = cn$grp),
                                  if (hasRES) list(xRE = cn$xRE)))
    invisible(conf)
  })
