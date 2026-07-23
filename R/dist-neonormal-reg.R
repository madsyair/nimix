## Neo-normal regression: a mixture of skewed/heavy-tailed regressions.
##
## Each neo-normal family (MSNBurr, SEP, GMSNBurr, ...) becomes a regression by
## replacing the per-component LOCATION mu_k with a linear predictor
## X beta_k, while its SHAPE parameters (sigma, alpha, nu, theta -- the set
## differs by family) stay per-component. The construction is identical across
## families; only the density kernel and the list of shape parameters change.
##
## Rather than writing eleven S4 methods per family, the boilerplate six --
## buildDataList, componentDensity, extractParamTraces, relabelComponents,
## simulateParams, customizeSamplers -- are generated once by .neoRegMethods()
## from a small spec: the family's density function and its shape-parameter
## names. Only the parts that genuinely differ (prior, model code, constants,
## validation) are written per family, and those are short.
##
## There is no conjugate sampler for the coefficients (the likelihood is not
## Gaussian in beta), so NIMBLE's defaults run throughout; like the heavy-tail
## and GLM regressions, budget more iterations than a Gaussian fit.

# ---------------------------------------------------------------------------
# Generic method generator.
#
# `class`      : the S4 class name (e.g. "MSNBurrRegSpec")
# `densR`      : the R-level density function (e.g. dmsnburr), signature
#                (x, mu, <shape...>, log)
# `shape`      : character vector of shape-parameter names in model order,
#                e.g. c("sigma", "alpha") or c("sigma", "alpha", "theta").
#                The corresponding node names are <shape>Tilde.
# ---------------------------------------------------------------------------
# Registry of shape-parameter names per neo-normal regression class, so the
# predictive can pull the right traces generically. Populated by
# .neoRegMethods() as each family is defined.
.neoShapeRegistry <- new.env(parent = emptyenv())

# Return the shape-parameter names for a spec, or NULL if it is not a
# neo-normal regression family (Normal/GLM/... skip the shape path).
.neoShapeNames <- function(spec) {
  cl <- class(spec)[1L]
  if (exists(cl, envir = .neoShapeRegistry, inherits = FALSE))
    get(cl, envir = .neoShapeRegistry) else NULL
}

.neoRegMethods <- function(class, densR, shape, rngR = NULL,
                           densName = NULL, priorLines = NULL) {
  nodeOf <- function(s) paste0(s, "Tilde")
  assign(class, shape, envir = .neoShapeRegistry)   # for the predictive path

  setMethod("buildDataList", class,
    function(spec, data, ...) list(y = as.numeric(data)))

  # componentDensity: mu plus the family's shape params, by name
  setMethod("componentDensity", class,
    function(spec, ...) {
      function(x, params) {
        args <- c(list(x, params[["mu"]]),
                  lapply(shape, function(s) params[[s]]))
        do.call(densR, args)
      }
    })

  # extractParamTraces: beta plus one trace per shape parameter
  setMethod("extractParamTraces", class,
    function(spec, samples, L, d = NULL, prior = NULL, ...) {
      p <- if (!is.null(prior$p)) prior$p else d
      if (is.null(p)) stop("extractParamTraces() needs p.", call. = FALSE)
      out <- list(beta = .nodeToArray(samples, "betaTilde", c(L, p)),
                  p = p,
                  coefNames = if (!is.null(prior$coefNames)) prior$coefNames
                              else paste0("b", seq_len(p)),
                  shape = shape)
      for (s in shape) out[[s]] <- .nodeToArray(samples, nodeOf(s), L)
      out
    })

  # relabelComponents: coefficients and every shape parameter follow the same
  # permutation; the summary carries beta plus <shape>_mean columns
  setMethod("relabelComponents", class,
    function(spec, paramTrace, idx, occList, perms, modalK, weights, ...) {
      p <- paramTrace$p; betaTr <- paramTrace$beta; coefNm <- paramTrace$coefNames
      m <- length(idx)
      betaRe <- array(NA_real_, dim = c(m, modalK, p))
      shRe <- lapply(shape, function(s) matrix(NA_real_, m, modalK))
      names(shRe) <- shape
      for (t in seq_len(m)) {
        r <- idx[t]; occ <- occList[[t]]; ord <- perms[t, ]
        for (k in seq_len(modalK)) {
          j <- occ[ord[k]]
          betaRe[t, k, ] <- betaTr[r, j, ]
          for (s in shape) shRe[[s]][t, k] <- paramTrace[[s]][r, j]
        }
      }
      betaMean <- apply(betaRe, c(2L, 3L), mean)
      summ <- data.frame(component = seq_len(modalK), weight = colMeans(weights))
      for (j in seq_len(p)) summ[[coefNm[j]]] <- betaMean[, j]
      for (s in shape) summ[[paste0(s, "_mean")]] <- colMeans(shRe[[s]])
      list(summary = summ,
           relabeled = c(list(beta = betaRe), shRe))
    })

  # simulateParams: beta from the g-prior; each shape param from its prior,
  # read generically from the prior list (aX/bX naming resolved by the family)
  setMethod("simulateParams", class,
    function(spec, prior, nClust, ...) {
      p <- prior$p
      beta <- matrix(NA_real_, nClust, p)
      for (j in seq_len(nClust)) beta[j, ] <- .rmvnorm1(prior$b0, prior$B0)
      out <- list(beta = beta)
      draws <- prior$shapeDraw            # family supplies a sampler closure
      sh <- draws(nClust)
      for (s in shape) out[[s]] <- sh[[s]]
      out
    })

  setMethod("customizeSamplers", class,
    function(spec, conf, model, ...) invisible(conf))

  # DPM engine: generated generically from densName + priorLines (the same
  # ingredients the fixed-K model code uses). Emission identical; allocations
  # follow a CRP. Only added when the family supplies densName/priorLines.
  if (!is.null(densName) && !is.null(priorLines)) {
    setMethod("buildModelCode", signature(class, "DPMEngine"),
      function(spec, engine, n, L, re = FALSE, reSlope = FALSE, ...) {
        code <- .neoRegDPMCode(densName, shape, priorLines,
                               re = isTRUE(re),
                               reSlope = isTRUE(re) && isTRUE(reSlope))
        list(code = code,
             monitors = c("betaTilde", paste0(shape, "Tilde"), "xi", "alpha",
                          if (isTRUE(re)) c("b", "tauRE"),
                          if (isTRUE(re) && isTRUE(reSlope)) c("sRE", "tauSlope")),
             paramNodes = c(beta = "betaTilde"), allocNode = "xi")
      })
  }

  # responseRng: draw a response from the family given the per-observation
  # linear predictor (= location, identity link) and its shape parameters.
  # Needed by posteriorPredictive: a link-aware family must supply both its
  # mean (epred) and its sampling law (predictive). Without this method a
  # skewed predictive falls back to the Gaussian default and loses the
  # skew/tails.
  if (!is.null(rngR)) {
    setMethod("responseRng", class,
      function(spec, eta, s2 = NULL, prior = NULL, shapeVals = NULL, ...) {
        # shapeVals: named list of per-draw shape vectors (sigma, alpha, ...)
        args <- c(list(length(eta), eta),
                  lapply(shape, function(sp) shapeVals[[sp]]))
        do.call(rngR, args)
      })
  }

  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Generic model-code builder. The likelihood is
#   y[i] ~ dFAM_k(inprod(X[i,], beta[z[i],]), shape1[z[i]], shape2[z[i]], ...)
# with a g-prior on beta and each shape parameter given its own prior line,
# supplied as a language object by the family.
# ---------------------------------------------------------------------------
.neoRegFixedKCode <- function(densName, shape, shapePriorLines, re = FALSE,
                              reSlope = FALSE) {
  shapeArgs <- paste0(vapply(shape, function(s) sprintf("%sTilde[z[i]]", s), ""),
                      collapse = ", ")
  # A random intercept enters the same location the coefficients do -- b[grp[i]]
  # added to the linear predictor -- with a sum-to-zero parameterisation (the
  # component intercepts absorb the group means, so b is a pure deviation and
  # mixes well; cf. the Gaussian RE code in engine-fixedk.R). A random slope
  # adds sRE[grp[i]] * xRE[i], also sum-to-zero. The constraint matters: freely
  # parameterised offsets are identified only jointly with the fixed intercept
  # and slope, which produces a correlated posterior ridge and poor effective
  # sample size even when a conjugate sampler is available. Neo-normal
  # emissions are location-scale, so location random effects mean the same
  # thing they do for a Gaussian, unlike a GLM where they enter via the link.
  loc <- "inprod(X[i, 1:p], betaObs[i, 1:p])"
  if (re) loc <- paste0(loc, " + b[grp[i]]")
  if (re && reSlope) loc <- paste0(loc, " + sRE[grp[i]] * xRE[i]")
  reBlock <- ""
  if (re && reSlope) {
    reBlock <- "
    for (g in 1:(G - 1)) {
      bf[g] ~ dnorm(0, sd = tauRE)
      sf[g] ~ dnorm(0, sd = tauSlope)
    }
    b[1:(G - 1)] <- bf[1:(G - 1)]
    b[G] <- -sum(bf[1:(G - 1)])
    sRE[1:(G - 1)] <- sf[1:(G - 1)]
    sRE[G] <- -sum(sf[1:(G - 1)])
    tauRE ~ dunif(tauMin, tauMax)
    tauSlope ~ dunif(tauMinSlope, tauMaxSlope)"
  } else if (re) {
    reBlock <- "
    for (g in 1:(G - 1)) bf[g] ~ dnorm(0, sd = tauRE)
    b[1:(G - 1)] <- bf[1:(G - 1)]
    b[G] <- -sum(bf[1:(G - 1)])
    tauRE ~ dunif(tauMin, tauMax)"
  }
  tmpl <- sprintf("{
    for (i in 1:n) {
      z[i] ~ dcat(weights[1:K])
      betaObs[i, 1:p] <- betaTilde[z[i], 1:p]
      muObs[i] <- %s
      y[i] ~ %s(muObs[i], %s)
    }
    weights[1:K] ~ ddirch(alphaVec[1:K])
    for (j in 1:K) {
      betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = B0[1:p, 1:p])
%s
    }%s
  }", loc, densName, shapeArgs, shapePriorLines, reBlock)
  str2lang(tmpl)
}

# Generic DPM model code: same emission as fixed-K but the allocations follow a
# CRP (Dirichlet process) instead of a fixed-K categorical. There is no
# conjugate cluster wrapper here (the neo-normal likelihood is not Gaussian in
# the coefficients), so NIMBLE samples the cluster parameters directly.
.neoRegDPMCode <- function(densName, shape, shapePriorLines, re = FALSE,
                           reSlope = FALSE) {
  shapeArgs <- paste0(vapply(shape, function(s) sprintf("%sTilde[xi[i]]", s), ""),
                      collapse = ", ")
  # RE enters the location exactly as in the fixed-K code: the DPM allocates
  # observations to clusters (iid, exchangeable), so an external grouping factor
  # and its random offsets sit cleanly on top of the linear predictor. The
  # random number of clusters does not interact with the fixed grouping.
  loc <- "inprod(X[i, 1:p], betaObs[i, 1:p])"
  if (re) loc <- paste0(loc, " + b[grp[i]]")
  if (re && reSlope) loc <- paste0(loc, " + sRE[grp[i]] * xRE[i]")
  reBlock <- ""
  if (re && reSlope) {
    reBlock <- "
    for (g in 1:(G - 1)) {
      bf[g] ~ dnorm(0, sd = tauRE)
      sf[g] ~ dnorm(0, sd = tauSlope)
    }
    b[1:(G - 1)] <- bf[1:(G - 1)]
    b[G] <- -sum(bf[1:(G - 1)])
    sRE[1:(G - 1)] <- sf[1:(G - 1)]
    sRE[G] <- -sum(sf[1:(G - 1)])
    tauRE ~ dunif(tauMin, tauMax)
    tauSlope ~ dunif(tauMinSlope, tauMaxSlope)"
  } else if (re) {
    reBlock <- "
    for (g in 1:(G - 1)) bf[g] ~ dnorm(0, sd = tauRE)
    b[1:(G - 1)] <- bf[1:(G - 1)]
    b[G] <- -sum(bf[1:(G - 1)])
    tauRE ~ dunif(tauMin, tauMax)"
  }
  tmpl <- sprintf("{
    for (i in 1:n) {
      betaObs[i, 1:p] <- betaTilde[xi[i], 1:p]
      muObs[i] <- %s
      y[i] ~ %s(muObs[i], %s)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
    for (j in 1:L) {
      betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = B0[1:p, 1:p])
%s
    }
    alpha ~ dgamma(shape = aAlpha, rate = bAlpha)%s
  }", loc, densName, shapeArgs, shapePriorLines, reBlock)
  str2lang(tmpl)
}

# Generic g-prior on the coefficients plus a family-supplied shape prior block.
# `shapePrior` is a function(y, control) returning a named list of hyper-
# parameters; `shapeDraw` a function(prior) returning function(nClust) that
# samples the shape parameters (used by simulateParams).
.neoRegPrior <- function(data, control, shapePrior) {
  y <- as.numeric(data)
  X <- control$X
  if (is.null(X) || !is.matrix(X))
    stop("neo-normal regression defaultPrior needs the design matrix in ",
         "control$X.", call. = FALSE)
  n <- length(y); p <- ncol(X)
  g <- if (!is.null(control$g)) control$g else n
  XtX <- crossprod(X)
  ridge <- 1e-8 * mean(diag(XtX)) + 1e-10
  B0 <- g * solve(XtX + diag(ridge, p))
  c(list(p = p, X = X, b0 = rep(0, p), B0 = B0), shapePrior(y, control))
}

# Append the random-effect constants (group index, group count, tau bounds, and
# the slope covariate if present) to a family's constant list. Every neo-normal
# family's buildConstants can share this -- the RE part is identical; only the
# shape hyperparameters differ. Keeps F10 support uniform across the framework.
.neoRegREConstants <- function(out, prior) {
  if (isTRUE(prior$hasRE)) {
    out$grp <- prior$reGrp; out$G <- prior$reG
    out$tauMin <- prior$tauMin; out$tauMax <- prior$tauMax
    if (isTRUE(prior$hasRESlope)) {
      out$xRE <- prior$reSlopeX
      out$tauMinSlope <- prior$tauMinSlope
      out$tauMaxSlope <- prior$tauMaxSlope
    }
  }
  out
}

# ===========================================================================
# MSNBurr regression (first family through the framework)
# ===========================================================================

#' @rdname MSNBurrRegSpec-class
#' @export
setClass("MSNBurrRegSpec", contains = "DistributionSpec")

#' MSNBurr regression specification
#'
#' A mixture of MSNBurr-IIa regressions: each component is a skewed regression
#' with its own coefficient vector, scale, and shape.
#' @name MSNBurrRegSpec-class
#' @return An object used internally by \code{\link{nimixReg}}.
#' @export
MSNBurrRegSpec <- function() new("MSNBurrRegSpec")

#' @describeIn isRegressionSpec MSNBurr regression is a regression spec.
setMethod("isRegressionSpec", "MSNBurrRegSpec", function(spec) TRUE)

.msnburrRegShapePrior <- function(y, control) {
  sy <- stats::sd(y); if (!is.finite(sy) || sy <= 0) sy <- 1
  aSig <- if (!is.null(control$aSig)) control$aSig else 3
  bSig <- if (!is.null(control$bSig)) control$bSig else (aSig - 1) * sy
  aA <- if (!is.null(control$aAlphaShape)) control$aAlphaShape else 2
  bA <- if (!is.null(control$bAlphaRate))  control$bAlphaRate  else 2
  list(aSig = aSig, bSig = bSig, aA = aA, bA = bA,
       shapeDraw = function(nClust) list(
         sigma = 1 / stats::rgamma(nClust, shape = aSig, rate = 1 / bSig),
         alpha = stats::rgamma(nClust, shape = aA, rate = bA)))
}

#' @describeIn defaultPrior MSNBurr regression prior.
#' @export
setMethod("defaultPrior", "MSNBurrRegSpec",
  function(spec, data, control = list(), ...)
    .neoRegPrior(data, control, .msnburrRegShapePrior))

#' @describeIn buildConstants MSNBurr regression constants.
#' @export
setMethod("buildConstants", "MSNBurrRegSpec",
  function(spec, prior, n, ...) {
    out <- list(n = n, p = prior$p, X = prior$X, b0 = prior$b0, B0 = prior$B0,
                aSig = prior$aSig, bSig = prior$bSig, aA = prior$aA, bA = prior$bA)
    if (isTRUE(prior$hasRE)) {
      out$grp <- prior$reGrp; out$G <- prior$reG
      out$tauMin <- prior$tauMin; out$tauMax <- prior$tauMax
      if (isTRUE(prior$hasRESlope)) {
        out$xRE <- prior$reSlopeX
        out$tauMinSlope <- prior$tauMinSlope
        out$tauMaxSlope <- prior$tauMaxSlope
      }
    }
    out
  })

#' @describeIn validateParams MSNBurr regression prior validation.
#' @export
setMethod("validateParams", "MSNBurrRegSpec",
  function(spec, params, ...) {
    req <- c("b0", "B0", "aSig", "bSig", "aA", "bA", "p")
    miss <- setdiff(req, names(params))
    if (length(miss))
      stop("MSNBurrRegSpec prior is missing: ", paste(miss, collapse = ", "),
           call. = FALSE)
    if (params$aSig <= 1 || params$bSig <= 0 || params$aA <= 0 || params$bA <= 0)
      stop("MSNBurr scale/shape hyperparameters must be positive (aSig > 1).",
           call. = FALSE)
    invisible(TRUE)
  })

#' @describeIn buildModelCode MSNBurr regression finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("MSNBurrRegSpec", "FixedKEngine"),
  function(spec, engine, n, L, re = FALSE, reSlope = FALSE, ...) {
    code <- .neoRegFixedKCode("dMSNBurr_k", c("sigma", "alpha"),
      paste("      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)",
            "      alphaTilde[j] ~ dgamma(shape = aA, rate = bA)", sep = "\n"),
      re = isTRUE(re), reSlope = isTRUE(re) && isTRUE(reSlope))
    list(code = code,
         monitors = c("betaTilde", "sigmaTilde", "alphaTilde", "weights", "z",
                      if (isTRUE(re)) c("b", "tauRE"),
                      if (isTRUE(re) && isTRUE(reSlope)) c("sRE", "tauSlope")),
         paramNodes = c(beta = "betaTilde"), allocNode = "z")
  })

#' @describeIn linkInv MSNBurr regression identity link.
#' @export
setMethod("linkInv", "MSNBurrRegSpec",
  function(spec, eta, prior = NULL, ...) eta)

#' @describeIn componentInits MSNBurr regression initial values.
#' @export
setMethod("componentInits", "MSNBurrRegSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    y <- as.numeric(data); X <- prior$X; p <- prior$p
    bhat <- tryCatch(as.numeric(solve(crossprod(X) + diag(1e-6, p),
                                      crossprod(X, y))),
                     error = function(e) rep(0, p))
    beta <- matrix(rep(bhat, each = count), nrow = count)
    sig0 <- prior$bSig / (prior$aSig - 1)
    xiInit <- .initClusters(y, min(count, 2L), initMethod)
    if (is.null(xiInit)) xiInit <- rep(1L, length(y))
    list(alloc = xiInit,
         params = list(betaTilde = beta, sigmaTilde = rep(sig0, count),
                       alphaTilde = rep(1, count)))
  })

# generate the boilerplate six
.neoRegMethods("MSNBurrRegSpec", dmsnburr, c("sigma", "alpha"), rmsnburr,
               densName = "dMSNBurr_k", priorLines = paste("      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)",
            "      alphaTilde[j] ~ dgamma(shape = aA, rate = bA)", sep = "\n"))

# ===========================================================================
# SEP regression (second family: different shape set -- sigma, nu)
# ===========================================================================

#' @rdname SEPRegSpec-class
#' @export
setClass("SEPRegSpec", contains = "DistributionSpec")

#' SEP regression specification
#' @name SEPRegSpec-class
#' @return An object used internally by \code{\link{nimixReg}}.
#' @export
SEPRegSpec <- function() new("SEPRegSpec")

#' @describeIn isRegressionSpec SEP regression is a regression spec.
setMethod("isRegressionSpec", "SEPRegSpec", function(spec) TRUE)

.sepRegShapePrior <- function(y, control) {
  sy <- stats::sd(y); if (!is.finite(sy) || sy <= 0) sy <- 1
  aSig <- if (!is.null(control$aSig)) control$aSig else 3
  bSig <- if (!is.null(control$bSig)) control$bSig else (aSig - 1) * sy
  aNu <- if (!is.null(control$aNu)) control$aNu else 2
  bNu <- if (!is.null(control$bNu)) control$bNu else 1
  list(aSig = aSig, bSig = bSig, aNu = aNu, bNu = bNu,
       shapeDraw = function(nClust) list(
         sigma = 1 / stats::rgamma(nClust, shape = aSig, rate = 1 / bSig),
         nu = stats::rgamma(nClust, shape = aNu, rate = bNu)))
}

#' @describeIn defaultPrior SEP regression prior.
#' @export
setMethod("defaultPrior", "SEPRegSpec",
  function(spec, data, control = list(), ...)
    .neoRegPrior(data, control, .sepRegShapePrior))

#' @describeIn buildConstants SEP regression constants.
#' @export
setMethod("buildConstants", "SEPRegSpec",
  function(spec, prior, n, ...)
    .neoRegREConstants(list(n = n, p = prior$p, X = prior$X, b0 = prior$b0, B0 = prior$B0,
         aSig = prior$aSig, bSig = prior$bSig, aNu = prior$aNu, bNu = prior$bNu), prior))

#' @describeIn validateParams SEP regression validation.
#' @export
setMethod("validateParams", "SEPRegSpec",
  function(spec, params, ...) {
    req <- c("b0", "B0", "aSig", "bSig", "aNu", "bNu", "p")
    miss <- setdiff(req, names(params))
    if (length(miss))
      stop("SEPRegSpec prior is missing: ", paste(miss, collapse = ", "),
           call. = FALSE)
    invisible(TRUE)
  })

#' @describeIn buildModelCode SEP regression finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("SEPRegSpec", "FixedKEngine"),
  function(spec, engine, n, L, re = FALSE, reSlope = FALSE, ...) {
    code <- .neoRegFixedKCode("dSEP_k", c("sigma", "nu"),
      paste("      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)",
            "      nuTilde[j] ~ dgamma(shape = aNu, rate = bNu)", sep = "\n"),
      re = isTRUE(re), reSlope = isTRUE(re) && isTRUE(reSlope))
    list(code = code,
         monitors = c("betaTilde", "sigmaTilde", "nuTilde", "weights", "z",
                      if (isTRUE(re)) c("b", "tauRE"),
                      if (isTRUE(re) && isTRUE(reSlope)) c("sRE", "tauSlope")),
         paramNodes = c(beta = "betaTilde"), allocNode = "z")
  })

#' @describeIn linkInv SEP regression identity link.
#' @export
setMethod("linkInv", "SEPRegSpec", function(spec, eta, prior = NULL, ...) eta)

#' @describeIn componentInits SEP regression initial values.
#' @export
setMethod("componentInits", "SEPRegSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    y <- as.numeric(data); X <- prior$X; p <- prior$p
    bhat <- tryCatch(as.numeric(solve(crossprod(X) + diag(1e-6, p),
                                      crossprod(X, y))),
                     error = function(e) rep(0, p))
    beta <- matrix(rep(bhat, each = count), nrow = count)
    sig0 <- prior$bSig / (prior$aSig - 1)
    xiInit <- .initClusters(y, min(count, 2L), initMethod)
    if (is.null(xiInit)) xiInit <- rep(1L, length(y))
    list(alloc = xiInit,
         params = list(betaTilde = beta, sigmaTilde = rep(sig0, count),
                       nuTilde = rep(2, count)))
  })

.neoRegMethods("SEPRegSpec", dsep, c("sigma", "nu"), rsep,
               densName = "dSEP_k", priorLines = paste("      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)",
            "      nuTilde[j] ~ dgamma(shape = aNu, rate = bNu)", sep = "\n"))

# ===========================================================================
# MSNBurr-IIa regression (same shape set as MSNBurr: sigma, alpha; only the
# density kernel differs -- dMSNBurr2a_k / dmsnburr2a)
# ===========================================================================

#' @rdname MSNBurr2aRegSpec-class
#' @export
setClass("MSNBurr2aRegSpec", contains = "DistributionSpec")

#' MSNBurr-IIa regression specification
#' @name MSNBurr2aRegSpec-class
#' @return An object used internally by \code{\link{nimixReg}}.
#' @export
MSNBurr2aRegSpec <- function() new("MSNBurr2aRegSpec")

#' @describeIn isRegressionSpec MSNBurr-IIa regression is a regression spec.
setMethod("isRegressionSpec", "MSNBurr2aRegSpec", function(spec) TRUE)

#' @describeIn defaultPrior MSNBurr-IIa regression prior (shares the MSNBurr
#'   scale/shape prior).
#' @export
setMethod("defaultPrior", "MSNBurr2aRegSpec",
  function(spec, data, control = list(), ...)
    .neoRegPrior(data, control, .msnburrRegShapePrior))

#' @describeIn buildConstants MSNBurr-IIa regression constants.
#' @export
setMethod("buildConstants", "MSNBurr2aRegSpec",
  function(spec, prior, n, ...)
    .neoRegREConstants(list(n = n, p = prior$p, X = prior$X, b0 = prior$b0, B0 = prior$B0,
         aSig = prior$aSig, bSig = prior$bSig, aA = prior$aA, bA = prior$bA), prior))

#' @describeIn validateParams MSNBurr-IIa regression validation.
#' @export
setMethod("validateParams", "MSNBurr2aRegSpec",
  function(spec, params, ...) {
    req <- c("b0", "B0", "aSig", "bSig", "aA", "bA", "p")
    miss <- setdiff(req, names(params))
    if (length(miss))
      stop("MSNBurr2aRegSpec prior is missing: ", paste(miss, collapse = ", "),
           call. = FALSE)
    if (params$aSig <= 1 || params$bSig <= 0 || params$aA <= 0 || params$bA <= 0)
      stop("MSNBurr-IIa scale/shape hyperparameters must be positive.",
           call. = FALSE)
    invisible(TRUE)
  })

#' @describeIn buildModelCode MSNBurr-IIa regression finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("MSNBurr2aRegSpec", "FixedKEngine"),
  function(spec, engine, n, L, re = FALSE, reSlope = FALSE, ...) {
    code <- .neoRegFixedKCode("dMSNBurr2a_k", c("sigma", "alpha"),
      paste("      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)",
            "      alphaTilde[j] ~ dgamma(shape = aA, rate = bA)", sep = "\n"),
      re = isTRUE(re), reSlope = isTRUE(re) && isTRUE(reSlope))
    list(code = code,
         monitors = c("betaTilde", "sigmaTilde", "alphaTilde", "weights", "z",
                      if (isTRUE(re)) c("b", "tauRE"),
                      if (isTRUE(re) && isTRUE(reSlope)) c("sRE", "tauSlope")),
         paramNodes = c(beta = "betaTilde"), allocNode = "z")
  })

#' @describeIn linkInv MSNBurr-IIa regression identity link.
#' @export
setMethod("linkInv", "MSNBurr2aRegSpec",
  function(spec, eta, prior = NULL, ...) eta)

#' @describeIn componentInits MSNBurr-IIa regression initial values.
#' @export
setMethod("componentInits", "MSNBurr2aRegSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    y <- as.numeric(data); X <- prior$X; p <- prior$p
    bhat <- tryCatch(as.numeric(solve(crossprod(X) + diag(1e-6, p),
                                      crossprod(X, y))),
                     error = function(e) rep(0, p))
    beta <- matrix(rep(bhat, each = count), nrow = count)
    sig0 <- prior$bSig / (prior$aSig - 1)
    xiInit <- .initClusters(y, min(count, 2L), initMethod)
    if (is.null(xiInit)) xiInit <- rep(1L, length(y))
    list(alloc = xiInit,
         params = list(betaTilde = beta, sigmaTilde = rep(sig0, count),
                       alphaTilde = rep(1, count)))
  })

.neoRegMethods("MSNBurr2aRegSpec", dmsnburr2a, c("sigma", "alpha"), rmsnburr2a,
               densName = "dMSNBurr2a_k", priorLines = paste("      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)",
            "      alphaTilde[j] ~ dgamma(shape = aA, rate = bA)", sep = "\n"))

# ===========================================================================
# FSSN regression (sigma, alpha -- but alpha has a LOG-NORMAL prior, not the
# Gamma used by MSNBurr; shows the framework handling a different shape prior)
# ===========================================================================

#' @rdname FSSNRegSpec-class
#' @export
setClass("FSSNRegSpec", contains = "DistributionSpec")

#' FSSN regression specification
#' @name FSSNRegSpec-class
#' @return An object used internally by \code{\link{nimixReg}}.
#' @export
FSSNRegSpec <- function() new("FSSNRegSpec")

#' @describeIn isRegressionSpec FSSN regression is a regression spec.
setMethod("isRegressionSpec", "FSSNRegSpec", function(spec) TRUE)

.fssnRegShapePrior <- function(y, control) {
  sy <- stats::sd(y); if (!is.finite(sy) || sy <= 0) sy <- 1
  aSig <- if (!is.null(control$aSig)) control$aSig else 3
  bSig <- if (!is.null(control$bSig)) control$bSig else (aSig - 1) * sy
  aScale <- if (!is.null(control$alphaLogSd)) control$alphaLogSd else 1
  list(aSig = aSig, bSig = bSig, aScale = aScale,
       shapeDraw = function(nClust) list(
         sigma = 1 / stats::rgamma(nClust, shape = aSig, rate = 1 / bSig),
         alpha = stats::rlnorm(nClust, 0, aScale)))
}

#' @describeIn defaultPrior FSSN regression prior (log-normal skewness).
#' @export
setMethod("defaultPrior", "FSSNRegSpec",
  function(spec, data, control = list(), ...)
    .neoRegPrior(data, control, .fssnRegShapePrior))

#' @describeIn buildConstants FSSN regression constants.
#' @export
setMethod("buildConstants", "FSSNRegSpec",
  function(spec, prior, n, ...)
    .neoRegREConstants(list(n = n, p = prior$p, X = prior$X, b0 = prior$b0, B0 = prior$B0,
         aSig = prior$aSig, bSig = prior$bSig, aScale = prior$aScale), prior))

#' @describeIn validateParams FSSN regression validation.
#' @export
setMethod("validateParams", "FSSNRegSpec",
  function(spec, params, ...) {
    req <- c("b0", "B0", "aSig", "bSig", "aScale", "p")
    miss <- setdiff(req, names(params))
    if (length(miss))
      stop("FSSNRegSpec prior is missing: ", paste(miss, collapse = ", "),
           call. = FALSE)
    if (params$aSig <= 1 || params$bSig <= 0 || params$aScale <= 0)
      stop("FSSN scale/shape hyperparameters must be positive (aSig > 1).",
           call. = FALSE)
    invisible(TRUE)
  })

#' @describeIn buildModelCode FSSN regression finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("FSSNRegSpec", "FixedKEngine"),
  function(spec, engine, n, L, re = FALSE, reSlope = FALSE, ...) {
    code <- .neoRegFixedKCode("dFSSN_k", c("sigma", "alpha"),
      paste("      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)",
            "      alphaTilde[j] ~ dlnorm(0, sd = aScale)", sep = "\n"),
      re = isTRUE(re), reSlope = isTRUE(re) && isTRUE(reSlope))
    list(code = code,
         monitors = c("betaTilde", "sigmaTilde", "alphaTilde", "weights", "z",
                      if (isTRUE(re)) c("b", "tauRE"),
                      if (isTRUE(re) && isTRUE(reSlope)) c("sRE", "tauSlope")),
         paramNodes = c(beta = "betaTilde"), allocNode = "z")
  })

#' @describeIn linkInv FSSN regression identity link.
#' @export
setMethod("linkInv", "FSSNRegSpec", function(spec, eta, prior = NULL, ...) eta)

#' @describeIn componentInits FSSN regression initial values.
#' @export
setMethod("componentInits", "FSSNRegSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    y <- as.numeric(data); X <- prior$X; p <- prior$p
    bhat <- tryCatch(as.numeric(solve(crossprod(X) + diag(1e-6, p),
                                      crossprod(X, y))),
                     error = function(e) rep(0, p))
    beta <- matrix(rep(bhat, each = count), nrow = count)
    sig0 <- prior$bSig / (prior$aSig - 1)
    xiInit <- .initClusters(y, min(count, 2L), initMethod)
    if (is.null(xiInit)) xiInit <- rep(1L, length(y))
    list(alloc = xiInit,
         params = list(betaTilde = beta, sigmaTilde = rep(sig0, count),
                       alphaTilde = rep(1, count)))
  })

.neoRegMethods("FSSNRegSpec", dfssn, c("sigma", "alpha"), rfssn,
               densName = "dFSSN_k", priorLines = paste("      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)",
            "      alphaTilde[j] ~ dlnorm(0, sd = aScale)", sep = "\n"))

# ===========================================================================
# GMSNBurr regression (THREE shape parameters: sigma, alpha, theta -- the
# first family to exercise the framework beyond two shape parameters)
# ===========================================================================

#' @rdname GMSNBurrRegSpec-class
#' @export
setClass("GMSNBurrRegSpec", contains = "DistributionSpec")

#' GMSNBurr regression specification
#' @name GMSNBurrRegSpec-class
#' @return An object used internally by \code{\link{nimixReg}}.
#' @export
GMSNBurrRegSpec <- function() new("GMSNBurrRegSpec")

#' @describeIn isRegressionSpec GMSNBurr regression is a regression spec.
setMethod("isRegressionSpec", "GMSNBurrRegSpec", function(spec) TRUE)

.gmsnburrRegShapePrior <- function(y, control) {
  sy <- stats::sd(y); if (!is.finite(sy) || sy <= 0) sy <- 1
  aSig <- if (!is.null(control$aSig)) control$aSig else 3
  bSig <- if (!is.null(control$bSig)) control$bSig else (aSig - 1) * sy
  aA <- if (!is.null(control$aAlphaShape)) control$aAlphaShape else 2
  bA <- if (!is.null(control$bAlphaRate))  control$bAlphaRate  else 2
  aT <- if (!is.null(control$aThetaShape)) control$aThetaShape else 2
  bT <- if (!is.null(control$bThetaRate))  control$bThetaRate  else 2
  list(aSig = aSig, bSig = bSig, aA = aA, bA = bA, aT = aT, bT = bT,
       shapeDraw = function(nClust) list(
         sigma = 1 / stats::rgamma(nClust, shape = aSig, rate = 1 / bSig),
         alpha = stats::rgamma(nClust, shape = aA, rate = bA),
         theta = stats::rgamma(nClust, shape = aT, rate = bT)))
}

#' @describeIn defaultPrior GMSNBurr regression prior.
#' @export
setMethod("defaultPrior", "GMSNBurrRegSpec",
  function(spec, data, control = list(), ...)
    .neoRegPrior(data, control, .gmsnburrRegShapePrior))

#' @describeIn buildConstants GMSNBurr regression constants.
#' @export
setMethod("buildConstants", "GMSNBurrRegSpec",
  function(spec, prior, n, ...)
    .neoRegREConstants(list(n = n, p = prior$p, X = prior$X, b0 = prior$b0, B0 = prior$B0,
         aSig = prior$aSig, bSig = prior$bSig, aA = prior$aA, bA = prior$bA,
         aT = prior$aT, bT = prior$bT), prior))

#' @describeIn validateParams GMSNBurr regression validation.
#' @export
setMethod("validateParams", "GMSNBurrRegSpec",
  function(spec, params, ...) {
    req <- c("b0", "B0", "aSig", "bSig", "aA", "bA", "aT", "bT", "p")
    miss <- setdiff(req, names(params))
    if (length(miss))
      stop("GMSNBurrRegSpec prior is missing: ", paste(miss, collapse = ", "),
           call. = FALSE)
    if (params$aSig <= 1 || params$bSig <= 0 || params$aA <= 0 ||
        params$bA <= 0 || params$aT <= 0 || params$bT <= 0)
      stop("GMSNBurr scale/shape hyperparameters must be positive (aSig > 1).",
           call. = FALSE)
    invisible(TRUE)
  })

#' @describeIn buildModelCode GMSNBurr regression finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("GMSNBurrRegSpec", "FixedKEngine"),
  function(spec, engine, n, L, re = FALSE, reSlope = FALSE, ...) {
    code <- .neoRegFixedKCode("dGMSNBurr_k", c("sigma", "alpha", "theta"),
      paste("      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)",
            "      alphaTilde[j] ~ dgamma(shape = aA, rate = bA)",
            "      thetaTilde[j] ~ dgamma(shape = aT, rate = bT)", sep = "\n"),
      re = isTRUE(re), reSlope = isTRUE(re) && isTRUE(reSlope))
    list(code = code,
         monitors = c("betaTilde", "sigmaTilde", "alphaTilde", "thetaTilde", "weights", "z",
                      if (isTRUE(re)) c("b", "tauRE"),
                      if (isTRUE(re) && isTRUE(reSlope)) c("sRE", "tauSlope")),
         paramNodes = c(beta = "betaTilde"), allocNode = "z")
  })

#' @describeIn linkInv GMSNBurr regression identity link.
#' @export
setMethod("linkInv", "GMSNBurrRegSpec",
  function(spec, eta, prior = NULL, ...) eta)

#' @describeIn componentInits GMSNBurr regression initial values.
#' @export
setMethod("componentInits", "GMSNBurrRegSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    y <- as.numeric(data); X <- prior$X; p <- prior$p
    bhat <- tryCatch(as.numeric(solve(crossprod(X) + diag(1e-6, p),
                                      crossprod(X, y))),
                     error = function(e) rep(0, p))
    beta <- matrix(rep(bhat, each = count), nrow = count)
    sig0 <- prior$bSig / (prior$aSig - 1)
    xiInit <- .initClusters(y, min(count, 2L), initMethod)
    if (is.null(xiInit)) xiInit <- rep(1L, length(y))
    list(alloc = xiInit,
         params = list(betaTilde = beta, sigmaTilde = rep(sig0, count),
                       alphaTilde = rep(1, count), thetaTilde = rep(1, count)))
  })

.neoRegMethods("GMSNBurrRegSpec", dgmsnburr, c("sigma", "alpha", "theta"), rgmsnburr,
               densName = "dGMSNBurr_k", priorLines = paste("      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)",
            "      alphaTilde[j] ~ dgamma(shape = aA, rate = bA)",
            "      thetaTilde[j] ~ dgamma(shape = aT, rate = bT)", sep = "\n"))

# ===========================================================================
# LEP regression (sigma, nu -- like SEP, shares the .epDefaultPrior scaling)
# ===========================================================================

#' @rdname LEPRegSpec-class
#' @export
setClass("LEPRegSpec", contains = "DistributionSpec")

#' LEP regression specification
#' @name LEPRegSpec-class
#' @return An object used internally by \code{\link{nimixReg}}.
#' @export
LEPRegSpec <- function() new("LEPRegSpec")

#' @describeIn isRegressionSpec LEP regression is a regression spec.
setMethod("isRegressionSpec", "LEPRegSpec", function(spec) TRUE)

.lepRegShapePrior <- function(y, control) {
  sy <- stats::sd(y); if (!is.finite(sy) || sy <= 0) sy <- 1
  aSig <- if (!is.null(control$aSig)) control$aSig else 3
  bSig <- if (!is.null(control$bSig)) control$bSig else (aSig - 1) * sy
  aNu <- if (!is.null(control$aNuShape)) control$aNuShape else 4
  bNu <- if (!is.null(control$bNuRate))  control$bNuRate  else 2
  list(aSig = aSig, bSig = bSig, aNu = aNu, bNu = bNu,
       shapeDraw = function(nClust) list(
         sigma = 1 / stats::rgamma(nClust, shape = aSig, rate = 1 / bSig),
         nu = stats::rgamma(nClust, shape = aNu, rate = bNu)))
}

#' @describeIn defaultPrior LEP regression prior.
#' @export
setMethod("defaultPrior", "LEPRegSpec",
  function(spec, data, control = list(), ...)
    .neoRegPrior(data, control, .lepRegShapePrior))

#' @describeIn buildConstants LEP regression constants.
#' @export
setMethod("buildConstants", "LEPRegSpec",
  function(spec, prior, n, ...)
    .neoRegREConstants(list(n = n, p = prior$p, X = prior$X, b0 = prior$b0, B0 = prior$B0,
         aSig = prior$aSig, bSig = prior$bSig, aNu = prior$aNu, bNu = prior$bNu), prior))

#' @describeIn validateParams LEP regression validation.
#' @export
setMethod("validateParams", "LEPRegSpec",
  function(spec, params, ...) {
    req <- c("b0", "B0", "aSig", "bSig", "aNu", "bNu", "p")
    miss <- setdiff(req, names(params))
    if (length(miss))
      stop("LEPRegSpec prior is missing: ", paste(miss, collapse = ", "),
           call. = FALSE)
    invisible(TRUE)
  })

#' @describeIn buildModelCode LEP regression finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("LEPRegSpec", "FixedKEngine"),
  function(spec, engine, n, L, re = FALSE, reSlope = FALSE, ...) {
    code <- .neoRegFixedKCode("dLEP_k", c("sigma", "nu"),
      paste("      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)",
            "      nuTilde[j] ~ dgamma(shape = aNu, rate = bNu)", sep = "\n"),
      re = isTRUE(re), reSlope = isTRUE(re) && isTRUE(reSlope))
    list(code = code,
         monitors = c("betaTilde", "sigmaTilde", "nuTilde", "weights", "z",
                      if (isTRUE(re)) c("b", "tauRE"),
                      if (isTRUE(re) && isTRUE(reSlope)) c("sRE", "tauSlope")),
         paramNodes = c(beta = "betaTilde"), allocNode = "z")
  })

#' @describeIn linkInv LEP regression identity link.
#' @export
setMethod("linkInv", "LEPRegSpec", function(spec, eta, prior = NULL, ...) eta)

#' @describeIn componentInits LEP regression initial values.
#' @export
setMethod("componentInits", "LEPRegSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    y <- as.numeric(data); X <- prior$X; p <- prior$p
    bhat <- tryCatch(as.numeric(solve(crossprod(X) + diag(1e-6, p),
                                      crossprod(X, y))),
                     error = function(e) rep(0, p))
    beta <- matrix(rep(bhat, each = count), nrow = count)
    sig0 <- prior$bSig / (prior$aSig - 1)
    xiInit <- .initClusters(y, min(count, 2L), initMethod)
    if (is.null(xiInit)) xiInit <- rep(1L, length(y))
    list(alloc = xiInit,
         params = list(betaTilde = beta, sigmaTilde = rep(sig0, count),
                       nuTilde = rep(2, count)))
  })

.neoRegMethods("LEPRegSpec", dlep, c("sigma", "nu"), rlep,
               densName = "dLEP_k", priorLines = paste("      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)",
            "      nuTilde[j] ~ dgamma(shape = aNu, rate = bNu)", sep = "\n"))

# ===========================================================================
# FSST regression (sigma, alpha, nu -- alpha log-normal, nu truncated-Gamma>2)
# ===========================================================================

#' @rdname FSSTRegSpec-class
#' @export
setClass("FSSTRegSpec", contains = "DistributionSpec")

#' FSST regression specification
#' @name FSSTRegSpec-class
#' @return An object used internally by \code{\link{nimixReg}}.
#' @export
FSSTRegSpec <- function() new("FSSTRegSpec")

#' @describeIn isRegressionSpec FSST regression is a regression spec.
setMethod("isRegressionSpec", "FSSTRegSpec", function(spec) TRUE)

.fsstRegShapePrior <- function(y, control) {
  sy <- stats::sd(y); if (!is.finite(sy) || sy <= 0) sy <- 1
  aSig <- if (!is.null(control$aSig)) control$aSig else 3
  bSig <- if (!is.null(control$bSig)) control$bSig else (aSig - 1) * sy
  aScale <- if (!is.null(control$alphaLogSd)) control$alphaLogSd else 1
  aNu <- if (!is.null(control$aNuShape)) control$aNuShape else 2
  bNu <- if (!is.null(control$bNuRate))  control$bNuRate  else 0.15
  list(aSig = aSig, bSig = bSig, aScale = aScale, aNu = aNu, bNu = bNu,
       shapeDraw = function(nClust) list(
         sigma = 1 / stats::rgamma(nClust, shape = aSig, rate = 1 / bSig),
         alpha = stats::rlnorm(nClust, 0, aScale),
         nu = 2 + stats::rgamma(nClust, shape = aNu, rate = bNu)))
}

#' @describeIn defaultPrior FSST regression prior.
#' @export
setMethod("defaultPrior", "FSSTRegSpec",
  function(spec, data, control = list(), ...)
    .neoRegPrior(data, control, .fsstRegShapePrior))

#' @describeIn buildConstants FSST regression constants.
#' @export
setMethod("buildConstants", "FSSTRegSpec",
  function(spec, prior, n, ...)
    .neoRegREConstants(list(n = n, p = prior$p, X = prior$X, b0 = prior$b0, B0 = prior$B0,
         aSig = prior$aSig, bSig = prior$bSig, aScale = prior$aScale,
         aNu = prior$aNu, bNu = prior$bNu), prior))

#' @describeIn validateParams FSST regression validation.
#' @export
setMethod("validateParams", "FSSTRegSpec",
  function(spec, params, ...) {
    req <- c("b0", "B0", "aSig", "bSig", "aScale", "aNu", "bNu", "p")
    miss <- setdiff(req, names(params))
    if (length(miss))
      stop("FSSTRegSpec prior is missing: ", paste(miss, collapse = ", "),
           call. = FALSE)
    invisible(TRUE)
  })

#' @describeIn buildModelCode FSST regression finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("FSSTRegSpec", "FixedKEngine"),
  function(spec, engine, n, L, re = FALSE, reSlope = FALSE, ...) {
    code <- .neoRegFixedKCode("dFSST_k", c("sigma", "alpha", "nu"),
      paste("      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)",
            "      alphaTilde[j] ~ dlnorm(0, sd = aScale)",
            "      nuTilde[j] ~ T(dgamma(shape = aNu, rate = bNu), 2, )",
            sep = "\n"),
      re = isTRUE(re), reSlope = isTRUE(re) && isTRUE(reSlope))
    list(code = code,
         monitors = c("betaTilde", "sigmaTilde", "alphaTilde", "nuTilde", "weights", "z",
                      if (isTRUE(re)) c("b", "tauRE"),
                      if (isTRUE(re) && isTRUE(reSlope)) c("sRE", "tauSlope")),
         paramNodes = c(beta = "betaTilde"), allocNode = "z")
  })

#' @describeIn linkInv FSST regression identity link.
#' @export
setMethod("linkInv", "FSSTRegSpec", function(spec, eta, prior = NULL, ...) eta)

#' @describeIn componentInits FSST regression initial values.
#' @export
setMethod("componentInits", "FSSTRegSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    y <- as.numeric(data); X <- prior$X; p <- prior$p
    bhat <- tryCatch(as.numeric(solve(crossprod(X) + diag(1e-6, p),
                                      crossprod(X, y))),
                     error = function(e) rep(0, p))
    beta <- matrix(rep(bhat, each = count), nrow = count)
    sig0 <- prior$bSig / (prior$aSig - 1)
    xiInit <- .initClusters(y, min(count, 2L), initMethod)
    if (is.null(xiInit)) xiInit <- rep(1L, length(y))
    list(alloc = xiInit,
         params = list(betaTilde = beta, sigmaTilde = rep(sig0, count),
                       alphaTilde = rep(1, count), nuTilde = rep(5, count)))
  })

.neoRegMethods("FSSTRegSpec", dfsst, c("sigma", "alpha", "nu"), rfsst,
               densName = "dFSST_k", priorLines = paste("      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)",
            "      alphaTilde[j] ~ dlnorm(0, sd = aScale)",
            "      nuTilde[j] ~ T(dgamma(shape = aNu, rate = bNu), 2, )",
            sep = "\n"))

# ===========================================================================
# FOSSEP regression (sigma, alpha, theta -- alpha log-normal, theta Gamma)
# ===========================================================================

#' @rdname FOSSEPRegSpec-class
#' @export
setClass("FOSSEPRegSpec", contains = "DistributionSpec")

#' FOSSEP regression specification
#' @name FOSSEPRegSpec-class
#' @return An object used internally by \code{\link{nimixReg}}.
#' @export
FOSSEPRegSpec <- function() new("FOSSEPRegSpec")

#' @describeIn isRegressionSpec FOSSEP regression is a regression spec.
setMethod("isRegressionSpec", "FOSSEPRegSpec", function(spec) TRUE)

.fossepRegShapePrior <- function(y, control) {
  sy <- stats::sd(y); if (!is.finite(sy) || sy <= 0) sy <- 1
  aSig <- if (!is.null(control$aSig)) control$aSig else 3
  bSig <- if (!is.null(control$bSig)) control$bSig else (aSig - 1) * sy
  aScale <- if (!is.null(control$alphaLogSd)) control$alphaLogSd else 1
  aTheta <- if (!is.null(control$aThetaShape)) control$aThetaShape else 4
  bTheta <- if (!is.null(control$bThetaRate))  control$bThetaRate  else 2
  list(aSig = aSig, bSig = bSig, aScale = aScale, aTheta = aTheta,
       bTheta = bTheta,
       shapeDraw = function(nClust) list(
         sigma = 1 / stats::rgamma(nClust, shape = aSig, rate = 1 / bSig),
         alpha = stats::rlnorm(nClust, 0, aScale),
         theta = stats::rgamma(nClust, shape = aTheta, rate = bTheta)))
}

#' @describeIn defaultPrior FOSSEP regression prior.
#' @export
setMethod("defaultPrior", "FOSSEPRegSpec",
  function(spec, data, control = list(), ...)
    .neoRegPrior(data, control, .fossepRegShapePrior))

#' @describeIn buildConstants FOSSEP regression constants.
#' @export
setMethod("buildConstants", "FOSSEPRegSpec",
  function(spec, prior, n, ...)
    .neoRegREConstants(list(n = n, p = prior$p, X = prior$X, b0 = prior$b0, B0 = prior$B0,
         aSig = prior$aSig, bSig = prior$bSig, aScale = prior$aScale,
         aTheta = prior$aTheta, bTheta = prior$bTheta), prior))

#' @describeIn validateParams FOSSEP regression validation.
#' @export
setMethod("validateParams", "FOSSEPRegSpec",
  function(spec, params, ...) {
    req <- c("b0", "B0", "aSig", "bSig", "aScale", "aTheta", "bTheta", "p")
    miss <- setdiff(req, names(params))
    if (length(miss))
      stop("FOSSEPRegSpec prior is missing: ", paste(miss, collapse = ", "),
           call. = FALSE)
    invisible(TRUE)
  })

#' @describeIn buildModelCode FOSSEP regression finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("FOSSEPRegSpec", "FixedKEngine"),
  function(spec, engine, n, L, re = FALSE, reSlope = FALSE, ...) {
    code <- .neoRegFixedKCode("dFOSSEP_k", c("sigma", "alpha", "theta"),
      paste("      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)",
            "      alphaTilde[j] ~ dlnorm(0, sd = aScale)",
            "      thetaTilde[j] ~ dgamma(shape = aTheta, rate = bTheta)",
            sep = "\n"),
      re = isTRUE(re), reSlope = isTRUE(re) && isTRUE(reSlope))
    list(code = code,
         monitors = c("betaTilde", "sigmaTilde", "alphaTilde", "thetaTilde", "weights", "z",
                      if (isTRUE(re)) c("b", "tauRE"),
                      if (isTRUE(re) && isTRUE(reSlope)) c("sRE", "tauSlope")),
         paramNodes = c(beta = "betaTilde"), allocNode = "z")
  })

#' @describeIn linkInv FOSSEP regression identity link.
#' @export
setMethod("linkInv", "FOSSEPRegSpec", function(spec, eta, prior = NULL, ...) eta)

#' @describeIn componentInits FOSSEP regression initial values.
#' @export
setMethod("componentInits", "FOSSEPRegSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    y <- as.numeric(data); X <- prior$X; p <- prior$p
    bhat <- tryCatch(as.numeric(solve(crossprod(X) + diag(1e-6, p),
                                      crossprod(X, y))),
                     error = function(e) rep(0, p))
    beta <- matrix(rep(bhat, each = count), nrow = count)
    sig0 <- prior$bSig / (prior$aSig - 1)
    xiInit <- .initClusters(y, min(count, 2L), initMethod)
    if (is.null(xiInit)) xiInit <- rep(1L, length(y))
    list(alloc = xiInit,
         params = list(betaTilde = beta, sigmaTilde = rep(sig0, count),
                       alphaTilde = rep(1, count), thetaTilde = rep(2, count)))
  })

.neoRegMethods("FOSSEPRegSpec", dfossep, c("sigma", "alpha", "theta"), rfossep,
               densName = "dFOSSEP_k", priorLines = paste("      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)",
            "      alphaTilde[j] ~ dlnorm(0, sd = aScale)",
            "      thetaTilde[j] ~ dgamma(shape = aTheta, rate = bTheta)",
            sep = "\n"))

# ===========================================================================
# JFST regression (sigma, alpha, theta -- both alpha and theta Gamma, aSh/bSh)
# ===========================================================================

#' @rdname JFSTRegSpec-class
#' @export
setClass("JFSTRegSpec", contains = "DistributionSpec")

#' JFST regression specification
#' @name JFSTRegSpec-class
#' @return An object used internally by \code{\link{nimixReg}}.
#' @export
JFSTRegSpec <- function() new("JFSTRegSpec")

#' @describeIn isRegressionSpec JFST regression is a regression spec.
setMethod("isRegressionSpec", "JFSTRegSpec", function(spec) TRUE)

.jfstRegShapePrior <- function(y, control) {
  sy <- stats::sd(y); if (!is.finite(sy) || sy <= 0) sy <- 1
  aSig <- if (!is.null(control$aSig)) control$aSig else 3
  bSig <- if (!is.null(control$bSig)) control$bSig else (aSig - 1) * sy
  aSh <- if (!is.null(control$aShapeShape)) control$aShapeShape else 4
  bSh <- if (!is.null(control$bShapeRate))  control$bShapeRate  else 1
  list(aSig = aSig, bSig = bSig, aSh = aSh, bSh = bSh,
       shapeDraw = function(nClust) list(
         sigma = 1 / stats::rgamma(nClust, shape = aSig, rate = 1 / bSig),
         alpha = stats::rgamma(nClust, shape = aSh, rate = bSh),
         theta = stats::rgamma(nClust, shape = aSh, rate = bSh)))
}

#' @describeIn defaultPrior JFST regression prior.
#' @export
setMethod("defaultPrior", "JFSTRegSpec",
  function(spec, data, control = list(), ...)
    .neoRegPrior(data, control, .jfstRegShapePrior))

#' @describeIn buildConstants JFST regression constants.
#' @export
setMethod("buildConstants", "JFSTRegSpec",
  function(spec, prior, n, ...)
    .neoRegREConstants(list(n = n, p = prior$p, X = prior$X, b0 = prior$b0, B0 = prior$B0,
         aSig = prior$aSig, bSig = prior$bSig, aSh = prior$aSh, bSh = prior$bSh), prior))

#' @describeIn validateParams JFST regression validation.
#' @export
setMethod("validateParams", "JFSTRegSpec",
  function(spec, params, ...) {
    req <- c("b0", "B0", "aSig", "bSig", "aSh", "bSh", "p")
    miss <- setdiff(req, names(params))
    if (length(miss))
      stop("JFSTRegSpec prior is missing: ", paste(miss, collapse = ", "),
           call. = FALSE)
    invisible(TRUE)
  })

#' @describeIn buildModelCode JFST regression finite mixture (fixed K).
#' @export
setMethod("buildModelCode", signature("JFSTRegSpec", "FixedKEngine"),
  function(spec, engine, n, L, re = FALSE, reSlope = FALSE, ...) {
    code <- .neoRegFixedKCode("dJFST_k", c("sigma", "alpha", "theta"),
      paste("      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)",
            "      alphaTilde[j] ~ dgamma(shape = aSh, rate = bSh)",
            "      thetaTilde[j] ~ dgamma(shape = aSh, rate = bSh)", sep = "\n"),
      re = isTRUE(re), reSlope = isTRUE(re) && isTRUE(reSlope))
    list(code = code,
         monitors = c("betaTilde", "sigmaTilde", "alphaTilde", "thetaTilde", "weights", "z",
                      if (isTRUE(re)) c("b", "tauRE"),
                      if (isTRUE(re) && isTRUE(reSlope)) c("sRE", "tauSlope")),
         paramNodes = c(beta = "betaTilde"), allocNode = "z")
  })

#' @describeIn linkInv JFST regression identity link.
#' @export
setMethod("linkInv", "JFSTRegSpec", function(spec, eta, prior = NULL, ...) eta)

#' @describeIn componentInits JFST regression initial values.
#' @export
setMethod("componentInits", "JFSTRegSpec",
  function(spec, prior, data, count, initMethod = "kmeans", ...) {
    y <- as.numeric(data); X <- prior$X; p <- prior$p
    bhat <- tryCatch(as.numeric(solve(crossprod(X) + diag(1e-6, p),
                                      crossprod(X, y))),
                     error = function(e) rep(0, p))
    beta <- matrix(rep(bhat, each = count), nrow = count)
    sig0 <- prior$bSig / (prior$aSig - 1)
    xiInit <- .initClusters(y, min(count, 2L), initMethod)
    if (is.null(xiInit)) xiInit <- rep(1L, length(y))
    list(alloc = xiInit,
         params = list(betaTilde = beta, sigmaTilde = rep(sig0, count),
                       alphaTilde = rep(1, count), thetaTilde = rep(3, count)))
  })

.neoRegMethods("JFSTRegSpec", djfst, c("sigma", "alpha", "theta"), rjfst,
               densName = "dJFST_k", priorLines = paste("      sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)",
            "      alphaTilde[j] ~ dgamma(shape = aSh, rate = bSh)",
            "      thetaTilde[j] ~ dgamma(shape = aSh, rate = bSh)", sep = "\n"))
