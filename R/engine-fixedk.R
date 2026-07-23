#' @include class-EngineConfig.R
NULL

## ---------------------------------------------------------------------------
## engine-fixedk.R
##
## Finite-mixture engine with a fixed, known K. The model uses a symmetric
## Dirichlet prior on the mixing weights and a categorical allocation per
## observation; there is no Chinese Restaurant Process, so the truncation
## reminder of the DPM does not arise and NIMBLE assigns conjugate samplers to
## the weights and the component parameters.
##
## The component-parameter blocks are the same conjugate priors used by the DPM
## models (Normal-Inverse-Gamma for the univariate and regression kernels,
## Normal-Inverse-Wishart for the multivariate kernel); only the allocation
## mechanism differs. The same deterministic-node patterns are reused: a
## multivariate distribution parameter may not be an expression, and a
## dynamically indexed cluster parameter is resolved through a per-observation
## node before being passed to the kernel.
## ---------------------------------------------------------------------------

#' @describeIn buildModelCode Univariate Gaussian finite-mixture code (fixed K).
#' @export
setMethod("buildModelCode", signature("NormalUvSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        z[i] ~ dcat(weights[1:K])
        y[i] ~ dnorm(muTilde[z[i]], var = s2Tilde[z[i]])
      }
      weights[1:K] ~ ddirch(alphaVec[1:K])
      for (j in 1:K) {
        s2Tilde[j] ~ dinvgamma(shape = nu0, scale = s0)
        muTilde[j] ~ dnorm(mu0, var = s2Tilde[j] / kappa0)
      }
    })
    list(code = code,
         monitors  = c("z", "muTilde", "s2Tilde", "weights"),
         paramNodes = c(mu = "muTilde", s2 = "s2Tilde"),
         allocNode  = "z")
  }
)

#' @describeIn buildModelCode Multivariate Gaussian finite-mixture code (fixed K).
#' @export
setMethod("buildModelCode", signature("NormalMvSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        z[i] ~ dcat(weights[1:K])
        muObs[i, 1:d]       <- muTilde[z[i], 1:d]
        covObs[i, 1:d, 1:d] <- covTilde[z[i], 1:d, 1:d]
        y[i, 1:d] ~ dmnorm(muObs[i, 1:d], cov = covObs[i, 1:d, 1:d])
      }
      weights[1:K] ~ ddirch(alphaVec[1:K])
      for (j in 1:K) {
        covTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
        covMu[j, 1:d, 1:d] <- covTilde[j, 1:d, 1:d] / kappa0
        muTilde[j, 1:d] ~ dmnorm(mu0[1:d], cov = covMu[j, 1:d, 1:d])
      }
    })
    list(code = code,
         monitors  = c("z", "muTilde", "covTilde", "weights"),
         paramNodes = c(mu = "muTilde", cov = "covTilde"),
         allocNode  = "z")
  }
)

#' @describeIn buildModelCode Mixture-of-linear-regressions finite-mixture code
#'   (fixed K).
#' @export
setMethod("buildModelCode", signature("NormalRegSpec", "FixedKEngine"),
  function(spec, engine, n, L, re = FALSE, reSlope = FALSE, ...) {
    code <- if (re && reSlope) nimble::nimbleCode({
      # Random intercept AND random slope. Both use the sum-to-zero
      # parameterisation. Freely parameterised offsets are identified only
      # jointly with the fixed intercept and slope, which produces correlated
      # posterior ridges (measured here: cor(beta1, mean(s)) = -0.929,
      # cor(beta0, mean(b)) = -0.953) and a minimum ESS of 52 against 226 under
      # the constraint -- and the free version had conjugate samplers and still
      # mixed worse, so the geometry, not the sampler, is what matters. The offsets are independent by design: with
      # correlated truth (rho = 0.92) the independent model still recovered
      # both (cor 0.97) and the correlation itself came back empirically.
      for (i in 1:n) {
        z[i] ~ dcat(weights[1:K])
        betaObs[i, 1:p] <- betaTilde[z[i], 1:p]
        mu[i] <- inprod(X[i, 1:p], betaObs[i, 1:p]) + b[grp[i]] +
                 sRE[grp[i]] * xRE[i]
        s2Obs[i] <- s2Tilde[z[i]]
        y[i] ~ dnorm(mu[i], var = s2Obs[i])
      }
      weights[1:K] ~ ddirch(alphaVec[1:K])
      for (j in 1:K) {
        s2Tilde[j] ~ dinvgamma(shape = nu0, scale = s0)
        covBeta[j, 1:p, 1:p] <- s2Tilde[j] * B0[1:p, 1:p]
        betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = covBeta[j, 1:p, 1:p])
      }
      for (g in 1:(G - 1)) {
        bf[g] ~ dnorm(0, sd = tauRE)
        sf[g] ~ dnorm(0, sd = tauSlope)
      }
      b[1:(G - 1)] <- bf[1:(G - 1)]
      b[G] <- -sum(bf[1:(G - 1)])
      sRE[1:(G - 1)] <- sf[1:(G - 1)]
      sRE[G] <- -sum(sf[1:(G - 1)])
      tauRE ~ dunif(tauMin, tauMax)
      tauSlope ~ dunif(tauMinSlope, tauMaxSlope)
    }) else if (!re) nimble::nimbleCode({
      for (i in 1:n) {
        z[i] ~ dcat(weights[1:K])
        betaObs[i, 1:p] <- betaTilde[z[i], 1:p]
        mu[i] <- inprod(X[i, 1:p], betaObs[i, 1:p])
        s2Obs[i] <- s2Tilde[z[i]]
        y[i] ~ dnorm(mu[i], var = s2Obs[i])
      }
      weights[1:K] ~ ddirch(alphaVec[1:K])
      for (j in 1:K) {
        s2Tilde[j] ~ dinvgamma(shape = nu0, scale = s0)
        covBeta[j, 1:p, 1:p] <- s2Tilde[j] * B0[1:p, 1:p]
        betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = covBeta[j, 1:p, 1:p])
      }
    }) else nimble::nimbleCode({
      # Random-intercept variant. Sum-to-zero parameterisation of b: the F4
      # gate measured a pure translation ridge cor(beta0, mean(b)) = -0.979
      # under free b (min ESS 25/2500); the constraint restored min ESS 205
      # with recovery intact, and keeps user-facing traces healthy.
      for (i in 1:n) {
        z[i] ~ dcat(weights[1:K])
        betaObs[i, 1:p] <- betaTilde[z[i], 1:p]
        mu[i] <- inprod(X[i, 1:p], betaObs[i, 1:p]) + b[grp[i]]
        s2Obs[i] <- s2Tilde[z[i]]
        y[i] ~ dnorm(mu[i], var = s2Obs[i])
      }
      weights[1:K] ~ ddirch(alphaVec[1:K])
      for (j in 1:K) {
        s2Tilde[j] ~ dinvgamma(shape = nu0, scale = s0)
        covBeta[j, 1:p, 1:p] <- s2Tilde[j] * B0[1:p, 1:p]
        betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = covBeta[j, 1:p, 1:p])
      }
      for (g in 1:(G - 1)) bf[g] ~ dnorm(0, sd = tauRE)
      b[1:(G - 1)] <- bf[1:(G - 1)]
      b[G] <- -sum(bf[1:(G - 1)])
      tauRE ~ dunif(tauMin, tauMax)
    })
    list(code = code,
         monitors  = c("z", "betaTilde", "s2Tilde", "weights",
                       if (re) c("b", "tauRE"),
                       if (re && reSlope) c("sRE", "tauSlope")),
         paramNodes = c(beta = "betaTilde", s2 = "s2Tilde"),
         allocNode  = "z")
  }
)

#' @describeIn buildModelCode Univariate Student-t finite-mixture code (fixed K).
#' @export
setMethod("buildModelCode", signature("StudentTUvSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        z[i] ~ dcat(weights[1:K])
        y[i] ~ dt(mu = muTilde[z[i]], tau = tauTilde[z[i]], df = df)
      }
      weights[1:K] ~ ddirch(alphaVec[1:K])
      for (j in 1:K) {
        tauTilde[j] ~ dgamma(shape = aTau, rate = bTau)
        muTilde[j]  ~ dnorm(mu0, sd = muSd)
      }
    })
    list(code = code,
         monitors  = c("z", "muTilde", "tauTilde", "weights"),
         paramNodes = c(mu = "muTilde", tau = "tauTilde"),
         allocNode  = "z")
  }
)

#' @describeIn buildModelCode Univariate Normal-Gamma (scale-mixture Student-t)
#'   finite-mixture code (fixed K).
#' @export
setMethod("buildModelCode", signature("NormalGammaUvSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        z[i] ~ dcat(weights[1:K])
        omega[i] ~ dgamma(shape = df / 2, rate = df / 2)
        y[i] ~ dnorm(muTilde[z[i]], var = s2Tilde[z[i]] / omega[i])
      }
      weights[1:K] ~ ddirch(alphaVec[1:K])
      for (j in 1:K) {
        s2Tilde[j] ~ dinvgamma(shape = nu0, scale = s0)
        muTilde[j] ~ dnorm(mu0, var = s2Tilde[j] / kappa0)
      }
    })
    list(code = code,
         monitors  = c("z", "muTilde", "s2Tilde", "weights"),
         paramNodes = c(mu = "muTilde", s2 = "s2Tilde"),
         allocNode  = "z")
  }
)

#' @describeIn runEngine Finite-mixture run with fixed K (Dirichlet weights +
#'   categorical allocation).
setMethod("runEngine", "FixedKEngine",
  function(engine, model, mcmcControl = list(), initMethod = "kmeans",
           seed = 1L, verbose = TRUE, ...) {
    spec  <- model@distSpec
    data  <- model@data
    prior <- model@prior
    K     <- model@Kmax           # for the finite mixture this slot holds K
    n     <- .nObs(data)
    d     <- .dataDimOf(data)

    hasRE <- isTRUE(prior$hasRE)
    hasRES <- isTRUE(prior$hasRESlope)
    mc <- buildModelCode(spec, engine, n = n, L = K, d = d, re = hasRE,
                         reSlope = hasRES)
    dataList <- buildDataList(spec, data)
    initRatio <- .resolveInitRatio(mcmcControl)
    initsFn <- function(s) {
      ci <- .withSeed(s, function() componentInits(spec, prior, data, K,
                      initMethod = initMethod, initRatio = initRatio))
      # Inits must sit inside the data-scaled tau prior: a fixed tauRE = 1
      # would fall outside dunif(tauMin, tauMax) whenever the response scale
      # is small.
      reInits <- if (hasRE)
        c(list(bf = rep(0, prior$reG - 1L), tauRE = prior$tauMax / 10),
          if (hasRES) list(sf = rep(0, prior$reG - 1L),
                           tauSlope = prior$tauMaxSlope / 10))
        else NULL
      if (K == 1L) c(reInits, ci$params)
      else c(list(z = ci$alloc, weights = rep(1 / K, K)), reInits, ci$params)
    }
    baseConst <- buildConstants(spec, prior, n)

    if (K == 1L) {
      # K = 1 is a single component: a one-category Dirichlet/categorical is
      # degenerate and will not build, so drop the allocation layer and pin the
      # per-observation lookups to the single cluster (z == 1).
      mc$code     <- .stripMixtureLayer(mc$code)
      mc$monitors <- setdiff(mc$monitors, c("z", "weights"))
      constants   <- c(baseConst, list(K = K, z = rep(1L, n)))
    } else {
      constants <- c(baseConst,
                     list(K = K, alphaVec = rep(engine@dirichletConc, K)))
    }

    paramDim <- if (!is.null(prior$p)) prior$p else d
    .runNimbleMixture(spec, mc, constants, dataList, initsFn,
                      n = n, count = K, paramDim = paramDim, prior = prior,
                      mcmcControl = mcmcControl, seed = seed, verbose = verbose)
  }
)

#' @describeIn buildModelCode Multivariate Normal-Gamma scale-mixture
#'   finite-mixture code (fixed K).
#' @export
setMethod("buildModelCode", signature("NormalGammaMvSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        z[i] ~ dcat(weights[1:K])
        omega[i] ~ dgamma(shape = df / 2, rate = df / 2)
        muObs[i, 1:d]       <- muTilde[z[i], 1:d]
        covObs[i, 1:d, 1:d] <- covTilde[z[i], 1:d, 1:d] / omega[i]
        y[i, 1:d] ~ dmnorm(muObs[i, 1:d], cov = covObs[i, 1:d, 1:d])
      }
      weights[1:K] ~ ddirch(alphaVec[1:K])
      for (j in 1:K) {
        covTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
        covMu[j, 1:d, 1:d] <- covTilde[j, 1:d, 1:d] / kappa0
        muTilde[j, 1:d] ~ dmnorm(mu0[1:d], cov = covMu[j, 1:d, 1:d])
      }
    })
    list(code = code,
         monitors  = c("z", "muTilde", "covTilde", "weights"),
         paramNodes = c(mu = "muTilde", cov = "covTilde"),
         allocNode  = "z")
  }
)

#' @describeIn buildModelCode Multivariate Student-t finite-mixture code
#'   (fixed K) using the user-defined \code{dmvt_nimix} kernel.
#' @export
setMethod("buildModelCode", signature("StudentTMvSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        z[i] ~ dcat(weights[1:K])
        muObs[i, 1:d]       <- muTilde[z[i], 1:d]
        covObs[i, 1:d, 1:d] <- covTilde[z[i], 1:d, 1:d]
        y[i, 1:d] ~ dmvt_nimix(muObs[i, 1:d], covObs[i, 1:d, 1:d], df)
      }
      weights[1:K] ~ ddirch(alphaVec[1:K])
      for (j in 1:K) {
        covTilde[j, 1:d, 1:d] ~ dinvwish(S = S0[1:d, 1:d], df = df0)
        covMu[j, 1:d, 1:d] <- covTilde[j, 1:d, 1:d] / kappa0
        muTilde[j, 1:d] ~ dmnorm(mu0[1:d], cov = covMu[j, 1:d, 1:d])
      }
    })
    list(code = code,
         monitors  = c("z", "muTilde", "covTilde", "weights"),
         paramNodes = c(mu = "muTilde", cov = "covTilde"),
         allocNode  = "z")
  }
)

#' @describeIn buildModelCode Poisson finite-mixture code (fixed K).
#' @export
setMethod("buildModelCode", signature("PoissonSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        z[i] ~ dcat(weights[1:K])
        y[i] ~ dpois(lambda[z[i]])
      }
      weights[1:K] ~ ddirch(alphaVec[1:K])
      for (j in 1:K) lambda[j] ~ dgamma(shape = a0, rate = b0)
    })
    list(code = code, monitors = c("z", "lambda", "weights"),
         paramNodes = c(lambda = "lambda"), allocNode = "z")
  }
)

#' @describeIn buildModelCode Binomial finite-mixture code (fixed K).
#' @export
setMethod("buildModelCode", signature("BinomialSpec", "FixedKEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      for (i in 1:n) {
        z[i] ~ dcat(weights[1:K])
        y[i] ~ dbin(prob[z[i]], size)
      }
      weights[1:K] ~ ddirch(alphaVec[1:K])
      for (j in 1:K) prob[j] ~ dbeta(a0, b0)
    })
    list(code = code, monitors = c("z", "prob", "weights"),
         paramNodes = c(prob = "prob"), allocNode = "z")
  }
)
