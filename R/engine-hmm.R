## engine-hmm.R -----------------------------------------------------------------
## Hidden-Markov mixture engine: component labels follow a first-order Markov
## chain in the DATA ORDER (which must therefore be time order), and the
## discrete state path is MARGINALISED out of the likelihood by the forward
## algorithm, so the MCMC only ever sees the continuous parameters.
##
## Design decisions (locked by the F2 gate prototype; numbers from this
## container, T = 300, S = 2):
##  * Forward marginalisation, not latent-state sampling: exact vs a pure-R
##    reference (difference 0 compiled and uncompiled), 4000 iterations in
##    ~1 s, min ESS/sec 456 vs 144 for the naive z-sampling model (x3.2 on an
##    easy, well-separated setting; the gap grows as emissions overlap), and
##    it removes T discrete nodes from the model graph.
##  * nimbleEcology::dHMM was evaluated and is categorical-emission only
##    (probObs is an S x O matrix); continuous emissions need this manual
##    forward, which the gate showed compiles exactly with nimix kernels too.
##  * Allocation draws are recovered POST-HOC by forward-filter
##    backward-sampling (FFBS) per retained draw, so every downstream tool
##    (relabel, psm, binderPartition, plots, ppc) works unchanged.
##  * Per-step scaling in the forward pass (not full log-space): exact and
##    stable at T in the hundreds; revisit if T >> 1e4.

#' @include class-DistributionSpec.R
#' @include class-EngineConfig.R
#' @include dist-normal-uv.R
#' @include dist-student-t.R
#' @include dist-msnburr.R
#' @include dist-poisson-binomial.R
NULL

# --- lazy kernel (globalenv; own guard -- see knowledge Sec 9.12) --------------
# The forward-likelihood must be a nimbleFunction registered as a distribution
# BEFORE nimbleModel() sees the code. Like every nimix kernel it is defined
# lazily in globalenv(): namespace-frame definitions fail during NIMBLE
# codegen once the package is installed. This guard is separate from
# .nimixEnsureMSNBurr()'s so neither needs to know about the other's newest
# kernel.
.nimixEnsureHMM <- function() {
  # Guard on the NEWEST kernel (Sec 9.12): adding a kernel means updating this.
  if (exists("dRegimeHMMGMSNB_k", envir = globalenv(), inherits = FALSE))
    return(invisible())

  # The skewed forward kernels call the primitive neo-normal kernels
  # (dMSNBurr_k, dMSNBurr2a_k, dGMSNBurr_k), which are defined lazily by a
  # SEPARATE guard. Ensure them first, or codegen fails to resolve them.
  .nimixEnsureMSNBurr()

  # One forward kernel per emission family. Each marginalises the discrete
  # state path by a scaled forward pass; the paired r-kernel exists only to
  # satisfy registerDistributions (the marginalised model never simulates
  # from it, so a plain rnorm placeholder is fine and avoids codegen clashes
  # with R-level r* functions such as rmsnburr).
  ## --- Normal --------------------------------------------------------------
  assign("dRegimeHMMNorm_k", nimble::nimbleFunction(
    run = function(x = double(1), mu = double(1), s2 = double(1),
                   P = double(2), init = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0)); T <- length(x); S <- length(mu)
      a <- numeric(S); an <- numeric(S); ll <- 0
      for (s in 1:S) a[s] <- init[s] * dnorm(x[1], mu[s], sqrt(s2[s]))
      c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
      if (T >= 2) for (t in 2:T) {
        for (sp in 1:S) { acc <- 0
          for (s in 1:S) acc <- acc + a[s] * P[s, sp]
          an[sp] <- acc * dnorm(x[t], mu[sp], sqrt(s2[sp])) }
        ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
        ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
      if (log) return(ll) else return(exp(ll))
    }), envir = globalenv())
  assign("rRegimeHMMNorm_k", nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(1), s2 = double(1),
                   P = double(2), init = double(1)) {
      returnType(double(1)); Tlen <- dim(P)[1]; out <- numeric(Tlen)
      z <- rcat(1, init)
      for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sqrt(s2[z])); z <- rcat(1, P[z, ]) }
      return(out) }), envir = globalenv())

  ## --- Student-t -----------------------------------------------------------
  assign("dRegimeHMMT_k", nimble::nimbleFunction(
    run = function(x = double(1), mu = double(1), tau = double(1),
                   df = double(0), P = double(2), init = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0)); T <- length(x); S <- length(mu)
      a <- numeric(S); an <- numeric(S); ll <- 0
      for (s in 1:S) a[s] <- init[s] * dt_nonstandard(x[1], df, mu[s], 1/sqrt(tau[s]))
      c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
      if (T >= 2) for (t in 2:T) {
        for (sp in 1:S) { acc <- 0
          for (s in 1:S) acc <- acc + a[s] * P[s, sp]
          an[sp] <- acc * dt_nonstandard(x[t], df, mu[sp], 1/sqrt(tau[sp])) }
        ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
        ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
      if (log) return(ll) else return(exp(ll))
    }), envir = globalenv())
  assign("rRegimeHMMT_k", nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(1), tau = double(1),
                   df = double(0), P = double(2), init = double(1)) {
      returnType(double(1)); Tlen <- dim(P)[1]; out <- numeric(Tlen)
      z <- rcat(1, init)
      for (t in 1:Tlen) { out[t] <- rt_nonstandard(1, df, mu[z], 1/sqrt(tau[z])); z <- rcat(1, P[z, ]) }
      return(out) }), envir = globalenv())

  ## --- Poisson (count regimes; not location-scale) -------------------------
  assign("dRegimeHMMPois_k", nimble::nimbleFunction(
    run = function(x = double(1), lambda = double(1),
                   P = double(2), init = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0)); T <- length(x); S <- length(lambda)
      a <- numeric(S); an <- numeric(S); ll <- 0
      for (s in 1:S) a[s] <- init[s] * dpois(x[1], lambda[s])
      c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
      if (T >= 2) for (t in 2:T) {
        for (sp in 1:S) { acc <- 0
          for (s in 1:S) acc <- acc + a[s] * P[s, sp]
          an[sp] <- acc * dpois(x[t], lambda[sp]) }
        ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
        ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
      if (log) return(ll) else return(exp(ll))
    }), envir = globalenv())
  assign("rRegimeHMMPois_k", nimble::nimbleFunction(
    run = function(n = integer(0), lambda = double(1),
                   P = double(2), init = double(1)) {
      returnType(double(1)); Tlen <- dim(P)[1]; out <- numeric(Tlen)
      z <- rcat(1, init)
      for (t in 1:Tlen) { out[t] <- rpois(1, lambda[z]); z <- rcat(1, P[z, ]) }
      return(out) }), envir = globalenv())

  ## --- MSNBurr (skewed regimes) --------------------------------------------
  assign("dRegimeHMMMSNB_k", nimble::nimbleFunction(
    run = function(x = double(1), mu = double(1), sigma = double(1),
                   alpha = double(1), P = double(2), init = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0)); T <- length(x); S <- length(mu)
      a <- numeric(S); an <- numeric(S); ll <- 0
      for (s in 1:S) a[s] <- init[s] * exp(dMSNBurr_k(x[1], mu[s], sigma[s], alpha[s], 1))
      c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
      if (T >= 2) for (t in 2:T) {
        for (sp in 1:S) { acc <- 0
          for (s in 1:S) acc <- acc + a[s] * P[s, sp]
          an[sp] <- acc * exp(dMSNBurr_k(x[t], mu[sp], sigma[sp], alpha[sp], 1)) }
        ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
        ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
      if (log) return(ll) else return(exp(ll))
    }), envir = globalenv())
  assign("rRegimeHMMMSNB_k", nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(1), sigma = double(1),
                   alpha = double(1), P = double(2), init = double(1)) {
      returnType(double(1)); Tlen <- dim(P)[1]; out <- numeric(Tlen)
      z <- rcat(1, init)
      for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sigma[z]); z <- rcat(1, P[z, ]) }
      return(out) }), envir = globalenv())

  ## --- MSNBurr-IIa (mirror parameterisation, distinct kernel) --------------
  assign("dRegimeHMMMSNB2a_k", nimble::nimbleFunction(
    run = function(x = double(1), mu = double(1), sigma = double(1),
                   alpha = double(1), P = double(2), init = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0)); T <- length(x); S <- length(mu)
      a <- numeric(S); an <- numeric(S); ll <- 0
      for (s in 1:S) a[s] <- init[s] * exp(dMSNBurr2a_k(x[1], mu[s], sigma[s], alpha[s], 1))
      c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
      if (T >= 2) for (t in 2:T) {
        for (sp in 1:S) { acc <- 0
          for (s in 1:S) acc <- acc + a[s] * P[s, sp]
          an[sp] <- acc * exp(dMSNBurr2a_k(x[t], mu[sp], sigma[sp], alpha[sp], 1)) }
        ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
        ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
      if (log) return(ll) else return(exp(ll))
    }), envir = globalenv())
  assign("rRegimeHMMMSNB2a_k", nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(1), sigma = double(1),
                   alpha = double(1), P = double(2), init = double(1)) {
      returnType(double(1)); Tlen <- dim(P)[1]; out <- numeric(Tlen)
      z <- rcat(1, init)
      for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sigma[z]); z <- rcat(1, P[z, ]) }
      return(out) }), envir = globalenv())

  ## --- GMSNBurr (four-parameter skewed regimes) ----------------------------
  assign("dRegimeHMMGMSNB_k", nimble::nimbleFunction(
    run = function(x = double(1), mu = double(1), sigma = double(1),
                   alpha = double(1), theta = double(1),
                   P = double(2), init = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0)); T <- length(x); S <- length(mu)
      a <- numeric(S); an <- numeric(S); ll <- 0
      for (s in 1:S) a[s] <- init[s] * exp(dGMSNBurr_k(x[1], mu[s], sigma[s], alpha[s], theta[s], 1))
      c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
      if (T >= 2) for (t in 2:T) {
        for (sp in 1:S) { acc <- 0
          for (s in 1:S) acc <- acc + a[s] * P[s, sp]
          an[sp] <- acc * exp(dGMSNBurr_k(x[t], mu[sp], sigma[sp], alpha[sp], theta[sp], 1)) }
        ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
        ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
      if (log) return(ll) else return(exp(ll))
    }), envir = globalenv())
  assign("rRegimeHMMGMSNB_k", nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(1), sigma = double(1),
                   alpha = double(1), theta = double(1),
                   P = double(2), init = double(1)) {
      returnType(double(1)); Tlen <- dim(P)[1]; out <- numeric(Tlen)
      z <- rcat(1, init)
      for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sigma[z]); z <- rcat(1, P[z, ]) }
      return(out) }), envir = globalenv())

  nimble::registerDistributions(list(
    dRegimeHMMNorm_k = list(
      BUGSdist = "dRegimeHMMNorm_k(mu, s2, P, init)",
      types = c("value = double(1)", "mu = double(1)", "s2 = double(1)",
                "P = double(2)", "init = double(1)"), discrete = FALSE),
    dRegimeHMMT_k = list(
      BUGSdist = "dRegimeHMMT_k(mu, tau, df, P, init)",
      types = c("value = double(1)", "mu = double(1)", "tau = double(1)",
                "df = double(0)", "P = double(2)", "init = double(1)"), discrete = FALSE),
    dRegimeHMMPois_k = list(
      BUGSdist = "dRegimeHMMPois_k(lambda, P, init)",
      types = c("value = double(1)", "lambda = double(1)",
                "P = double(2)", "init = double(1)"), discrete = TRUE),
    dRegimeHMMMSNB_k = list(
      BUGSdist = "dRegimeHMMMSNB_k(mu, sigma, alpha, P, init)",
      types = c("value = double(1)", "mu = double(1)", "sigma = double(1)",
                "alpha = double(1)", "P = double(2)", "init = double(1)"), discrete = FALSE),
    dRegimeHMMMSNB2a_k = list(
      BUGSdist = "dRegimeHMMMSNB2a_k(mu, sigma, alpha, P, init)",
      types = c("value = double(1)", "mu = double(1)", "sigma = double(1)",
                "alpha = double(1)", "P = double(2)", "init = double(1)"), discrete = FALSE),
    dRegimeHMMGMSNB_k = list(
      BUGSdist = "dRegimeHMMGMSNB_k(mu, sigma, alpha, theta, P, init)",
      types = c("value = double(1)", "mu = double(1)", "sigma = double(1)",
                "alpha = double(1)", "theta = double(1)",
                "P = double(2)", "init = double(1)"), discrete = FALSE)),
    verbose = FALSE)
  invisible()
}

# --- engine class --------------------------------------------------------------

#' Hidden-Markov mixture engine (regime switching in time)
#'
#' Component labels follow a first-order Markov chain in the order of the
#' observations (which must therefore be time order), instead of being
#' independent: \eqn{z_t \mid z_{t-1} \sim \mathrm{Cat}(P[z_{t-1}, ])}. The
#' state path is marginalised out of the likelihood by the forward algorithm,
#' so the MCMC samples only the continuous parameters; allocation draws are
#' then recovered exactly by forward-filter backward-sampling (FFBS) per
#' retained draw, which is what makes every downstream tool (\code{relabel},
#' \code{psm}, \code{binderPartition}, plots) work unchanged.
#'
#' @slot transConc Positive scalar: symmetric Dirichlet concentration of the
#'   prior on each row of the transition matrix. Values above 1 favour
#'   persistence-agnostic rows; the default 1 is uniform on each row simplex.
#' @keywords internal
#' @export
setClass(
  "HMMEngine",
  contains = "EngineConfig",
  representation(transConc = "numeric"),
  prototype = prototype(name = "hmm", transConc = 1)
)

#' Construct a hidden-Markov mixture engine configuration
#'
#' @param transConc Positive scalar concentration for the symmetric Dirichlet
#'   prior on each transition-matrix row. Default \code{1}.
#' @return An \code{\linkS4class{HMMEngine}} object.
#' @examples
#' eng <- HMMEngine()
#' eng
#' @export
HMMEngine <- function(transConc = 1) {
  if (length(transConc) != 1L || !is.finite(transConc) || transConc <= 0)
    stop("transConc must be a positive scalar.", call. = FALSE)
  new("HMMEngine", name = "hmm", transConc = as.numeric(transConc))
}

# --- model code (Normal univariate emission; further families are follow-up
# increments, one at a time, per the gated plan) --------------------------------

#' @describeIn buildModelCode Univariate Gaussian regime-switching HMM code
#'   (states marginalised by the forward algorithm; no allocation node).
#' @export
setMethod("buildModelCode", signature("NormalUvSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMNorm_k(muTilde[1:K], s2Tilde[1:K],
                                P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        s2Tilde[j] ~ dinvgamma(shape = nu0, scale = s0)
        muTilde[j] ~ dnorm(mu0, var = s2Tilde[j] / kappa0)
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code,
         monitors  = c("muTilde", "s2Tilde", "P"),
         paramNodes = c(mu = "muTilde", s2 = "s2Tilde"),
         allocNode  = "zFFBS")   # not in the model; filled post-hoc by FFBS
  }
)

#' @describeIn buildModelCode Univariate Student-t regime-switching HMM code
#'   (fixed df; states marginalised by the forward algorithm).
#' @export
setMethod("buildModelCode", signature("StudentTUvSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMT_k(muTilde[1:K], tauTilde[1:K], df,
                             P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        tauTilde[j] ~ dgamma(shape = aTau, rate = bTau)
        muTilde[j]  ~ dnorm(mu0, sd = muSd)
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code,
         monitors  = c("muTilde", "tauTilde", "P"),
         paramNodes = c(mu = "muTilde", tau = "tauTilde"),
         allocNode  = "zFFBS")
  }
)

#' @describeIn buildModelCode Poisson regime-switching HMM code (count data;
#'   states marginalised by the forward algorithm).
#' @export
setMethod("buildModelCode", signature("PoissonSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMPois_k(lambda[1:K], P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        lambda[j] ~ dgamma(shape = a0, rate = b0)
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code, monitors = c("lambda", "P"),
         paramNodes = c(lambda = "lambda"), allocNode = "zFFBS")
  }
)

#' @describeIn buildModelCode MSNBurr (neo-normal, skewed) regime-switching
#'   HMM code (states marginalised by the forward algorithm).
#' @export
setMethod("buildModelCode", signature("MSNBurrUvSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMMSNB_k(muTilde[1:K], sigmaTilde[1:K],
                                alphaTilde[1:K], P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        muTilde[j]    ~ dnorm(mu0, sd = muSd)
        sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
        alphaTilde[j] ~ dgamma(shape = aA, rate = bA)
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code, monitors = c("muTilde", "sigmaTilde", "alphaTilde", "P"),
         paramNodes = c(mu = "muTilde", sigma = "sigmaTilde",
                        alpha = "alphaTilde"),
         allocNode = "zFFBS")
  }
)

#' @describeIn buildModelCode MSNBurr-IIa regime-switching HMM code.
#' @export
setMethod("buildModelCode", signature("MSNBurr2aUvSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMMSNB2a_k(muTilde[1:K], sigmaTilde[1:K],
                                  alphaTilde[1:K], P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        muTilde[j]    ~ dnorm(mu0, sd = muSd)
        sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
        alphaTilde[j] ~ dgamma(shape = aA, rate = bA)
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code, monitors = c("muTilde", "sigmaTilde", "alphaTilde", "P"),
         paramNodes = c(mu = "muTilde", sigma = "sigmaTilde",
                        alpha = "alphaTilde"),
         allocNode = "zFFBS")
  }
)

#' @describeIn buildModelCode GMSNBurr (four-parameter, skewed) regime-switching
#'   HMM code.
#' @export
setMethod("buildModelCode", signature("GMSNBurrUvSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMGMSNB_k(muTilde[1:K], sigmaTilde[1:K],
                                 alphaTilde[1:K], thetaTilde[1:K],
                                 P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        muTilde[j]    ~ dnorm(mu0, sd = muSd)
        sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
        alphaTilde[j] ~ dgamma(shape = aA, rate = bA)
        thetaTilde[j] ~ dgamma(shape = aT, rate = bT)
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code,
         monitors = c("muTilde", "sigmaTilde", "alphaTilde", "thetaTilde", "P"),
         paramNodes = c(mu = "muTilde", sigma = "sigmaTilde",
                        alpha = "alphaTilde", theta = "thetaTilde"),
         allocNode = "zFFBS")
  }
)

# --- per-family emission densities for FFBS / Viterbi ---------------------------
# T x K matrix of emission densities; `draw = g` for a single retained draw
# (FFBS), `draw = NULL` at the posterior means (Viterbi). New emission
# families implement this one method and a forward kernel; the engine,
# FFBS, and viterbiPath() need no changes.
setGeneric(".hmmEmisDens",
           function(spec, samples, y, K, draw = NULL, prior = NULL)
             standardGeneric(".hmmEmisDens"))

setMethod(".hmmEmisDens", "NormalUvSpec",
  function(spec, samples, y, K, draw = NULL, prior = NULL) {
    muM <- .nodeToArray(samples, "muTilde", K)
    s2M <- .nodeToArray(samples, "s2Tilde", K)
    mu <- if (is.null(draw)) colMeans(muM) else muM[draw, ]
    sg <- sqrt(if (is.null(draw)) colMeans(s2M) else s2M[draw, ])
    vapply(seq_len(K), function(s) stats::dnorm(y, mu[s], sg[s]),
           numeric(length(y)))
  })

setMethod(".hmmEmisDens", "PoissonSpec",
  function(spec, samples, y, K, draw = NULL, prior = NULL) {
    lamM <- .nodeToArray(samples, "lambda", K)
    lam <- if (is.null(draw)) colMeans(lamM) else lamM[draw, ]
    vapply(seq_len(K), function(s) stats::dpois(round(y), lam[s]),
           numeric(length(y)))
  })

setMethod(".hmmEmisDens", "MSNBurrUvSpec",
  function(spec, samples, y, K, draw = NULL, prior = NULL) {
    muM <- .nodeToArray(samples, "muTilde", K)
    sgM <- .nodeToArray(samples, "sigmaTilde", K)
    alM <- .nodeToArray(samples, "alphaTilde", K)
    mu <- if (is.null(draw)) colMeans(muM) else muM[draw, ]
    sg <- if (is.null(draw)) colMeans(sgM) else sgM[draw, ]
    al <- if (is.null(draw)) colMeans(alM) else alM[draw, ]
    vapply(seq_len(K), function(s) dmsnburr(y, mu[s], sg[s], al[s]),
           numeric(length(y)))
  })

setMethod(".hmmEmisDens", "MSNBurr2aUvSpec",
  function(spec, samples, y, K, draw = NULL, prior = NULL) {
    muM <- .nodeToArray(samples, "muTilde", K)
    sgM <- .nodeToArray(samples, "sigmaTilde", K)
    alM <- .nodeToArray(samples, "alphaTilde", K)
    mu <- if (is.null(draw)) colMeans(muM) else muM[draw, ]
    sg <- if (is.null(draw)) colMeans(sgM) else sgM[draw, ]
    al <- if (is.null(draw)) colMeans(alM) else alM[draw, ]
    vapply(seq_len(K), function(s) dmsnburr2a(y, mu[s], sg[s], al[s]),
           numeric(length(y)))
  })

setMethod(".hmmEmisDens", "GMSNBurrUvSpec",
  function(spec, samples, y, K, draw = NULL, prior = NULL) {
    muM <- .nodeToArray(samples, "muTilde", K)
    sgM <- .nodeToArray(samples, "sigmaTilde", K)
    alM <- .nodeToArray(samples, "alphaTilde", K)
    thM <- .nodeToArray(samples, "thetaTilde", K)
    mu <- if (is.null(draw)) colMeans(muM) else muM[draw, ]
    sg <- if (is.null(draw)) colMeans(sgM) else sgM[draw, ]
    al <- if (is.null(draw)) colMeans(alM) else alM[draw, ]
    th <- if (is.null(draw)) colMeans(thM) else thM[draw, ]
    vapply(seq_len(K), function(s) dgmsnburr(y, mu[s], sg[s], al[s], th[s]),
           numeric(length(y)))
  })

setMethod(".hmmEmisDens", "StudentTUvSpec",
  function(spec, samples, y, K, draw = NULL, prior = NULL) {
    muM  <- .nodeToArray(samples, "muTilde", K)
    tauM <- .nodeToArray(samples, "tauTilde", K)
    df <- if (!is.null(prior$df)) prior$df else 4
    mu <- if (is.null(draw)) colMeans(muM) else muM[draw, ]
    sg <- 1 / sqrt(if (is.null(draw)) colMeans(tauM) else tauM[draw, ])
    vapply(seq_len(K), function(s) stats::dt((y - mu[s]) / sg[s], df) / sg[s],
           numeric(length(y)))
  })

#' @describeIn buildModelCode MSNBurr regime-switching HMM code (skewed
#'   regimes; states marginalised by the forward algorithm).
#' @export
setMethod("buildModelCode", signature("MSNBurrUvSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMMSNB_k(muTilde[1:K], sigmaTilde[1:K],
                                alphaTilde[1:K], P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        muTilde[j]    ~ dnorm(mu0, sd = muSd)
        sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
        alphaTilde[j] ~ dgamma(shape = aA, rate = bA)
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code,
         monitors  = c("muTilde", "sigmaTilde", "alphaTilde", "P"),
         paramNodes = c(mu = "muTilde", sigma = "sigmaTilde",
                        alpha = "alphaTilde"),
         allocNode  = "zFFBS")
  }
)

#' @describeIn buildModelCode Poisson regime-switching HMM code (count
#'   regimes; states marginalised by the forward algorithm).
#' @export
setMethod("buildModelCode", signature("PoissonSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMPois_k(lambda[1:K], P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        lambda[j] ~ dgamma(shape = a0, rate = b0)
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code,
         monitors  = c("lambda", "P"),
         paramNodes = c(lambda = "lambda"),
         allocNode  = "zFFBS")
  }
)

setMethod(".hmmEmisDens", "MSNBurrUvSpec",
  function(spec, samples, y, K, draw = NULL, prior = NULL) {
    muM <- .nodeToArray(samples, "muTilde", K)
    sgM <- .nodeToArray(samples, "sigmaTilde", K)
    alM <- .nodeToArray(samples, "alphaTilde", K)
    mu <- if (is.null(draw)) colMeans(muM) else muM[draw, ]
    sg <- if (is.null(draw)) colMeans(sgM) else sgM[draw, ]
    al <- if (is.null(draw)) colMeans(alM) else alM[draw, ]
    vapply(seq_len(K), function(s) dmsnburr(y, mu[s], sg[s], al[s]),
           numeric(length(y)))
  })

setMethod(".hmmEmisDens", "PoissonSpec",
  function(spec, samples, y, K, draw = NULL, prior = NULL) {
    lM <- .nodeToArray(samples, "lambda", K)
    lam <- if (is.null(draw)) colMeans(lM) else lM[draw, ]
    vapply(seq_len(K), function(s) stats::dpois(y, lam[s]),
           numeric(length(y)))
  })

# --- FFBS: exact allocation draws given each retained parameter draw -----------
# Forward pass with per-step scaling, then backward sampling
# z_T ~ alpha_T, z_t | z_{t+1} proportional to alpha_t(s) P[s, z_{t+1}].
# Pure R; cost O(draws x T x S^2), measured well under a second per thousand
# draws at T = 300, S <= 4.
.hmmFFBS <- function(spec, samples, y, K, Parr, init, prior, seed = 1L) {
  ndraw <- nrow(samples); T <- length(y); S <- K
  out <- matrix(1L, ndraw, T)
  # Numerical-stability helper: emission densities can underflow to 0 for a
  # thin-tailed family (e.g. an outlying y under every state's current draw),
  # which would make a/sum(a) = 0/0 = NaN and crash sample.int with
  # "NA in probability vector". When the whole weight vector underflows we
  # fall back to the uniform over states -- the forward pass here mirrors the
  # per-step scaling the compiled kernel already does, so this only bites the
  # post-hoc R decode, and only in the degenerate all-zero case.
  safeNorm <- function(v) {
    s <- sum(v)
    if (!is.finite(s) || s <= 0) rep(1 / length(v), length(v)) else v / s
  }
  set.seed(seed)
  for (g in seq_len(ndraw)) {
    D <- .hmmEmisDens(spec, samples, y, K, draw = g, prior = prior)  # T x S
    D[!is.finite(D)] <- 0
    P <- Parr[g, , ]
    A <- matrix(0, T, S)
    A[1, ] <- safeNorm(init * D[1L, ])
    for (t in 2:T) {
      a <- as.vector(crossprod(P, A[t - 1L, ])) * D[t, ]
      A[t, ] <- safeNorm(a)
    }
    z <- integer(T)
    z[T] <- sample.int(S, 1L, prob = A[T, ])
    for (t in (T - 1L):1L) {
      w <- safeNorm(A[t, ] * P[, z[t + 1L]])
      z[t] <- sample.int(S, 1L, prob = w)
    }
    out[g, ] <- z
  }
  out
}

# --- runEngine ------------------------------------------------------------------

#' @describeIn runEngine Hidden-Markov mixture run (forward-marginalised
#'   likelihood; allocations recovered post-hoc by FFBS).
setMethod("runEngine", "HMMEngine",
  function(engine, model, mcmcControl = list(), initMethod = "kmeans",
           seed = 1L, verbose = TRUE, ...) {
    .nimixEnsureHMM()
    spec  <- model@distSpec
    data  <- model@data
    prior <- model@prior
    K     <- model@Kmax
    n     <- .nObs(data)
    d     <- .dataDimOf(data)

    okSpec <- is(spec, "NormalUvSpec") || is(spec, "StudentTUvSpec") ||
              is(spec, "MSNBurrUvSpec") || is(spec, "PoissonSpec") ||
              is(spec, "MSNBurr2aUvSpec") || is(spec, "GMSNBurrUvSpec")
    if (!okSpec)
      stop("method = 'hmm' currently supports distribution = \"normal\", ",
           "\"student-t\", \"poisson\", \"msnburr\", \"msnburr2a\", or ",
           "\"gmsnburr\" (univariate); further emission families follow ",
           "the gated plan.", call. = FALSE)
    if (K < 2L)
      stop("method = 'hmm' needs K >= 2 states.", call. = FALSE)

    mc <- buildModelCode(spec, engine, n = n, L = K, d = d)
    constants <- c(buildConstants(spec, prior, n),
                   list(K = K, alphaP = rep(engine@transConc, K),
                        init = rep(1 / K, K)))
    dataList <- buildDataList(spec, data)

    initRatio <- .resolveInitRatio(mcmcControl)
    Pinit <- matrix((1 - 0.8) / (K - 1), K, K); diag(Pinit) <- 0.8
    initsFn <- function(s) {
      ci <- .withSeed(s, function() componentInits(spec, prior, data, K,
                      initMethod = initMethod, initRatio = initRatio))
      c(list(P = Pinit), ci$params)     # no allocation node to initialise
    }

    out <- .runNimbleMixture(spec, mc, constants, dataList, initsFn,
                      n = n, count = K, paramDim = d, prior = prior,
                      mcmcControl = mcmcControl, seed = seed, verbose = verbose,
                      cacheExtra = list(engine = "HMMEngine"))

    # The runner found no allocation columns (marginalised model), so its
    # alloc/Kposterior/entropy are placeholders. Recover exact allocation
    # draws by FFBS, then rebuild the derived quantities and the per-chain
    # diagnostics from them (chainId marks the pooled rows per chain).
    S <- out$mcmcSamples
    y <- as.numeric(data)
    Parr <- array(0, dim = c(nrow(S), K, K))
    for (r in seq_len(K)) for (cc in seq_len(K))
      Parr[, r, cc] <- S[, sprintf("P[%d, %d]", r, cc)]
    alloc <- .hmmFFBS(spec, S, y, K, Parr, init = rep(1 / K, K),
                      prior = prior, seed = seed)

    out$clusterAllocation <- alloc
    out$Kposterior <- .rowDistinct(alloc, K)
    chainId <- out$diagnostics$chainId
    chains  <- split(seq_len(nrow(alloc)), chainId)
    chainK  <- lapply(chains, function(ix) .rowDistinct(alloc[ix, , drop = FALSE], K))
    chainEnt <- lapply(chains, function(ix) .allocEntropy(alloc[ix, , drop = FALSE], K))
    diag2 <- .multiChainDiag(chainK, NULL, NULL, chainEnt)
    diag2$chainId <- chainId
    out$diagnostics <- diag2
    out
  }
)

# --- Viterbi --------------------------------------------------------------------

#' Most probable state path of a hidden-Markov mixture fit
#'
#' Computes the Viterbi (maximum a posteriori joint) state sequence of a
#' \code{method = "hmm"} fit, by default at the posterior means of the state
#' parameters and transition matrix.
#'
#' Note the difference from \code{binderPartition()}: Viterbi gives the single
#' jointly most probable \emph{path} under the Markov prior, while the Binder
#' partition summarises marginal co-clustering across all FFBS draws. They
#' usually agree on well-separated regimes and differ exactly where the state
#' is genuinely uncertain.
#'
#' @param fit A \code{FitResult} from \code{nimixClust(..., method = "hmm")}.
#' @return Integer vector of length \eqn{n}: the decoded state per time point.
#' @export
viterbiPath <- function(fit) {
  stopifnot(is(fit, "FitResult"))
  if (!identical(fit@engineUsed, "hmm"))
    stop("viterbiPath() is for method = 'hmm' fits.", call. = FALSE)
  S <- fit@mcmcSamples
  K <- fit@Kmax
  y <- as.numeric(fit@data)
  P <- matrix(0, K, K)
  for (r in seq_len(K)) for (cc in seq_len(K))
    P[r, cc] <- mean(S[, sprintf("P[%d, %d]", r, cc)])
  Tn <- length(y)
  D <- .hmmEmisDens(fit@distSpec, S, y, K, draw = NULL, prior = fit@prior)
  # Floor underflowed/zero densities so log() does not produce -Inf that would
  # make an entire time point unreachable in the Viterbi recursion. The floor
  # is far below any real density and only affects genuine underflow.
  D[!is.finite(D) | D <= 0] <- .Machine$double.xmin
  ld <- log(D)
  V <- matrix(-Inf, Tn, K); B <- matrix(0L, Tn, K)
  V[1, ] <- log(rep(1 / K, K)) + ld[1, ]
  lP <- log(P)
  for (t in 2:Tn) for (s in seq_len(K)) {
    cand <- V[t - 1L, ] + lP[, s]
    B[t, s] <- which.max(cand)
    V[t, s] <- max(cand) + ld[t, s]
  }
  z <- integer(Tn)
  z[Tn] <- which.max(V[Tn, ])
  for (t in (Tn - 1L):1L) z[t] <- B[t + 1L, z[t + 1L]]
  z
}
