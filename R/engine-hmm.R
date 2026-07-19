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
# BUGSdist registration entries for the neo-normal regression HMM kernels,
# generated from each family's kernel base name and shape list -- the same
# facts .makeNeoHMMRegKernel() uses. Keeps the registration in step with the
# kernels without hand-maintaining nine near-identical list entries.
.neoHMMRegBUGS <- function() {
  specs <- list(
    list("MSNBReg",   c("sigma", "alpha")),
    list("MSNB2aReg", c("sigma", "alpha")),
    list("FSSNReg",   c("sigma", "alpha")),
    list("SEPReg",    c("sigma", "nu")),
    list("LEPReg",    c("sigma", "nu")),
    list("GMSNBReg",  c("sigma", "alpha", "theta")),
    list("FSSTReg",   c("sigma", "alpha", "nu")),
    list("FOSSEPReg", c("sigma", "alpha", "theta")),
    list("JFSTReg",   c("sigma", "alpha", "theta")))
  out <- list()
  for (sp in specs) {
    base <- sp[[1]]; sh <- sp[[2]]
    nm <- paste0("dRegimeHMM", base, "_k")
    shArgs <- paste(sh, collapse = ", ")
    shTypes <- paste(sprintf('"%s = double(1)"', sh), collapse = ", ")
    out[[nm]] <- list(
      BUGSdist = sprintf("%s(X, beta, %s, P, init)", nm, shArgs),
      types = eval(parse(text = sprintf(
        'c("value = double(1)", "X = double(2)", "beta = double(2)", %s, "P = double(2)", "init = double(1)")',
        shTypes))),
      discrete = FALSE)
  }
  out
}

.nimixEnsureHMM <- function() {
  # Guard on the NEWEST kernel (Sec 9.12): adding a kernel means updating this.
  if (exists("dRegimeHMMStudentTReg_k", envir = globalenv(), inherits = FALSE))
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

  ## --- FSSN (Ferreira-Steel skew normal) ------------------------------------
  assign("dRegimeHMMFSSN_k", nimble::nimbleFunction(
    run = function(x = double(1), mu = double(1), sigma = double(1),
                   alpha = double(1), P = double(2), init = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0)); T <- length(x); S <- length(mu)
      a <- numeric(S); an <- numeric(S); ll <- 0
      for (s in 1:S) a[s] <- init[s] * exp(dFSSN_k(x[1], mu[s], sigma[s], alpha[s], 1))
      c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
      if (T >= 2) for (t in 2:T) {
        for (sp in 1:S) { acc <- 0
          for (s in 1:S) acc <- acc + a[s] * P[s, sp]
          an[sp] <- acc * exp(dFSSN_k(x[t], mu[sp], sigma[sp], alpha[sp], 1)) }
        ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
        ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
      if (log) return(ll) else return(exp(ll))
    }), envir = globalenv())
  assign("rRegimeHMMFSSN_k", nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(1), sigma = double(1),
                   alpha = double(1), P = double(2), init = double(1)) {
      returnType(double(1)); Tlen <- dim(P)[1]; out <- numeric(Tlen)
      z <- rcat(1, init)
      for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sigma[z]); z <- rcat(1, P[z, ]) }
      return(out) }), envir = globalenv())

  ## --- FSST (Ferreira-Steel skew t: heavy-tailed AND skewed) ----------------
  assign("dRegimeHMMFSST_k", nimble::nimbleFunction(
    run = function(x = double(1), mu = double(1), sigma = double(1),
                   alpha = double(1), nu = double(1),
                   P = double(2), init = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0)); T <- length(x); S <- length(mu)
      a <- numeric(S); an <- numeric(S); ll <- 0
      for (s in 1:S) a[s] <- init[s] * exp(dFSST_k(x[1], mu[s], sigma[s], alpha[s], nu[s], 1))
      c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
      if (T >= 2) for (t in 2:T) {
        for (sp in 1:S) { acc <- 0
          for (s in 1:S) acc <- acc + a[s] * P[s, sp]
          an[sp] <- acc * exp(dFSST_k(x[t], mu[sp], sigma[sp], alpha[sp], nu[sp], 1)) }
        ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
        ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
      if (log) return(ll) else return(exp(ll))
    }), envir = globalenv())
  assign("rRegimeHMMFSST_k", nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(1), sigma = double(1),
                   alpha = double(1), nu = double(1),
                   P = double(2), init = double(1)) {
      returnType(double(1)); Tlen <- dim(P)[1]; out <- numeric(Tlen)
      z <- rcat(1, init)
      for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sigma[z]); z <- rcat(1, P[z, ]) }
      return(out) }), envir = globalenv())

  ## --- SEP ------------------------------------------------------------
  assign("dRegimeHMMSEP_k", nimble::nimbleFunction(
    run = function(x = double(1), mu = double(1), sigma = double(1),
                   nu = double(1), P = double(2), init = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0)); T <- length(x); S <- length(mu)
      a <- numeric(S); an <- numeric(S); ll <- 0
      for (s in 1:S) a[s] <- init[s] * exp(dSEP_k(x[1], mu[s], sigma[s], nu[s], 1))
      c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
      if (T >= 2) for (t in 2:T) {
        for (sp in 1:S) { acc <- 0
          for (s in 1:S) acc <- acc + a[s] * P[s, sp]
          an[sp] <- acc * exp(dSEP_k(x[t], mu[sp], sigma[sp], nu[sp], 1)) }
        ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
        ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
      if (log) return(ll) else return(exp(ll))
    }), envir = globalenv())
  assign("rRegimeHMMSEP_k", nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(1), sigma = double(1),
                   nu = double(1), P = double(2), init = double(1)) {
      returnType(double(1)); Tlen <- dim(P)[1]; out <- numeric(Tlen)
      z <- rcat(1, init)
      for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sigma[z]); z <- rcat(1, P[z, ]) }
      return(out) }), envir = globalenv())

  ## --- LEP ------------------------------------------------------------
  assign("dRegimeHMMLEP_k", nimble::nimbleFunction(
    run = function(x = double(1), mu = double(1), sigma = double(1),
                   nu = double(1), P = double(2), init = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0)); T <- length(x); S <- length(mu)
      a <- numeric(S); an <- numeric(S); ll <- 0
      for (s in 1:S) a[s] <- init[s] * exp(dLEP_k(x[1], mu[s], sigma[s], nu[s], 1))
      c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
      if (T >= 2) for (t in 2:T) {
        for (sp in 1:S) { acc <- 0
          for (s in 1:S) acc <- acc + a[s] * P[s, sp]
          an[sp] <- acc * exp(dLEP_k(x[t], mu[sp], sigma[sp], nu[sp], 1)) }
        ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
        ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
      if (log) return(ll) else return(exp(ll))
    }), envir = globalenv())
  assign("rRegimeHMMLEP_k", nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(1), sigma = double(1),
                   nu = double(1), P = double(2), init = double(1)) {
      returnType(double(1)); Tlen <- dim(P)[1]; out <- numeric(Tlen)
      z <- rcat(1, init)
      for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sigma[z]); z <- rcat(1, P[z, ]) }
      return(out) }), envir = globalenv())

  ## --- FOSSEP ------------------------------------------------------------
  assign("dRegimeHMMFOSSEP_k", nimble::nimbleFunction(
    run = function(x = double(1), mu = double(1), sigma = double(1),
                   alpha = double(1), theta = double(1),
                   P = double(2), init = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0)); T <- length(x); S <- length(mu)
      a <- numeric(S); an <- numeric(S); ll <- 0
      for (s in 1:S) a[s] <- init[s] * exp(dFOSSEP_k(x[1], mu[s], sigma[s], alpha[s], theta[s], 1))
      c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
      if (T >= 2) for (t in 2:T) {
        for (sp in 1:S) { acc <- 0
          for (s in 1:S) acc <- acc + a[s] * P[s, sp]
          an[sp] <- acc * exp(dFOSSEP_k(x[t], mu[sp], sigma[sp], alpha[sp], theta[sp], 1)) }
        ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
        ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
      if (log) return(ll) else return(exp(ll))
    }), envir = globalenv())
  assign("rRegimeHMMFOSSEP_k", nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(1), sigma = double(1),
                   alpha = double(1), theta = double(1),
                   P = double(2), init = double(1)) {
      returnType(double(1)); Tlen <- dim(P)[1]; out <- numeric(Tlen)
      z <- rcat(1, init)
      for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sigma[z]); z <- rcat(1, P[z, ]) }
      return(out) }), envir = globalenv())

  ## --- JFST ------------------------------------------------------------
  assign("dRegimeHMMJFST_k", nimble::nimbleFunction(
    run = function(x = double(1), mu = double(1), sigma = double(1),
                   alpha = double(1), theta = double(1),
                   P = double(2), init = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0)); T <- length(x); S <- length(mu)
      a <- numeric(S); an <- numeric(S); ll <- 0
      for (s in 1:S) a[s] <- init[s] * exp(dJFST_k(x[1], mu[s], sigma[s], alpha[s], theta[s], 1))
      c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
      if (T >= 2) for (t in 2:T) {
        for (sp in 1:S) { acc <- 0
          for (s in 1:S) acc <- acc + a[s] * P[s, sp]
          an[sp] <- acc * exp(dJFST_k(x[t], mu[sp], sigma[sp], alpha[sp], theta[sp], 1)) }
        ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
        ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
      if (log) return(ll) else return(exp(ll))
    }), envir = globalenv())
  assign("rRegimeHMMJFST_k", nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(1), sigma = double(1),
                   alpha = double(1), theta = double(1),
                   P = double(2), init = double(1)) {
      returnType(double(1)); Tlen <- dim(P)[1]; out <- numeric(Tlen)
      z <- rcat(1, init)
      for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sigma[z]); z <- rcat(1, P[z, ]) }
      return(out) }), envir = globalenv())

  ## --- Binomial (regime-switching proportions; size known) ------------------
  assign("dRegimeHMMBinom_k", nimble::nimbleFunction(
    run = function(x = double(1), prob = double(1), size = double(0),
                   P = double(2), init = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0)); T <- length(x); S <- length(prob)
      a <- numeric(S); an <- numeric(S); ll <- 0
      for (s in 1:S) a[s] <- init[s] * exp(dbinom(x[1], size, prob[s], 1))
      c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
      if (T >= 2) for (t in 2:T) {
        for (sp in 1:S) { acc <- 0
          for (s in 1:S) acc <- acc + a[s] * P[s, sp]
          an[sp] <- acc * exp(dbinom(x[t], size, prob[sp], 1)) }
        ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
        ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
      if (log) return(ll) else return(exp(ll))
    }), envir = globalenv())
  assign("rRegimeHMMBinom_k", nimble::nimbleFunction(
    run = function(n = integer(0), prob = double(1), size = double(0),
                   P = double(2), init = double(1)) {
      returnType(double(1)); Tlen <- dim(P)[1]; out <- numeric(Tlen)
      z <- rcat(1, init)
      for (t in 1:Tlen) { out[t] <- rbinom(1, size, prob[z]); z <- rcat(1, P[z, ]) }
      return(out) }), envir = globalenv())

  ## --- Student-t regression (Markov-switching heavy-tail regression) --------
  ## Same shape as the Gaussian regression kernel, but the emission is a
  ## non-standard Student-t: location X beta_s, scale sqrt(s2_s), df fixed.
  ## Covers both the direct Student-t and Normal-Gamma reg specs (same t
  ## marginal); the guard refuses the Gaussian inheritance so this OWN kernel
  ## is always the one that runs.
  assign("dRegimeHMMStudentTReg_k", nimble::nimbleFunction(
    run = function(x = double(1), X = double(2), beta = double(2),
                   s2 = double(1), df = double(0), P = double(2),
                   init = double(1), log = integer(0, default = 0)) {
      returnType(double(0))
      T <- length(x); S <- length(s2); p <- dim(X)[2]
      a <- numeric(S); an <- numeric(S); ll <- 0
      for (s in 1:S)
        a[s] <- init[s] * dt_nonstandard(x[1], df,
                   sum(X[1, 1:p] * beta[s, 1:p]), sqrt(s2[s]))
      c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
      if (T >= 2) for (t in 2:T) {
        for (sp in 1:S) { acc <- 0
          for (s in 1:S) acc <- acc + a[s] * P[s, sp]
          an[sp] <- acc * dt_nonstandard(x[t], df,
                      sum(X[t, 1:p] * beta[sp, 1:p]), sqrt(s2[sp])) }
        ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
        ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
      if (log) return(ll) else return(exp(ll))
    }), envir = globalenv())
  assign("rRegimeHMMStudentTReg_k", nimble::nimbleFunction(
    run = function(n = integer(0), X = double(2), beta = double(2),
                   s2 = double(1), df = double(0), P = double(2),
                   init = double(1)) {
      returnType(double(1)); Tlen <- dim(X)[1]; out <- numeric(Tlen)
      z <- rcat(1, init)
      for (t in 1:Tlen) { out[t] <- rnorm(1, 0, 1); z <- rcat(1, P[z, ]) }
      return(out) }), envir = globalenv())
  ## --- Neo-normal regression HMM kernels (generated) ------------------------
  ## Each family's forward kernel is the same FFBS marginalisation, differing
  ## only in the primitive density name and the list of per-state shape
  ## parameters. .makeNeoHMMRegKernel() builds the d- and r-kernels from those
  ## two facts, so a family needs no hand-written nimbleFunction. The shape
  ## arguments are woven into the run() signature and the density call; the
  ## r-kernel is a placeholder (rnorm) because the marginalised model never
  ## simulates from the emission (cf. the other -Reg kernels).
  .makeNeoHMMRegKernel <- function(kernelBase, densName, shapeNames) {
    shapeSig <- paste(sprintf("%s = double(1)", shapeNames), collapse = ", ")
    shapeArg1 <- paste(sprintf("%s[s]", shapeNames), collapse = ", ")
    shapeArgT <- paste(sprintf("%s[sp]", shapeNames), collapse = ", ")
    firstShape <- shapeNames[1]
    dSrc <- sprintf('
      nimble::nimbleFunction(run = function(x = double(1), X = double(2),
          beta = double(2), %s, P = double(2), init = double(1),
          log = integer(0, default = 0)) {
        returnType(double(0))
        T <- length(x); S <- length(%s); p <- dim(X)[2]
        a <- numeric(S); an <- numeric(S); ll <- 0
        for (s in 1:S)
          a[s] <- init[s] * exp(%s(x[1],
                    sum(X[1, 1:p] * beta[s, 1:p]), %s, 1))
        c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
        ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
        if (T >= 2) for (t in 2:T) {
          for (sp in 1:S) { acc <- 0
            for (s in 1:S) acc <- acc + a[s] * P[s, sp]
            an[sp] <- acc * exp(%s(x[t],
                        sum(X[t, 1:p] * beta[sp, 1:p]), %s, 1)) }
          ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
          ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
        if (log) return(ll) else return(exp(ll))
      })',
      shapeSig, firstShape, densName, shapeArg1, densName, shapeArgT)
    assign(paste0("dRegimeHMM", kernelBase, "_k"),
           eval(parse(text = dSrc)), envir = globalenv())
    rSrc <- sprintf('
      nimble::nimbleFunction(run = function(n = integer(0), X = double(2),
          beta = double(2), %s, P = double(2), init = double(1)) {
        returnType(double(1)); Tlen <- dim(X)[1]; out <- numeric(Tlen)
        z <- rcat(1, init)
        for (t in 1:Tlen) { out[t] <- rnorm(1, 0, 1); z <- rcat(1, P[z, ]) }
        return(out) })', shapeSig)
    assign(paste0("rRegimeHMM", kernelBase, "_k"),
           eval(parse(text = rSrc)), envir = globalenv())
  }

  .makeNeoHMMRegKernel("MSNBReg", "dMSNBurr_k", c("sigma", "alpha"))
  .makeNeoHMMRegKernel("MSNB2aReg", "dMSNBurr2a_k", c("sigma", "alpha"))
  .makeNeoHMMRegKernel("FSSNReg", "dFSSN_k", c("sigma", "alpha"))
  .makeNeoHMMRegKernel("SEPReg", "dSEP_k", c("sigma", "nu"))
  .makeNeoHMMRegKernel("LEPReg", "dLEP_k", c("sigma", "nu"))
  .makeNeoHMMRegKernel("GMSNBReg", "dGMSNBurr_k", c("sigma", "alpha", "theta"))
  .makeNeoHMMRegKernel("FSSTReg", "dFSST_k", c("sigma", "alpha", "nu"))
  .makeNeoHMMRegKernel("FOSSEPReg", "dFOSSEP_k", c("sigma", "alpha", "theta"))
  .makeNeoHMMRegKernel("JFSTReg", "dJFST_k", c("sigma", "alpha", "theta"))
  ## --- Binomial regression (Markov-switching proportion regression) ---------
  ## Logit link: p = plogis(X beta_s), known size. Same shape as the Poisson
  ## regression kernel with dbinom in place of dpois.
  assign("dRegimeHMMBinomReg_k", nimble::nimbleFunction(
    run = function(x = double(1), X = double(2), beta = double(2),
                   size = double(0), P = double(2), init = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      T <- length(x); S <- dim(beta)[1]; p <- dim(X)[2]
      a <- numeric(S); an <- numeric(S); ll <- 0
      for (s in 1:S)
        a[s] <- init[s] * dbinom(x[1], size,
                                 ilogit(sum(X[1, 1:p] * beta[s, 1:p])))
      c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
      if (T >= 2) for (t in 2:T) {
        for (sp in 1:S) { acc <- 0
          for (s in 1:S) acc <- acc + a[s] * P[s, sp]
          an[sp] <- acc * dbinom(x[t], size,
                                 ilogit(sum(X[t, 1:p] * beta[sp, 1:p]))) }
        ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
        ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
      if (log) return(ll) else return(exp(ll))
    }), envir = globalenv())
  assign("rRegimeHMMBinomReg_k", nimble::nimbleFunction(
    run = function(n = integer(0), X = double(2), beta = double(2),
                   size = double(0), P = double(2), init = double(1)) {
      returnType(double(1)); Tlen <- dim(X)[1]; out <- numeric(Tlen)
      z <- rcat(1, init)
      for (t in 1:Tlen) { out[t] <- rbinom(1, size, 0.5); z <- rcat(1, P[z, ]) }
      return(out) }), envir = globalenv())

  ## --- Poisson regression (Markov-switching count regression) ---------------
  ## Same shape as the Gaussian regression kernel, but the emission is a
  ## Poisson with a log link: mu = exp(X beta_s). No error-variance parameter.
  assign("dRegimeHMMPoisReg_k", nimble::nimbleFunction(
    run = function(x = double(1), X = double(2), beta = double(2),
                   P = double(2), init = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      T <- length(x); S <- dim(beta)[1]; p <- dim(X)[2]
      a <- numeric(S); an <- numeric(S); ll <- 0
      for (s in 1:S)
        a[s] <- init[s] * dpois(x[1], exp(sum(X[1, 1:p] * beta[s, 1:p])))
      c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
      if (T >= 2) for (t in 2:T) {
        for (sp in 1:S) { acc <- 0
          for (s in 1:S) acc <- acc + a[s] * P[s, sp]
          an[sp] <- acc * dpois(x[t], exp(sum(X[t, 1:p] * beta[sp, 1:p]))) }
        ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
        ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
      if (log) return(ll) else return(exp(ll))
    }), envir = globalenv())
  assign("rRegimeHMMPoisReg_k", nimble::nimbleFunction(
    run = function(n = integer(0), X = double(2), beta = double(2),
                   P = double(2), init = double(1)) {
      returnType(double(1)); Tlen <- dim(X)[1]; out <- numeric(Tlen)
      z <- rcat(1, init)
      for (t in 1:Tlen) { out[t] <- rpois(1, 1); z <- rcat(1, P[z, ]) }
      return(out) }), envir = globalenv())

  ## --- Normal regression (Markov-switching regression, Hamilton 1989) -------
  ## The emission mean is X[t, ] %*% beta[s, ], so unlike every other kernel
  ## here the density at t depends on t through the design matrix, not only
  ## through y[t]. Everything else -- the scaled forward pass, the state
  ## marginalisation -- is unchanged.
  assign("dRegimeHMMNormReg_k", nimble::nimbleFunction(
    run = function(x = double(1), X = double(2), beta = double(2),
                   s2 = double(1), P = double(2), init = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      T <- length(x); S <- length(s2); p <- dim(X)[2]
      a <- numeric(S); an <- numeric(S); ll <- 0
      for (s in 1:S)
        a[s] <- init[s] * dnorm(x[1], sum(X[1, 1:p] * beta[s, 1:p]),
                                sqrt(s2[s]))
      c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
      if (T >= 2) for (t in 2:T) {
        for (sp in 1:S) { acc <- 0
          for (s in 1:S) acc <- acc + a[s] * P[s, sp]
          an[sp] <- acc * dnorm(x[t], sum(X[t, 1:p] * beta[sp, 1:p]),
                                sqrt(s2[sp])) }
        ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
        ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
      if (log) return(ll) else return(exp(ll))
    }), envir = globalenv())
  assign("rRegimeHMMNormReg_k", nimble::nimbleFunction(
    run = function(n = integer(0), X = double(2), beta = double(2),
                   s2 = double(1), P = double(2), init = double(1)) {
      returnType(double(1)); Tlen <- dim(X)[1]; out <- numeric(Tlen)
      z <- rcat(1, init)
      for (t in 1:Tlen) { out[t] <- rnorm(1, 0, 1); z <- rcat(1, P[z, ]) }
      return(out) }), envir = globalenv())

  hmmDists <- list(
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
                "P = double(2)", "init = double(1)"), discrete = FALSE),
    dRegimeHMMFSSN_k = list(
      BUGSdist = "dRegimeHMMFSSN_k(mu, sigma, alpha, P, init)",
      types = c("value = double(1)", "mu = double(1)", "sigma = double(1)",
                "alpha = double(1)", "P = double(2)", "init = double(1)"),
      discrete = FALSE),
    dRegimeHMMFSST_k = list(
      BUGSdist = "dRegimeHMMFSST_k(mu, sigma, alpha, nu, P, init)",
      types = c("value = double(1)", "mu = double(1)", "sigma = double(1)",
                "alpha = double(1)", "nu = double(1)",
                "P = double(2)", "init = double(1)"), discrete = FALSE),
    dRegimeHMMSEP_k = list(
      BUGSdist = "dRegimeHMMSEP_k(mu, sigma, nu, P, init)",
      types = c("value = double(1)", "mu = double(1)", "sigma = double(1)",
                "nu = double(1)", "P = double(2)", "init = double(1)"),
      discrete = FALSE),
    dRegimeHMMLEP_k = list(
      BUGSdist = "dRegimeHMMLEP_k(mu, sigma, nu, P, init)",
      types = c("value = double(1)", "mu = double(1)", "sigma = double(1)",
                "nu = double(1)", "P = double(2)", "init = double(1)"),
      discrete = FALSE),
    dRegimeHMMFOSSEP_k = list(
      BUGSdist = "dRegimeHMMFOSSEP_k(mu, sigma, alpha, theta, P, init)",
      types = c("value = double(1)", "mu = double(1)", "sigma = double(1)",
                "alpha = double(1)", "theta = double(1)",
                "P = double(2)", "init = double(1)"), discrete = FALSE),
    dRegimeHMMJFST_k = list(
      BUGSdist = "dRegimeHMMJFST_k(mu, sigma, alpha, theta, P, init)",
      types = c("value = double(1)", "mu = double(1)", "sigma = double(1)",
                "alpha = double(1)", "theta = double(1)",
                "P = double(2)", "init = double(1)"), discrete = FALSE),
    dRegimeHMMBinom_k = list(
      BUGSdist = "dRegimeHMMBinom_k(prob, size, P, init)",
      types = c("value = double(1)", "prob = double(1)", "size = double(0)",
                "P = double(2)", "init = double(1)"), discrete = TRUE),
    dRegimeHMMNormReg_k = list(
      BUGSdist = "dRegimeHMMNormReg_k(X, beta, s2, P, init)",
      types = c("value = double(1)", "X = double(2)", "beta = double(2)",
                "s2 = double(1)", "P = double(2)", "init = double(1)"),
      discrete = FALSE),
    dRegimeHMMPoisReg_k = list(
      BUGSdist = "dRegimeHMMPoisReg_k(X, beta, P, init)",
      types = c("value = double(1)", "X = double(2)", "beta = double(2)",
                "P = double(2)", "init = double(1)"), discrete = TRUE),
    dRegimeHMMBinomReg_k = list(
      BUGSdist = "dRegimeHMMBinomReg_k(X, beta, size, P, init)",
      types = c("value = double(1)", "X = double(2)", "beta = double(2)",
                "size = double(0)", "P = double(2)", "init = double(1)"),
      discrete = TRUE),
    dRegimeHMMStudentTReg_k = list(
      BUGSdist = "dRegimeHMMStudentTReg_k(X, beta, s2, df, P, init)",
      types = c("value = double(1)", "X = double(2)", "beta = double(2)",
                "s2 = double(1)", "df = double(0)", "P = double(2)",
                "init = double(1)"), discrete = FALSE))
  nimble::registerDistributions(c(hmmDists, .neoHMMRegBUGS()),
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

#' @describeIn buildModelCode FSSN (Ferreira-Steel skew normal)
#'   regime-switching HMM code.
#' @export
setMethod("buildModelCode", signature("FSSNUvSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    # Prior node names mirror .fssnFixedKCode() exactly and were checked
    # against buildConstants: the alpha prior is dlnorm(0, sd = aScale), NOT
    # a gamma with aA/bA as the MSNBurr families use. Getting this wrong does
    # not error -- it leaves the prior undefined and the chain diverges.
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMFSSN_k(muTilde[1:K], sigmaTilde[1:K],
                                alphaTilde[1:K], P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        muTilde[j]    ~ dnorm(mu0, sd = muSd)
        sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
        alphaTilde[j] ~ dlnorm(0, sd = aScale)
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code, monitors = c("muTilde", "sigmaTilde", "alphaTilde", "P"),
         paramNodes = c(mu = "muTilde", sigma = "sigmaTilde",
                        alpha = "alphaTilde"),
         allocNode = "zFFBS")
  }
)

#' @describeIn buildModelCode FSST (Ferreira-Steel skew t) regime-switching
#'   HMM code -- heavy-tailed and skewed regimes at once.
#' @export
setMethod("buildModelCode", signature("FSSTUvSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMFSST_k(muTilde[1:K], sigmaTilde[1:K],
                                alphaTilde[1:K], nuTilde[1:K],
                                P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        muTilde[j]    ~ dnorm(mu0, sd = muSd)
        sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
        alphaTilde[j] ~ dlnorm(0, sd = aScale)
        nuTilde[j]    ~ T(dgamma(shape = aNu, rate = bNu), 2, )
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code,
         monitors = c("muTilde", "sigmaTilde", "alphaTilde", "nuTilde", "P"),
         paramNodes = c(mu = "muTilde", sigma = "sigmaTilde",
                        alpha = "alphaTilde", nu = "nuTilde"),
         allocNode = "zFFBS")
  }
)

#' @describeIn buildModelCode SEP (skew exponential power) regime-switching HMM code.
#' @export
setMethod("buildModelCode", signature("SEPUvSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    # Prior node names mirror the family's fixedk template exactly (checked
    # against buildConstants): nu ~ dgamma(aNu, bNu), untruncated.
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMSEP_k(muTilde[1:K], sigmaTilde[1:K],
                                nuTilde[1:K], P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        muTilde[j]    ~ dnorm(mu0, sd = muSd)
        sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
        nuTilde[j]    ~ dgamma(shape = aNu, rate = bNu)
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code, monitors = c("muTilde", "sigmaTilde", "nuTilde", "P"),
         paramNodes = c(mu = "muTilde", sigma = "sigmaTilde", nu = "nuTilde"),
         allocNode = "zFFBS")
  }
)

#' @describeIn buildModelCode LEP regime-switching HMM code.
#' @export
setMethod("buildModelCode", signature("LEPUvSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    # Prior node names mirror the family's fixedk template exactly (checked
    # against buildConstants): nu ~ dgamma(aNu, bNu), untruncated.
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMLEP_k(muTilde[1:K], sigmaTilde[1:K],
                                nuTilde[1:K], P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        muTilde[j]    ~ dnorm(mu0, sd = muSd)
        sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
        nuTilde[j]    ~ dgamma(shape = aNu, rate = bNu)
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code, monitors = c("muTilde", "sigmaTilde", "nuTilde", "P"),
         paramNodes = c(mu = "muTilde", sigma = "sigmaTilde", nu = "nuTilde"),
         allocNode = "zFFBS")
  }
)

#' @describeIn buildModelCode FOSSEP regime-switching HMM code.
#' @export
setMethod("buildModelCode", signature("FOSSEPUvSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMFOSSEP_k(muTilde[1:K], sigmaTilde[1:K],
                                alphaTilde[1:K], thetaTilde[1:K],
                                P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        muTilde[j]    ~ dnorm(mu0, sd = muSd)
        sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
        alphaTilde[j] ~ dlnorm(0, sd = aScale)
        thetaTilde[j] ~ dgamma(shape = aTheta, rate = bTheta)
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code,
         monitors = c("muTilde", "sigmaTilde", "alphaTilde", "thetaTilde",
                      "P"),
         paramNodes = c(mu = "muTilde", sigma = "sigmaTilde",
                        alpha = "alphaTilde", theta = "thetaTilde"),
         allocNode = "zFFBS")
  }
)

#' @describeIn buildModelCode JFST (Jones-Faddy skew t) regime-switching HMM code.
#' @export
setMethod("buildModelCode", signature("JFSTUvSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMJFST_k(muTilde[1:K], sigmaTilde[1:K],
                                alphaTilde[1:K], thetaTilde[1:K],
                                P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        muTilde[j]    ~ dnorm(mu0, sd = muSd)
        sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)
        alphaTilde[j] ~ dgamma(shape = aSh, rate = bSh)
        thetaTilde[j] ~ dgamma(shape = aSh, rate = bSh)
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code,
         monitors = c("muTilde", "sigmaTilde", "alphaTilde", "thetaTilde",
                      "P"),
         paramNodes = c(mu = "muTilde", sigma = "sigmaTilde",
                        alpha = "alphaTilde", theta = "thetaTilde"),
         allocNode = "zFFBS")
  }
)

#' @describeIn buildModelCode Binomial regime-switching HMM code
#'   (regime-specific success probabilities, known \code{size}).
#' @export
setMethod("buildModelCode", signature("BinomialSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    # Prior node mirrors the fixedk/DPM template: prob ~ dbeta(a0, b0); the
    # number of trials `size` is a known constant (prior$size), exactly as in
    # buildConstants. The Poisson case already proved the engine is generic
    # over non-location-scale, discrete emissions.
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMBinom_k(prob[1:K], size, P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        prob[j] ~ dbeta(a0, b0)
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code, monitors = c("prob", "P"),
         paramNodes = c(prob = "prob"), allocNode = "zFFBS")
  }
)

#' @describeIn buildModelCode Markov-switching Student-t (heavy-tail)
#'   regression HMM code: location coefficients and scale switch with the
#'   regime, df a fixed hyperparameter.
#' @export
setMethod("buildModelCode", signature("StudentTRegSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMStudentTReg_k(X[1:n, 1:p], betaTilde[1:K, 1:p],
                                       s2Tilde[1:K], df, P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        s2Tilde[j] ~ dinvgamma(shape = nu0, scale = s0)
        covBeta[j, 1:p, 1:p] <- s2Tilde[j] * B0[1:p, 1:p]
        betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = covBeta[j, 1:p, 1:p])
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code, monitors = c("betaTilde", "s2Tilde", "P"),
         paramNodes = c(beta = "betaTilde", s2 = "s2Tilde"),
         allocNode = "zFFBS")
  }
)

#' @describeIn buildModelCode Markov-switching Normal-Gamma regression HMM
#'   code: same Student-t marginal as the direct parameterisation, so the same
#'   kernel. Defined explicitly because the two heavy-tail specs are siblings
#'   under NormalRegSpec, not parent and child (9.29/9.41).
#' @export
setMethod("buildModelCode", signature("NormalGammaRegSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMStudentTReg_k(X[1:n, 1:p], betaTilde[1:K, 1:p],
                                       s2Tilde[1:K], df, P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        s2Tilde[j] ~ dinvgamma(shape = nu0, scale = s0)
        covBeta[j, 1:p, 1:p] <- s2Tilde[j] * B0[1:p, 1:p]
        betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = covBeta[j, 1:p, 1:p])
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code, monitors = c("betaTilde", "s2Tilde", "P"),
         paramNodes = c(beta = "betaTilde", s2 = "s2Tilde"),
         allocNode = "zFFBS")
  }
)

#' @describeIn buildModelCode Markov-switching Binomial (proportion)
#'   regression HMM code: logit-link coefficients switch with the regime,
#'   known number of trials \code{size}.
#' @export
setMethod("buildModelCode", signature("BinomialRegSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMBinomReg_k(X[1:n, 1:p], betaTilde[1:K, 1:p], size,
                                    P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = B0[1:p, 1:p])
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code, monitors = c("betaTilde", "P"),
         paramNodes = c(beta = "betaTilde"), allocNode = "zFFBS")
  }
)

#' @describeIn buildModelCode Markov-switching Poisson (count) regression HMM
#'   code: log-link coefficients switch with a latent Markov regime, no
#'   error-variance parameter.
#' @export
setMethod("buildModelCode", signature("PoissonRegSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    # Prior mirrors the Poisson fixedk template: betaTilde ~ dmnorm(b0, B0).
    # As with the Gaussian regression HMM there is no conjugate sampler here
    # (allocations marginalised), and the guard below refuses the specs that
    # merely inherit this code.
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMPoisReg_k(X[1:n, 1:p], betaTilde[1:K, 1:p],
                                   P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = B0[1:p, 1:p])
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code, monitors = c("betaTilde", "P"),
         paramNodes = c(beta = "betaTilde"), allocNode = "zFFBS")
  }
)

#' @describeIn buildModelCode Markov-switching regression HMM code: the
#'   regression coefficients and error variance switch with a latent
#'   first-order Markov regime (Hamilton 1989).
#' @export
setMethod("buildModelCode", signature("NormalRegSpec", "HMMEngine"),
  function(spec, engine, n, L, ...) {
    # Prior nodes mirror the fixedk regression template exactly: the
    # Normal-Inverse-Gamma pair with the s2-scaled coefficient covariance.
    #
    # Note what is NOT here: the conjugate NIG sampler. customizeSamplers for
    # this spec returns early unless a "z" node exists, and under the HMM the
    # allocations are marginalised away -- so NIMBLE's defaults are used. That
    # is not a shortcut but a requirement: the conjugate update conditions on
    # the allocations, and with them integrated out the coefficient full
    # conditional is no longer Normal-Inverse-Gamma. Installing it anyway is
    # exactly the silent-wrong-sampler failure of the Student-t regression bug.
    code <- nimble::nimbleCode({
      y[1:n] ~ dRegimeHMMNormReg_k(X[1:n, 1:p], betaTilde[1:K, 1:p],
                                   s2Tilde[1:K], P[1:K, 1:K], init[1:K])
      for (j in 1:K) {
        s2Tilde[j] ~ dinvgamma(shape = nu0, scale = s0)
        covBeta[j, 1:p, 1:p] <- s2Tilde[j] * B0[1:p, 1:p]
        betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = covBeta[j, 1:p, 1:p])
        P[j, 1:K] ~ ddirch(alphaP[1:K])
      }
    })
    list(code = code, monitors = c("betaTilde", "s2Tilde", "P"),
         paramNodes = c(beta = "betaTilde", s2 = "s2Tilde"),
         allocNode = "zFFBS")
  }
)

# --- per-family emission RNG for forecasting -----------------------------------
#
# The forecast mirror of .hmmEmisDens: given a posterior draw index and a
# regime, return ONE observation drawn from that regime's emission. Node names
# are taken from the matching .hmmEmisDens method rather than remembered -- the
# families disagree (s2Tilde vs sigmaTilde vs tauTilde vs lambda vs prob) and
# guessing has cost this package a debugging session before.
#
# Vectorised over (draw, state) pairs: both arrive as equal-length vectors, so
# a whole horizon step for every posterior draw is one call.
setGeneric(".hmmEmisRng",
  function(spec, samples, K, draws, states, prior = NULL, xrow = NULL)
    standardGeneric(".hmmEmisRng"))

.rngPick <- function(samples, node, K, draws, states) {
  A <- .nodeToArray(samples, node, K)
  A[cbind(draws, states)]
}

setMethod(".hmmEmisRng", "NormalUvSpec",
  function(spec, samples, K, draws, states, prior = NULL, xrow = NULL)
    stats::rnorm(length(draws), .rngPick(samples, "muTilde", K, draws, states),
                 sqrt(.rngPick(samples, "s2Tilde", K, draws, states))))

setMethod(".hmmEmisRng", "StudentTUvSpec",
  function(spec, samples, K, draws, states, prior = NULL, xrow = NULL) {
    mu <- .rngPick(samples, "muTilde", K, draws, states)
    sg <- 1 / sqrt(.rngPick(samples, "tauTilde", K, draws, states))
    mu + sg * stats::rt(length(draws), df = prior$df)
  })

setMethod(".hmmEmisRng", "PoissonSpec",
  function(spec, samples, K, draws, states, prior = NULL, xrow = NULL)
    stats::rpois(length(draws), .rngPick(samples, "lambda", K, draws, states)))

setMethod(".hmmEmisRng", "BinomialSpec",
  function(spec, samples, K, draws, states, prior = NULL, xrow = NULL)
    stats::rbinom(length(draws), prior$size,
                  .rngPick(samples, "prob", K, draws, states)))

setMethod(".hmmEmisRng", "NormalRegSpec",
  function(spec, samples, K, draws, states, prior = NULL, xrow = NULL) {
    p <- prior$p
    bA <- .nodeToArray(samples, "betaTilde", c(K, p))
    # xrow is either ONE row shared by every draw (fixed covariates) or a
    # draws x p matrix (autoregressive: each posterior path carries its own
    # trajectory, so its own lag).
    X <- if (is.matrix(xrow)) xrow else
      matrix(xrow, nrow = length(draws), ncol = p, byrow = TRUE)
    mu <- vapply(seq_along(draws),
                 function(i) sum(X[i, ] * bA[draws[i], states[i], ]),
                 numeric(1))
    stats::rnorm(length(draws), mu,
                 sqrt(.rngPick(samples, "s2Tilde", K, draws, states)))
  })

setMethod(".hmmEmisRng", "PoissonRegSpec",
  function(spec, samples, K, draws, states, prior = NULL, xrow = NULL) {
    p <- prior$p
    bA <- .nodeToArray(samples, "betaTilde", c(K, p))
    X <- if (is.matrix(xrow)) xrow else
      matrix(xrow, nrow = length(draws), ncol = p, byrow = TRUE)
    mu <- vapply(seq_along(draws),
                 function(i) exp(sum(X[i, ] * bA[draws[i], states[i], ])),
                 numeric(1))
    stats::rpois(length(draws), mu)
  })

setMethod(".hmmEmisRng", "BinomialRegSpec",
  function(spec, samples, K, draws, states, prior = NULL, xrow = NULL) {
    p <- prior$p; sz <- prior$size
    bA <- .nodeToArray(samples, "betaTilde", c(K, p))
    X <- if (is.matrix(xrow)) xrow else
      matrix(xrow, nrow = length(draws), ncol = p, byrow = TRUE)
    pr <- vapply(seq_along(draws),
                 function(i) stats::plogis(sum(X[i, ] * bA[draws[i], states[i], ])),
                 numeric(1))
    stats::rbinom(length(draws), sz, pr)
  })

.hmmEmisRngStudentTReg <- function(spec, samples, K, draws, states, prior = NULL, xrow = NULL) {
  p <- prior$p; df <- prior$df
  bA <- .nodeToArray(samples, "betaTilde", c(K, p))
  s2A <- .nodeToArray(samples, "s2Tilde", K)
  X <- if (is.matrix(xrow)) xrow else
    matrix(xrow, nrow = length(draws), ncol = p, byrow = TRUE)
  mu <- vapply(seq_along(draws),
               function(i) sum(X[i, ] * bA[draws[i], states[i], ]), numeric(1))
  sg <- sqrt(vapply(seq_along(draws), function(i) s2A[draws[i], states[i]], numeric(1)))
  mu + sg * stats::rt(length(draws), df = df)
}
setMethod(".hmmEmisRng", "StudentTRegSpec", .hmmEmisRngStudentTReg)
setMethod(".hmmEmisRng", "NormalGammaRegSpec", .hmmEmisRngStudentTReg)

.mkRng3 <- function(rfun, n3) function(spec, samples, K, draws, states,
                                       prior = NULL, xrow = NULL) {
  rfun(length(draws),
       .rngPick(samples, "muTilde", K, draws, states),
       .rngPick(samples, "sigmaTilde", K, draws, states),
       .rngPick(samples, n3, K, draws, states))
}
.mkRng4 <- function(rfun) function(spec, samples, K, draws, states,
                                   prior = NULL, xrow = NULL) {
  rfun(length(draws),
       .rngPick(samples, "muTilde", K, draws, states),
       .rngPick(samples, "sigmaTilde", K, draws, states),
       .rngPick(samples, "alphaTilde", K, draws, states),
       .rngPick(samples, "thetaTilde", K, draws, states))
}

setMethod(".hmmEmisRng", "MSNBurrUvSpec",   .mkRng3(rmsnburr,   "alphaTilde"))
setMethod(".hmmEmisRng", "MSNBurr2aUvSpec", .mkRng3(rmsnburr2a, "alphaTilde"))
setMethod(".hmmEmisRng", "FSSNUvSpec",      .mkRng3(rfssn,      "alphaTilde"))
setMethod(".hmmEmisRng", "SEPUvSpec",       .mkRng3(rsep,       "nuTilde"))
setMethod(".hmmEmisRng", "LEPUvSpec",       .mkRng3(rlep,       "nuTilde"))
setMethod(".hmmEmisRng", "GMSNBurrUvSpec",  .mkRng4(rgmsnburr))
setMethod(".hmmEmisRng", "FOSSEPUvSpec",    .mkRng4(rfossep))
setMethod(".hmmEmisRng", "JFSTUvSpec",      .mkRng4(rjfst))

setMethod(".hmmEmisRng", "FSSTUvSpec",
  function(spec, samples, K, draws, states, prior = NULL, xrow = NULL)
    rfsst(length(draws),
          .rngPick(samples, "muTilde", K, draws, states),
          .rngPick(samples, "sigmaTilde", K, draws, states),
          .rngPick(samples, "alphaTilde", K, draws, states),
          .rngPick(samples, "nuTilde", K, draws, states)))

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

setMethod(".hmmEmisDens", "FSSNUvSpec",
  function(spec, samples, y, K, draw = NULL, prior = NULL) {
    muM <- .nodeToArray(samples, "muTilde", K)
    sgM <- .nodeToArray(samples, "sigmaTilde", K)
    alM <- .nodeToArray(samples, "alphaTilde", K)
    mu <- if (is.null(draw)) colMeans(muM) else muM[draw, ]
    sg <- if (is.null(draw)) colMeans(sgM) else sgM[draw, ]
    al <- if (is.null(draw)) colMeans(alM) else alM[draw, ]
    vapply(seq_len(K), function(s) dfssn(y, mu[s], sg[s], al[s]),
           numeric(length(y)))
  })

setMethod(".hmmEmisDens", "FSSTUvSpec",
  function(spec, samples, y, K, draw = NULL, prior = NULL) {
    muM <- .nodeToArray(samples, "muTilde", K)
    sgM <- .nodeToArray(samples, "sigmaTilde", K)
    alM <- .nodeToArray(samples, "alphaTilde", K)
    nuM <- .nodeToArray(samples, "nuTilde", K)
    mu <- if (is.null(draw)) colMeans(muM) else muM[draw, ]
    sg <- if (is.null(draw)) colMeans(sgM) else sgM[draw, ]
    al <- if (is.null(draw)) colMeans(alM) else alM[draw, ]
    nu <- if (is.null(draw)) colMeans(nuM) else nuM[draw, ]
    vapply(seq_len(K), function(s) dfsst(y, mu[s], sg[s], al[s], nu[s]),
           numeric(length(y)))
  })

setMethod(".hmmEmisDens", "SEPUvSpec",
  function(spec, samples, y, K, draw = NULL, prior = NULL) {
    muM <- .nodeToArray(samples, "muTilde", K)
    sgM <- .nodeToArray(samples, "sigmaTilde", K)
    nuM <- .nodeToArray(samples, "nuTilde", K)
    mu <- if (is.null(draw)) colMeans(muM) else muM[draw, ]
    sg <- if (is.null(draw)) colMeans(sgM) else sgM[draw, ]
    nu <- if (is.null(draw)) colMeans(nuM) else nuM[draw, ]
    vapply(seq_len(K), function(s) dsep(y, mu[s], sg[s], nu[s]),
           numeric(length(y)))
  })

setMethod(".hmmEmisDens", "LEPUvSpec",
  function(spec, samples, y, K, draw = NULL, prior = NULL) {
    muM <- .nodeToArray(samples, "muTilde", K)
    sgM <- .nodeToArray(samples, "sigmaTilde", K)
    nuM <- .nodeToArray(samples, "nuTilde", K)
    mu <- if (is.null(draw)) colMeans(muM) else muM[draw, ]
    sg <- if (is.null(draw)) colMeans(sgM) else sgM[draw, ]
    nu <- if (is.null(draw)) colMeans(nuM) else nuM[draw, ]
    vapply(seq_len(K), function(s) dlep(y, mu[s], sg[s], nu[s]),
           numeric(length(y)))
  })

setMethod(".hmmEmisDens", "FOSSEPUvSpec",
  function(spec, samples, y, K, draw = NULL, prior = NULL) {
    muM <- .nodeToArray(samples, "muTilde", K)
    sgM <- .nodeToArray(samples, "sigmaTilde", K)
    alM <- .nodeToArray(samples, "alphaTilde", K)
    thM <- .nodeToArray(samples, "thetaTilde", K)
    mu <- if (is.null(draw)) colMeans(muM) else muM[draw, ]
    sg <- if (is.null(draw)) colMeans(sgM) else sgM[draw, ]
    al <- if (is.null(draw)) colMeans(alM) else alM[draw, ]
    th <- if (is.null(draw)) colMeans(thM) else thM[draw, ]
    vapply(seq_len(K), function(s) dfossep(y, mu[s], sg[s], al[s], th[s]),
           numeric(length(y)))
  })

setMethod(".hmmEmisDens", "JFSTUvSpec",
  function(spec, samples, y, K, draw = NULL, prior = NULL) {
    muM <- .nodeToArray(samples, "muTilde", K)
    sgM <- .nodeToArray(samples, "sigmaTilde", K)
    alM <- .nodeToArray(samples, "alphaTilde", K)
    thM <- .nodeToArray(samples, "thetaTilde", K)
    mu <- if (is.null(draw)) colMeans(muM) else muM[draw, ]
    sg <- if (is.null(draw)) colMeans(sgM) else sgM[draw, ]
    al <- if (is.null(draw)) colMeans(alM) else alM[draw, ]
    th <- if (is.null(draw)) colMeans(thM) else thM[draw, ]
    vapply(seq_len(K), function(s) djfst(y, mu[s], sg[s], al[s], th[s]),
           numeric(length(y)))
  })

setMethod(".hmmEmisDens", "BinomialSpec",
  function(spec, samples, y, K, draw = NULL, prior = NULL) {
    pM <- .nodeToArray(samples, "prob", K)
    p <- if (is.null(draw)) colMeans(pM) else pM[draw, ]
    sz <- prior$size
    vapply(seq_len(K), function(s) stats::dbinom(round(y), sz, p[s]),
           numeric(length(y)))
  })

.hmmEmisDensStudentTReg <- function(spec, samples, y, K, draw = NULL, prior = NULL) {
  X <- prior$X; p <- prior$p; df <- prior$df
  bA <- .nodeToArray(samples, "betaTilde", c(K, p))
  s2A <- .nodeToArray(samples, "s2Tilde", K)
  B <- if (is.null(draw)) apply(bA, c(2L, 3L), mean) else matrix(bA[draw, , ], K, p)
  s2 <- if (is.null(draw)) colMeans(s2A) else s2A[draw, ]
  vapply(seq_len(K), function(s)
    stats::dt((y - as.numeric(X %*% B[s, ])) / sqrt(s2[s]), df = df) / sqrt(s2[s]),
    numeric(length(y)))
}
setMethod(".hmmEmisDens", "StudentTRegSpec", .hmmEmisDensStudentTReg)
setMethod(".hmmEmisDens", "NormalGammaRegSpec", .hmmEmisDensStudentTReg)

# ---------------------------------------------------------------------------
# Neo-normal regression HMM methods (generated). For each family, three methods
# depend only on: the kernel base name, the shape parameter names, the shape
# prior lines (for the model code), and the R-level density/RNG (for the FFBS
# allocation density and the forecast RNG). .makeNeoHMMRegMethods() builds all
# three, so a family under the HMM engine needs no hand-written S4 method.
# ---------------------------------------------------------------------------
.makeNeoHMMRegMethods <- function(class, kernelBase, shapeNames, priorLines,
                                  densR, rngR) {
  kernel <- paste0("dRegimeHMM", kernelBase, "_k")
  shapeNodes <- paste0(shapeNames, "Tilde")

  # buildModelCode: y ~ kernel(X, beta, <shapeTilde[1:K]...>, P, init)
  setMethod("buildModelCode", signature(class, "HMMEngine"),
    function(spec, engine, n, L, ...) {
      shapeVec <- paste(sprintf("%s[1:K]", shapeNodes), collapse = ", ")
      codeSrc <- sprintf("nimble::nimbleCode({
        y[1:n] ~ %s(X[1:n, 1:p], betaTilde[1:K, 1:p], %s, P[1:K, 1:K], init[1:K])
        for (j in 1:K) {
%s
          betaTilde[j, 1:p] ~ dmnorm(b0[1:p], cov = B0[1:p, 1:p])
          P[j, 1:K] ~ ddirch(alphaP[1:K])
        }
      })", kernel, shapeVec, priorLines)
      code <- eval(parse(text = codeSrc))
      list(code = code,
           monitors = c("betaTilde", shapeNodes, "P"),
           paramNodes = c(beta = "betaTilde"), allocNode = "zFFBS")
    })

  # .hmmEmisDens: per-state emission density at location X beta_s
  setMethod(".hmmEmisDens", class,
    function(spec, samples, y, K, draw = NULL, prior = NULL) {
      X <- prior$X; p <- prior$p
      bA <- .nodeToArray(samples, "betaTilde", c(K, p))
      shA <- lapply(shapeNodes, function(nd) .nodeToArray(samples, nd, K))
      B <- if (is.null(draw)) apply(bA, c(2L, 3L), mean) else matrix(bA[draw, , ], K, p)
      shv <- lapply(shA, function(A) if (is.null(draw)) colMeans(A) else A[draw, ])
      vapply(seq_len(K), function(s) {
        args <- c(list(y, as.numeric(X %*% B[s, ])),
                  lapply(shv, function(v) v[s]))
        do.call(densR, args)
      }, numeric(length(y)))
    })

  # .hmmEmisRng: draw a response at X beta_s with per-state shape
  setMethod(".hmmEmisRng", class,
    function(spec, samples, K, draws, states, prior = NULL, xrow = NULL) {
      p <- prior$p
      bA <- .nodeToArray(samples, "betaTilde", c(K, p))
      shA <- lapply(shapeNodes, function(nd) .nodeToArray(samples, nd, K))
      X <- if (is.matrix(xrow)) xrow else
        matrix(xrow, nrow = length(draws), ncol = p, byrow = TRUE)
      mu <- vapply(seq_along(draws),
                   function(i) sum(X[i, ] * bA[draws[i], states[i], ]), numeric(1))
      vapply(seq_along(draws), function(i) {
        args <- c(list(1, mu[i]),
                  lapply(shA, function(A) A[draws[i], states[i]]))
        do.call(rngR, args)
      }, numeric(1))
    })
  invisible(TRUE)
}

# Wire up all nine neo-normal families under the HMM engine.
.makeNeoHMMRegMethods("MSNBurrRegSpec", "MSNBReg", c("sigma", "alpha"),
  "          sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)\n          alphaTilde[j] ~ dgamma(shape = aA, rate = bA)",
  dmsnburr, rmsnburr)
.makeNeoHMMRegMethods("MSNBurr2aRegSpec", "MSNB2aReg", c("sigma", "alpha"),
  "          sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)\n          alphaTilde[j] ~ dgamma(shape = aA, rate = bA)",
  dmsnburr2a, rmsnburr2a)
.makeNeoHMMRegMethods("FSSNRegSpec", "FSSNReg", c("sigma", "alpha"),
  "          sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)\n          alphaTilde[j] ~ dlnorm(0, sd = aScale)",
  dfssn, rfssn)
.makeNeoHMMRegMethods("SEPRegSpec", "SEPReg", c("sigma", "nu"),
  "          sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)\n          nuTilde[j] ~ dgamma(shape = aNu, rate = bNu)",
  dsep, rsep)
.makeNeoHMMRegMethods("LEPRegSpec", "LEPReg", c("sigma", "nu"),
  "          sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)\n          nuTilde[j] ~ dgamma(shape = aNu, rate = bNu)",
  dlep, rlep)
.makeNeoHMMRegMethods("GMSNBurrRegSpec", "GMSNBReg", c("sigma", "alpha", "theta"),
  "          sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)\n          alphaTilde[j] ~ dgamma(shape = aA, rate = bA)\n          thetaTilde[j] ~ dgamma(shape = aT, rate = bT)",
  dgmsnburr, rgmsnburr)
.makeNeoHMMRegMethods("FSSTRegSpec", "FSSTReg", c("sigma", "alpha", "nu"),
  "          sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)\n          alphaTilde[j] ~ dlnorm(0, sd = aScale)\n          nuTilde[j] ~ T(dgamma(shape = aNu, rate = bNu), 2, )",
  dfsst, rfsst)
.makeNeoHMMRegMethods("FOSSEPRegSpec", "FOSSEPReg", c("sigma", "alpha", "theta"),
  "          sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)\n          alphaTilde[j] ~ dlnorm(0, sd = aScale)\n          thetaTilde[j] ~ dgamma(shape = aTheta, rate = bTheta)",
  dfossep, rfossep)
.makeNeoHMMRegMethods("JFSTRegSpec", "JFSTReg", c("sigma", "alpha", "theta"),
  "          sigmaTilde[j] ~ dinvgamma(shape = aSig, scale = bSig)\n          alphaTilde[j] ~ dgamma(shape = aSh, rate = bSh)\n          thetaTilde[j] ~ dgamma(shape = aSh, rate = bSh)",
  djfst, rjfst)


setMethod(".hmmEmisDens", "BinomialRegSpec",
  function(spec, samples, y, K, draw = NULL, prior = NULL) {
    X <- prior$X; p <- prior$p; sz <- prior$size
    bA <- .nodeToArray(samples, "betaTilde", c(K, p))
    B <- if (is.null(draw)) apply(bA, c(2L, 3L), mean) else
      matrix(bA[draw, , ], nrow = K, ncol = p)
    vapply(seq_len(K),
           function(s) stats::dbinom(round(y), sz,
                                     stats::plogis(as.numeric(X %*% B[s, ]))),
           numeric(length(y)))
  })

setMethod(".hmmEmisDens", "PoissonRegSpec",
  function(spec, samples, y, K, draw = NULL, prior = NULL) {
    X <- prior$X; p <- prior$p
    bA <- .nodeToArray(samples, "betaTilde", c(K, p))   # draws x K x p
    B <- if (is.null(draw)) apply(bA, c(2L, 3L), mean) else
      matrix(bA[draw, , ], nrow = K, ncol = p)
    vapply(seq_len(K),
           function(s) stats::dpois(round(y), exp(as.numeric(X %*% B[s, ]))),
           numeric(length(y)))
  })

setMethod(".hmmEmisDens", "NormalRegSpec",
  function(spec, samples, y, K, draw = NULL, prior = NULL) {
    # The only emission whose density varies with t through a covariate: the
    # design matrix rides along in the prior, as buildConstants already
    # assumes.
    X <- prior$X; p <- prior$p
    bA <- .nodeToArray(samples, "betaTilde", c(K, p))   # draws x K x p
    s2M <- .nodeToArray(samples, "s2Tilde", K)
    if (is.null(draw)) {
      B <- apply(bA, c(2L, 3L), mean)
      s2 <- colMeans(s2M)
    } else {
      B <- matrix(bA[draw, , ], nrow = K, ncol = p)
      s2 <- s2M[draw, ]
    }
    vapply(seq_len(K),
           function(s) stats::dnorm(y, as.numeric(X %*% B[s, ]), sqrt(s2[s])),
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

    # EXACT class match, not is(): NormalGammaUvSpec `contains` NormalUvSpec,
    # so the old is()-based guard silently accepted
    # distribution = "normal-gamma" and dispatch then fell through to the
    # inherited (NormalUvSpec, HMMEngine) method -- the user asked for heavy
    # tails and got a plain Gaussian HMM with no error (the 9.29 inheritance
    # class again; the permanent guard test caught it). normal-gamma stays
    # excluded by design: its augmented representation is what the
    # marginalised forward kernel exists to avoid (9.13), and direct
    # student-t serves the heavy-tail case.
    hmmFamilies <- c("NormalUvSpec", "StudentTUvSpec", "PoissonSpec",
                     "BinomialSpec", "NormalRegSpec", "PoissonRegSpec",
                     "BinomialRegSpec", "StudentTRegSpec", "NormalGammaRegSpec",
                     "MSNBurrRegSpec", "MSNBurr2aRegSpec", "FSSNRegSpec",
                     "SEPRegSpec", "LEPRegSpec", "GMSNBurrRegSpec",
                     "FSSTRegSpec", "FOSSEPRegSpec", "JFSTRegSpec",
                     "MSNBurrUvSpec", "MSNBurr2aUvSpec", "GMSNBurrUvSpec",
                     "FSSNUvSpec", "FSSTUvSpec", "SEPUvSpec", "LEPUvSpec",
                     "FOSSEPUvSpec", "JFSTUvSpec")
    okSpec <- class(spec)[1L] %in% hmmFamilies
    if (!okSpec)
      stop("method = 'hmm' currently supports distribution = \"normal\", ",
           "\"student-t\", \"poisson\", \"msnburr\", \"msnburr2a\", ",
           "\"binomial\", \"gmsnburr\", \"fssn\", \"fsst\", \"sep\", ",
           "\"lep\", \"fossep\", or \"jfst\" (univariate); further ",
           "emission families follow the gated plan.", call. = FALSE)
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

# --- forecasting ---------------------------------------------------------------

#' Forecast ahead from a regime-switching (HMM) fit
#'
#' Draws from the posterior predictive distribution \code{h} steps beyond the
#' end of the series, integrating over both parameter uncertainty (one path
#' per posterior draw) and regime uncertainty (the regime is sampled, not
#' fixed).
#'
#' The recipe is the forward algorithm run one step further: filter the regime
#' distribution to the last observation, push it through the transition matrix
#' \eqn{h} times, sample a regime from it, and draw from that regime's
#' emission.
#'
#' @section What an HMM can and cannot forecast:
#' The point forecast decays toward the stationary mixture as \code{h} grows,
#' at rate \eqn{|\lambda_2|^h} in the second eigenvalue of the transition
#' matrix. That is the model working, not failing: an HMM forecasts the
#' \emph{regime}, and once the regime is unknowable the best guess for
#' \code{y} is the long-run average over regimes. Measured on a two-regime
#' series at +/-2 with persistence 0.95/0.90, the median forecast fell from
#' 1.8 at \code{h = 1} to -1.1 by \code{h = 12}, while the regime
#' probabilities went from 0.11/0.89 to 0.59/0.41 -- converging on the
#' stationary 0.67/0.33. Reversion, not signal. Interval coverage over the
#' same horizon was 0.917 against a nominal 0.90.
#'
#' \strong{Everything hinges on whether the regime persists.} On a
#' Markov-switching regression whose regime held across the forecast
#' boundary, the median forecast beat a constant-mean benchmark fivefold
#' (RMSE 0.50 against 2.65); on a series that happened to switch at exactly
#' the first forecast step -- a 5\% event under the fitted persistence -- it
#' lost (4.42 against 2.45). That is not a defect to tune away: a switch one
#' step ahead is unpredictable by construction, and the model correctly
#' assigned it low probability beforehand. The interval still covered.
#'
#' \strong{Beware the point summary when regimes disagree.} With well
#' separated regimes the predictive is \emph{bimodal}, and its median sits
#' in the trough between the modes -- where the model puts almost no mass.
#' \code{$summary} is convenient, not faithful; \code{$draws} is the honest
#' object, and \code{$regime} usually answers the question actually being
#' asked ("which regime will we be in"), degrading far more gracefully than
#' any forecast of \code{y}.
#'
#' For a trend or a seasonal pattern this is the wrong model; nimix does not
#' fit those.
#'
#' @section Autoregression:
#' A Markov-switching autoregression is an ordinary Markov-switching
#' regression with the response lagged into the design matrix
#' (\code{nimixReg(y ~ ylag, df, K = 2, method = "hmm")}), and forecasting
#' one needs \code{lags} rather than \code{newdata}: nobody knows the future
#' \eqn{y}, so the lag has to come from the forecast itself. Each posterior
#' draw carries its own trajectory and feeds it back as its own lag, which is
#' what makes the interval widen properly -- measured on an MS-AR(1), from
#' 2.15 at \code{h = 1} to 5.12 at \code{h = 10}, with RMSE 0.54 against
#' 3.97 for a constant-mean benchmark. With exogenous predictors as well,
#' pass those in \code{newdata} and the lags in \code{lags} (measured on an
#' MS-ARX: RMSE 0.83 against 5.29).
#'
#' \strong{Do not fake a lag through \code{newdata}.} It will be taken as a
#' known covariate, held fixed, and the forecast will inherit whatever you
#' invented -- confidently and without complaint. \code{lags} exists so that
#' you do not have to.
#'
#' @param object A \code{\linkS4class{FitResult}} fitted with
#'   \code{method = "hmm"}.
#' @param h Integer horizon (>= 1).
#' @param newdata For a Markov-switching regression, a data frame of
#'   \code{h} future rows supplying the \emph{exogenous} predictors. Ignored
#'   otherwise.
#' @param lags Named integer vector marking predictors that are the response
#'   in disguise, e.g. \code{lags = c(ylag = 1)} for a column holding
#'   \eqn{y_{t-1}}. Those columns are \strong{generated by the forecast},
#'   not supplied to it: each posterior draw feeds its own trajectory back as
#'   its own lag, so the uncertainty compounds the way an autoregression's
#'   should. Give a predictor in \code{lags} or in \code{newdata}, never
#'   both.
#' @param draws Maximum posterior draws to use (thinned evenly).
#' @param level Central credible level for the interval.
#' @return A list with \code{summary} (a data frame of \code{h} rows:
#'   \code{mean}, \code{median}, \code{lower}, \code{upper}), \code{regime}
#'   (an \code{h} x \code{K} matrix of regime probabilities), and
#'   \code{draws} (the raw \code{draws} x \code{h} predictive sample).
#' @examples
#' \dontrun{
#' f  <- nimixClust(y, K = 2, method = "hmm")
#' fc <- nimixForecast(f, h = 6)
#' fc$summary
#' fc$regime          # usually the more informative half
#'
#' # a Markov-switching AR(1): the lag is generated, not supplied
#' g  <- nimixReg(y ~ ylag, df, K = 2, method = "hmm")
#' nimixForecast(g, h = 10, lags = c(ylag = 1))
#' }
#' @seealso \code{\link{viterbiPath}}
#' @export
nimixForecast <- function(object, h = 1L, newdata = NULL, lags = NULL,
                          draws = 500L, level = 0.9) {
  if (!identical(object@engineUsed, "hmm"))
    stop("nimixForecast() needs a fit from method = 'hmm'; this one used '",
         object@engineUsed, "'. Forecasting means pushing a regime forward ",
         "in time, which the other engines do not model.", call. = FALSE)
  h <- as.integer(h)
  if (length(h) != 1L || is.na(h) || h < 1L)
    stop("h must be a single integer >= 1.", call. = FALSE)
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must be a single number strictly between 0 and 1.",
         call. = FALSE)

  spec <- object@distSpec
  prior <- object@prior
  K <- object@Kmax
  S <- object@mcmcSamples
  y <- object@data
  # Any regression emission (Gaussian, Poisson, ...) carries a design matrix
  # and needs future covariates; test the shared trait, not one class. Using
  # is(spec, "NormalRegSpec") would have missed PoissonRegSpec, a sibling
  # rather than a subclass.
  isReg <- isRegressionSpec(spec)

  Xf <- NULL; lagIdx <- NULL
  if (isReg) {
    tt <- stats::delete.response(prior$terms)
    need <- all.vars(tt)
    # `lags` names the predictors that ARE the response, shifted. They cannot
    # come from newdata -- nobody knows the future y -- so they are fed back
    # from the forecast itself, one trajectory per posterior draw.
    if (!is.null(lags)) {
      if (is.null(names(lags)) || anyDuplicated(names(lags)) ||
          any(!is.finite(lags)) || any(lags < 1) ||
          any(lags != as.integer(lags)))
        stop("`lags` must be a named vector of positive integers, e.g. ",
             "lags = c(ylag = 1) meaning the column `ylag` holds y at t - 1.",
             call. = FALSE)
      bad <- setdiff(names(lags), need)
      if (length(bad))
        stop("`lags` names predictor(s) the model does not use: ",
             paste(bad, collapse = ", "), ".", call. = FALSE)
      if (any(names(lags) %in% names(newdata)))
        stop("Give a predictor in `lags` or in `newdata`, not both: ",
             paste(intersect(names(lags), names(newdata)), collapse = ", "),
             ". A lagged response is generated by the forecast, not supplied ",
             "to it.", call. = FALSE)
      # as.integer() drops names; keep them, they are the whole contract.
      lags <- stats::setNames(as.integer(lags), names(lags))
    }
    exo <- setdiff(need, names(lags))
    if (length(exo) && is.null(newdata))
      stop("A Markov-switching regression needs `newdata` for its ",
           "non-lagged predictor(s): ", paste(exo, collapse = ", "),
           ". The regime can be projected forward, exogenous covariates ",
           "cannot.", call. = FALSE)
    if (is.null(lags)) {
      Xf <- stats::model.matrix(tt, newdata)
      if (nrow(Xf) != h)
        stop("`newdata` has ", nrow(Xf), " row(s) but h = ", h, ".",
             call. = FALSE)
      if (ncol(Xf) != prior$p)
        stop("`newdata` gives ", ncol(Xf), " predictor column(s); the fit ",
             "used ", prior$p, ".", call. = FALSE)
    } else {
      if (length(exo)) {
        if (nrow(newdata) != h)
          stop("`newdata` has ", nrow(newdata), " row(s) but h = ", h, ".",
               call. = FALSE)
      }
      # Build one template row so the column layout (and factor contrasts)
      # match the fit exactly; the lag columns are overwritten per step.
      tmpl <- as.data.frame(lapply(need, function(v) {
        if (v %in% names(lags)) y[length(y)] else newdata[[v]][1L]
      }))
      names(tmpl) <- need
      Xtmpl <- stats::model.matrix(tt, tmpl)
      if (ncol(Xtmpl) != prior$p)
        stop("The lag/newdata combination gives ", ncol(Xtmpl),
             " predictor column(s); the fit used ", prior$p, ".",
             call. = FALSE)
      lagIdx <- match(names(lags), colnames(Xtmpl))
      if (anyNA(lagIdx))
        stop("Lagged predictor(s) ", paste(names(lags)[is.na(lagIdx)],
             collapse = ", "), " do not appear as plain columns of the design ",
             "matrix; nimixForecast() cannot feed a lag back through a ",
             "transformation.", call. = FALSE)
      if (max(lags) > length(y))
        stop("lags reach back ", max(lags), " step(s) but the series has only ",
             length(y), " observation(s).", call. = FALSE)
    }
  }

  Parr <- .nodeToArray(S, "P", c(K, K))
  m <- nrow(S)
  use <- if (m > draws)
    as.integer(round(seq(1, m, length.out = draws))) else seq_len(m)
  nu <- length(use)

  # Filter each draw to the end of the series. This repeats the forward pass
  # that the compiled kernel already does internally -- the kernel returns only
  # the likelihood, so alpha_T is not recoverable from the fit and has to be
  # recomputed here.
  Tn <- length(y)
  alpha <- matrix(0, nu, K)
  for (i in seq_len(nu)) {
    g <- use[i]
    D <- .hmmEmisDens(spec, S, y, K, draw = g, prior = prior)
    D[!is.finite(D)] <- 0
    P <- matrix(Parr[g, , ], K, K)
    a <- D[1L, ] / K
    s <- sum(a); a <- if (s > 0) a / s else rep(1 / K, K)
    if (Tn >= 2L) for (t in 2:Tn) {
      a <- as.numeric(crossprod(P, a)) * D[t, ]
      s <- sum(a); a <- if (s > 0) a / s else rep(1 / K, K)
    }
    alpha[i, ] <- a
  }

  fc <- matrix(NA_real_, nu, h)
  regime <- matrix(0, h, K)
  # For the autoregressive case each draw needs its own design row, because
  # each draw has its own forecast trajectory to lag. Xrec is that row set,
  # rebuilt every step; the exogenous columns (if any) are shared.
  Xrec <- NULL
  if (!is.null(lagIdx)) {
    # Start from the template row: that fixes the column layout AND carries
    # the columns that are neither lagged nor exogenous -- the intercept,
    # which is 1 forever and belongs to neither category.
    Xrec <- matrix(Xtmpl[1L, ], nu, prior$p, byrow = TRUE)
    hasExo <- length(exo) > 0L
  }
  for (step in seq_len(h)) {
    for (i in seq_len(nu))
      alpha[i, ] <- as.numeric(crossprod(matrix(Parr[use[i], , ], K, K),
                                         alpha[i, ]))
    regime[step, ] <- colMeans(alpha)
    st <- vapply(seq_len(nu),
                 function(i) sample.int(K, 1L, prob = alpha[i, ]), integer(1))
    xr <- NULL
    if (isReg && is.null(lagIdx)) {
      xr <- Xf[step, ]
    } else if (isReg) {
      # exogenous part: this step's row of newdata, shared by every draw.
      # Rebuilt through model.matrix() so factor contrasts stay the fit's.
      if (hasExo) {
        rowDf <- newdata[step, , drop = FALSE]
        for (v in names(lags)) rowDf[[v]] <- y[length(y)]   # placeholder
        Xstep <- stats::model.matrix(tt, rowDf)
        Xrec[] <- matrix(Xstep[1L, ], nu, prior$p, byrow = TRUE)
      }
      # lagged part: y at t - lag, which is an OBSERVED value while the lag
      # still reaches into the data and this draw's own forecast once it does
      # not. That is the whole point -- the uncertainty compounds because each
      # path feeds on itself.
      for (j in seq_along(lags)) {
        back <- step - lags[[j]]
        Xrec[, lagIdx[j]] <- if (back >= 1L) fc[, back] else
          y[length(y) + back]
      }
      xr <- Xrec
    }
    fc[, step] <- .hmmEmisRng(spec, S, K, use, st, prior = prior, xrow = xr)
  }

  a <- (1 - level) / 2
  qs <- apply(fc, 2L, stats::quantile, probs = c(a, 0.5, 1 - a),
              names = FALSE, na.rm = TRUE)
  colnames(regime) <- paste0("regime", seq_len(K))
  list(
    summary = data.frame(h = seq_len(h), mean = colMeans(fc),
                         median = qs[2, ], lower = qs[1, ], upper = qs[3, ]),
    regime  = regime,
    draws   = fc
  )
}
