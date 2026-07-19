#' @include dist-msnburr.R
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
#' @include dist-poisson-binomial.R
#' @include dist-mv-reg.R
#' @include dist-heavytail-reg.R
#' @include dist-glm-reg.R
#' @include dist-normal-gamma.R
#' @include dist-student-t.R
#' @include dist-student-t-mv.R
#' @include dist-normal-gamma-mv.R
NULL

## ---------------------------------------------------------------------------
## registerDistribution.R
##
## A tiny registry so advanced users can register their own DistributionSpec
## subclasses and refer to them by name in nimixClust(distribution = ...).
## Built-in distributions are registered on package load (see .onLoad).
## ---------------------------------------------------------------------------

.distRegistry <- new.env(parent = emptyenv())

#' Register a component distribution
#'
#' Adds a \code{\linkS4class{DistributionSpec}} to the registry under its
#' \code{name} slot so it can be selected by name. New built-in distributions
#' (Student-t, Poisson/Binomial) are planned for v0.4.0.
#'
#' @param spec A \code{\linkS4class{DistributionSpec}} instance.
#' @param overwrite Logical; overwrite an existing entry of the same name?
#' @return Invisibly, the registered name.
#' @examples
#' registerDistribution(NormalUvSpec(), overwrite = TRUE)
#' listDistributions()
#' @export
registerDistribution <- function(spec, overwrite = FALSE) {
  if (!methods::is(spec, "DistributionSpec"))
    stop("spec must inherit from DistributionSpec.", call. = FALSE)
  nm <- spec@name
  if (exists(nm, envir = .distRegistry, inherits = FALSE) && !overwrite)
    stop("Distribution '", nm, "' already registered; use overwrite = TRUE.",
         call. = FALSE)
  assign(nm, spec, envir = .distRegistry)
  invisible(nm)
}

#' Retrieve a registered distribution by name
#' @param name Character scalar.
#' @return A \code{\linkS4class{DistributionSpec}}.
#' @export
getDistribution <- function(name) {
  if (!exists(name, envir = .distRegistry, inherits = FALSE))
    stop("Unknown distribution '", name, "'. Registered: ",
         paste(listDistributions(), collapse = ", "), call. = FALSE)
  get(name, envir = .distRegistry, inherits = FALSE)
}

#' List registered distribution names
#' @return Character vector of registered names.
#' @export
listDistributions <- function() sort(ls(envir = .distRegistry))

#' @keywords internal
.nimixDefineMSNBurr <- function() {
  if (exists("dSkewMvITOG_k", envir = globalenv(), inherits = FALSE)) return(invisible())
  ge <- globalenv()
  softlomega <- quote(if (alpha < 1e-300) {
    lomega <- -(alpha + 1) * log(alpha) - 0.9189385332046727
  } else {
    lomega <- (alpha + 1) * log1p(1 / alpha) - 0.9189385332046727
  })
  # Build each density with nimbleFunction() evaluated so its enclosure is the
  # global environment, then assign it there explicitly. (A namespace-frame
  # enclosure fails NIMBLE C++ codegen for these scalar densities.)
  makeIn <- function(expr) eval(expr, envir = ge)
  assign("dMSNBurr_k", makeIn(quote(nimble::nimbleFunction(
    run = function(x = double(0), mu = double(0), sigma = double(0),
                   alpha = double(0), log = integer(0, default = 0)) {
      returnType(double(0))
      if (alpha < 1e-300) {
        lomega <- -(alpha + 1) * log(alpha) - 0.9189385332046727
      } else {
        lomega <- (alpha + 1) * log1p(1 / alpha) - 0.9189385332046727
      }
      omega <- exp(lomega); zo <- -omega * ((x - mu) / sigma)
      u <- zo - log(alpha); sp <- max(u, 0) + log1p(exp(-abs(u)))
      lp <- lomega - log(sigma) + zo - (alpha + 1) * sp
      if (log) return(lp) else return(exp(lp))
    }))), envir = ge)
  assign("rMSNBurr_k", makeIn(quote(nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(0), sigma = double(0),
                   alpha = double(0)) {
      returnType(double(0))
      if (alpha < 1e-300) {
        lomega <- -(alpha + 1) * log(alpha) - 0.9189385332046727
      } else {
        lomega <- (alpha + 1) * log1p(1 / alpha) - 0.9189385332046727
      }
      omega <- exp(lomega); p <- runif(1)
      lt <- log(exp(-log(p) / alpha) - 1)
      return(mu - (sigma / omega) * (log(alpha) + lt))
    }))), envir = ge)
  assign("dMSNBurr2a_k", makeIn(quote(nimble::nimbleFunction(
    run = function(x = double(0), mu = double(0), sigma = double(0),
                   alpha = double(0), log = integer(0, default = 0)) {
      returnType(double(0))
      if (alpha < 1e-300) {
        lomega <- -(alpha + 1) * log(alpha) - 0.9189385332046727
      } else {
        lomega <- (alpha + 1) * log1p(1 / alpha) - 0.9189385332046727
      }
      omega <- exp(lomega); zt <- omega * ((x - mu) / sigma)
      u <- log(alpha) - zt; sp <- max(u, 0) + log1p(exp(-abs(u)))
      lp <- lomega - log(sigma) + (alpha + 1) * log(alpha) - alpha * zt -
        (alpha + 1) * sp
      if (log) return(lp) else return(exp(lp))
    }))), envir = ge)
  assign("dGMSNBurr_k", makeIn(quote(nimble::nimbleFunction(
    run = function(x = double(0), mu = double(0), sigma = double(0),
                   alpha = double(0), theta = double(0),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      lbeta_at <- lgamma(alpha) + lgamma(theta) - lgamma(alpha + theta)
      if (alpha < 1e-300) {
        lomega <- -0.9189385332046727 + lbeta_at +
          alpha * (log(theta) - log(alpha))
      } else {
        lomega <- -0.9189385332046727 + lbeta_at -
          theta * (log(theta) - log(alpha)) +
          (alpha + theta) * log1p(theta / alpha)
      }
      omega <- exp(lomega)
      zo <- -omega * ((x - mu) / sigma)
      zoa <- zo + log(theta) - log(alpha)
      sp <- max(zoa, 0) + log1p(exp(-abs(zoa)))
      lp <- lomega - log(sigma) + theta * (log(theta) - log(alpha)) +
        theta * zo - (alpha + theta) * sp - lbeta_at
      if (log) return(lp) else return(exp(lp))
    }))), envir = ge)
  assign("rGMSNBurr_k", makeIn(quote(nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(0), sigma = double(0),
                   alpha = double(0), theta = double(0)) {
      returnType(double(0))
      lbeta_at <- lgamma(alpha) + lgamma(theta) - lgamma(alpha + theta)
      if (alpha < 1e-300) {
        lomega <- -0.9189385332046727 + lbeta_at +
          alpha * (log(theta) - log(alpha))
      } else {
        lomega <- -0.9189385332046727 + lbeta_at -
          theta * (log(theta) - log(alpha)) +
          (alpha + theta) * log1p(theta / alpha)
      }
      omega <- exp(lomega)
      Xg <- rgamma(1, shape = alpha, rate = 1)
      Yg <- rgamma(1, shape = theta, rate = 1)
      log_ratio <- log(Yg) - log(Xg) + log(alpha) - log(theta)
      return(mu - (sigma / omega) * log_ratio)
    }))), envir = ge)
  # --- Batch B: SEP (symmetric exponential power) ---
  assign("dSEP_k", makeIn(quote(nimble::nimbleFunction(
    run = function(x = double(0), mu = double(0), sigma = double(0),
                   nu = double(0), log = integer(0, default = 0)) {
      returnType(double(0))
      z <- abs(x - mu) / sigma
      lp <- -log(2) - (1 / nu) * log(2) - lgamma(1 + 1 / nu) -
        log(sigma) - 0.5 * z^nu
      if (log) return(lp) else return(exp(lp))
    }))), envir = ge)
  assign("rSEP_k", makeIn(quote(nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(0), sigma = double(0),
                   nu = double(0)) {
      returnType(double(0))
      W <- rgamma(1, shape = 1 / nu, rate = 0.5)
      za <- W^(1 / nu)
      s <- 2 * (runif(1, 0, 1) < 0.5) - 1
      return(mu + sigma * s * za)
    }))), envir = ge)
  # --- Batch B: LEP (exponential power, alternative parameterisation) ---
  assign("dLEP_k", makeIn(quote(nimble::nimbleFunction(
    run = function(x = double(0), mu = double(0), sigma = double(0),
                   nu = double(0), log = integer(0, default = 0)) {
      returnType(double(0))
      z <- abs(x - mu) / sigma
      lp <- -log(2) - (1 / nu) * log(nu) - lgamma(1 + 1 / nu) -
        log(sigma) - z^nu / nu
      if (log) return(lp) else return(exp(lp))
    }))), envir = ge)
  assign("rLEP_k", makeIn(quote(nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(0), sigma = double(0),
                   nu = double(0)) {
      returnType(double(0))
      W <- rgamma(1, shape = 1 / nu, rate = 1 / nu)
      za <- W^(1 / nu)
      s <- 2 * (runif(1, 0, 1) < 0.5) - 1
      return(mu + sigma * s * za)
    }))), envir = ge)
  # --- Batch B: FSSN (Fernandez-Steel skew Normal) ---
  assign("dFSSN_k", makeIn(quote(nimble::nimbleFunction(
    run = function(x = double(0), mu = double(0), sigma = double(0),
                   alpha = double(0), log = integer(0, default = 0)) {
      returnType(double(0))
      z <- (x - mu) / sigma
      lc <- log(2) - log(sigma) - log(alpha + 1 / alpha) - 0.9189385332046727
      # FS convention: alpha == gamma, so alpha > 1 skews right.
      scale <- 1 / alpha
      if (z < 0) scale <- alpha
      lp <- lc - 0.5 * (z * scale)^2
      if (log) return(lp) else return(exp(lp))
    }))), envir = ge)
  assign("rFSSN_k", makeIn(quote(nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(0), sigma = double(0),
                   alpha = double(0)) {
      returnType(double(0))
      a2 <- alpha * alpha
      W <- abs(rnorm(1, 0, 1))
      if (runif(1, 0, 1) < a2 / (1 + a2)) {
        z <- W * alpha        # positive side, sd alpha
      } else {
        z <- -W / alpha       # negative side, sd 1/alpha
      }
      return(mu + sigma * z)
    }))), envir = ge)
  # --- Batch B: FOSSEP (Fernandez-Steel skew exponential power) ---
  assign("dFOSSEP_k", makeIn(quote(nimble::nimbleFunction(
    run = function(x = double(0), mu = double(0), sigma = double(0),
                   alpha = double(0), theta = double(0),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      z <- (x - mu) / sigma
      az <- abs(z)
      if (z < 0) {
        base <- -0.5 * (alpha * az)^theta
      } else {
        base <- -0.5 * (az / alpha)^theta
      }
      lp <- base - log(sigma) + log(alpha) - log1p(alpha^2) -
        (1 / theta) * log(2) - lgamma(1 + 1 / theta)
      if (log) return(lp) else return(exp(lp))
    }))), envir = ge)
  assign("rFOSSEP_k", makeIn(quote(nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(0), sigma = double(0),
                   alpha = double(0), theta = double(0)) {
      returnType(double(0))
      a2 <- alpha * alpha
      W <- rgamma(1, shape = 1 / theta, rate = 1)
      mag <- (2 * W)^(1 / theta)
      if (runif(1, 0, 1) < a2 / (1 + a2)) {
        z <- alpha * mag
      } else {
        z <- -mag / alpha
      }
      return(mu + sigma * z)
    }))), envir = ge)
  # --- Batch B: FSST (Fernandez-Steel skew Student-t; t-kernel inlined) ---
  assign("dFSST_k", makeIn(quote(nimble::nimbleFunction(
    run = function(x = double(0), mu = double(0), sigma = double(0),
                   alpha = double(0), nu = double(0),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      z <- (x - mu) / sigma
      # FS convention: alpha == gamma, so alpha > 1 skews right.
      tscale <- 1 / alpha
      if (z < 0) tscale <- alpha
      t <- z * tscale
      tlp <- lgamma((nu + 1) / 2) - lgamma(nu / 2) -
        0.5 * (log(nu) + 1.1447298858494002) -
        ((nu + 1) / 2) * log1p(t * t / nu)
      lp <- log(2) - log(sigma) - log(alpha + 1 / alpha) + tlp
      if (log) return(lp) else return(exp(lp))
    }))), envir = ge)
  assign("rFSST_k", makeIn(quote(nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(0), sigma = double(0),
                   alpha = double(0), nu = double(0)) {
      returnType(double(0))
      a2 <- alpha * alpha
      Zt <- rnorm(1, 0, 1)
      Wt <- rgamma(1, shape = nu / 2, rate = 0.5)
      Tt <- abs(Zt * sqrt(nu / Wt))
      if (runif(1, 0, 1) < a2 / (1 + a2)) {
        z <- Tt * alpha
      } else {
        z <- -Tt / alpha
      }
      return(mu + sigma * z)
    }))), envir = ge)
  # --- Batch B: JFST (Jones-Faddy skew-t; branch-free rz) ---
  assign("dJFST_k", makeIn(quote(nimble::nimbleFunction(
    run = function(x = double(0), mu = double(0), sigma = double(0),
                   alpha = double(0), theta = double(0),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      z <- (x - mu) / sigma
      rz <- z / sqrt(alpha + theta + z * z)
      lbeta_at <- lgamma(alpha) + lgamma(theta) - lgamma(alpha + theta)
      lp <- (alpha + 0.5) * log1p(rz) + (theta + 0.5) * log1p(-rz) -
        (alpha + theta - 1) * log(2) - 0.5 * log(alpha + theta) -
        lbeta_at - log(sigma)
      if (log) return(lp) else return(exp(lp))
    }))), envir = ge)
  assign("rJFST_k", makeIn(quote(nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(0), sigma = double(0),
                   alpha = double(0), theta = double(0)) {
      returnType(double(0))
      B <- rbeta(1, alpha, theta)
      rz <- 2 * B - 1
      z <- rz * sqrt((alpha + theta) / (1 - rz * rz))
      return(mu + sigma * z)
    }))), envir = ge)
  # --- Batch C: Ferreira-Steel skew multivariate Normal (A = chol(Sigma)) ---
  assign("dSkewMvN_k", makeIn(quote(nimble::nimbleFunction(
    run = function(x = double(1), mu = double(1), Sigma = double(2),
                   gam = double(1), log = integer(0, default = 0)) {
      returnType(double(0))
      m <- length(x)
      U <- chol(Sigma)
      Ui <- inverse(U)
      d <- x - mu
      lp <- 0
      for (j in 1:m) {
        s <- 0
        for (i in 1:m) s <- s + d[i] * Ui[i, j]
        sc <- 1 / gam[j]
        if (s < 0) sc <- gam[j]
        lp <- lp + log(2) - log(gam[j] + 1 / gam[j]) - 0.9189385332046727 -
          0.5 * (s * sc)^2 - log(U[j, j])
      }
      if (log) return(lp) else return(exp(lp))
    }))), envir = ge)
  assign("rSkewMvN_k", makeIn(quote(nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(1), Sigma = double(2),
                   gam = double(1)) {
      returnType(double(1))
      m <- length(mu)
      U <- chol(Sigma)
      eps <- numeric(m)
      for (j in 1:m) {
        a2 <- gam[j] * gam[j]
        W <- abs(rnorm(1, 0, 1))
        if (runif(1, 0, 1) < a2 / (1 + a2)) eps[j] <- W * gam[j]
        else eps[j] <- -W / gam[j]
      }
      out <- numeric(m)
      for (c in 1:m) {
        s <- 0
        for (r in 1:m) s <- s + U[r, c] * eps[r]
        out[c] <- mu[c] + s
      }
      return(out)
    }))), envir = ge)
  # --- Batch C: FS skew mv Normal with estimated orthogonal factor O (m = 2) ---
  # A = O U, U = chol(Sigma), O = I - 2 v v' the Householder reflection with
  # v = (sin theta, cos theta) (FS 2007 Appendix A). Note |O| = -1 always, so
  # O = I is NOT in the FS restricted set O_2; theta = 0 gives O = diag(1, -1),
  # which equals the O = I family with gamma_2 replaced by 1/gamma_2.
  assign("dSkewMvNO_k", makeIn(quote(nimble::nimbleFunction(
    run = function(x = double(1), mu = double(1), Sigma = double(2),
                   gam = double(1), theta = double(0),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      m <- length(x)
      U <- chol(Sigma)
      Ui <- inverse(U)
      dv <- x - mu
      w <- numeric(m)
      for (j in 1:m) {
        s <- 0
        for (i in 1:m) s <- s + dv[i] * Ui[i, j]
        w[j] <- s
      }
      s1 <- sin(theta)
      c1 <- cos(theta)
      vw <- s1 * w[1] + c1 * w[2]
      eps <- numeric(2)
      eps[1] <- w[1] - 2 * s1 * vw
      eps[2] <- w[2] - 2 * c1 * vw
      lp <- 0
      for (j in 1:m) {
        sc <- 1 / gam[j]
        if (eps[j] < 0) sc <- gam[j]
        lp <- lp + log(2) - log(gam[j] + 1 / gam[j]) - 0.9189385332046727 -
          0.5 * (eps[j] * sc)^2 - log(U[j, j])
      }
      if (log) return(lp) else return(exp(lp))
    }))), envir = ge)
  assign("rSkewMvNO_k", makeIn(quote(nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(1), Sigma = double(2),
                   gam = double(1), theta = double(0)) {
      returnType(double(1))
      m <- length(mu)
      U <- chol(Sigma)
      eps <- numeric(m)
      for (j in 1:m) {
        a2 <- gam[j] * gam[j]
        W <- abs(rnorm(1, 0, 1))
        if (runif(1, 0, 1) < a2 / (1 + a2)) eps[j] <- W * gam[j]
        else eps[j] <- -W / gam[j]
      }
      # eta = A' eps + mu = U' O' eps + mu ; O symmetric so O' = O
      s1 <- sin(theta)
      c1 <- cos(theta)
      ve <- s1 * eps[1] + c1 * eps[2]
      oe <- numeric(2)
      oe[1] <- eps[1] - 2 * s1 * ve
      oe[2] <- eps[2] - 2 * c1 * ve
      out <- numeric(m)
      for (c in 1:m) {
        s <- 0
        for (r in 1:m) s <- s + U[r, c] * oe[r]
        out[c] <- mu[c] + s
      }
      return(out)
    }))), envir = ge)
  # --- Batch C: FS skew mv Normal with estimated O, general m ---
  # O = O_{th^m} ... O_{th^2} applied to w = (U')^{-1}(x - mu); the blocks are
  # applied j = 2, ..., m, which reproduces the matrix product (verified against
  # the R reference and against the m = 2 kernel).
  assign("dSkewMvNOG_k", makeIn(quote(nimble::nimbleFunction(
    run = function(x = double(1), mu = double(1), Sigma = double(2),
                   gam = double(1), theta = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      m <- length(x)
      U <- chol(Sigma)
      Ui <- inverse(U)
      dv <- x - mu
      w <- numeric(m)
      for (j in 1:m) {
        s <- 0
        for (i in 1:m) s <- s + dv[i] * Ui[i, j]
        w[j] <- s
      }
      eps <- numeric(m)
      for (i in 1:m) eps[i] <- w[i]
      for (j in 2:m) {
        start <- 1 + (j - 2) * (j - 1) / 2
        v <- numeric(j)
        cp <- 1
        v[1] <- sin(theta[start])
        if (j > 2) {
          for (i in 2:(j - 1)) {
            cp <- cp * cos(theta[start + i - 2])
            v[i] <- cp * sin(theta[start + i - 1])
          }
          cp <- cp * cos(theta[start + j - 2])
        } else {
          cp <- cos(theta[start])
        }
        v[j] <- cp
        off <- m - j
        ve <- 0
        for (i in 1:j) ve <- ve + v[i] * eps[off + i]
        vv <- 0
        for (i in 1:j) vv <- vv + v[i] * v[i]
        for (i in 1:j) eps[off + i] <- eps[off + i] - 2 * v[i] * ve / vv
      }
      lp <- 0
      for (j in 1:m) {
        sc <- 1 / gam[j]
        if (eps[j] < 0) sc <- gam[j]
        lp <- lp + log(2) - log(gam[j] + 1 / gam[j]) - 0.9189385332046727 -
          0.5 * (eps[j] * sc)^2 - log(U[j, j])
      }
      if (log) return(lp) else return(exp(lp))
    }))), envir = ge)
  # --- Batch C: FS skew mv independent-Student, estimated O, general m ---
  assign("dSkewMvITOG_k", makeIn(quote(nimble::nimbleFunction(
    run = function(x = double(1), mu = double(1), Sigma = double(2),
                   gam = double(1), nu = double(1), theta = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      m <- length(x)
      U <- chol(Sigma)
      Ui <- inverse(U)
      dv <- x - mu
      w <- numeric(m)
      for (j in 1:m) {
        s <- 0
        for (i in 1:m) s <- s + dv[i] * Ui[i, j]
        w[j] <- s
      }
      eps <- numeric(m)
      for (i in 1:m) eps[i] <- w[i]
      for (j in 2:m) {
        start <- 1 + (j - 2) * (j - 1) / 2
        v <- numeric(j)
        cp <- 1
        v[1] <- sin(theta[start])
        if (j > 2) {
          for (i in 2:(j - 1)) {
            cp <- cp * cos(theta[start + i - 2])
            v[i] <- cp * sin(theta[start + i - 1])
          }
          cp <- cp * cos(theta[start + j - 2])
        } else {
          cp <- cos(theta[start])
        }
        v[j] <- cp
        off <- m - j
        ve <- 0
        for (i in 1:j) ve <- ve + v[i] * eps[off + i]
        vv <- 0
        for (i in 1:j) vv <- vv + v[i] * v[i]
        for (i in 1:j) eps[off + i] <- eps[off + i] - 2 * v[i] * ve / vv
      }
      lp <- 0
      for (j in 1:m) {
        sc <- 1 / gam[j]
        if (eps[j] < 0) sc <- gam[j]
        t <- eps[j] * sc
        tlp <- lgamma((nu[j] + 1) / 2) - lgamma(nu[j] / 2) -
          0.5 * (log(nu[j]) + 1.1447298858494002) -
          ((nu[j] + 1) / 2) * log1p(t * t / nu[j])
        lp <- lp + log(2) - log(gam[j] + 1 / gam[j]) + tlp - log(U[j, j])
      }
      if (log) return(lp) else return(exp(lp))
    }))), envir = ge)
  # --- Batch C: FS skew mv independent-Student with estimated O (m = 2) ---
  assign("dSkewMvITO_k", makeIn(quote(nimble::nimbleFunction(
    run = function(x = double(1), mu = double(1), Sigma = double(2),
                   gam = double(1), nu = double(1), theta = double(0),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      m <- length(x)
      U <- chol(Sigma)
      Ui <- inverse(U)
      dv <- x - mu
      w <- numeric(m)
      for (j in 1:m) {
        s <- 0
        for (i in 1:m) s <- s + dv[i] * Ui[i, j]
        w[j] <- s
      }
      s1 <- sin(theta)
      c1 <- cos(theta)
      vw <- s1 * w[1] + c1 * w[2]
      eps <- numeric(2)
      eps[1] <- w[1] - 2 * s1 * vw
      eps[2] <- w[2] - 2 * c1 * vw
      lp <- 0
      for (j in 1:m) {
        sc <- 1 / gam[j]
        if (eps[j] < 0) sc <- gam[j]
        t <- eps[j] * sc
        tlp <- lgamma((nu[j] + 1) / 2) - lgamma(nu[j] / 2) -
          0.5 * (log(nu[j]) + 1.1447298858494002) -
          ((nu[j] + 1) / 2) * log1p(t * t / nu[j])
        lp <- lp + log(2) - log(gam[j] + 1 / gam[j]) + tlp - log(U[j, j])
      }
      if (log) return(lp) else return(exp(lp))
    }))), envir = ge)
  assign("rSkewMvITO_k", makeIn(quote(nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(1), Sigma = double(2),
                   gam = double(1), nu = double(1), theta = double(0)) {
      returnType(double(1))
      m <- length(mu)
      U <- chol(Sigma)
      eps <- numeric(m)
      for (j in 1:m) {
        a2 <- gam[j] * gam[j]
        Zt <- rnorm(1, 0, 1)
        Wt <- rgamma(1, shape = nu[j] / 2, rate = 0.5)
        Tt <- abs(Zt * sqrt(nu[j] / Wt))
        if (runif(1, 0, 1) < a2 / (1 + a2)) eps[j] <- Tt * gam[j]
        else eps[j] <- -Tt / gam[j]
      }
      s1 <- sin(theta)
      c1 <- cos(theta)
      ve <- s1 * eps[1] + c1 * eps[2]
      oe <- numeric(2)
      oe[1] <- eps[1] - 2 * s1 * ve
      oe[2] <- eps[2] - 2 * c1 * ve
      out <- numeric(m)
      for (c in 1:m) {
        s <- 0
        for (r in 1:m) s <- s + U[r, c] * oe[r]
        out[c] <- mu[c] + s
      }
      return(out)
    }))), envir = ge)
  # --- Batch C: FS skew multivariate independent-Student (t-kernel inlined) ---
  assign("dSkewMvIT_k", makeIn(quote(nimble::nimbleFunction(
    run = function(x = double(1), mu = double(1), Sigma = double(2),
                   gam = double(1), nu = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      m <- length(x)
      U <- chol(Sigma)
      Ui <- inverse(U)
      d <- x - mu
      lp <- 0
      for (j in 1:m) {
        s <- 0
        for (i in 1:m) s <- s + d[i] * Ui[i, j]
        sc <- 1 / gam[j]
        if (s < 0) sc <- gam[j]
        t <- s * sc
        tlp <- lgamma((nu[j] + 1) / 2) - lgamma(nu[j] / 2) -
          0.5 * (log(nu[j]) + 1.1447298858494002) -
          ((nu[j] + 1) / 2) * log1p(t * t / nu[j])
        lp <- lp + log(2) - log(gam[j] + 1 / gam[j]) + tlp - log(U[j, j])
      }
      if (log) return(lp) else return(exp(lp))
    }))), envir = ge)
  assign("rSkewMvIT_k", makeIn(quote(nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(1), Sigma = double(2),
                   gam = double(1), nu = double(1)) {
      returnType(double(1))
      m <- length(mu)
      U <- chol(Sigma)
      eps <- numeric(m)
      for (j in 1:m) {
        a2 <- gam[j] * gam[j]
        Zt <- rnorm(1, 0, 1)
        Wt <- rgamma(1, shape = nu[j] / 2, rate = 0.5)
        Tt <- abs(Zt * sqrt(nu[j] / Wt))
        if (runif(1, 0, 1) < a2 / (1 + a2)) eps[j] <- Tt * gam[j]
        else eps[j] <- -Tt / gam[j]
      }
      out <- numeric(m)
      for (c in 1:m) {
        s <- 0
        for (r in 1:m) s <- s + U[r, c] * eps[r]
        out[c] <- mu[c] + s
      }
      return(out)
    }))), envir = ge)
  assign("rMSNBurr2a_k", makeIn(quote(nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(0), sigma = double(0),
                   alpha = double(0)) {
      returnType(double(0))
      if (alpha < 1e-300) {
        lomega <- -(alpha + 1) * log(alpha) - 0.9189385332046727
      } else {
        lomega <- (alpha + 1) * log1p(1 / alpha) - 0.9189385332046727
      }
      omega <- exp(lomega); q <- 1 - runif(1)
      lt <- log(exp(-log(q) / alpha) - 1)
      return(mu + (sigma / omega) * (log(alpha) + lt))
    }))), envir = ge)
  invisible()
}

.onLoad <- function(libname, pkgname) {
  # The MSNBurr / MSNBurr-IIa NIMBLE densities MUST be created in the global
  # environment, not as namespace objects. A nimbleFunction created inside a
  # non-global frame (a package namespace, or any enclosing function) fails
  # NIMBLE C++ code generation for scalar user-defined distributions with
  # "argument is of length zero"; the same body created at top level compiles.
  # We therefore build them once here via eval() in globalenv(). Branch-free
  # softplus for stability. Iriawan (2000); Choir (2020).
  # Register built-ins. nimixClust() resolves "normal" to the univariate or
  # multivariate spec by data shape (see .selectClusterSpec); the "normal"
  # alias below is the univariate default for direct getDistribution() calls.
  assign("normal-uv", NormalUvSpec(), envir = .distRegistry)
  assign("normal-mv", NormalMvSpec(), envir = .distRegistry)
  assign("normal-reg", NormalRegSpec(), envir = .distRegistry)
  assign("student-t", StudentTUvSpec(), envir = .distRegistry)
  assign("normal-gamma", NormalGammaUvSpec(), envir = .distRegistry)
  assign("student-t-mv", StudentTMvSpec(), envir = .distRegistry)
  assign("normal-gamma-mv", NormalGammaMvSpec(), envir = .distRegistry)
  assign("poisson", PoissonSpec(), envir = .distRegistry)
  assign("binomial", BinomialSpec(), envir = .distRegistry)
  assign("poisson-reg", PoissonRegSpec(), envir = .distRegistry)
  assign("binomial-reg", BinomialRegSpec(), envir = .distRegistry)
  assign("msnburr-reg", MSNBurrRegSpec(), envir = .distRegistry)
  assign("sep-reg", SEPRegSpec(), envir = .distRegistry)
  assign("msnburr2a-reg", MSNBurr2aRegSpec(), envir = .distRegistry)
  assign("fssn-reg", FSSNRegSpec(), envir = .distRegistry)
  assign("gmsnburr-reg", GMSNBurrRegSpec(), envir = .distRegistry)
  assign("lep-reg", LEPRegSpec(), envir = .distRegistry)
  assign("fsst-reg", FSSTRegSpec(), envir = .distRegistry)
  assign("fossep-reg", FOSSEPRegSpec(), envir = .distRegistry)
  assign("jfst-reg", JFSTRegSpec(), envir = .distRegistry)
  assign("student-t-reg", StudentTRegSpec(), envir = .distRegistry)
  assign("normal-gamma-reg", NormalGammaRegSpec(), envir = .distRegistry)
  assign("normal-mv-reg", NormalMvRegSpec(), envir = .distRegistry)
  assign("student-t-mv-reg", StudentTMvRegSpec(), envir = .distRegistry)
  assign("normal-gamma-mv-reg", NormalGammaMvRegSpec(), envir = .distRegistry)
  assign("msnburr", MSNBurrUvSpec(), envir = .distRegistry)
  assign("msnburr2a", MSNBurr2aUvSpec(), envir = .distRegistry)
  assign("gmsnburr", GMSNBurrUvSpec(), envir = .distRegistry)
  assign("sep", SEPUvSpec(), envir = .distRegistry)
  assign("lep", LEPUvSpec(), envir = .distRegistry)
  assign("fssn", FSSNUvSpec(), envir = .distRegistry)
  assign("fossep", FOSSEPUvSpec(), envir = .distRegistry)
  assign("fsst", FSSTUvSpec(), envir = .distRegistry)
  assign("jfst", JFSTUvSpec(), envir = .distRegistry)
  assign("skewnormal-mv", SkewNormalMvSpec(), envir = .distRegistry)
  assign("skewistudent-mv", SkewIStudentMvSpec(), envir = .distRegistry)
  assign("skewnormal-mv-o", SkewNormalMvOSpec(), envir = .distRegistry)
  assign("skewistudent-mv-o", SkewIStudentMvOSpec(), envir = .distRegistry)
  assign("skewnormal-mv-og", SkewNormalMvOGenSpec(), envir = .distRegistry)
  assign("skewistudent-mv-og", SkewIStudentMvOGenSpec(), envir = .distRegistry)
  assign("normal", NormalUvSpec(), envir = .distRegistry)
  # Register the user-defined multivariate-t density with NIMBLE so the
  # StudentTMvSpec kernel resolves at model-build time.
  suppressMessages(suppressWarnings(try(
    nimble::registerDistributions(list(
      dmvt_nimix = list(
        BUGSdist = "dmvt_nimix(mu, cov, df)",
        types = c("value = double(1)", "mu = double(1)",
                  "cov = double(2)", "df = double(0)")))),
    silent = TRUE)))
  # The unnormalised Potts prior for the MRF engine is built and registered
  # lazily in globalenv by .nimixDefinePotts()/.nimixEnsureMSNBurr(); doing it
  # here (namespace frame) makes NIMBLE fail to find rPottsNimix at code-gen.
  invisible(NULL)
}

# Lazily ensure the MSNBurr densities exist in the global environment AND are
# registered with NIMBLE. Registration is deferred to first use (not .onLoad)
# because a registration performed while the objects are being (re)built during
# package load can bind a distribution name to a namespace-frame object that
# fails C++ code generation; binding the name once, here, to the global-frame
# objects avoids that. Idempotent and cheap.
.nimixEnsureMSNBurr <- function() {
  .nimixDefineMSNBurr()
  .nimixDefinePotts()
  if (isTRUE(.nimixState$msnburrRegistered)) return(invisible())
  eval(quote(suppressMessages(suppressWarnings(try(nimble::registerDistributions(list(
    dMSNBurr_k = list(
      BUGSdist = "dMSNBurr_k(mu, sigma, alpha)",
      types = c("value = double(0)", "mu = double(0)",
                "sigma = double(0)", "alpha = double(0)"),
      discrete = FALSE),
    dMSNBurr2a_k = list(
      BUGSdist = "dMSNBurr2a_k(mu, sigma, alpha)",
      types = c("value = double(0)", "mu = double(0)",
                "sigma = double(0)", "alpha = double(0)"),
      discrete = FALSE),
    dGMSNBurr_k = list(
      BUGSdist = "dGMSNBurr_k(mu, sigma, alpha, theta)",
      types = c("value = double(0)", "mu = double(0)", "sigma = double(0)",
                "alpha = double(0)", "theta = double(0)"),
      discrete = FALSE),
    dPottsNimix = list(
      BUGSdist = "dPottsNimix(beta, e1, e2)",
      types = c("value = double(1)", "beta = double(0)",
                "e1 = double(1)", "e2 = double(1)"),
      discrete = TRUE, mixedSizes = TRUE),
    dSEP_k = list(
      BUGSdist = "dSEP_k(mu, sigma, nu)",
      types = c("value = double(0)", "mu = double(0)", "sigma = double(0)",
                "nu = double(0)"), discrete = FALSE),
    dLEP_k = list(
      BUGSdist = "dLEP_k(mu, sigma, nu)",
      types = c("value = double(0)", "mu = double(0)", "sigma = double(0)",
                "nu = double(0)"), discrete = FALSE),
    dFSSN_k = list(
      BUGSdist = "dFSSN_k(mu, sigma, alpha)",
      types = c("value = double(0)", "mu = double(0)", "sigma = double(0)",
                "alpha = double(0)"), discrete = FALSE),
    dFOSSEP_k = list(
      BUGSdist = "dFOSSEP_k(mu, sigma, alpha, theta)",
      types = c("value = double(0)", "mu = double(0)", "sigma = double(0)",
                "alpha = double(0)", "theta = double(0)"), discrete = FALSE),
    dFSST_k = list(
      BUGSdist = "dFSST_k(mu, sigma, alpha, nu)",
      types = c("value = double(0)", "mu = double(0)", "sigma = double(0)",
                "alpha = double(0)", "nu = double(0)"), discrete = FALSE),
    dJFST_k = list(
      BUGSdist = "dJFST_k(mu, sigma, alpha, theta)",
      types = c("value = double(0)", "mu = double(0)", "sigma = double(0)",
                "alpha = double(0)", "theta = double(0)"),
      discrete = FALSE),
    dSkewMvN_k = list(
      BUGSdist = "dSkewMvN_k(mu, Sigma, gam)",
      types = c("value = double(1)", "mu = double(1)", "Sigma = double(2)",
                "gam = double(1)"),
      discrete = FALSE),
    dSkewMvIT_k = list(
      BUGSdist = "dSkewMvIT_k(mu, Sigma, gam, nu)",
      types = c("value = double(1)", "mu = double(1)", "Sigma = double(2)",
                "gam = double(1)", "nu = double(1)"),
      discrete = FALSE),
    dSkewMvNO_k = list(
      BUGSdist = "dSkewMvNO_k(mu, Sigma, gam, theta)",
      types = c("value = double(1)", "mu = double(1)", "Sigma = double(2)",
                "gam = double(1)", "theta = double(0)"),
      discrete = FALSE),
    dSkewMvITO_k = list(
      BUGSdist = "dSkewMvITO_k(mu, Sigma, gam, nu, theta)",
      types = c("value = double(1)", "mu = double(1)", "Sigma = double(2)",
                "gam = double(1)", "nu = double(1)", "theta = double(0)"),
      discrete = FALSE),
    dSkewMvNOG_k = list(
      BUGSdist = "dSkewMvNOG_k(mu, Sigma, gam, theta)",
      types = c("value = double(1)", "mu = double(1)", "Sigma = double(2)",
                "gam = double(1)", "theta = double(1)"),
      discrete = FALSE),
    dSkewMvITOG_k = list(
      BUGSdist = "dSkewMvITOG_k(mu, Sigma, gam, nu, theta)",
      types = c("value = double(1)", "mu = double(1)", "Sigma = double(2)",
                "gam = double(1)", "nu = double(1)", "theta = double(1)"),
      discrete = FALSE))), silent = TRUE)))),
    envir = globalenv())
  .nimixState$msnburrRegistered <- TRUE
  invisible()
}

# Build the Potts prior's d/r functions in the GLOBAL environment. Registering
# them from a namespace frame makes NIMBLE fail to find rPottsNimix during code
# generation for the latent label node (it is invoked to simulate z), the same
# class of failure that affects the scalar neo-normal densities. Building here,
# in globalenv, resolves it.
.nimixDefinePotts <- function() {
  if (exists("rPottsNimix", envir = globalenv(), inherits = FALSE))
    return(invisible())
  ge <- globalenv()
  assign("dPottsNimix", eval(quote(nimble::nimbleFunction(
    run = function(x = double(1), beta = double(0),
                   e1 = double(1), e2 = double(1),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      s <- 0
      nE <- length(e1)
      for (m in 1:nE) if (x[e1[m]] == x[e2[m]]) s <- s + 1
      lp <- beta * s
      if (log) return(lp) else return(exp(lp))
    })), envir = ge), envir = ge)
  assign("rPottsNimix", eval(quote(nimble::nimbleFunction(
    run = function(n = integer(0), beta = double(0),
                   e1 = double(1), e2 = double(1)) {
      returnType(double(1))
      ## Labels are always supplied as inits and updated by the Gibbs sweep;
      ## exact Potts simulation is not needed for inference.
      out <- numeric(length = 1)
      return(out)
    })), envir = ge), envir = ge)
  invisible()
}

.nimixState <- new.env(parent = emptyenv())
