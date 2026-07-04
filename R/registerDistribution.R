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
  if (exists("dMSNBurr_k", envir = globalenv(), inherits = FALSE)) return(invisible())
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
  assign("student-t-reg", StudentTRegSpec(), envir = .distRegistry)
  assign("normal-gamma-reg", NormalGammaRegSpec(), envir = .distRegistry)
  assign("normal-mv-reg", NormalMvRegSpec(), envir = .distRegistry)
  assign("student-t-mv-reg", StudentTMvRegSpec(), envir = .distRegistry)
  assign("normal-gamma-mv-reg", NormalGammaMvRegSpec(), envir = .distRegistry)
  assign("msnburr", MSNBurrUvSpec(), envir = .distRegistry)
  assign("msnburr2a", MSNBurr2aUvSpec(), envir = .distRegistry)
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
  # Register the unnormalised Potts prior for the MRF engine (valid for MCMC
  # because beta is fixed; see engine-mrf.R header).
  suppressMessages(suppressWarnings(try(
    nimble::registerDistributions(list(
      dPottsNimix = list(
        BUGSdist = "dPottsNimix(beta, e1, e2)",
        types = c("value = double(1)", "beta = double(0)",
                  "e1 = double(1)", "e2 = double(1)"),
        discrete = TRUE, mixedSizes = TRUE))),
    silent = TRUE)))
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
      discrete = FALSE))), silent = TRUE)))), envir = globalenv())
  .nimixState$msnburrRegistered <- TRUE
  invisible()
}

.nimixState <- new.env(parent = emptyenv())
