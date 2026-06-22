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

.onLoad <- function(libname, pkgname) {
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
  invisible(NULL)
}
