## ---------------------------------------------------------------------------
## nimixReg.R
##
## User-facing entry point for a DPM mixture of linear regressions (v0.3.0).
## Baseline gating is CONSTANT (mixing weights from the CRP; not covariate
## dependent) --. Covariate-dependent
## (concomitant) gating is an explicit future opt-in and currently errors with a
## pointer, rather than silently enabling a harder-to-identify model.
## ---------------------------------------------------------------------------

#' Bayesian mixture of linear regressions
#'
#' Fit a mixture of Gaussian linear regressions. Each component is
#' \eqn{y \sim N(x^\top \beta_k, \sigma^2_k)} with a conjugate
#' Normal-Inverse-Gamma cluster prior. The number of components can be inferred
#' with a Dirichlet Process Mixture (\code{method = "dpm"}, using \code{K_max})
#' or fixed (\code{method = "fixedk"}, using \code{K}). Gating is constant: the
#' mixing weights do not depend on the covariates.
#'
#' @param formula A model formula, e.g. \code{y ~ x1 + x2}.
#' @param data A data frame containing the formula variables.
#' @param K Integer number of components for the finite mixture
#'   (\code{method = "fixedk"}). Required for that method.
#' @param K_max Integer truncation level for the DPM (\code{method = "dpm"});
#'   it must sit comfortably above the expected number of clusters, since the
#'   dCRP sampler errors if the occupied-cluster count ever needs to exceed it.
#'   A generous data-aware default (with headroom) is used when missing.
#' @param distribution Component distribution. Currently \code{"normal"}
#'   (Gaussian linear component). Other GLM families are planned.
#' @param method Engine: \code{"dpm"} (default; estimate the number of
#'   components), \code{"fixedk"} (fixed \code{K}), or \code{"mrf"}
#'   (spatially constrained: fixed \code{K}, Potts smoothing of the labels on
#'   a \code{spatialWeights} graph; Gaussian response, fixed
#'   \code{prior$beta}, default 0.8).
#' @param spatialWeights Optional \code{\linkS4class{SpatialWeightSpec}} (one
#'   region per observation). Required by, and only used with,
#'   \code{method = "mrf"}.
#' @param gating Mixing-weight model: \code{"constant"} (default; weights do not
#'   depend on covariates). \code{"covariate"} (concomitant gating) is a planned
#'   opt-in and currently errors.
#' @param prior A named list of prior overrides passed to
#'   \code{\link{defaultPrior}} (e.g. \code{g} for the g-prior factor,
#'   \code{nu0} for the InvGamma shape) plus, for the DPM, optional
#'   \code{concPrior = c(shape, rate)}, or, for the finite mixture,
#'   \code{dirichletConc}.
#' @param mcmcControl A named list of MCMC controls: \code{niter},
#'   \code{nburnin}, \code{thin}, and the optional \code{initRatio} -- the
#'   fraction of the truncation / component cap (\code{K_max} or \code{K})
#'   seeded by the dispersed cluster initialisation (default 0.8; must lie in
#'   (0, 1)). Lower it to leave more headroom below the truncation; raising it
#'   to 0.95 or above is allowed but warns, as it leaves little headroom.
#' @param initMethod Initialisation: \code{"kmeans"} (default) or
#'   \code{"single"}.
#' @param seed Integer RNG seed.
#' @param verbose Logical; print NIMBLE's configuration and progress output.
#'   Defaults to \code{FALSE} (quiet): NIMBLE's compilation notes and the benign
#'   dCRP truncation note are silenced, while nimix's own diagnostics (e.g. a
#'   censored-posterior warning) and any error still surface. Set \code{TRUE} to
#'   see NIMBLE's configuration and a progress bar.
#'
#' @return A \code{\linkS4class{FitResult}}. \code{summary()} reports relabelled
#'   per-component regression coefficients and residual variances;
#'   \code{predict(fit, newdata)} returns the posterior predictive mean;
#'   \code{plot(fit, type = "fitted")} shows observed vs fitted.
#'
#' @references
#' Hurn, M., Justel, A., & Robert, C.P. (2003). Estimating mixtures of
#' regressions. \emph{JCGS}, 12(1), 55--79. \doi{10.1198/1061860031329}
#'
#' Grün, B., & Leisch, F. (2008). FlexMix version 2. \emph{JSS}, 28(4), 1--35.
#' \doi{10.18637/jss.v028.i04}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' x <- runif(200, -3, 3)
#' grp <- rep(1:2, each = 100)
#' y <- ifelse(grp == 1, 2 * x, -2 * x) + rnorm(200, 0, 0.7)
#' df <- data.frame(y = y, x = x)
#'
#' ## number of regimes estimated (DPM)
#' fit <- nimixReg(y ~ x, df, K_max = 8,
#'                 mcmcControl = list(niter = 2000, nburnin = 1000),
#'                 verbose = FALSE)
#' summary(fit)
#' predict(fit, newdata = data.frame(x = c(-2, 0, 2)))
#'
#' ## fixed number of regimes (finite mixture)
#' fit2 <- nimixReg(y ~ x, df, K = 2, method = "fixedk",
#'                  mcmcControl = list(niter = 2000, nburnin = 1000),
#'                  verbose = FALSE)
#' summary(fit2)
#' }
#' @export
nimixReg <- function(formula, data,
                     K = NULL,
                     K_max = NULL,
                     distribution = "normal",
                     method = c("dpm", "fixedk", "mrf"),
                     gating = c("constant", "covariate"),
                     prior = list(),
                     mcmcControl = list(),
                     initMethod = c("kmeans", "single"),
                     seed = 1L,
                     verbose = FALSE,
                     spatialWeights = NULL) {
  cl <- match.call()
  method <- match.arg(method)
  gating <- match.arg(gating)
  initMethod <- match.arg(initMethod)

  if (!is.null(spatialWeights) && !methods::is(spatialWeights, "SpatialWeightSpec"))
    stop("spatialWeights must be a SpatialWeightSpec (see ?spatialWeights).",
         call. = FALSE)
  if (method == "mrf" && is.null(spatialWeights))
    stop("method = 'mrf' needs a spatialWeights neighbourhood ",
         "(see ?spatialWeights, ?gridAdjacency).", call. = FALSE)
  if (method != "mrf" && !is.null(spatialWeights))
    stop("spatialWeights is only used by method = 'mrf'; with method = '",
         method, "' leave it NULL.", call. = FALSE)

  # --- scope guards --------------------------------------------------------
  if (gating == "covariate")
    stop("gating = 'covariate' (concomitant / covariate-dependent gating) is a ",
         "planned opt-in and is not yet implemented. The default constant ",
         "gating avoids the extra identifiability risk of letting the mixing ",
         "weights depend on the same covariates as the response.",
         call. = FALSE)
  distribution <- tolower(distribution)

  # K is for the finite mixture (fixed components); K_max is the DPM truncation.
  if (method %in% c("fixedk", "mrf")) {
    if (is.null(K))
      stop("method = '", method, "' needs the number of components K.", call. = FALSE)
    if (!is.null(K_max))
      stop("Use K (not K_max) with method = '", method, "'.", call. = FALSE)
  } else if (!is.null(K)) {
    stop("Use K_max (not K) with method = 'dpm'.", call. = FALSE)
  }

  # --- parse the formula / build the design matrix -------------------------
  mf <- stats::model.frame(formula, data)
  y  <- stats::model.response(mf)
  if (is.null(y)) stop("The formula must have a response (left-hand side).",
                       call. = FALSE)
  isMvResp <- is.matrix(y) && ncol(y) > 1L
  if (!isMvResp) y <- as.numeric(y)
  tt <- stats::terms(mf)
  X  <- stats::model.matrix(tt, mf)
  if (anyNA(y) || anyNA(X))
    stop("Missing values are not supported; please remove NAs.", call. = FALSE)

  regName <- if (isMvResp) switch(distribution,
      normal = "normal-mv-reg",
      studentt = , "student-t" = , t = "student-t-mv-reg",
      normalgamma = , "normal-gamma" = "normal-gamma-mv-reg",
      NA_character_)
    else switch(distribution,
      normal = "normal-reg", poisson = "poisson-reg", binomial = "binomial-reg",
      studentt = , "student-t" = , t = "student-t-reg",
      normalgamma = , "normal-gamma" = "normal-gamma-reg",
      NA_character_)
  if (is.na(regName))
    stop("distribution = '", distribution, "' is not available for nimixReg",
         if (isMvResp) " with a multivariate response" else "", "; use ",
         if (isMvResp) "'normal', 'studentt', or 'normalgamma'."
         else "'normal', 'studentt'/'normalgamma', 'poisson', or 'binomial'.",
         call. = FALSE)

  # Count families need a non-negative integer response; reject continuous y
  # early with a clear message instead of letting it surface as a numerical
  # underflow deep inside the CRP/categorical sampler.
  if (regName %in% c("poisson-reg", "binomial-reg")) {
    yv <- if (isMvResp) as.numeric(y) else y
    if (any(yv < 0) || any(abs(yv - round(yv)) > 1e-8))
      stop("distribution = '", distribution, "' needs a non-negative integer ",
           "count response.", call. = FALSE)
  }

  n <- if (isMvResp) nrow(y) else length(y); p <- ncol(X)
  if (n < 2L) stop("Need at least 2 observations.", call. = FALSE)
  if (n <= p) stop("Need more observations than coefficients (n = ", n,
                   ", p = ", p, ").", call. = FALSE)

  spec <- getDistribution(regName)

  if (method %in% c("fixedk", "mrf")) {
    nComp <- as.integer(K)
    if (nComp < 1L) stop("K must be >= 1.", call. = FALSE)
  } else {
    if (is.null(K_max)) K_max <- .defaultTruncation(n)
    nComp <- as.integer(K_max)
    if (nComp < 2L) stop("K_max must be >= 2.", call. = FALSE)
  }
  if (nComp > n) stop("The number of components cannot exceed n.", call. = FALSE)

  # --- prior + engine ------------------------------------------------------
  priorList <- defaultPrior(spec, y, control = c(prior, list(X = X)))
  validateParams(spec, priorList)
  # carry items needed for labelling and prediction
  priorList$coefNames <- colnames(X)
  priorList$terms     <- tt
  priorList$formula   <- formula
  if (isMvResp && is.null(priorList$respNames))
    priorList$respNames <- colnames(y)

  if (method == "fixedk") {
    dConc <- if (!is.null(prior$dirichletConc)) prior$dirichletConc else 1
    engine <- FixedKEngine(dirichletConc = dConc)
  } else if (method == "mrf") {
    if (nrow(getAdjacency(spatialWeights)) != n)
      stop("spatialWeights has ", nrow(getAdjacency(spatialWeights)),
           " regions but the data has ", n, " observations; they must match.",
           call. = FALSE)
    mrfBeta <- if (!is.null(prior$beta)) prior$beta else 0.8
    engine <- MRFEngine(beta = mrfBeta, spatial = spatialWeights,
                        estimateBeta = isTRUE(prior$estimateBeta),
                        betaMax = if (!is.null(prior$betaMax)) prior$betaMax else 2)
  } else {
    concPrior <- if (!is.null(prior$concPrior)) prior$concPrior else c(2, 4)
    engine <- DPMEngine(concPrior = concPrior)
  }

  model <- RegressionMixModel(data = y, X = X, formula = formula,
                              distSpec = spec, engine = engine,
                              Kmax = nComp, prior = priorList)

  # --- run -----------------------------------------------------------------
  raw <- runEngine(engine, model, mcmcControl = mcmcControl,
                   initMethod = initMethod, seed = seed, verbose = verbose)

  new("FitResult",
      mcmcSamples       = raw$mcmcSamples,
      Kposterior        = raw$Kposterior,
      clusterAllocation = raw$clusterAllocation,
      paramTrace        = raw$paramTrace,
      engineUsed        = engine@name,
      distSpec          = spec,
      data              = y,
      Kmax              = nComp,
      prior             = priorList,
      relabeled         = list(),
      mcmcControl       = raw$mcmcControl,
      diagnostics       = raw$diagnostics,
      call              = cl)
}
