## ---------------------------------------------------------------------------
## nimixClust.R
##
## User-facing entry point for mixture clustering.
##   v0.1.0: univariate Gaussian, DPM.
##   v0.2.0: ADDS multivariate Gaussian (NormalMvSpec), DPM. A numeric vector
##           routes to NormalUvSpec; a numeric matrix (one row per observation,
##           >= 2 columns) routes to NormalMvSpec.
## Out-of-scope requests (rjmcmc; non-normal distributions) fail with an
## informative pointer to the roadmap version that will deliver them (project
## knowledge: never silently absorb a mismatch).
## ---------------------------------------------------------------------------

# Resolve the component DistributionSpec from the requested distribution name
# and the observed data shape. Keeps nimixClust() free of inline class logic.
.selectClusterSpec <- function(distribution, isMv, d) {
  distribution <- tolower(distribution)
  if (distribution %in% c("normal-mv", "mvnormal", "multivariate-normal")) {
    if (!isMv)
      stop("distribution = 'normal-mv' requires a numeric matrix with >= 2 ",
           "columns.", call. = FALSE)
    return(getDistribution("normal-mv"))
  }
  if (distribution == "normal-uv") {
    if (isMv)
      stop("distribution = 'normal-uv' requires a univariate numeric vector.",
           call. = FALSE)
    return(getDistribution("normal-uv"))
  }
  if (distribution %in% c("studentt", "student-t", "t"))
    return(getDistribution(if (isMv) "student-t-mv" else "student-t"))
  if (distribution %in% c("normalgamma", "normal-gamma"))
    return(getDistribution(if (isMv) "normal-gamma-mv" else "normal-gamma"))
  if (distribution == "student-t-mv") return(getDistribution("student-t-mv"))
  if (distribution == "normal-gamma-mv")
    return(getDistribution("normal-gamma-mv"))
  if (distribution == "poisson") {
    if (isMv) stop("distribution = 'poisson' is univariate (counts).",
                   call. = FALSE)
    return(getDistribution("poisson"))
  }
  if (distribution == "binomial") {
    if (isMv) stop("distribution = 'binomial' is univariate (counts).",
                   call. = FALSE)
    return(getDistribution("binomial"))
  }
  if (distribution == "normal")
    return(getDistribution(if (isMv) "normal-mv" else "normal-uv"))
  stop("distribution = '", distribution, "' is not available. Supported ",
       "clustering families: 'normal' (univariate/multivariate), 'studentt', ",
       "'normalgamma', 'poisson', 'binomial'.", call. = FALSE)
}

#' Bayesian mixture clustering
#'
#' Fit a Bayesian Gaussian mixture model for clustering. Two engines are
#' available: a finite mixture with a fixed, known number of components
#' (\code{method = "fixedk"}, using the argument \code{K}) and a Dirichlet
#' Process Mixture that infers the number of components (\code{method = "dpm"},
#' using the truncation level \code{K_max}). The component family (univariate or
#' multivariate Gaussian) is chosen from the shape of \code{data} and the
#' \code{distribution} argument.
#'
#' @param data A numeric vector (univariate) or a numeric matrix with one row
#'   per observation and one column per dimension (multivariate). A single-column
#'   matrix is treated as univariate.
#' @param K Integer number of components for the finite mixture
#'   (\code{method = "fixedk"}). Required for that method; must not be given for
#'   \code{method = "dpm"} (use \code{K_max} there).
#' @param K_max Integer truncation level for the Dirichlet Process Mixture
#'   (\code{method = "dpm"}); the number of components is estimated up to this
#'   bound. Because the dCRP sampler errors if the occupied-cluster count ever
#'   needs to exceed it, \code{K_max} should sit comfortably above the expected
#'   number of clusters; a generous data-aware default (giving headroom above
#'   that count) is used when missing. Must not be given for
#'   \code{method = "fixedk"} (use \code{K}).
#' @param distribution Component distribution. \code{"normal"} (default) picks
#'   the univariate or multivariate Gaussian automatically from the data shape;
#'   \code{"normal-uv"} / \code{"normal-mv"} force a specific one. Student-t /
#'   Poisson / Binomial are planned for v0.4.0.
#' @param method Engine: \code{"dpm"} (default; estimate the number of
#'   components), \code{"fixedk"} (finite mixture with known \code{K}), or
#'   \code{"rjmcmc"} (planned for v0.5.0, currently errors).
#' @param prior A named list of prior overrides passed to
#'   \code{\link{defaultPrior}} (univariate: \code{cLoc}, \code{nu0};
#'   multivariate: \code{cLoc}, \code{df0}) plus, for the DPM, optional
#'   \code{concPrior = c(shape, rate)} for the concentration, or, for the finite
#'   mixture, \code{dirichletConc} for the Dirichlet weight prior.
#' @param mcmcControl A named list with \code{niter}, \code{nburnin},
#'   \code{thin}.
#' @param initMethod Initialisation for the cluster allocation: \code{"kmeans"}
#'   (default, dispersed start) or \code{"single"}.
#' @param seed Integer RNG seed for reproducibility.
#' @param verbose Logical; print NIMBLE's configuration and progress output.
#'   Defaults to \code{FALSE} (quiet): NIMBLE's compilation notes and the benign
#'   dCRP truncation note are silenced, while nimix's own diagnostics (e.g. a
#'   censored-posterior warning) and any error still surface. Set \code{TRUE} to
#'   see NIMBLE's configuration and a progress bar.
#'
#' @return A \code{\linkS4class{FitResult}}. Call \code{summary()} for
#'   relabelled estimates, \code{plot()} for diagnostics, and \code{predict()}
#'   for the posterior predictive density.
#'
#' @references
#' de Valpine, P., et al. (2017). Programming with models ... with NIMBLE.
#' \emph{JCGS}, 26(2), 403--413. \doi{10.1080/10618600.2016.1172487}
#'
#' Neal, R.M. (2000). Markov chain sampling methods for Dirichlet process
#' mixture models. \emph{JCGS}, 9(2), 249--265.
#' \doi{10.1080/10618600.2000.10474879}
#'
#' McLachlan, G.J., & Peel, D. (2000). \emph{Finite Mixture Models}. Wiley.
#' \doi{10.1002/0471721182}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#'
#' ## Univariate, number of clusters estimated (DPM)
#' y <- c(rnorm(100, -3, 1), rnorm(100, 3, 1))
#' fit <- nimixClust(y, K_max = 8,
#'                   mcmcControl = list(niter = 2000, nburnin = 1000),
#'                   verbose = FALSE)
#' summary(fit)
#' plot(fit, type = "K")
#'
#' ## Univariate, fixed number of components (finite mixture)
#' fit2 <- nimixClust(y, K = 2, method = "fixedk",
#'                    mcmcControl = list(niter = 2000, nburnin = 1000),
#'                    verbose = FALSE)
#' summary(fit2)
#'
#' ## Multivariate (2-D), DPM
#' Y <- rbind(matrix(rnorm(200, -2), ncol = 2),
#'            matrix(rnorm(200,  2), ncol = 2))
#' fitMv <- nimixClust(Y, K_max = 8,
#'                     mcmcControl = list(niter = 2000, nburnin = 1000),
#'                     verbose = FALSE)
#' summary(fitMv)
#' plot(fitMv, type = "cluster")
#' }
#' @export
nimixClust <- function(data,
                       K = NULL,
                       K_max = NULL,
                       distribution = "normal",
                       method = c("dpm", "fixedk", "rjmcmc"),
                       prior = list(),
                       mcmcControl = list(),
                       initMethod = c("kmeans", "single"),
                       seed = 1L,
                       verbose = FALSE) {
  cl <- match.call()
  method <- match.arg(method)
  initMethod <- match.arg(initMethod)

  # --- roadmap / scope guards ----------------------------------------------
  if (method == "rjmcmc")
    stop("method = 'rjmcmc' is planned for v0.5.0 and is not yet implemented. ",
         "Use method = 'dpm' (estimate K) or method = 'fixedk' (known K).",
         call. = FALSE)

  # K is for the finite mixture (fixed, known number of components); K_max is
  # the truncation level for the DPM (K is estimated). Guard against swapping
  # them so the error is informative rather than a downstream surprise.
  if (method == "fixedk") {
    if (is.null(K))
      stop("method = 'fixedk' needs the number of components K (a known/assumed ",
           "value), e.g. nimixClust(data, K = 3, method = 'fixedk').",
           call. = FALSE)
    if (!is.null(K_max))
      stop("Use K (not K_max) with method = 'fixedk'. K_max is for the DPM, ",
           "which estimates the number of components.", call. = FALSE)
  } else {  # dpm
    if (!is.null(K))
      stop("Use K_max (not K) with method = 'dpm'. The DPM estimates the ",
           "number of components; K (a fixed value) is for method = 'fixedk'.",
           call. = FALSE)
  }

  # --- resolve data shape --------------------------------------------------
  if (is.matrix(data)) {
    if (!is.numeric(data)) stop("data matrix must be numeric.", call. = FALSE)
    if (ncol(data) == 1L) data <- as.numeric(data)   # 1 column == univariate
  } else {
    data <- as.numeric(data)
  }
  if (anyNA(data))
    stop("Missing values are not supported; please remove NAs.", call. = FALSE)

  isMv <- is.matrix(data)
  n <- if (isMv) nrow(data) else length(data)
  d <- if (isMv) ncol(data) else 1L
  if (n < 2L) stop("Need at least 2 observations.", call. = FALSE)
  if (isMv && n <= d)
    stop("Need more observations than dimensions (n = ", n, ", d = ", d, ").",
         call. = FALSE)

  spec <- .selectClusterSpec(distribution, isMv, d)

  # Count families need a non-negative integer response; reject continuous data
  # early rather than failing later inside the sampler.
  if (spec@name %in% c("poisson", "binomial")) {
    if (any(data < 0) || any(abs(data - round(data)) > 1e-8))
      stop("distribution = '", distribution, "' needs non-negative integer ",
           "counts.", call. = FALSE)
  }

  # Resolve the component count: fixed K, or a data-aware default truncation.
  if (method == "fixedk") {
    nComp <- as.integer(K)
    if (nComp < 1L) stop("K must be >= 1.", call. = FALSE)
  } else {
    if (is.null(K_max)) K_max <- .defaultTruncation(n)
    nComp <- as.integer(K_max)
    if (nComp < 2L) stop("K_max must be >= 2.", call. = FALSE)
  }
  if (nComp > n) stop("The number of components cannot exceed the number of ",
                      "observations.", call. = FALSE)

  # --- prior + engine ------------------------------------------------------
  priorList <- defaultPrior(spec, data, control = prior)
  validateParams(spec, priorList)

  if (method == "fixedk") {
    dConc <- if (!is.null(prior$dirichletConc)) prior$dirichletConc else 1
    engine <- FixedKEngine(dirichletConc = dConc)
  } else {
    concPrior <- if (!is.null(prior$concPrior)) prior$concPrior else c(2, 4)
    engine <- DPMEngine(concPrior = concPrior)
  }

  model <- ClusterModel(data = data, distSpec = spec, engine = engine,
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
      data              = data,
      Kmax              = nComp,
      prior             = priorList,
      relabeled         = list(),
      mcmcControl       = raw$mcmcControl,
      call              = cl)
}
