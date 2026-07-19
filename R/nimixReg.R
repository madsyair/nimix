#' @include registerDistribution.R
NULL

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
#' @param random Optional one-sided formula for group-level random effects.
#'   \code{random = ~ region} gives a random intercept: every component's
#'   linear predictor gains a shared group offset
#'   \eqn{b_{g(i)} \sim N(0, \tau^2)}. \code{random = ~ x | region}
#'   additionally gives a random slope for \code{x} (which must be a term of
#'   \code{formula}); as in \code{lme4}'s \code{(x|g)}, the intercept comes
#'   along. Both offsets use a sum-to-zero constraint -- the parameterisation
#'   that mixes well -- so the component coefficients absorb the group means
#'   and the reported \code{b} / \code{sRE} are centred. The two offsets are
#'   given independent priors by default. \code{tauRE} and \code{tauSlope}
#'   estimate the spread of \emph{your} groups, which with few groups is
#'   itself variable. Currently supported with \code{method = "fixedk"} and
#'   \code{distribution = "normal"} or \code{"studentt"} (heavy-tailed
#'   residuals).
#' @param K Integer number of components for the finite mixture
#'   (\code{method = "fixedk"}). Required for that method.
#' @param K_max Integer truncation level for the DPM (\code{method = "dpm"});
#'   it must sit comfortably above the expected number of clusters, since the
#'   dCRP sampler errors if the occupied-cluster count ever needs to exceed it.
#'   A generous data-aware default (with headroom) is used when missing.
#' @param distribution Component distribution. Gaussian (\code{"normal"}),
#'   heavy-tailed (\code{"studentt"}, \code{"normalgamma"}), GLM
#'   (\code{"poisson"}, \code{"binomial"}), and the nine neo-normal skew
#'   families (\code{"msnburr"}, \code{"msnburr2a"}, \code{"gmsnburr"},
#'   \code{"fssn"}, \code{"fsst"}, \code{"sep"}, \code{"lep"}, \code{"fossep"},
#'   \code{"jfst"}) are available; multivariate variants exist for the Gaussian
#'   and heavy-tailed families.
#'
#'   \strong{Non-Gaussian families cost effective sample size.} Only the
#'   Gaussian regression has a conjugate coefficient update; every other family
#'   (neo-normal, heavy-tailed, GLM) is sampled by NIMBLE's defaults, which mix
#'   more slowly. Measured on \code{"fixedk"} with two components, the Gaussian
#'   reaches about 10.7 ESS/s against 1.9 for \code{"msnburr"} and 1.4 for the
#'   three-shape \code{"gmsnburr"} -- roughly a 5--7x factor. Budget
#'   proportionally more iterations for a skewed or heavy-tailed fit; the
#'   three-shape families (gmsnburr, fsst, fossep, jfst) are the slowest.
#' @param method Fitting method: \code{"dpm"}, \code{"fixedk"},
#'   \code{"mrf"}, or \code{"hmm"}.
#'
#'   \code{"hmm"} fits a \strong{Markov-switching regression} (Hamilton
#'   1989): the coefficients and error variance switch with a latent
#'   first-order Markov regime, so the rows of \code{data} are a
#'   \emph{time series} rather than an exchangeable sample -- their order
#'   carries meaning here as it does nowhere else in \code{nimixReg}. Give
#'   \code{K} for the number of regimes, as for \code{"fixedk"}; the regime
#'   path is marginalised out of the likelihood and decoded afterwards, so
#'   \code{\link{viterbiPath}} gives the most probable regime sequence.
#'   Currently \code{distribution = "normal"} (Gaussian) or \code{"poisson"} (log-link counts), or \code{"studentt"} / \code{"normalgamma"} (heavy-tailed), or \code{"binomial"} (logit-link
#'   proportions, with the number of trials in \code{prior = list(size = )}).
#'
#'   \strong{Budget more iterations than for \code{"fixedk"}.} Marginalising
#'   the regime path rules out the conjugate Normal-Inverse-Gamma update, so
#'   the coefficients are sampled by NIMBLE's defaults instead. Measured on a
#'   two-regime series, that costs roughly four times the effective sample
#'   size per second (ESS/s 2.3 against 8.9 for the conjugate \code{"fixedk"}
#'   sampler), even though the marginalised chain has a shorter wall time. A
#'   light run gives usably-decoded regimes but wide coefficient intervals;
#'   raise \code{niter} until the intervals settle.
#' @param spatialWeights Optional \code{\linkS4class{SpatialWeightSpec}} (one
#'   region per observation). Required by, and only used with,
#'   \code{method = "mrf"}.
#' @param gating Mixing-weight model: \code{"constant"} (default; weights do not
#'   depend on covariates). \code{"covariate"} (concomitant gating) is a planned
#'   opt-in and currently errors.
#' @param prior A named list of prior overrides passed to
#'   \code{\link{defaultPrior}} (e.g. \code{g} for the g-prior factor,
#'   \code{nu0} for the InvGamma shape, and \code{s2Guess} for its scale --
#'   see \sQuote{Reading the error variance} below; with a multivariate
#'   response, \code{sigmaGuess} plays the same role for the residual
#'   covariance) plus,
#'   for the DPM, optional \code{concPrior = c(shape, rate)}, or, for the
#'   finite mixture, \code{dirichletConc}.
#'
#' @section Reading the error variance:
#' The prior on each component's error variance \code{s2} is centred on the
#' residual variance of a \emph{global} OLS fit, which ignores the mixture.
#' For separated components that quantity measures the spread \emph{between}
#' components as much as the spread within one, so the prior is deliberately
#' conservative and \code{s2} is biased upward. The bias is largest exactly
#' where mixtures are most useful -- well-separated components at moderate
#' sample size -- and it disappears as the per-component sample size grows:
#' with a prior/truth scale ratio of about 60, the measured bias runs near
#' 9x at 25 observations per component, 2.5x at 150, 1.2x at 1000, and 1.0x
#' by 5000. It is a bias in the safe direction (wider predictive intervals),
#' and the coefficients are unaffected.
#'
#' The conservatism is load-bearing rather than incidental: it also
#' regularises against splitting when \code{K} is over-specified, so it
#' stays the default.
#'
#' If you know the within-component scale -- from a pilot study, the
#' literature, or domain knowledge -- say so directly:
#' \code{prior = list(s2Guess = 0.3)} makes 0.3 the prior mean of \code{s2}.
#' (\code{s0} sets the raw InvGamma scale instead, for callers who think in
#' those terms; give one or the other.) On a simulated benchmark with a true
#' \code{s2} of 0.25 that moved the estimate from 0.78 (3.1x) to 0.36
#' (1.4x), slopes untouched. The override is deliberately absolute: a
#' multiplier on the automatic scale would only ever mean "a fraction of a
#' quantity that measures the wrong thing", which is a knob, not a
#' statement.
#'
#' \strong{Where a tighter prior is safe.} A tenfold tighter scale left the
#' DPM's recovery of \code{K} untouched (modal \code{K} = 2 throughout a
#' benchmark with two true components), and is likewise safe for
#' \code{fixedk} when \code{K} is correct, where it halved the error in
#' \code{s2}. It is \strong{not} safe for \code{fixedk} with an
#' over-specified \code{K}: on the same data with \code{K = 4} against two
#' true components, the default occupied two components while a tenfold
#' tighter prior occupied three. Wide components make two suffice; narrow
#' ones make the model want more. If you are unsure whether \code{K} is
#' over-specified, leave the default alone or use the DPM.
#'
#' \strong{Recipe: two-stage empirical Bayes.} If you need the scale right
#' and have no external knowledge, fit once with the default, read the
#' within-component residual variance off that fit, and refit with it as
#' \code{s2Guess}. It costs a second compile-and-run, which is why it is not
#' the default, but it needs no guesswork:
#'
#' \preformatted{
#' fit1 <- nimixReg(y ~ x, df, K = 2, method = "fixedk")
#' fit1 <- relabel(fit1)
#' # residuals of each point under its own MAP component
#' z    <- binderPartition(fit1)$partition
#' cf   <- fit1@relabeled$summary
#' Xm   <- model.matrix(y ~ x, df)
#' res  <- df$y - rowSums(Xm * as.matrix(cf[z, c("(Intercept)", "x")]))
#' s2hat <- sum(res^2) / (nrow(df) - 2 * ncol(Xm))
#'
#' fit2 <- nimixReg(y ~ x, df, K = 2, method = "fixedk",
#'                  prior = list(s2Guess = s2hat))
#' }
#'
#' The first fit's \emph{allocation} is reliable even where its \code{s2} is
#' not -- that is what makes the recipe work.
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
                     random = NULL,
                     K = NULL,
                     K_max = NULL,
                     distribution = "normal",
                     method = c("dpm", "fixedk", "mrf", "hmm"),
                     gating = c("constant", "covariate"),
                     prior = list(),
                     mcmcControl = list(),
                     initMethod = c("kmeans", "single", "spread"),
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

  # K is for the finite mixture (fixed components) and for the HMM's regimes;
  # K_max is the DPM truncation.
  if (method %in% c("fixedk", "mrf", "hmm")) {
    if (is.null(K))
      stop("method = '", method, "' needs the number of ",
           if (method == "hmm") "regimes K." else "components K.",
           call. = FALSE)
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
      msnburr = "msnburr-reg", msnburr2a = "msnburr2a-reg", sep = "sep-reg", fssn = "fssn-reg", gmsnburr = "gmsnburr-reg", lep = "lep-reg", fsst = "fsst-reg", fossep = "fossep-reg", jfst = "jfst-reg",
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

  if (method %in% c("fixedk", "mrf", "hmm")) {
    # "hmm" belongs here, not in the truncation branch: its K is the fixed
    # number of regimes, exactly like fixedk's components. Omitting it sent K
    # to .defaultTruncation(n) and silently fitted an 8-regime chain.
    nComp <- as.integer(K)
    if (nComp < 1L) stop("K must be >= 1.", call. = FALSE)
    if (method == "hmm" && nComp < 2L)
      stop("method = 'hmm' needs K >= 2 regimes.", call. = FALSE)
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
  if (!is.null(random)) {
    if (!inherits(random, "formula") || length(random) != 2L)
      stop("`random` must be a one-sided formula: ~ group for a random ",
           "intercept, or ~ x | group to add a random slope for x.",
           call. = FALSE)
    # ~ group          -> random intercept
    # ~ x | group      -> random intercept AND random slope for x (the
    #                     intercept always comes along, as in lme4's (x|g))
    rexpr <- random[[2L]]
    if (is.call(rexpr) && identical(as.character(rexpr[[1L]]), "|")) {
      slopeVars <- all.vars(rexpr[[2L]])
      gv <- all.vars(rexpr[[3L]])
      if (length(slopeVars) != 1L)
        stop("`random` supports exactly one random-slope variable, ",
             "e.g. random = ~ x | region.", call. = FALSE)
    } else {
      slopeVars <- character(0)
      gv <- all.vars(random)
    }
    if (length(gv) != 1L)
      stop("`random` currently supports exactly one grouping factor.",
           call. = FALSE)
    if (!gv %in% names(data))
      stop("Grouping variable '", gv, "' not found in `data`.", call. = FALSE)
    neoRE <- c("msnburr", "msnburr2a", "gmsnburr", "fssn", "fsst",
               "sep", "lep", "fossep", "jfst")
    glmRE <- c("poisson", "binomial")
    reDist <- c("normal", "studentt", neoRE, glmRE)
    # fixed-K takes any RE-capable family. DPM adds a random-cluster allocation
    # on top, which sits cleanly with a fixed external grouping (both act on iid
    # exchangeable observations), so it is allowed for the neo-normal families.
    # HMM is excluded: its single ordered time series has no exchangeable
    # grouping, so a random effect is unidentified and confounds with the
    # regime transitions -- see the RE-in-HMM design note. Panel HMM would be
    # the prerequisite, and it does not exist yet.
    okRE <- (method == "fixedk" && distribution %in% reDist) ||
            (method == "dpm" && distribution %in% c("normal", neoRE))
    if (!okRE)
      stop("`random` is supported for method = 'fixedk' with a Gaussian, ",
           "Student-t, neo-normal, or GLM distribution, and for method = ",
           "'dpm' with a Gaussian or neo-normal distribution. HMM is excluded ",
           "(a single time series has no exchangeable grouping for a random ",
           "effect); panel-HMM support would be the prerequisite.",
           call. = FALSE)
    grp <- as.integer(factor(data[[gv]]))
    G <- max(grp)
    if (G < 3L)
      stop("Random effects need at least 3 groups; '", gv, "' has ", G,
           ".", call. = FALSE)
    priorList$hasRE  <- TRUE
    priorList$reGrp  <- grp
    priorList$reG    <- G
    priorList$reVar  <- gv
    # Data-scaled bounds for the random-effect SDs, following the same
    # principle as the rest of nimix's priors. Fixed bounds silently broke
    # the offsets whenever the response was on a large scale: with y x1000
    # the needed tauRE was 771 against a hard ceiling of 5, and
    # cor(b_hat, truth) collapsed from 0.992 to 0.091. tauRE has the units of
    # y, so its scale is sd(y); a factor of 5 is generous (offsets cannot
    # plausibly spread several times wider than the response itself).
    # tauRE has the units of the linear predictor. For a Gaussian/neo-normal
    # (identity link) that is the scale of y itself; for a GLM the random
    # effect lives on the LINK scale (log or logit), where sd(y) is the wrong
    # ruler -- a Poisson count of thousands has huge sd(y) but log-scale
    # offsets of order 1. So scale by sd(y) for identity-link families and by
    # a fixed link-scale bound for GLMs.
    isGLM <- distribution %in% c("poisson", "binomial")
    if (isGLM) {
      sdRE <- 1
      if (is.null(priorList$tauMax)) priorList$tauMax <- 5
      if (is.null(priorList$tauMin)) priorList$tauMin <- 1e-4
    } else {
      sdY <- stats::sd(as.numeric(y))
      if (!is.finite(sdY) || sdY <= 0) sdY <- 1
      sdRE <- sdY
      if (is.null(priorList$tauMax)) priorList$tauMax <- 5 * sdY
      if (is.null(priorList$tauMin)) priorList$tauMin <- 1e-4 * sdY
    }
    if (length(slopeVars) == 1L) {
      sv <- slopeVars
      if (!sv %in% colnames(X))
        stop("Random-slope variable '", sv, "' must be a term of the fixed ",
             "effects formula (columns: ",
             paste(setdiff(colnames(X), "(Intercept)"), collapse = ", "),
             ").", call. = FALSE)
      priorList$hasRESlope <- TRUE
      priorList$reSlopeX   <- as.numeric(X[, sv])
      priorList$reSlopeVar <- sv
      # tauSlope has the units of y/x, so its scale is sd(y)/sd(x).
      sdXs <- stats::sd(priorList$reSlopeX)  # slope tau: y/x units
      if (!is.finite(sdXs) || sdXs <= 0)
        stop("Random-slope variable '", sv, "' is constant; it carries no ",
             "group-varying slope information.", call. = FALSE)
      scS <- sdRE / sdXs
      if (is.null(priorList$tauMaxSlope)) priorList$tauMaxSlope <- 5 * scS
      if (is.null(priorList$tauMinSlope)) priorList$tauMinSlope <- 1e-4 * scS
    }
  }
  priorList$formula   <- formula
  if (isMvResp && is.null(priorList$respNames))
    priorList$respNames <- colnames(y)

  if (method == "fixedk") {
    dConc <- if (!is.null(prior$dirichletConc)) prior$dirichletConc else 1
    engine <- FixedKEngine(dirichletConc = dConc)
  } else if (method == "hmm") {
    # Markov-switching regression (Hamilton 1989): the coefficients and error
    # variance switch with a latent first-order Markov regime, so the
    # observation ORDER carries meaning here as it does nowhere else in
    # nimixReg -- rows are a time series, not an exchangeable sample.
    tConc <- if (!is.null(prior$transConc)) prior$transConc else 1
    engine <- HMMEngine(transConc = tConc)
  } else if (method == "mrf") {
    if (nRegions(spatialWeights) != n)
      stop("spatialWeights has ", nRegions(spatialWeights),
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
