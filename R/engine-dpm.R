#' @include class-EngineConfig.R
#' @include class-DistributionSpec.R
NULL

## ---------------------------------------------------------------------------
## engine-dpm.R
##
## Shared MCMC orchestration plus the DPM engine. This layer is
## dimension-agnostic: the DistributionSpec builds the constants, data, inits,
## and parses the parameter traces, so adding a distribution needs no engine
## edits. The DPM path uses NIMBLE's default configuration so its optimised CRP
## samplers and conjugacy recognition are used, including the collapsed sampler
## that does not update empty-component parameters.
## ---------------------------------------------------------------------------

# Muffle ONLY NIMBLE's known-benign configuration chatter, by pattern, and let
# every other condition through. Used in BOTH verbose modes. The dCRP truncation
# reminder and the not-fully-initialized note are emitted as R messages and are
# purely cosmetic: the actual breach of the truncation is a separate hard error
# (translated by the caller below), so dropping these notes -- even under
# verbose = TRUE -- hides nothing actionable. Genuine warnings (anything that
# could signal invalid draws) are never matched here and always propagate;
# errors are never caught here.
.benignNimbleChatter <- paste(
  "less than the number of potential clusters",
  "not strictly valid if it ever proposes",
  "model is not fully initialized",
  "No samplers assigned",   # deliberate: fixed-beta MRF leaves beta unsampled
  sep = "|")

.muffleBenignChatter <- function(expr) {
  ## Robust to zero-length or multi-line condition messages (NIMBLE emits
  ## both): any() collapses, isTRUE() guards NA/empty.
  withCallingHandlers(
    expr,
    message = function(m) {
      if (isTRUE(any(grepl(.benignNimbleChatter, conditionMessage(m)))))
        invokeRestart("muffleMessage")
    },
    warning = function(w) {
      if (isTRUE(any(grepl(.benignNimbleChatter, conditionMessage(w)))))
        invokeRestart("muffleWarning")
      # otherwise: do nothing, so the warning propagates to the user
    })
}

# verbose = FALSE path: in addition to muffling the benign chatter above, also
# drop NIMBLE's cat-based progress notes (capture.output) and turn off the
# progress bar. This is still NOT a blanket warning suppressor -- non-benign
# warnings continue to reach the user even in quiet mode.
.quietly <- function(expr) {
  oldVerbose <- nimble::getNimbleOption("verbose")
  oldBar     <- nimble::getNimbleOption("MCMCprogressBar")
  nimble::nimbleOptions(verbose = FALSE, MCMCprogressBar = FALSE)
  on.exit(nimble::nimbleOptions(verbose = oldVerbose,
                                MCMCprogressBar = oldBar), add = TRUE)
  res <- NULL
  utils::capture.output(
    res <- .muffleBenignChatter(force(expr)),
    file = nullfile())
  res
}

# Build, compile, run, and extract a NIMBLE mixture model. Engine-specific code
# --- compiled-model cache (compile-once-reuse, v0.5.0) ----------------------
# Compiling a NIMBLE model and its MCMC is the dominant cost of a fit. When a
# later fit requests an identical model STRUCTURE -- same generated code, same
# constants/hyperparameters, same monitors, same component spec class -- the
# compiled model and MCMC are reused: only the data and initial values are reset
# and the chain re-run. This reproduces a fresh compile exactly while skipping
# recompilation. Data and inits are values, not structure, so they are never
# part of the cache key.
.nimixModelCache <- new.env(parent = emptyenv())
.nimixModelCache$entries <- list()
.nimixModelCache$builds  <- 0L           # number of genuine (re)compilations
.NIMIX_CACHE_MAX <- 6L

#' Clear the compiled-model cache
#'
#' nimix reuses compiled NIMBLE models across fits that share an identical model
#' structure (see the \code{reuse} entry of \code{mcmcControl} in
#' \code{\link{nimixClust}}). Compiled models are large; this empties the cache
#' and releases them.
#'
#' @return Invisibly, the number of cached compiled models that were removed.
#' @examples
#' nimixClearCache()
#' @export
nimixClearCache <- function() {
  k <- length(.nimixModelCache$entries)
  .nimixModelCache$entries <- list()
  invisible(k)
}

# Structural key: data and inits are deliberately excluded.
.cacheKey <- function(mc, constants, spec, extra = NULL) {
  list(code      = deparse(mc$code),
       constants = constants,
       monitors  = sort(as.character(mc$monitors)),
       specClass = as.character(class(spec)),
       extra     = extra)
}

.cacheGet <- function(key) {
  ents <- .nimixModelCache$entries
  for (i in seq_along(ents)) {
    if (identical(ents[[i]]$key, key)) {
      hit <- ents[[i]]                                  # LRU: move to front
      .nimixModelCache$entries <- c(list(hit), ents[-i])
      return(hit$compiled)
    }
  }
  NULL
}

.cachePut <- function(key, compiled) {
  ents <- c(list(list(key = key, compiled = compiled)),
            .nimixModelCache$entries)
  if (length(ents) > .NIMIX_CACHE_MAX) ents <- ents[seq_len(.NIMIX_CACHE_MAX)]
  .nimixModelCache$entries <- ents
  invisible(NULL)
}

# --- shared NIMBLE runner ---------------------------------------------------
# Every engine (DPM, fixed-K) funnels through this single runner: a spec only
# has to assemble the model code, constants, inits, monitors, and the name
# of the allocation node (xi for the CRP, z for the finite mixture); everything
# downstream is shared.
.runNimbleMixture <- function(spec, mc, constants, dataList, initsFn,
                              n, count, paramDim, prior,
                              mcmcControl, seed, verbose,
                              configureHook = NULL, cacheExtra = NULL) {
  .nimixEnsureMSNBurr()
  ctrl <- utils::modifyList(
    list(niter = 11000L, nburnin = 1000L, thin = 1L, reuse = TRUE,
         nchains = 1L), mcmcControl)
  reuse   <- isTRUE(ctrl$reuse)
  nchains <- max(1L, as.integer(ctrl$nchains))
  # Parallel chains (opt-in). Each worker builds and compiles its OWN model, so
  # no compiled C++ objects are shared across processes -- the only NIMBLE-safe
  # way to parallelise. Uses forking (parallel::mclapply), so it has no effect
  # on Windows (which cannot fork) and falls back to sequential there.
  parallel <- isTRUE(ctrl$parallel) && nchains > 1L
  if (parallel && (.Platform$OS.type == "windows" ||
                   !requireNamespace("parallel", quietly = TRUE))) {
    warning("parallel = TRUE needs a forking platform (not Windows) and the ",
            "'parallel' package; running chains sequentially.", call. = FALSE)
    parallel <- FALSE
  }
  ncores <- if (!is.null(ctrl$ncores)) as.integer(ctrl$ncores) else
    if (parallel) max(1L, parallel::detectCores()) else 1L
  buildInits <- initsFn(seed)          # inits used only to build the structure

  # Expensive, cacheable: build + compile the model and its MCMC.
  buildCompiled <- function(dirName = NULL) {
    .nimixModelCache$builds <- .nimixModelCache$builds + 1L
    rmodel <- nimble::nimbleModel(code = mc$code, constants = constants,
                                  data = dataList, inits = buildInits,
                                  calculate = FALSE)
    # Stochastic variable names (index-stripped) let runWith reset any node not
    # supplied by inits to NA, so the MCMC's initializeModel re-simulates it
    # deterministically under setSeed -- exactly as a fresh compile does --
    # rather than inheriting a previous run's end state.
    stochNodes <- rmodel$getNodeNames(stochOnly = TRUE, includeData = FALSE)
    stochVars  <- unique(sub("\\[.*$", "", stochNodes))
    cmodel <- nimble::compileNimble(rmodel, showCompilerOutput = FALSE,
                                    dirName = dirName)
    conf <- nimble::configureMCMC(rmodel, monitors = mc$monitors,
                                  print = verbose)
    customizeSamplers(spec, conf, rmodel)
    if (!is.null(configureHook)) configureHook(conf, rmodel)
    cmcmc <- nimble::compileNimble(nimble::buildMCMC(conf), project = rmodel,
                                   showCompilerOutput = FALSE)
    list(cmodel = cmodel, cmcmc = cmcmc, stochVars = stochVars)
  }

  # Cheap: reset data + inits on the (possibly reused) compiled model, then run.
  runWith <- function(compiled, chainInits, chainSeed) {
    cm <- compiled$cmodel
    cm$setData(dataList)
    for (nm in names(chainInits))
      tryCatch(cm[[nm]] <- chainInits[[nm]], error = function(e) NULL)
    for (nm in setdiff(compiled$stochVars, names(chainInits)))
      tryCatch(cm[[nm]][] <- NA, error = function(e) NULL)
    nimble::runMCMC(compiled$cmcmc, niter = ctrl$niter, nburnin = ctrl$nburnin,
                    thin = ctrl$thin, setSeed = chainSeed, progressBar = verbose)
  }

  getCompiled <- function() {
    if (!reuse) return(buildCompiled())
    key <- .cacheKey(mc, constants, spec, cacheExtra)
    compiled <- .cacheGet(key)
    if (is.null(compiled)) {
      compiled <- buildCompiled()
      .cachePut(key, compiled)
    }
    compiled
  }

  # Run a single chain, muffling benign chatter and translating the dCRP
  # truncation error into actionable guidance naming K_max.
  runChain <- function(chainInits, chainSeed) {
    go <- function() runWith(getCompiled(), chainInits, chainSeed)
    runOnce <- function() if (verbose) .muffleBenignChatter(go()) else .quietly(go())
    as.matrix(tryCatch(
      runOnce(),
      error = function(e) {
        msg <- conditionMessage(e)
        if (grepl("proper model|more components than|cluster parameters", msg,
                  ignore.case = TRUE))
          stop("The DPM sampler tried to use more components than the ",
               "truncation level K_max = ", count, ". Increase K_max (e.g. ",
               "K_max = ", 2L * count, ") and re-run; a larger truncation only ",
               "adds headroom and does not change the posterior on the number ",
               "of clusters.", call. = FALSE)
        stop(e)
      }))
  }

  parseAlloc <- function(samples) {
    if (is.null(.nodeColInfo(samples, mc$allocNode)))
      matrix(1L, nrow = nrow(samples), ncol = n)
    else
      matrix(as.integer(round(.nodeToArray(samples, mc$allocNode, n))),
             nrow = nrow(samples))
  }

  # --- run chains (chains 2..M reuse the compiled model via the cache) -------
  chainSamples <- vector("list", nchains)
  chainAlloc   <- vector("list", nchains)
  chainK       <- vector("list", nchains)
  chainAlpha   <- vector("list", nchains)
  chainBeta    <- vector("list", nchains)
  chainEntropy <- vector("list", nchains)
  hasAlpha <- FALSE
  hasBeta  <- FALSE
  # Produce the raw samples matrix for one chain. In parallel mode each worker
  # builds & compiles its own model (safe under forking); in sequential mode
  # chains reuse the cached compiled model (compile-once-reuse).
  chainSamplesFor <- function(ch) {
    chainSeed <- seed + (ch - 1L)          # distinct seed + dispersed inits
    if (parallel) {
      # unique compile dir per worker so forked NIMBLE builds never collide
      compiled <- buildCompiled(dirName = tempfile("nimix_chain_"))
      as.matrix(runWith(compiled, initsFn(chainSeed), chainSeed))
    } else {
      runChain(initsFn(chainSeed), chainSeed)
    }
  }

  if (parallel) {
    chainSamples <- parallel::mclapply(
      seq_len(nchains), chainSamplesFor,
      mc.cores = min(ncores, nchains), mc.set.seed = FALSE)
    bad <- which(vapply(chainSamples,
                        function(x) inherits(x, "try-error") || is.null(x),
                        logical(1)))
    if (length(bad))
      stop("Parallel chain ", bad[1], " failed: ",
           conditionMessage(attr(chainSamples[[bad[1]]], "condition")),
           call. = FALSE)
  } else {
    chainSamples <- lapply(seq_len(nchains), chainSamplesFor)
  }

  # Derived per-chain quantities (cheap; always sequential).
  for (ch in seq_len(nchains)) {
    S <- chainSamples[[ch]]
    a <- parseAlloc(S)
    chainAlloc[[ch]]   <- a
    chainK[[ch]]       <- .rowDistinct(a, count)
    chainEntropy[[ch]] <- .allocEntropy(a, count)
    if ("alpha" %in% colnames(S)) {
      hasAlpha <- TRUE
      chainAlpha[[ch]] <- as.numeric(S[, "alpha"])
    }
    if ("beta" %in% colnames(S) && stats::var(as.numeric(S[, "beta"])) > 0) {
      hasBeta <- TRUE
      chainBeta[[ch]] <- as.numeric(S[, "beta"])
    }
  }
  samples <- do.call(rbind, chainSamples)      # pooled draws (all chains)
  alloc   <- do.call(rbind, chainAlloc)
  Kpost   <- unlist(chainK, use.names = FALSE)

  # For the DPM (allocation node "xi") a posterior sitting at the truncation
  # level is censored. Warn only on sustained censoring across the pooled draws.
  # The FixedK path (node "z") has K fixed, so this never applies.
  if (identical(mc$allocNode, "xi")) {
    atCeiling <- mean(Kpost >= count)
    if (atCeiling > 0.02)
      warning("The DPM used all K_max = ", count, " truncation slots in ",
              round(100 * atCeiling), "% of iterations: the posterior on the ",
              "number of clusters is likely censored. Re-run with a larger ",
              "K_max (e.g. ", 2L * count, ").", call. = FALSE)
  }

  diagnostics <- .multiChainDiag(chainK, if (hasAlpha) chainAlpha else NULL,
                                 if (hasBeta) chainBeta else NULL,
                                 chainEntropy)
  # Preserve chain identity for the pooled draws. Stacking chains without a
  # marker makes per-chain diagnostics (post-hoc R-hat, per-chain traces, and
  # bayesplot's iter x chain x param draws array) impossible to reconstruct.
  # Stored in the diagnostics list, so the FitResult class is unchanged.
  diagnostics$chainId <- rep(seq_len(nchains),
                             vapply(chainSamples, nrow, integer(1)))
  paramTrace  <- extractParamTraces(spec, samples, count, d = paramDim,
                                    prior = prior)

  list(mcmcSamples = samples, Kposterior = as.integer(Kpost),
       clusterAllocation = alloc, paramTrace = paramTrace,
       diagnostics = diagnostics, mcmcControl = ctrl)
}

#' @describeIn runEngine Dirichlet Process Mixture run (NIMBLE dCRP).
setMethod("runEngine", "DPMEngine",
  function(engine, model, mcmcControl = list(), initMethod = "kmeans",
           seed = 1L, verbose = TRUE, ...) {
    spec  <- model@distSpec
    data  <- model@data
    prior <- model@prior
    L     <- model@Kmax
    n     <- .nObs(data)
    d     <- .dataDimOf(data)

    mc <- buildModelCode(spec, engine, n = n, L = L, d = d)
    constants <- c(buildConstants(spec, prior, n),
                   list(L = L,
                        aAlpha = engine@concPrior[1],
                        bAlpha = engine@concPrior[2]))
    dataList <- buildDataList(spec, data)

    initRatio <- .resolveInitRatio(mcmcControl)
    initsFn <- function(s) {
      ci <- .withSeed(s, function() componentInits(spec, prior, data, L,
                      initMethod = initMethod, initRatio = initRatio))
      c(list(xi = ci$alloc, alpha = 1), ci$params)
    }

    paramDim <- if (!is.null(prior$p)) prior$p else d
    .runNimbleMixture(spec, mc, constants, dataList, initsFn,
                      n = n, count = L, paramDim = paramDim, prior = prior,
                      mcmcControl = mcmcControl, seed = seed, verbose = verbose)
  }
)
