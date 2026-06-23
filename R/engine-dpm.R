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

# Run an expression with NIMBLE's *cosmetic* progress/diagnostic chatter
# silenced, used when verbose = FALSE. Importantly this is NOT a blanket warning
# suppressor: only NIMBLE's known-benign configuration chatter is muffled. Any
# other warning -- in particular anything that could signal the sampler produced
# invalid draws -- is allowed to propagate to the user even when verbose = FALSE,
# so quiet mode never hides a genuine validity problem. Errors are never caught
# here.
#
# .benignNimbleChatter holds patterns for NIMBLE messages/warnings that are
# purely informational (the dCRP truncation reminder is a heads-up that the
# truncated representation is in use; the actual breach of the truncation is a
# hard error, raised separately and translated by the caller, not a warning).
.benignNimbleChatter <- paste(
  "less than the number of potential clusters",
  "not strictly valid if it ever proposes",
  "model is not fully initialized",
  sep = "|")

.quietly <- function(expr) {
  oldVerbose <- nimble::getNimbleOption("verbose")
  oldBar     <- nimble::getNimbleOption("MCMCprogressBar")
  nimble::nimbleOptions(verbose = FALSE, MCMCprogressBar = FALSE)
  on.exit(nimble::nimbleOptions(verbose = oldVerbose,
                                MCMCprogressBar = oldBar), add = TRUE)
  res <- NULL
  # capture.output drops NIMBLE's cat-based console notes; the handlers below
  # muffle ONLY the benign config messages/warnings by pattern and let every
  # other condition through (a non-benign warning is re-signalled so the user
  # still sees it).
  withCallingHandlers(
    utils::capture.output(
      res <- force(expr),
      file = nullfile()),
    message = function(m) {
      if (grepl(.benignNimbleChatter, conditionMessage(m)))
        invokeRestart("muffleMessage")
    },
    warning = function(w) {
      if (grepl(.benignNimbleChatter, conditionMessage(w)))
        invokeRestart("muffleWarning")
      # otherwise: do nothing, so the warning propagates to the user
    })
  res
}

# Build, compile, run, and extract a NIMBLE mixture model. Engine-specific code
# only has to assemble the model code, constants, inits, monitors, and the name
# of the allocation node (xi for the CRP, z for the finite mixture); everything
# downstream is shared.
.runNimbleMixture <- function(spec, mc, constants, dataList, inits,
                              n, count, paramDim, prior,
                              mcmcControl, seed, verbose) {
  ctrl <- utils::modifyList(
    list(niter = 11000L, nburnin = 1000L, thin = 1L), mcmcControl)

  go <- function() {
    rmodel <- nimble::nimbleModel(code = mc$code, constants = constants,
                                  data = dataList, inits = inits,
                                  calculate = FALSE)
    cmodel <- nimble::compileNimble(rmodel, showCompilerOutput = FALSE)
    conf <- nimble::configureMCMC(rmodel, monitors = mc$monitors,
                                  print = verbose)
    customizeSamplers(spec, conf, rmodel)
    mcmc <- nimble::buildMCMC(conf)
    cmcmc <- nimble::compileNimble(mcmc, project = rmodel,
                                   showCompilerOutput = FALSE)
    nimble::runMCMC(cmcmc, niter = ctrl$niter, nburnin = ctrl$nburnin,
                    thin = ctrl$thin, setSeed = seed, progressBar = verbose)
  }

  runOnce <- function() if (verbose) go() else .quietly(go())
  samples <- tryCatch(
    runOnce(),
    error = function(e) {
      msg <- conditionMessage(e)
      # The dCRP sampler stops with this message when the chain needs more
      # occupied clusters than the truncation level provides. Translate NIMBLE's
      # raw text into actionable guidance naming K_max (= the truncation level).
      if (grepl("proper model|more components than|cluster parameters", msg,
                ignore.case = TRUE))
        stop("The DPM sampler tried to use more components than the truncation ",
             "level K_max = ", count, ". Increase K_max (e.g. K_max = ",
             2L * count, ") and re-run; a larger truncation only adds headroom ",
             "and does not change the posterior on the number of clusters.",
             call. = FALSE)
      stop(e)
    })
  samples <- as.matrix(samples)

  # When the allocation node is monitored we parse it; for a single-component
  # fit (K = 1) the allocation is fixed to 1 and is not in the samples.
  if (is.null(.nodeColInfo(samples, mc$allocNode))) {
    alloc <- matrix(1L, nrow = nrow(samples), ncol = n)
  } else {
    allocArr <- .nodeToArray(samples, mc$allocNode, n)
    alloc <- matrix(as.integer(round(allocArr)), nrow = nrow(samples))
  }
  Kpost <- apply(alloc, 1L, function(r) length(unique(r)))

  # For the DPM (allocation node "xi") a posterior that sits at the truncation
  # level is censored: the data may support more clusters than K_max allows. Warn
  # only on *sustained* censoring (a meaningful share of kept draws using every
  # slot), not a one-off excursion to the ceiling. The FixedK path (node "z") has
  # K fixed, so this never applies.
  if (identical(mc$allocNode, "xi")) {
    atCeiling <- mean(Kpost >= count)
    if (atCeiling > 0.02)
      warning("The DPM used all K_max = ", count, " truncation slots in ",
              round(100 * atCeiling), "% of iterations: the posterior on the ",
              "number of clusters is likely censored. Re-run with a larger ",
              "K_max (e.g. ", 2L * count, ").", call. = FALSE)
  }
  paramTrace <- extractParamTraces(spec, samples, count, d = paramDim,
                                   prior = prior)

  list(mcmcSamples = samples, Kposterior = as.integer(Kpost),
       clusterAllocation = alloc, paramTrace = paramTrace, mcmcControl = ctrl)
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

    ci <- componentInits(spec, prior, data, L, initMethod = initMethod)
    inits <- c(list(xi = ci$alloc, alpha = 1), ci$params)

    paramDim <- if (!is.null(prior$p)) prior$p else d
    .runNimbleMixture(spec, mc, constants, dataList, inits,
                      n = n, count = L, paramDim = paramDim, prior = prior,
                      mcmcControl = mcmcControl, seed = seed, verbose = verbose)
  }
)
