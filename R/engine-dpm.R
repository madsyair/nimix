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

# Run an expression with NIMBLE's progress/diagnostic chatter silenced. Used
# when verbose = FALSE so quiet runs (e.g. examples) do not print compilation
# notes or the CRP truncation reminder.
.quietly <- function(expr) {
  oldVerbose <- nimble::getNimbleOption("verbose")
  oldBar     <- nimble::getNimbleOption("MCMCprogressBar")
  nimble::nimbleOptions(verbose = FALSE, MCMCprogressBar = FALSE)
  on.exit(nimble::nimbleOptions(verbose = oldVerbose,
                                MCMCprogressBar = oldBar), add = TRUE)
  res <- NULL
  utils::capture.output(
    suppressWarnings(suppressMessages(res <- force(expr))),
    file = nullfile()
  )
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

  samples <- if (verbose) go() else .quietly(go())
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
