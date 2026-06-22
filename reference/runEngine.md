# Run a mixture engine on a model (internal generic)

Dispatches on the engine:
[`DPMEngine`](https://madsyair.github.io/nimix/reference/DPMEngine-class.md)
builds a CRP model,
[`FixedKEngine`](https://madsyair.github.io/nimix/reference/FixedKEngine-class.md)
builds a finite-mixture model. Returns the raw pieces used to construct
a
[`FitResult`](https://madsyair.github.io/nimix/reference/FitResult-class.md).

## Usage

``` r
runEngine(
  engine,
  model,
  mcmcControl = list(),
  initMethod = "kmeans",
  seed = 1L,
  verbose = TRUE,
  ...
)

# S4 method for class 'DPMEngine'
runEngine(
  engine,
  model,
  mcmcControl = list(),
  initMethod = "kmeans",
  seed = 1L,
  verbose = TRUE,
  ...
)

# S4 method for class 'FixedKEngine'
runEngine(
  engine,
  model,
  mcmcControl = list(),
  initMethod = "kmeans",
  seed = 1L,
  verbose = TRUE,
  ...
)
```

## Arguments

- engine:

  An
  [`EngineConfig`](https://madsyair.github.io/nimix/reference/EngineConfig-class.md).

- model:

  A
  [`MixtureModel`](https://madsyair.github.io/nimix/reference/MixtureModel-class.md).

- mcmcControl, initMethod, seed, verbose:

  Sampler controls.

- ...:

  Reserved.

## Value

A named list with the MCMC samples, the posterior of the number of
occupied components, the allocation matrix, the parsed parameter traces,
and the resolved MCMC control list.

## Functions

- `runEngine(DPMEngine)`: Dirichlet Process Mixture run (NIMBLE dCRP).

- `runEngine(FixedKEngine)`: Finite-mixture run with fixed K (Dirichlet
  weights + categorical allocation).
