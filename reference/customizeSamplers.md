# Customise MCMC samplers for a component spec

Hook called by the engine after
[`configureMCMC()`](https://rdrr.io/pkg/nimble/man/configureMCMC.html)
and before
[`buildMCMC()`](https://rdrr.io/pkg/nimble/man/buildMCMC.html), letting
a spec swap NIMBLE's default sampler on its own nodes. The default is a
no-op; the scale-mixture specs override it to put a slice sampler on the
latent precision multipliers, which mixes the partition markedly better
than the default random walk.

## Usage

``` r
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'DistributionSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'NormalGammaRegSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'NormalGammaMvRegSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'NormalGammaUvSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'NormalGammaMvSpec'
customizeSamplers(spec, conf, model, ...)
```

## Arguments

- spec:

  A
  [`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md).

- conf:

  An MCMC configuration object.

- model:

  The (uncompiled) NIMBLE model.

- ...:

  Unused.

## Value

The (possibly modified) `conf`, invisibly.

## Functions

- `customizeSamplers(DistributionSpec)`: Default: leave NIMBLE's
  samplers unchanged.

- `customizeSamplers(NormalGammaRegSpec)`: Slice-sample the latent
  precision multipliers.

- `customizeSamplers(NormalGammaMvRegSpec)`: Slice-sample the latent
  precision multipliers.

- `customizeSamplers(NormalGammaUvSpec)`: Slice-sample the latent
  precision multipliers, which mixes the partition better than the
  default random walk.

- `customizeSamplers(NormalGammaMvSpec)`: Slice-sample the latent
  precision multipliers.
