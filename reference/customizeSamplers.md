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

# S4 method for class 'NormalRegSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'StudentTRegSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'NormalGammaRegSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'NormalGammaMvRegSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'NormalGammaMvSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'NormalGammaUvSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'SkewNormalMvOSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'SkewIStudentMvOSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'SkewNormalMvOGenSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'SkewIStudentMvOGenSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'MSNBurrUvSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'MSNBurr2aUvSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'GMSNBurrUvSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'FSSNUvSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'FOSSEPUvSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'FSSTUvSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'JFSTUvSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'SEPUvSpec'
customizeSamplers(spec, conf, model, ...)

# S4 method for class 'LEPUvSpec'
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

- `customizeSamplers(NormalRegSpec)`: Replace RW samplers on `betaTilde`
  and `s2Tilde` with the exact Normal-Inverse-Gamma conditional (FixedK
  path only; the DPM path already receives the conjugate CRP wrapper).

- `customizeSamplers(StudentTRegSpec)`: Student-t regression keeps
  NIMBLE's default samplers: the inherited exact NIG Gibbs step is only
  valid for Gaussian errors.

- `customizeSamplers(NormalGammaRegSpec)`: Slice-sample the latent
  precision multipliers.

- `customizeSamplers(NormalGammaMvRegSpec)`: Slice-sample the latent
  precision multipliers.

- `customizeSamplers(NormalGammaMvSpec)`: Slice-sample the latent
  precision multipliers.

- `customizeSamplers(NormalGammaUvSpec)`: Slice-sample the latent
  precision multipliers, which mixes the partition better than the
  default random walk.

- `customizeSamplers(SkewNormalMvOSpec)`: Slice-sample the Householder
  angles: the FS likelihood in theta is bounded and can be multimodal
  near the edge of Theta^2, where an adaptive random walk mixes poorly.

- `customizeSamplers(SkewIStudentMvOSpec)`: Slice-sample the Householder
  angles.

- `customizeSamplers(SkewNormalMvOGenSpec)`: Slice-sample the
  Householder angles.

- `customizeSamplers(SkewIStudentMvOGenSpec)`: Slice-sample the
  Householder angles.

- `customizeSamplers(MSNBurrUvSpec)`: AF_slice block over (mu, sigma,
  alpha).

- `customizeSamplers(MSNBurr2aUvSpec)`: AF_slice block over (mu, sigma,
  alpha).

- `customizeSamplers(GMSNBurrUvSpec)`: AF_slice block over (mu, sigma,
  alpha, theta).

- `customizeSamplers(FSSNUvSpec)`: AF_slice block over (mu, sigma,
  alpha).

- `customizeSamplers(FOSSEPUvSpec)`: AF_slice block over (mu, sigma,
  alpha, theta).

- `customizeSamplers(FSSTUvSpec)`: AF_slice block over (mu, sigma,
  alpha, nu); the truncated nu node poses no problem for slice sampling.

- `customizeSamplers(JFSTUvSpec)`: AF_slice block over (mu, sigma,
  alpha, theta).

- `customizeSamplers(SEPUvSpec)`: AF_slice block over (mu, sigma, nu).

- `customizeSamplers(LEPUvSpec)`: AF_slice block over (mu, sigma, nu).
