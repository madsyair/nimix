# Build the NIMBLE model code for a (distribution, engine) pair

Dispatches on BOTH the component distribution and the engine. This is
the extensibility seam described in adding a new distribution to the DPM
engine means adding one method here, not editing the engine or the spec
base class.

## Usage

``` r
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'PoissonRegSpec,DPMEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'PoissonRegSpec,FixedKEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'BinomialRegSpec,DPMEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'BinomialRegSpec,FixedKEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'GMSNBurrUvSpec,FixedKEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'GMSNBurrUvSpec,DPMEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'NormalRegSpec,DPMEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'StudentTRegSpec,DPMEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'StudentTRegSpec,FixedKEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'NormalGammaRegSpec,DPMEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'NormalGammaRegSpec,FixedKEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'MSNBurrUvSpec,FixedKEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'MSNBurr2aUvSpec,FixedKEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'MSNBurrUvSpec,DPMEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'MSNBurr2aUvSpec,DPMEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'NormalMvRegSpec,DPMEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'NormalMvRegSpec,FixedKEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'StudentTMvRegSpec,DPMEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'StudentTMvRegSpec,FixedKEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'NormalGammaMvRegSpec,DPMEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'NormalGammaMvRegSpec,FixedKEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'NormalMvSpec,DPMEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'NormalGammaMvSpec,DPMEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'NormalUvSpec,DPMEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'NormalGammaUvSpec,DPMEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'PoissonSpec,DPMEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'BinomialSpec,DPMEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'StudentTMvSpec,DPMEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'StudentTUvSpec,DPMEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'NormalUvSpec,FixedKEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'NormalMvSpec,FixedKEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'NormalRegSpec,FixedKEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'StudentTUvSpec,FixedKEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'NormalGammaUvSpec,FixedKEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'NormalGammaMvSpec,FixedKEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'StudentTMvSpec,FixedKEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'PoissonSpec,FixedKEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'BinomialSpec,FixedKEngine'
buildModelCode(spec, engine, n, L, ...)

# S4 method for class 'DistributionSpec,MRFEngine'
buildModelCode(spec, engine, n, L, ...)
```

## Arguments

- spec:

  A
  [`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md).

- engine:

  An
  [`EngineConfig`](https://madsyair.github.io/nimix/reference/EngineConfig-class.md).

- n:

  Integer, number of observations.

- L:

  Integer, cluster-parameter truncation length (`= K_max`).

- ...:

  Reserved for methods.

## Value

A list with elements `code` (a `nimbleCode` object), `monitors`
(character vector), and `paramNodes` (named character vector mapping
logical parameter names to model nodes).

## Functions

- `buildModelCode(spec = PoissonRegSpec, engine = DPMEngine)`: Poisson
  GLM regression DPM code (log link).

- `buildModelCode(spec = PoissonRegSpec, engine = FixedKEngine)`:
  Poisson GLM regression fixed-K code (log link).

- `buildModelCode(spec = BinomialRegSpec, engine = DPMEngine)`: Binomial
  GLM regression DPM code (logit link).

- `buildModelCode(spec = BinomialRegSpec, engine = FixedKEngine)`:
  Binomial GLM regression fixed-K code (logit link).

- `buildModelCode(spec = GMSNBurrUvSpec, engine = FixedKEngine)`:
  GMSNBurr finite mixture (fixed K).

- `buildModelCode(spec = GMSNBurrUvSpec, engine = DPMEngine)`: GMSNBurr
  DPM mixture.

- `buildModelCode(spec = NormalRegSpec, engine = DPMEngine)`: DPM
  mixture-of-linear-regressions model code (dCRP) with a conjugate
  Normal-Inverse-Gamma cluster prior and constant (CRP) gating.

  Builds NIMBLE code for \$\$y_i \sim N(x_i^\top \beta\_{\xi_i},
  s^2\_{\xi_i}), \quad \xi\_{1:n} \sim CRP(\alpha, n),\$\$ with \\s^2_j
  \sim InvGamma(nu_0, s_0)\\ and \\\beta_j \mid s^2_j \sim N_p(b_0,
  s^2_j B_0)\\. The number of predictors \\p\\ and the design matrix
  \\X\\ are constants.

- `buildModelCode(spec = StudentTRegSpec, engine = DPMEngine)`:
  Student-t regression DPM code (direct t density; the scale enters as
  precision `tau = 1 / s2`).

- `buildModelCode(spec = StudentTRegSpec, engine = FixedKEngine)`:
  Student-t regression fixed-K code.

- `buildModelCode(spec = NormalGammaRegSpec, engine = DPMEngine)`:
  Normal-Gamma regression DPM code (scale mixture; conjugate
  coefficients conditional on `omega`).

- `buildModelCode(spec = NormalGammaRegSpec, engine = FixedKEngine)`:
  Normal-Gamma regression fixed-K code.

- `buildModelCode(spec = MSNBurrUvSpec, engine = FixedKEngine)`: MSNBurr
  finite mixture (fixed K).

- `buildModelCode(spec = MSNBurr2aUvSpec, engine = FixedKEngine)`:
  MSNBurr-IIa finite mixture (fixed K).

- `buildModelCode(spec = MSNBurrUvSpec, engine = DPMEngine)`: MSNBurr
  DPM mixture.

- `buildModelCode(spec = MSNBurr2aUvSpec, engine = DPMEngine)`:
  MSNBurr-IIa DPM mixture.

- `buildModelCode(spec = NormalMvRegSpec, engine = DPMEngine)`:
  Multivariate-response Normal regression DPM code.

- `buildModelCode(spec = NormalMvRegSpec, engine = FixedKEngine)`:
  Multivariate-response Normal regression fixed-K.

- `buildModelCode(spec = StudentTMvRegSpec, engine = DPMEngine)`:
  Multivariate-response Student-t regression DPM code (direct
  multivariate-t kernel).

- `buildModelCode(spec = StudentTMvRegSpec, engine = FixedKEngine)`:
  Multivariate-response Student-t regression fixed-K.

- `buildModelCode(spec = NormalGammaMvRegSpec, engine = DPMEngine)`:
  Multivariate-response Normal-Gamma regression DPM code (scale mixture;
  conjugate cluster updates).

- `buildModelCode(spec = NormalGammaMvRegSpec, engine = FixedKEngine)`:
  Multivariate-response Normal-Gamma regression fixed-K code.

- `buildModelCode(spec = NormalMvSpec, engine = DPMEngine)`:
  Multivariate Gaussian DPM model code (dCRP) with a conjugate
  Normal-Inverse-Wishart cluster base measure.

  Builds NIMBLE code for \$\$y_i \sim N_d(\mu\_{\xi_i},
  \Sigma\_{\xi_i}), \quad \xi\_{1:n} \sim CRP(\alpha, n),\$\$ with
  \\\Sigma_j \sim InvWishart(S_0, df_0)\\ and \\\mu_j \sim N_d(\mu_0,
  \Sigma_j / \kappa_0)\\. The dimension \\d\\ is a constant so the index
  ranges `1:d` expand at model-build time.

- `buildModelCode(spec = NormalGammaMvSpec, engine = DPMEngine)`:
  Multivariate Normal-Gamma scale-mixture DPM code with a
  per-observation latent precision multiplier `omega`.

- `buildModelCode(spec = NormalUvSpec, engine = DPMEngine)`: Univariate
  Gaussian DPM model code (dCRP).

  Builds the NIMBLE code for \$\$y_i \sim N(\mu\_{\xi_i}, s^2\_{\xi_i}),
  \quad \xi\_{1:n} \sim CRP(\alpha, n),\$\$ with a conjugate
  Normal-Inverse-Gamma cluster prior and a Gamma hyperprior on
  \\\alpha\\. Cluster-parameter vectors have length `L = K_max`
  (NIMBLE's exact truncation; the sampler stays proper as long as the
  number of occupied clusters is strictly below `L`).

- `buildModelCode(spec = NormalGammaUvSpec, engine = DPMEngine)`:
  Normal-Gamma scale-mixture DPM model code (dCRP) with a
  per-observation latent precision multiplier `omega`.

- `buildModelCode(spec = PoissonSpec, engine = DPMEngine)`: Poisson DPM
  model code (Gamma-Poisson conjugate).

- `buildModelCode(spec = BinomialSpec, engine = DPMEngine)`: Binomial
  DPM model code (Beta-Binomial conjugate).

- `buildModelCode(spec = StudentTMvSpec, engine = DPMEngine)`:
  Multivariate Student-t DPM code using the user-defined `dmvt_nimix`
  kernel with a Normal-Inverse-Wishart cluster prior.

- `buildModelCode(spec = StudentTUvSpec, engine = DPMEngine)`:
  Univariate Student-t DPM model code (dCRP). The t density is evaluated
  directly; df is a constant.

- `buildModelCode(spec = NormalUvSpec, engine = FixedKEngine)`:
  Univariate Gaussian finite-mixture code (fixed K).

- `buildModelCode(spec = NormalMvSpec, engine = FixedKEngine)`:
  Multivariate Gaussian finite-mixture code (fixed K).

- `buildModelCode(spec = NormalRegSpec, engine = FixedKEngine)`:
  Mixture-of-linear-regressions finite-mixture code (fixed K).

- `buildModelCode(spec = StudentTUvSpec, engine = FixedKEngine)`:
  Univariate Student-t finite-mixture code (fixed K).

- `buildModelCode(spec = NormalGammaUvSpec, engine = FixedKEngine)`:
  Univariate Normal-Gamma (scale-mixture Student-t) finite-mixture code
  (fixed K).

- `buildModelCode(spec = NormalGammaMvSpec, engine = FixedKEngine)`:
  Multivariate Normal-Gamma scale-mixture finite-mixture code (fixed K).

- `buildModelCode(spec = StudentTMvSpec, engine = FixedKEngine)`:
  Multivariate Student-t finite-mixture code (fixed K) using the
  user-defined `dmvt_nimix` kernel.

- `buildModelCode(spec = PoissonSpec, engine = FixedKEngine)`: Poisson
  finite-mixture code (fixed K).

- `buildModelCode(spec = BinomialSpec, engine = FixedKEngine)`: Binomial
  finite-mixture code (fixed K).

- `buildModelCode(spec = DistributionSpec, engine = MRFEngine)`: Default
  MRF kernel for any family: the family's fixed-K kernel with the
  Dirichlet-categorical label layer replaced by the joint Potts node
  (derived mechanically; see the spatial design notes).

## References

Hurn et al. (2003)
[doi:10.1198/1061860031329](https://doi.org/10.1198/1061860031329) ; de
Valpine et al. (2017)
[doi:10.1080/10618600.2016.1172487](https://doi.org/10.1080/10618600.2016.1172487)
.

Zhang et al. (2004)
[doi:10.1023/B:STCO.0000039481.32735.0c](https://doi.org/10.1023/B%3ASTCO.0000039481.32735.0c)
; de Valpine et al. (2017)
[doi:10.1080/10618600.2016.1172487](https://doi.org/10.1080/10618600.2016.1172487)
.

Neal, R.M. (2000)
[doi:10.1080/10618600.2000.10474879](https://doi.org/10.1080/10618600.2000.10474879)
; de Valpine et al. (2017)
[doi:10.1080/10618600.2016.1172487](https://doi.org/10.1080/10618600.2016.1172487)
.
