# Dirichlet Process Mixture engine (native NIMBLE dCRP)

Wraps NIMBLE's Chinese Restaurant Process distribution `dCRP` and the
specialised `CRP`, `CRP_cluster_wrapper`, and `CRP_concentration`
samplers. Empty components are handled natively: the collapsed sampler
only updates parameters of occupied clusters.

## Slots

- `name`:

  Engine identifier, fixed to `"dpm"`.

- `concPrior`:

  A length-2 numeric `c(shape, rate)` for the Gamma hyperprior on the
  concentration parameter \\\alpha\\. A Gamma prior is required for
  NIMBLE to assign the `CRP_concentration` sampler, so the data can
  inform the level of concentration rather than fixing it.

## References

Neal, R.M. (2000). Markov chain sampling methods for Dirichlet process
mixture models. *JCGS*, 9(2), 249–265.
[doi:10.1080/10618600.2000.10474879](https://doi.org/10.1080/10618600.2000.10474879)

Ferguson, T.S. (1973). A Bayesian analysis of some nonparametric
problems. *The Annals of Statistics*, 1(2), 209–230.
[doi:10.1214/aos/1176342360](https://doi.org/10.1214/aos/1176342360)

Escobar, M.D., & West, M. (1995). Bayesian density estimation and
inference using mixtures. *JASA*, 90(430), 577–588.
[doi:10.1080/01621459.1995.10476550](https://doi.org/10.1080/01621459.1995.10476550)
