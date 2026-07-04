# nimix: Bayesian Mixture Clustering and Regression with NIMBLE

Bayesian mixture modelling built on the nimble platform. The package
implements univariate and multivariate mixture clustering and
mixture-of-regressions through two inference engines: a Dirichlet
Process Mixture (DPM) engine based on the Chinese Restaurant Process
(which estimates the number of occupied components) and a fixed-K
finite-mixture engine. It is organised around an extensible S4
[`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md)
contract so that new component distributions and engines can be added
without rewriting existing code.

## Inference engines

- `method = "dpm"`: Dirichlet process / Chinese restaurant process; the
  number of occupied components is estimated from the data.

- `method = "fixedk"`: finite mixture with a known number of components
  `K`.

## Component distributions

Gaussian (univariate and multivariate), Student-t and Normal-Gamma
(heavy-tailed, univariate and multivariate), and Poisson / Binomial
counts, for both clustering
([`nimixClust`](https://madsyair.github.io/nimix/reference/nimixClust.md))
and regression
([`nimixReg`](https://madsyair.github.io/nimix/reference/nimixReg.md),
including multivariate responses).

## References

de Valpine, P., Turek, D., Paciorek, C.J., Anderson-Bergman, C., Temple
Lang, D., & Bodik, R. (2017). Programming with models: writing
statistical algorithms for general model structures with NIMBLE.
*Journal of Computational and Graphical Statistics*, 26(2), 403–413.
[doi:10.1080/10618600.2016.1172487](https://doi.org/10.1080/10618600.2016.1172487)

Neal, R.M. (2000). Markov chain sampling methods for Dirichlet process
mixture models. *Journal of Computational and Graphical Statistics*,
9(2), 249–265.
[doi:10.1080/10618600.2000.10474879](https://doi.org/10.1080/10618600.2000.10474879)

## See also

Useful links:

- <https://github.com/madsyair/nimix>

- <https://madsyair.github.io/nimix/>

- Report bugs at <https://github.com/madsyair/nimix/issues>

## Author

**Maintainer**: Achmad Syahrul Choir <madsyair@stis.ac.id>
([ORCID](https://orcid.org/0000-0001-7088-0646))

Authors:

- Achmad Syahrul Choir <madsyair@stis.ac.id>
  ([ORCID](https://orcid.org/0000-0001-7088-0646))
