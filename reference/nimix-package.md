# nimix: Bayesian Mixture Clustering and Regression with NIMBLE

Bayesian mixture modelling built on the nimble platform. The release
implements univariate and multivariate Gaussian mixture clustering
through a Dirichlet Process Mixture (DPM) engine based on the Chinese
Restaurant Process. The package is organised around an extensible S4
[`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md)
contract so that new component distributions and (later) a
reversible-jump engine can be added without rewriting existing code.

## Roadmap (this is v0.2.0)

- v0.1.0: S4 foundation,
  [`NormalUvSpec`](https://madsyair.github.io/nimix/reference/NormalUvSpec-class.md),
  univariate
  [`nimixClust`](https://madsyair.github.io/nimix/reference/nimixClust.md)
  on the DPM engine.

- v0.2.0 (this release): multivariate clustering
  ([`NormalMvSpec`](https://madsyair.github.io/nimix/reference/NormalMvSpec-class.md)).

- v0.3.0: mixture-of-regressions
  ([`nimixReg`](https://madsyair.github.io/nimix/reference/nimixReg.md)).

- v0.5.0+: reversible jump MCMC engine.

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

**Maintainer**: nimix Developers <nimix@example.org>

Authors:

- nimix Developers <nimix@example.org>
