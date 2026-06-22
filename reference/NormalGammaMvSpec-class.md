# Multivariate Normal-Gamma (scale-mixture multivariate-t) component

Conjugate scale-mixture representation of a multivariate Student-t
component: identical marginal to
[`StudentTMvSpec`](https://madsyair.github.io/nimix/reference/StudentTMvSpec-class.md),
with conjugate Normal-Inverse-Wishart cluster updates. `df` is a fixed
hyperparameter.

## Slots

- `name`:

  Fixed to `"normal-gamma-mv"`.

- `paramNames`:

  `c("mu", "Sigma")`.

## References

Andrews, D.F., & Mallows, C.L. (1974). Scale mixtures of normal
distributions. *JRSS-B*, 36(1), 99–102.
[doi:10.1111/j.2517-6161.1974.tb00989.x](https://doi.org/10.1111/j.2517-6161.1974.tb00989.x)

Backlund, E., & Hobert, J.P. (2020). \[Gibbs sampling for multivariate
linear regression with errors that are scale mixtures of normals under a
conjugate Normal-Inverse-Wishart prior.\]

## See also

[`StudentTMvSpec`](https://madsyair.github.io/nimix/reference/StudentTMvSpec.md)
for the direct (non-conjugate) path.
