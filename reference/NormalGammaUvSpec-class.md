# Univariate Normal-Gamma (scale-mixture Student-t) component specification

A conjugate scale-mixture representation of a univariate Student-t
component: identical marginal to
[`StudentTUvSpec`](https://madsyair.github.io/nimix/reference/StudentTUvSpec-class.md),
but with conjugate cluster updates because the kernel is Gaussian
conditional on a latent per-observation precision multiplier. The
degrees of freedom `df` are a fixed hyperparameter.

## Slots

- `name`:

  Fixed to `"normal-gamma"`.

- `paramNames`:

  `c("mu", "s2")`.

## References

Andrews, D.F., & Mallows, C.L. (1974). Scale mixtures of normal
distributions. *JRSS-B*, 36(1), 99–102.
[doi:10.1111/j.2517-6161.1974.tb00989.x](https://doi.org/10.1111/j.2517-6161.1974.tb00989.x)

West, M. (1987). On scale mixtures of normal distributions.
*Biometrika*, 74(3), 646–648.
[doi:10.1093/biomet/74.3.646](https://doi.org/10.1093/biomet/74.3.646)

Lange, K.L., Little, R.J.A., & Taylor, J.M.G. (1989). Robust statistical
modeling using the t distribution. *JASA*, 84(408), 881–896.
[doi:10.1080/01621459.1989.10478852](https://doi.org/10.1080/01621459.1989.10478852)

## See also

[`StudentTUvSpec`](https://madsyair.github.io/nimix/reference/StudentTUvSpec.md)
for the direct (non-conjugate) path to the same marginal.
