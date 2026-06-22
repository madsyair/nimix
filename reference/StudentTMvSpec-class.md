# Multivariate Student-t component specification (direct density)

Direct multivariate-t kernel via a user-defined NIMBLE distribution.
Same marginal as
[`NormalGammaMvSpec`](https://madsyair.github.io/nimix/reference/NormalGammaMvSpec-class.md)
but non-conjugate. `df` is a fixed hyperparameter.

## Slots

- `name`:

  Fixed to `"student-t-mv"`.

- `paramNames`:

  `c("mu", "Sigma")`.

## References

Lange, K.L., Little, R.J.A., & Taylor, J.M.G. (1989). Robust statistical
modeling using the t distribution. *JASA*, 84(408), 881–896.
[doi:10.1080/01621459.1989.10478852](https://doi.org/10.1080/01621459.1989.10478852)

## See also

[`NormalGammaMvSpec`](https://madsyair.github.io/nimix/reference/NormalGammaMvSpec.md)
for the conjugate path.
