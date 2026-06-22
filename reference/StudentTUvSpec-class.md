# Univariate Student-t component specification

Univariate Student-t component specification

## Slots

- `name`:

  Fixed to `"student-t"`.

- `paramNames`:

  `c("mu", "tau")` (location and precision; the reported scale is
  \\\sigma = \tau^{-1/2}\\).

- `dataDim`:

  `1L`.

## References

Lange, K.L., Little, R.J.A., & Taylor, J.M.G. (1989). Robust statistical
modeling using the t distribution. *JASA*, 84(408), 881–896.
[doi:10.1080/01621459.1989.10478852](https://doi.org/10.1080/01621459.1989.10478852)

## See also

[`NormalGammaUvSpec`](https://madsyair.github.io/nimix/reference/NormalGammaUvSpec.md)
for the conjugate scale-mixture path to the same marginal distribution.
