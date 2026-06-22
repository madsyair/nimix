# Multivariate Gaussian component specification

Multivariate Gaussian component specification

## Slots

- `name`:

  Fixed to `"normal-mv"`.

- `paramNames`:

  `c("mu", "Sigma")`.

- `dataDim`:

  `NA_integer_`: the actual dimension \\d\\ is taken from the data at
  fit time and carried in the prior list.

## References

Zhang, Z., Chan, K.L., Wu, Y., & Chen, C. (2004). Learning a
multivariate Gaussian mixture model with the reversible jump MCMC
algorithm. *Statistics and Computing*, 14, 343–355.
[doi:10.1023/B:STCO.0000039481.32735.0c](https://doi.org/10.1023/B%3ASTCO.0000039481.32735.0c)

Dellaportas, P., & Papageorgiou, I. (2006). Multivariate mixtures of
normals with unknown number of components. *Statistics and Computing*,
16, 57–68.
[doi:10.1007/s11222-006-5338-6](https://doi.org/10.1007/s11222-006-5338-6)

## See also

[`nimixClust`](https://madsyair.github.io/nimix/reference/nimixClust.md),
[`NormalUvSpec`](https://madsyair.github.io/nimix/reference/NormalUvSpec-class.md)
