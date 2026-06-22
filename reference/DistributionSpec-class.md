# Virtual base class for mixture component distributions

`DistributionSpec` is the abstract S4 class that every component
distribution in nimix extends. It is never instantiated directly; use a
concrete subclass such as
[`NormalUvSpec`](https://madsyair.github.io/nimix/reference/NormalUvSpec-class.md).

## Slots

- `name`:

  Character scalar, a short identifier (e.g. `"normal-uv"`).

- `paramNames`:

  Character vector of component parameter names.

- `priorSpec`:

  A named list of prior hyperparameters. Empty until
  [`defaultPrior`](https://madsyair.github.io/nimix/reference/defaultPrior.md)
  (or the user) fills it.

- `dataDim`:

  Integer, the data dimension the spec is meant for (1 for univariate).
  Used for early validation in
  [`nimixClust`](https://madsyair.github.io/nimix/reference/nimixClust.md).

## References

Frühwirth-Schnatter, S. (2006). *Finite Mixture and Markov Switching
Models*. Springer.
[doi:10.1007/978-0-387-35768-3](https://doi.org/10.1007/978-0-387-35768-3)

## See also

[`NormalUvSpec`](https://madsyair.github.io/nimix/reference/NormalUvSpec-class.md),
[`registerDistribution`](https://madsyair.github.io/nimix/reference/registerDistribution.md)
