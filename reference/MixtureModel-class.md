# Virtual base class for mixture models

Virtual base class for mixture models

## Slots

- `data`:

  Numeric data (vector for univariate clustering in v0.1.0).

- `distSpec`:

  A
  [`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md).

- `engine`:

  An
  [`EngineConfig`](https://madsyair.github.io/nimix/reference/EngineConfig-class.md).

- `Kmax`:

  Integer truncation level for the number of components.

- `prior`:

  A named list of prior hyperparameters.
