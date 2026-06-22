# Finite-mixture engine with fixed, known K

The simplest engine: the number of components K is fixed (not inferred).
Mixing weights have a symmetric Dirichlet prior and each observation has
a categorical allocation. Because there is no Chinese Restaurant
Process, the truncation considerations of the DPM do not apply, and
NIMBLE assigns conjugate samplers to the weights and component
parameters. Useful as a fast baseline when K is known or assumed, and
for classical model selection by comparing fits across several values of
K.

## Slots

- `name`:

  Engine identifier, fixed to `"fixedk"`.

- `dirichletConc`:

  Positive scalar concentration of the symmetric Dirichlet prior on the
  mixing weights (\\1\\ is uniform on the simplex; values below 1 favour
  sparser weight vectors).

## References

McLachlan, G.J., & Peel, D. (2000). *Finite Mixture Models*. Wiley.
[doi:10.1002/0471721182](https://doi.org/10.1002/0471721182)

Frühwirth-Schnatter, S. (2006). *Finite Mixture and Markov Switching
Models*. Springer.
[doi:10.1007/978-0-387-35768-3](https://doi.org/10.1007/978-0-387-35768-3)
