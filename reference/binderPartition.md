# Binder-loss point partition (Dahl's least-squares criterion)

Selects, among the partitions actually visited by the chain, the one
minimising the posterior expected Binder loss – equivalently, the draw
whose pairwise co-clustering matrix is closest (in squared error) to the
posterior similarity matrix (Dahl 2006). All draws inform the similarity
matrix; none are discarded.

## Usage

``` r
binderPartition(fit, S = NULL)
```

## Arguments

- fit:

  A `FitResult`.

- S:

  Optional precomputed
  [`psm`](https://madsyair.github.io/nimix/reference/psm.md) matrix.

## Value

A list with `partition` (integer vector, labels recoded to `1..K`), `K`,
`draw` (index of the selected iteration), `score` (the least-squares
criterion value), and `psm`.

## References

Binder, D. A. (1978), Biometrika 65, 31–38. Dahl, D. B. (2006), in
*Bayesian Inference for Gene Expression and Proteomics*, Cambridge
University Press, 201–218.

## See also

[`psm`](https://madsyair.github.io/nimix/reference/psm.md),
[`relabel`](https://madsyair.github.io/nimix/reference/relabel.md).
