# Posterior similarity matrix

Returns the \\n \times n\\ matrix \\S\_{ij} = \Pr(z_i = z_j \mid y)\\,
estimated as the fraction of posterior draws in which observations \\i\\
and \\j\\ share a cluster. The quantity is invariant to label
permutations and to the number of occupied clusters, so *every* draw
contributes – unlike
[`relabel`](https://madsyair.github.io/nimix/reference/relabel.md),
which must condition on the modal number of clusters before component
parameters can be aligned.

## Usage

``` r
psm(fit)
```

## Arguments

- fit:

  A `FitResult`.

## Value

A symmetric matrix with unit diagonal.

## References

Binder, D. A. (1978), Biometrika 65, 31–38.

## See also

[`binderPartition`](https://madsyair.github.io/nimix/reference/binderPartition.md)
for a point-estimate partition,
[`relabel`](https://madsyair.github.io/nimix/reference/relabel.md) for
component-parameter summaries.
