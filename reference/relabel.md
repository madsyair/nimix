# Correct label switching in a fitted mixture

Post-hoc relabelling of MCMC output so that per-component posterior
summaries are meaningful. Conditions on the modal number of occupied
clusters, then applies an algorithm from the label.switching package.

## Usage

``` r
relabel(fit, method = "ECR-ITERATIVE-1", ...)

# S4 method for class 'FitResult'
relabel(fit, method = "ECR-ITERATIVE-1", ...)
```

## Arguments

- fit:

  A
  [`FitResult`](https://madsyair.github.io/nimix/reference/FitResult-class.md).

- method:

  One of `"ECR-ITERATIVE-1"` (default) or `"ECR"`. The default is the
  pivot-free iterative ECR; `"ECR"` uses the highest-posterior
  allocation as the pivot.

- ...:

  Reserved.

## Value

The `fit` with its `relabeled` slot populated.

## Functions

- `relabel(FitResult)`: Relabelling for a fitted result.

## References

Papastamoulis, P., & Iliopoulos, G. (2010). An artificial allocations
based solution to the label switching problem. *JCGS*, 19(2), 313–331.
[doi:10.1198/jcgs.2010.09008](https://doi.org/10.1198/jcgs.2010.09008)

Papastamoulis, P. (2016). label.switching: An R package for dealing with
the label switching problem in MCMC outputs. *JSS, Code Snippets*,
69(1).
[doi:10.18637/jss.v069.c01](https://doi.org/10.18637/jss.v069.c01)
