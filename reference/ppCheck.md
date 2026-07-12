# Posterior predictive check

Simulates replicated data sets from the fitted mixture – conditionally
on each posterior draw's allocation, which makes the check
label-invariant and valid for every engine including the spatial MRF –
and compares test statistics of the replicates against the observed data
via posterior predictive p-values (Gelman, Meng & Stern 1996; Gelman et
al. 2020, Section 6). Extreme p-values (near 0 or 1) flag aspects of the
data the model fails to reproduce; values near 0.5 indicate no evidence
of misfit for that statistic.

## Usage

``` r
ppCheck(
  fit,
  nrep = 200,
  statistics = c("mean", "sd", "min", "max"),
  seed = 1L,
  store_yrep = FALSE
)
```

## Arguments

- fit:

  A
  [`FitResult`](https://madsyair.github.io/nimix/reference/FitResult-class.md)
  from
  [`nimixClust`](https://madsyair.github.io/nimix/reference/nimixClust.md).

- nrep:

  Number of replicated data sets (posterior draws used), thinned evenly
  from the retained draws. Default 200.

- statistics:

  Named list of test-statistic functions, or a character vector naming
  built-ins from `mean`, `sd`, `min`, `max`, `skew`. Multivariate data
  applies each statistic column-wise and reports per-column results.

- seed:

  RNG seed for the replicate simulation.

- store_yrep:

  If `TRUE`, attach the simulated replicates as attribute `"yrep"` (with
  `"y"` and `"draws"`) on the result – the inputs graphical PPC
  functions consume.

## Value

An object of class `nimixPPC`: a data frame of observed value, replicate
mean, and posterior predictive p-value per statistic, printed with
guidance.

## References

Gelman, A., Meng, X.-L., & Stern, H. (1996). Posterior predictive
assessment of model fitness via realized discrepancies. *Statistica
Sinica*, 6, 733–807.

Gelman, A., et al. (2020). Bayesian workflow. *arXiv:2011.01808*,
Section 6.
