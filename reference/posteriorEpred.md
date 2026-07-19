# Expected value of the posterior predictive

\\E\[Y \mid X\]\\ for each posterior draw and row of data, on the
**response scale**, averaging over components: \\\sum_k w_k \\
g^{-1}(x'\beta_k)\\, where \\g^{-1}\\ is the inverse link of the family.
For a Normal fit that is just \\\sum_k w_k x'\beta_k\\; for a Poisson it
is \\\sum_k w_k \exp(x'\beta_k)\\, and for a Binomial \\\sum_k w_k \\
\mathrm{size} \cdot \mathrm{plogis}(x'\beta_k)\\. The link is applied to
each component before the average, as in `brms::posterior_epred` –
averaging first and transforming after would give a different, wrong
answer.

## Usage

``` r
posteriorEpred(object, newdata = NULL, draws = 500L)
```

## Arguments

- object:

  A
  [`FitResult`](https://madsyair.github.io/nimix/reference/FitResult-class.md)
  from
  [`nimixReg`](https://madsyair.github.io/nimix/reference/nimixReg.md).

- newdata:

  Optional data frame; defaults to the fitted data.

- draws:

  Maximum posterior draws to use (thinned evenly).

## Value

A `draws` x `n` matrix.

## Read this before using it on a mixture

For a mixture the expectation can describe a distribution that no
component has. Two crossing lines with equal weight have a flat mixture
mean – measured on such a fit, the expectation came back at 0.009,
-0.017 and -0.043 for x = -1, 0, 1, for data whose components have
slopes +1.5 and -1.5. Nothing is wrong with the number; it is simply not
a summary anyone wants.
[`posteriorLinpred`](https://madsyair.github.io/nimix/reference/posteriorLinpred.md)
is usually the right question, and for a regime-switching fit
[`nimixForecast`](https://madsyair.github.io/nimix/reference/nimixForecast.md)'s
`$regime` is another.

The weights are the posterior allocation probabilities of the fitted
rows. With `newdata` they become the mixture weights instead, since a
new row's component is unknown – and for `method = "hmm"` that is
refused outright, because a regime weight is a function of time and a
future row has no decoded regime. Project it with
[`nimixForecast`](https://madsyair.github.io/nimix/reference/nimixForecast.md).

## See also

[`posteriorLinpred`](https://madsyair.github.io/nimix/reference/posteriorLinpred.md),
[`posteriorPredictive`](https://madsyair.github.io/nimix/reference/posteriorPredictive.md)
