# Draws from the posterior predictive distribution

Like
[`posteriorEpred`](https://madsyair.github.io/nimix/reference/posteriorEpred.md)
but with the residual noise added: what a new observation would actually
look like, not merely its expectation. For a mixture this is the
well-behaved one of the three – it is genuinely bimodal where the
components disagree, rather than collapsing to a mean that sits between
them.

## Usage

``` r
posteriorPredictive(object, newdata = NULL, draws = 500L)
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

## See also

[`posteriorLinpred`](https://madsyair.github.io/nimix/reference/posteriorLinpred.md),
[`posteriorEpred`](https://madsyair.github.io/nimix/reference/posteriorEpred.md),
[`nimixForecast`](https://madsyair.github.io/nimix/reference/nimixForecast.md)
for projecting a regime forward in time.
