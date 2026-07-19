# Per-component linear predictors

For each posterior draw and each row of data, the linear predictor
\\\eta_k = x'\beta_k\\ of *every* component – no mixing, no noise. This
is usually the object you want from a mixture of regressions: it is what
each regime or cluster actually predicts, which
[`posteriorEpred`](https://madsyair.github.io/nimix/reference/posteriorEpred.md)
averages away.

## Usage

``` r
posteriorLinpred(object, newdata = NULL, transform = FALSE, draws = 500L)
```

## Arguments

- object:

  A
  [`FitResult`](https://madsyair.github.io/nimix/reference/FitResult-class.md)
  from
  [`nimixReg`](https://madsyair.github.io/nimix/reference/nimixReg.md).

- newdata:

  Optional data frame; defaults to the fitted data.

- transform:

  If `TRUE`, apply the inverse-link so each slice is a component
  response mean rather than a linear predictor. Default `FALSE` (the
  linear-predictor scale, as the name says).

- draws:

  Maximum posterior draws to use (thinned evenly).

## Value

A `draws` x `n` x `K` array, with the component index last.

## Details

Like `brms::posterior_linpred`, this returns the *linear predictor* –
the scale on which the coefficients are linear – not the response mean.
For a Poisson (log link) or Binomial (logit link) fit that is \\\log
\mu\\ or the log-odds, not \\\mu\\. Set `transform = TRUE` to apply the
inverse link and get each component's response mean instead (the
counterpart of `brms`'s `transform` argument); for a Normal fit the two
coincide.

## See also

[`posteriorEpred`](https://madsyair.github.io/nimix/reference/posteriorEpred.md),
[`posteriorPredictive`](https://madsyair.github.io/nimix/reference/posteriorPredictive.md)

## Examples

``` r
if (FALSE) { # \dontrun{
lp <- posteriorLinpred(fit, newdata = data.frame(x = c(-1, 0, 1)))
apply(lp, c(2, 3), mean)     # each component's fitted line
} # }
```
