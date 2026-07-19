# Draw a response given linear predictors and error scale

The family-specific tail of
[`posteriorPredictive`](https://madsyair.github.io/nimix/reference/posteriorPredictive.md):
given the per-observation linear predictor(s) and (for a Normal) the
error variance, draw one response each. Defaults to the Gaussian \\y
\sim N(\eta, \sigma^2)\\; Poisson and Binomial override it with their
own link and sampler, so a predictive draw from a count model is a
count, not a Gaussian jitter around the mean.

## Usage

``` r
responseRng(spec, eta, s2 = NULL, prior = NULL, ...)

# S4 method for class 'DistributionSpec'
responseRng(spec, eta, s2 = NULL, prior = NULL, ...)

# S4 method for class 'PoissonRegSpec'
responseRng(spec, eta, s2 = NULL, prior = NULL, ...)

# S4 method for class 'BinomialRegSpec'
responseRng(spec, eta, s2 = NULL, prior = NULL, ...)

# S4 method for class 'StudentTRegSpec'
responseRng(spec, eta, s2 = NULL, prior = NULL, ...)

# S4 method for class 'NormalGammaRegSpec'
responseRng(spec, eta, s2 = NULL, prior = NULL, ...)
```

## Arguments

- spec:

  A distribution spec.

- eta:

  Numeric vector of linear predictors.

- s2:

  Numeric vector of error variances (Normal only; ignored elsewhere).

- prior:

  Optional prior list (e.g. Binomial `size`).

- ...:

  Unused.

## Value

A numeric vector of draws, one per element of `eta`.

## Methods (by class)

- `responseRng(DistributionSpec)`: Gaussian: identity link, Normal
  noise.

- `responseRng(PoissonRegSpec)`: Poisson: log link, Poisson draw.

- `responseRng(BinomialRegSpec)`: Binomial: logit link, Binomial draw.

- `responseRng(StudentTRegSpec)`: Student-t regression: identity link,
  scaled-t noise. Both the direct Student-t and the Normal-Gamma
  augmentation have the same Student-t marginal, so they share this.
  `s2` is the error VARIANCE, which for a t with `df` degrees of freedom
  is \\\sigma^2 \\ \mathrm{df}/(\mathrm{df}-2)\\; the draw is \\\eta +
  \sigma \\ t\_{\mathrm{df}}\\ with \\\sigma = \sqrt{s2
  (\mathrm{df}-2)/\mathrm{df}}\\.

- `responseRng(NormalGammaRegSpec)`: Normal-Gamma regression: same
  Student-t marginal as the direct parameterisation, so the same
  scaled-t draw. Defined explicitly because the two heavy-tail specs are
  siblings under `NormalRegSpec`, not parent and child – inheritance
  would give the Gaussian default and silently drop the tails.
