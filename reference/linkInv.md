# Inverse link for a regression component

Maps the linear predictor to the response mean. Default is the identity
link (Normal-linear); GLM specs override it (log, logit).

## Usage

``` r
linkInv(spec, eta, prior = NULL, ...)

# S4 method for class 'DistributionSpec'
linkInv(spec, eta, prior = NULL, ...)

# S4 method for class 'PoissonRegSpec'
linkInv(spec, eta, prior = NULL, ...)

# S4 method for class 'BinomialRegSpec'
linkInv(spec, eta, prior = NULL, ...)
```

## Arguments

- spec:

  A
  [`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md).

- eta:

  Linear predictor value(s).

- prior:

  Optional prior list (for e.g. the Binomial `size`).

- ...:

  Unused.

## Value

The response mean.

## Functions

- `linkInv(DistributionSpec)`: Identity link (Normal-linear).

- `linkInv(PoissonRegSpec)`: Log link inverse (`exp`).

- `linkInv(BinomialRegSpec)`: Logit link inverse (`size * plogis`).
