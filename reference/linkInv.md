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

# S4 method for class 'MSNBurrRegSpec'
linkInv(spec, eta, prior = NULL, ...)

# S4 method for class 'SEPRegSpec'
linkInv(spec, eta, prior = NULL, ...)

# S4 method for class 'MSNBurr2aRegSpec'
linkInv(spec, eta, prior = NULL, ...)

# S4 method for class 'FSSNRegSpec'
linkInv(spec, eta, prior = NULL, ...)

# S4 method for class 'GMSNBurrRegSpec'
linkInv(spec, eta, prior = NULL, ...)

# S4 method for class 'LEPRegSpec'
linkInv(spec, eta, prior = NULL, ...)

# S4 method for class 'FSSTRegSpec'
linkInv(spec, eta, prior = NULL, ...)

# S4 method for class 'FOSSEPRegSpec'
linkInv(spec, eta, prior = NULL, ...)

# S4 method for class 'JFSTRegSpec'
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

- `linkInv(MSNBurrRegSpec)`: MSNBurr regression identity link.

- `linkInv(SEPRegSpec)`: SEP regression identity link.

- `linkInv(MSNBurr2aRegSpec)`: MSNBurr-IIa regression identity link.

- `linkInv(FSSNRegSpec)`: FSSN regression identity link.

- `linkInv(GMSNBurrRegSpec)`: GMSNBurr regression identity link.

- `linkInv(LEPRegSpec)`: LEP regression identity link.

- `linkInv(FSSTRegSpec)`: FSST regression identity link.

- `linkInv(FOSSEPRegSpec)`: FOSSEP regression identity link.

- `linkInv(JFSTRegSpec)`: JFST regression identity link.
