# Shape the observed data into the NIMBLE `data` list

Shape the observed data into the NIMBLE `data` list

## Usage

``` r
buildDataList(spec, data, ...)

# S4 method for class 'PoissonRegSpec'
buildDataList(spec, data, ...)

# S4 method for class 'BinomialRegSpec'
buildDataList(spec, data, ...)

# S4 method for class 'GMSNBurrUvSpec'
buildDataList(spec, data, ...)

# S4 method for class 'NormalRegSpec'
buildDataList(spec, data, ...)

# S4 method for class 'MSNBurrUvSpec'
buildDataList(spec, data, ...)

# S4 method for class 'MSNBurr2aUvSpec'
buildDataList(spec, data, ...)

# S4 method for class 'NormalMvRegSpec'
buildDataList(spec, data, ...)

# S4 method for class 'NormalMvSpec'
buildDataList(spec, data, ...)

# S4 method for class 'NormalUvSpec'
buildDataList(spec, data, ...)

# S4 method for class 'PoissonSpec'
buildDataList(spec, data, ...)

# S4 method for class 'BinomialSpec'
buildDataList(spec, data, ...)

# S4 method for class 'StudentTUvSpec'
buildDataList(spec, data, ...)
```

## Arguments

- spec:

  A
  [`DistributionSpec`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md).

- data:

  The observed data (vector or matrix).

- ...:

  Reserved for methods.

## Value

A named list, typically `list(y = ...)`.

## Functions

- `buildDataList(PoissonRegSpec)`: Response and design matrix.

- `buildDataList(BinomialRegSpec)`: Response and design matrix.

- `buildDataList(GMSNBurrUvSpec)`: Univariate response.

- `buildDataList(NormalRegSpec)`: Response vector for the regression
  mixture.

- `buildDataList(MSNBurrUvSpec)`: Univariate response.

- `buildDataList(MSNBurr2aUvSpec)`: Univariate response.

- `buildDataList(NormalMvRegSpec)`: Matrix response.

- `buildDataList(NormalMvSpec)`: Multivariate data matrix (one row per
  observation).

- `buildDataList(NormalUvSpec)`: Univariate data vector.

- `buildDataList(PoissonSpec)`: Count data vector.

- `buildDataList(BinomialSpec)`: Count data vector.

- `buildDataList(StudentTUvSpec)`: Univariate data vector.
