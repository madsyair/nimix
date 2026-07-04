# Validate Distribution Parameters

Internal function to validate parameters common to all neo-normal
distributions.

## Usage

``` r
.validate_params(mu, sigma, alpha, theta)
```

## Arguments

- mu:

  Location parameter.

- sigma:

  Scale parameter (must be positive).

- alpha:

  Shape parameter (must be positive).

- theta:

  Shape parameter (must be positive).

## Value

Invisible NULL if valid, otherwise stops with error.
