# US state poverty and income, 2023 (SAIPE official estimates)

Median household income and the all-ages poverty rate for the 48
contiguous United States plus the District of Columbia (49 regions),
from the U.S. Census Bureau's Small Area Income and Poverty Estimates
(SAIPE) programme, 2023 vintage. Used together with
[`usStateAdj`](https://madsyair.github.io/nimix/reference/usStateAdj.md)
by the spatial-mixture vignette.

## Usage

``` r
usStates2023
```

## Format

A data frame with 49 rows and 5 columns:

- state:

  State name (plus the District of Columbia).

- postal:

  Two-letter postal abbreviation.

- fips:

  Two-digit state FIPS code.

- medianIncome:

  SAIPE median household income estimate, 2023 (US\$).

- povertyRate:

  SAIPE poverty rate estimate, all ages, 2023 (percent).

## Source

U.S. Census Bureau, SAIPE 2023,
`https://www2.census.gov/programs-surveys/saipe/datasets/2023/2023-state-and-county/`

## Details

Parsed on 2026-07-03 from the official fixed-layout file `est23all.txt`
(state rows). Alaska, Hawaii and the territories are excluded so that
the companion contiguity graph is connected.
