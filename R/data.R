## R/data.R -------------------------------------------------------------------
## Documentation for datasets shipped with nimix.

#' World Development Indicators, 2022 (country-level official statistics)
#'
#' A cross-section of four official development indicators for 207 countries
#' and territories in 2022, retrieved from the World Bank's World Development
#' Indicators (WDI) API. The dataset is used by the package vignettes as a
#' realistic official-statistics example for mixture clustering (countries
#' group into latent development regimes) and mixture-of-regressions (the
#' Preston-curve relationship between national income and life expectancy).
#'
#' @format A data frame with 207 rows and 9 columns:
#' \describe{
#'   \item{country}{Country or territory name (World Bank).}
#'   \item{iso3}{ISO 3166-1 alpha-3 code.}
#'   \item{region}{World Bank geographic region.}
#'   \item{income}{World Bank income classification.}
#'   \item{gdp_pc}{GDP per capita, constant 2015 US$ (indicator
#'     \code{NY.GDP.PCAP.KD}).}
#'   \item{life_exp}{Life expectancy at birth, total, years
#'     (\code{SP.DYN.LE00.IN}).}
#'   \item{fertility}{Total fertility rate, births per woman
#'     (\code{SP.DYN.TFRT.IN}).}
#'   \item{urban}{Urban population, percent of total (\code{SP.URB.TOTL.IN.ZS}).}
#'   \item{year}{Observation year (2022 for all rows).}
#' }
#'
#' @details
#' Retrieved 2026-07-02 from \code{https://api.worldbank.org/v2/} (year 2022;
#' aggregates such as regions and income groups excluded via the country
#' metadata endpoint; only countries with complete values on all four
#' indicators retained). One record (Central African Republic) was excluded
#' because its retrieved 2022 life-expectancy value (18.8 years) is an evident
#' source-data artifact.
#'
#' World Bank data are published under the Creative Commons Attribution 4.0
#' license (CC BY 4.0); see
#' \code{https://datacatalog.worldbank.org/public-licenses}.
#'
#' @source World Bank, World Development Indicators.
#'   \code{https://databank.worldbank.org/source/world-development-indicators}
#'
#' @examples
#' data(wdi2022)
#' summary(log(wdi2022$gdp_pc))
"wdi2022"

#' US state poverty and income, 2023 (SAIPE official estimates)
#'
#' Median household income and the all-ages poverty rate for the 48 contiguous
#' United States plus the District of Columbia (49 regions), from the U.S.
#' Census Bureau's Small Area Income and Poverty Estimates (SAIPE) programme,
#' 2023 vintage. Used together with \code{\link{usStateAdj}} by the
#' spatial-mixture vignette.
#'
#' @format A data frame with 49 rows and 5 columns:
#' \describe{
#'   \item{state}{State name (plus the District of Columbia).}
#'   \item{postal}{Two-letter postal abbreviation.}
#'   \item{fips}{Two-digit state FIPS code.}
#'   \item{medianIncome}{SAIPE median household income estimate, 2023 (US$).}
#'   \item{povertyRate}{SAIPE poverty rate estimate, all ages, 2023 (percent).}
#' }
#' @details Parsed on 2026-07-03 from the official fixed-layout file
#' \code{est23all.txt} (state rows). Alaska, Hawaii and the territories are
#' excluded so that the companion contiguity graph is connected.
#' @source U.S. Census Bureau, SAIPE 2023,
#'   \code{https://www2.census.gov/programs-surveys/saipe/datasets/2023/2023-state-and-county/}
"usStates2023"

#' Contiguity of the 48 contiguous US states + DC (official derivation)
#'
#' Symmetric binary adjacency matrix (49 x 49, postal-code dimnames) between
#' the contiguous United States and the District of Columbia, derived from the
#' U.S. Census Bureau's 2023 county adjacency file: two states are adjacent
#' when any pair of their counties is listed as adjacent across the state
#' line. The Census county file's conventions (including some water
#' adjacencies) are inherited as-is.
#'
#' @format A 49 x 49 numeric 0/1 matrix; rows/columns ordered by state FIPS
#'   and named by postal code. Use \code{spatialWeights(usStateAdj)} to obtain
#'   a \code{\linkS4class{SpatialWeightSpec}}.
#' @details Derived on 2026-07-03 from \code{county_adjacency2023.txt};
#' 112 undirected edges; e.g. Tennessee has the well-known 8 neighbours and
#' Maine exactly 1.
#' @source U.S. Census Bureau, county adjacency file,
#'   \code{https://www2.census.gov/geo/docs/reference/county_adjacency/}
"usStateAdj"
