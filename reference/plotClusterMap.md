# Map the clusters of a spatial mixture fit

Draws the fitted partition as a choropleth: each map feature (region) is
coloured by its cluster label. Optionally shades each region by
allocation uncertainty, so that regions whose membership the posterior
is unsure about stand out.

## Usage

``` r
plotClusterMap(
  fit,
  shp,
  partition = c("binder", "modal"),
  idCol = NULL,
  uncertainty = FALSE,
  palette = NULL,
  main = NULL,
  legendPos = "topright",
  ...
)
```

## Arguments

- fit:

  A clustering
  [`FitResult`](https://madsyair.github.io/nimix/reference/FitResult-class.md)
  (any engine; the MRF engine is the usual source).

- shp:

  The map: either a path to a shapefile (`.shp`, read with
  [`sf::st_read`](https://r-spatial.github.io/sf/reference/st_read.html))
  or an `sf` object already in memory. Features must correspond to
  observations.

- partition:

  How to summarise the posterior into one label per region: `"binder"`
  (default) uses
  [`binderPartition`](https://madsyair.github.io/nimix/reference/binderPartition.md),
  the draw minimising Binder loss against the posterior similarity
  matrix – a label-invariant summary; `"modal"` uses the per-region
  posterior mode of the allocation trace.

- idCol:

  Optional name of a column in `shp` giving the observation index
  (1-based) of each feature. Default `NULL` assumes the features are
  already in observation order.

- uncertainty:

  Logical; if `TRUE`, regions are shaded towards white in proportion to
  their allocation entropy (0 = certain, log K = maximally uncertain),
  so pale regions are the ones the posterior cannot place.

- palette:

  Colours, one per cluster. Default is
  `grDevices::hcl.colors(K, "Dark 3")`.

- main:

  Plot title.

- legendPos:

  Legend position keyword (see
  [`legend`](https://rdrr.io/r/graphics/legend.html)), or `NA` to
  suppress the legend.

- ...:

  Passed on to `plot` of the sf geometry (e.g. `border`, `lwd`).

## Value

Invisibly, a data.frame with one row per region: `cluster` and (if
requested) `entropy`. Called for its side effect, the plot.

## Examples

``` r
if (FALSE) { # \dontrun{
sw  <- spatialWeights(nb)                      # neighbourhood used in the fit
fit <- nimixClust(y, K = 3, method = "mrf", spatialWeights = sw)
plotClusterMap(fit, "regions.shp")             # from a shapefile on disk
plotClusterMap(fit, mysf, uncertainty = TRUE)  # from an sf object
} # }
```
