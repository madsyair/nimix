# Internal validity indices for a fitted clustering

Computes standard internal cluster-validity indices – silhouette width,
Dunn index, Calinski-Harabasz – for the point partition of a clustering
fit. By default the partition is \[binderPartition()\]'s least-squares
partition, so every posterior draw informs it and no relabelling is
required.

## Usage

``` r
clusterValidity(
  fit,
  metrics = c("silhouette", "dunn", "ch"),
  partition = NULL,
  dist = NULL
)
```

## Arguments

- fit:

  A \[FitResult\]\[FitResult-class\] from a \*clustering\* fit
  (\[nimixClust()\]). Regression fits are refused: distance-based
  indices have no meaning for mixtures of regressions.

- metrics:

  Character vector, any of \`"silhouette"\`, \`"dunn"\`, \`"ch"\`.

- partition:

  Optional integer vector of cluster labels (length \`n\`). Default:
  \`binderPartition(fit)\$partition\`.

- dist:

  Optional \`stats::dist\` object. Default: Euclidean distance on the
  fitted data. Supply your own for scaled or non-Euclidean analyses.

## Value

Named numeric vector with one entry per requested metric.

## Interpretation caveat

Internal indices reward \*geometric\* separation. A mixture model is
\*density\*-based, and a fit with genuinely overlapping components –
often the scientifically correct model – will score a low silhouette
even when \[ppCheck()\] says the model reproduces the data perfectly.
Treat these indices as a secondary lens (useful, e.g., for comparing
partitions from different \`K\` on equal footing), never as a
model-adequacy verdict; that job belongs to \[ppCheck()\] and \[psm()\].

For the full battery of indices beyond these three, call the backends
directly, e.g. \`fpc::cluster.stats(d, part)\` or
\`clusterCrit::intCriteria(X, part, "all")\` – both accept exactly the
\`(partition, distance/data)\` pair this function assembles.

## See also

\[binderPartition()\], \[psm()\], \[ppCheck()\]

## Examples

``` r
# \donttest{
y <- c(rnorm(60, -3), rnorm(60, 3))
fit <- nimixClust(y, K_max = 6,
                  mcmcControl = list(niter = 800, nburnin = 300),
                  verbose = FALSE)
clusterValidity(fit)
#>   silhouette         dunn           ch 
#>    0.7993590    0.1217575 1049.4718680 
# }
```
