# Build a space-time adjacency for spatio-temporal mixtures

Expands a spatial adjacency over `nTime` time points, so that node \\(i,
t)\\ neighbours its spatial neighbours at the same time and itself at
\\t \pm 1\\. The result is an ordinary
[`SpatialWeightSpec`](https://madsyair.github.io/nimix/reference/SpatialWeightSpec-class.md):
pass it to `nimixClust(..., method = "mrf", spatialWeights = )` and the
Potts prior couples allocations across space *and* time, with no other
change.

## Usage

``` r
spacetimeAdjacency(spaceWeights, nTime, spatial = TRUE, temporal = TRUE)
```

## Arguments

- spaceWeights:

  A
  [`SpatialWeightSpec`](https://madsyair.github.io/nimix/reference/SpatialWeightSpec-class.md)
  over \\n\_{loc}\\ locations, e.g. from
  [`gridAdjacency`](https://madsyair.github.io/nimix/reference/gridAdjacency.md).

- nTime:

  Integer, number of time points (\>= 2).

- spatial:

  Logical; keep the spatial edges within each time point. `FALSE` gives
  a temporal-only graph (independent chains per location).

- temporal:

  Logical; keep the edges linking consecutive time points.

## Value

A
[`SpatialWeightSpec`](https://madsyair.github.io/nimix/reference/SpatialWeightSpec-class.md)
over \\n\_{loc} \times n\_{time}\\ nodes.

## Details

Observations must be ordered with time varying slowest – node \\(i, t)\\
is row \\(t - 1) n\_{loc} + i\\ – which is what
[`as.vector()`](https://rdrr.io/r/base/vector.html) of a \\n\_{loc}
\times n\_{time}\\ matrix gives.

On a 5x5 grid over 8 time points with deliberately overlapping
components, a plain mixture failed outright (allocation accuracy 0.50,
one component recovered); the spatial graph reached 0.95; the space-time
graph reached 0.995. The temporal edges earn their place.

Note the coupling is isotropic: one `beta` governs spatial and temporal
edges alike, because the Potts prior reads the adjacency as unweighted.
Misspecifying it is mild rather than fatal – on regimes that were random
across space but perfectly persistent in time, imposing the spatial
edges anyway cost 2.5 percentage points (0.900 against 0.925 for a
temporal-only graph). If your structure is purely temporal, build the
temporal-only graph with `spatial = FALSE`.

For a pure time series with no spatial component, `method = "hmm"` is
the better tool: it marginalises the state path (better mixing) and
offers
[`viterbiPath`](https://madsyair.github.io/nimix/reference/viterbiPath.md).
This function is for when space matters too.

## Scale limit

The graph itself no longer limits the problem size: `SpatialWeightSpec`
stores an edge list, so a 10 000-node space-time graph builds in a
fraction of a second within a few tens of MB (it was OOM-killed outright
before v1.5.0, during construction), and 50 000 nodes take about a
second. The binding constraint is now NIMBLE's per-node model memory
during `nimbleModel`/compilation: measured in a 4 GB container, fits ran
at 5000 nodes (~1.6 GB) and 7000 nodes (~2.8 GB) but died near 10 000.
If you need more, more RAM buys it roughly linearly – the graph is no
longer the wall.

## See also

[`gridAdjacency`](https://madsyair.github.io/nimix/reference/gridAdjacency.md),
[`viterbiPath`](https://madsyair.github.io/nimix/reference/viterbiPath.md)

## Examples

``` r
W  <- gridAdjacency(3, 3)
ST <- spacetimeAdjacency(W, nTime = 4)
nrow(getAdjacency(ST))   # 36 = 9 locations x 4 times
#> [1] 36
```
