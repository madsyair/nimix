# Spatially constrained mixtures (MRF engine)

In a standard mixture the latent labels are independent across
observations. For regional data that is often unrealistic: neighbouring
regions tend to belong to the same regime. `method = "mrf"` replaces the
independent labels with a **Potts Markov random field** on a
neighbourhood graph (Potts 1952; Besag 1974; spatially variant finite
mixtures, Blekas et al. 2005): labels of adjacent regions attract each
other with interaction strength `beta`.

This vignette analyses **official statistics**: 2023 SAIPE state poverty
rates (U.S. Census Bureau) on the official state contiguity graph, both
shipped with the package
([`?usStates2023`](https://madsyair.github.io/nimix/reference/usStates2023.md),
[`?usStateAdj`](https://madsyair.github.io/nimix/reference/usStateAdj.md)).
The identical workflow applies to any regional official statistics –
e.g. BPS indicators on a kabupaten/kota contiguity graph.

> MCMC chunks use `eval = FALSE` (CRAN time limits); the printed results
> are from an actual run.

## The data and the graph

``` r

library(nimix)
#> Loading required package: nimble
#> nimble version 1.4.2 is loaded.
#> For more information on NIMBLE and a User Manual,
#> please visit https://R-nimble.org.
#> 
#> Attaching package: 'nimble'
#> The following object is masked from 'package:stats':
#> 
#>     simulate
#> The following object is masked from 'package:base':
#> 
#>     declare
data(usStates2023); data(usStateAdj)
sw <- spatialWeights(usStateAdj)
sw
#> SpatialWeightSpec: 49 regions, 112 undirected edges
#>   degree range: 1 - 8 | weights: binary (contiguity)
neighborsOf(sw, "TN")     # Tennessee's famous 8 neighbours
#> [1] "AL" "AR" "GA" "KY" "MS" "MO" "NC" "VA"
summary(usStates2023$povertyRate)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    7.30   10.30   11.90   12.27   13.70   18.90
```

## Honest statistical note, up front

The Potts normalising constant is intractable. With a **fixed** `beta`
(default `prior$beta = 0.8`) the constant is absorbed and the MCMC over
the labels and component parameters is exact. **Estimating** `beta`
(`prior$estimateBeta = TRUE`) uses the classical **pseudo-likelihood**
Metropolis update (Besag 1975) – an approximation, slightly biased for
strong interactions; an exchange-algorithm refinement is a possible
future upgrade.

## Two poverty regimes with spatial cohesion

``` r

y <- usStates2023$povertyRate
fit <- nimixClust(y, K = 2, method = "mrf", spatialWeights = sw,
                  prior = list(estimateBeta = TRUE),
                  mcmcControl = list(niter = 4000, nburnin = 1500),
                  seed = 1)
fit <- relabel(fit)
summary(fit)
```

Actual results on the 2023 data:

    #> regimes: 11.4% vs 14.5% poverty (weights 0.73 / 0.27)
    #> beta   : posterior mean 1.22, P(beta > 0.2) = 1  (clear positive interaction)
    #> high-poverty regime (13 regions):
    #>   AL AR DC FL GA LA MS NC NM OK SC TN TX

The high-poverty regime is the spatially **contiguous Southern belt**
(plus New Mexico, with DC an enclave inside the belt) – exactly the kind
of structure an exchangeable mixture cannot represent. A `fixedk` fit on
the same data finds similar regime means but agrees with the MRF
allocation for only 82% of states: the spatial smoothing snaps boundary
states to their neighbours’ regime.

``` r

zMap <- apply(fit@clusterAllocation, 2L,
              function(v) as.integer(names(which.max(table(v)))))
split(usStates2023$postal, zMap)
```

## Choosing and diagnosing `beta`

- `beta = 0` removes the smoothing (labels independent given the
  components); larger `beta` increases cohesion. When it is estimated,
  [`summary()`](https://rdrr.io/r/base/summary.html) reports split-Rhat
  / ESS and the posterior mean of `beta`.
- The engine warns when the posterior of `beta` piles up near its upper
  prior bound `prior$betaMax` (default 2): a (near-)saturated field –
  raise `betaMax` or interpret the smoothing as effectively maximal.

## Spatially clustered regressions

The same field applies to
[`nimixReg()`](https://madsyair.github.io/nimix/reference/nimixReg.md):
regression coefficients that cluster across adjacent regions
(e.g. growth patterns across neighbouring areas).

``` r

## regimes in how income relates to poverty across neighbouring states
df <- data.frame(pov = usStates2023$povertyRate,
                 inc = log(usStates2023$medianIncome))
fr <- nimixReg(pov ~ inc, df, K = 2, method = "mrf", spatialWeights = sw,
               prior = list(estimateBeta = TRUE),
               mcmcControl = list(niter = 4000, nburnin = 1500), seed = 1)
summary(relabel(fr))
```

On this state-level pair the fit collapses to a single national regime
(`pov ~ 144.6 - 11.8 log(income)`, the strong income-poverty gradient):
with only 49 regions and one dominant relationship, the model does not
invent a second regime – the same honesty the Preston-curve vignette
shows for the exchangeable case.

## Limitations

Estimated `beta` is a pseudo-likelihood approximation (see above); the
MRF engine currently supports Gaussian components (univariate,
multivariate, and univariate-response regression); heavy-tailed
components under the spatial engine are planned. Synthetic-lattice
recovery benchmarks for all three model types ship in
`inst/harness/run_recovery_suite.R` (groups `mrf`, `mrfmv`, `mrfreg`).
