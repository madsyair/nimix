# Diagnostics, partitions, and bayesplot interoperability

Fitting a mixture is half the job; interpreting the posterior is the
other half, and mixtures make it subtle. The component labels are
exchangeable, so “component 1” means nothing across draws until you
condition on something, and under a Dirichlet process the *number* of
components moves too. This article shows the tools `nimix` provides for
that, and how to hand a fit to
[`bayesplot`](https://mc-stan.org/bayesplot/) without taking it on as a
dependency.

> MCMC chunks use `eval = FALSE` (CRAN time limits); the printed results
> are from an actual run.

``` r

library(nimix)

set.seed(1)
y <- c(rnorm(120, -3, 0.8), rnorm(90, 0, 0.6), rnorm(120, 3, 0.8))
fit <- nimixClust(y, method = "dpm", K_max = 10,
                  mcmcControl = list(niter = 2000, nburnin = 800, nchains = 2),
                  seed = 1)
```

## Two questions, two different tools

There are two questions one usually asks of a clustering posterior, and
they need different summaries.

**“What are the component parameters?”** This needs labels aligned
across draws, which
[`relabel()`](https://madsyair.github.io/nimix/reference/relabel.md)
does by conditioning on the modal number of occupied clusters and then
applying a label-switching algorithm. Conditioning is unavoidable – a
draw with four clusters has no “third component” to align with a draw
that has three – but it means some draws are set aside:

``` r

Kpost <- apply(fit@clusterAllocation, 1, function(z) length(unique(z)))
table(Kpost)
#> Kpost
#>    3    4    5    6    7    8
#> 1215  748  317   93   21    6
```

Here only about half the draws sit at the modal `K = 3`; the rest inform
`K` and the partition but not the per-component summaries.

**“Which observations belong together?”** This question does not need
labels at all, because co-clustering is label-invariant: whether
observations *i* and *j* share a component does not depend on what that
component is called.
[`psm()`](https://madsyair.github.io/nimix/reference/psm.md) estimates
the posterior similarity matrix $`S_{ij} = \Pr(z_i = z_j \mid y)`$ from
**every** draw, and
[`binderPartition()`](https://madsyair.github.io/nimix/reference/binderPartition.md)
returns the point partition minimising the expected Binder (1978) loss –
equivalently, Dahl’s (2006) least-squares criterion, searched over the
partitions the chain actually visited.

``` r

S  <- psm(fit)
bp <- binderPartition(fit, S)
bp$K
#> [1] 3
table(bp$partition)
#>
#>   1   2   3
#> 118  98 114
```

No [`relabel()`](https://madsyair.github.io/nimix/reference/relabel.md)
call was needed, and no draws were discarded. The two tools are
complements: use
[`relabel()`](https://madsyair.github.io/nimix/reference/relabel.md) for
component parameters,
[`psm()`](https://madsyair.github.io/nimix/reference/psm.md)/
[`binderPartition()`](https://madsyair.github.io/nimix/reference/binderPartition.md)
for the partition. On overlapping clusters the similarity matrix is
especially informative, because it reports genuine allocation
uncertainty (entries near 0.5) instead of forcing a hard assignment.

## Convergence diagnostics that are actually meaningful

Reaching for R-hat on `muTilde[1]` is a trap: under label switching that
trace refers to different components in different chains, so the number
looks plausible and means nothing. `nimix` steers you to
**label-invariant functionals** – the number of occupied clusters, the
allocation entropy, and (for DPM fits) the concentration `alpha` – which
are meaningful on raw draws and retain the per-chain structure that
R-hat needs.

[`drawsArray()`](https://madsyair.github.io/nimix/reference/drawsArray.md)
returns these in the `iterations x chains x parameters` layout that
`bayesplot`’s `mcmc_*` functions accept natively:

``` r

da <- drawsArray(fit)                 # params = "invariant" by default
dim(da)
#> [1] 1200    2    3
dimnames(da)[[3]]
#> [1] "K"       "entropy" "alpha"
```

Ask for component parameters before relabelling and the adaptor refuses,
on purpose:

``` r

drawsArray(fit, "components")
#> Error: drawsArray(params = "components") needs relabel() first.
#>   Component draws are not identified before relabelling: ...
```

After `fit <- relabel(fit)` the component view is available; because
relabelling conditions on the modal `K`, chains no longer have equal
lengths, so the chain dimension is honestly collapsed to one – fine for
posterior density and interval plots, not for cross-chain R-hat.

## bayesplot, without the dependency

`bayesplot` sits in `Suggests`. The adaptors return base R objects, so
nothing is pulled into your namespace until you load `bayesplot`
yourself.

``` r

library(bayesplot)

# convergence on the label-invariant functionals
mcmc_trace(drawsArray(fit))
mcmc_rhat_hist(rhat(drawsArray(fit)))

# posterior predictive overlay
pd <- ppcData(fit, ndraws = 50)       # list(y, yrep)
ppc_dens_overlay(pd$y, pd$yrep)
```

[`ppcData()`](https://madsyair.github.io/nimix/reference/ppcData.md)
wraps
[`posteriorPredict()`](https://madsyair.github.io/nimix/reference/posteriorPredict.md),
which returns the replicated data sets themselves (an `ndraws x n`
matrix, or `ndraws x n x d` for multivariate fits). For a numeric
summary rather than a plot,
[`ppCheck()`](https://madsyair.github.io/nimix/reference/ppCheck.md)
gives tail-area posterior predictive p-values for a set of statistics:

``` r

ppCheck(fit, nrep = 200)
#>   statistic  observed  repMean   ppp
#> 1      mean    ...      ...      0.490
#> 2        sd    ...      ...      0.710
#> 3       min    ...      ...      0.130
#> 4       max    ...      ...      0.735
```

Values away from 0 or 1 indicate the replicated data reproduce that
statistic of the observed data – here all four are comfortably interior.

## Every plot hands back its data

The built-in [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
methods draw with base graphics but return, invisibly, the tidy data
frame they plotted – so you can reproduce or restyle any of them with
ggplot2, lattice, or plotly without `nimix` depending on those packages.

``` r

d <- plot(fit, type = "density")      # draws the base plot
head(d)                               # ... and returns x / density
library(ggplot2)
ggplot(d, aes(x, density)) + geom_line()
```

## A secondary lens: internal validity indices

[`clusterValidity()`](https://madsyair.github.io/nimix/reference/clusterValidity.md)
bridges the Binder partition to the classical internal indices –
silhouette width (`cluster`), Dunn and Calinski-Harabasz (`fpc`), both
Suggests-only:

``` r

clusterValidity(fit)
#> silhouette       dunn         ch
#>      0.903      1.616   4299.717      # well-separated example
```

One caveat matters enough to repeat. These indices reward *geometric*
separation; a mixture is *density*-based. A fit with genuinely
overlapping components – often the scientifically correct model – scores
low:

``` r

# two components centred at -0.7 and +0.7: a legitimate model of this data
clusterValidity(fit_overlap)
#> silhouette       dunn         ch
#>      0.522      0.000    224.123
```

Low silhouette here is not evidence against the model –
[`ppCheck()`](https://madsyair.github.io/nimix/reference/ppCheck.md) on
the same fit is clean. Use these indices to compare partitions (say,
across candidate `K`) on an equal footing, and keep model adequacy
questions with
[`ppCheck()`](https://madsyair.github.io/nimix/reference/ppCheck.md) and
[`psm()`](https://madsyair.github.io/nimix/reference/psm.md). For the
full battery beyond these three, both backends accept the pair this
function assembles:
`fpc::cluster.stats(dist(fit@data), binderPartition(fit)$partition)`.

## References

Binder, D. A. (1978). Bayesian cluster analysis. *Biometrika* 65, 31–38.

Dahl, D. B. (2006). Model-based clustering for expression data via a
Dirichlet process mixture model. In *Bayesian Inference for Gene
Expression and Proteomics*, Cambridge University Press, 201–218.

Gabry, J., Simpson, D., Vehtari, A., Betancourt, M. & Gelman, A. (2019).
Visualization in Bayesian workflow. *JRSS A* 182, 389–402.
