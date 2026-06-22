# Heavy-tailed components: Student-t vs Normal-Gamma

`nimix` offers two ways to put **heavy-tailed (Student-t) components**
in a mixture. They have the *same* marginal distribution; they differ
only in how the sampler reaches it.

- `distribution = "studentt"` evaluates the **t density directly**. The
  cluster parameters are updated non-conjugately. For multivariate data
  the t density is supplied by `nimix` as a user-defined NIMBLE
  distribution, since NIMBLE has no built-in multivariate-t.
- `distribution = "normalgamma"` uses the **conjugate scale-mixture**
  representation: a per-observation latent precision multiplier `omega`
  with `y ~ N(mu, Sigma / omega)`, `omega ~ Gamma(df/2, df/2)`.
  Conditional on `omega` the kernel is Gaussian, so the conjugate
  cluster updates are kept, which is cheaper per iteration. The latent
  `omega` also act as robustness weights.

The univariate or multivariate variant is chosen automatically from the
data shape. `df` is a fixed hyperparameter (must exceed 2).

> These are the **same distribution, two routes** – not two different
> distributions. They are also unrelated to the “Normal-Gamma”
> *shrinkage prior* on coefficients, which is a different idea.

## Same marginal

``` r

mu <- 1.5; s2 <- 2.3; df <- 5
mix <- function(y) integrate(function(w)
  dnorm(y, mu, sqrt(s2 / w)) * dgamma(w, df / 2, df / 2), 0, Inf)$value
tdn <- function(y) dt((y - mu) / sqrt(s2), df) / sqrt(s2)
ys <- c(-2, 0, 1.5, 4)
data.frame(y = ys, scale_mixture = sapply(ys, mix), student_t = sapply(ys, tdn))
# the two columns agree to numerical-integration error
```

## Fitting both on the same data

``` r

library(nimix)
set.seed(1)
y <- c(rt(120, df = 4) - 5, rt(120, df = 4) + 5)   # two heavy-tailed clusters

fit_t <- nimixClust(y, K_max = 8, distribution = "studentt", prior = list(df = 4),
                    mcmcControl = list(niter = 4000, nburnin = 1000),
                    verbose = FALSE)
fit_ng <- nimixClust(y, K_max = 8, distribution = "normalgamma", prior = list(df = 4),
                     mcmcControl = list(niter = 4000, nburnin = 1000),
                     verbose = FALSE)
summary(fit_t)
summary(fit_ng)
```

Both recover the same two clusters. To compare sampling efficiency on
your own problem, time each fit and divide by the effective sample size
of the parameters of interest (effective samples per second); the
Normal-Gamma route is usually cheaper per iteration, while the direct-t
route stores no latent `omega`.

## Which to choose – measured, not assumed

The honest comparison is by effective samples per second of the number
of clusters `K` (a label-invariant quantity). The script
`inst/harness/run_benchmark_heavytail.R` runs this; representative
numbers on one machine (R 4.3.3, NIMBLE 1.4.2, 6000 iterations) were:

| task                      | route        | ESS(K)/s | wall time |
|---------------------------|--------------|---------:|----------:|
| univariate **clustering** | Student-t    |      ~41 |     ~39 s |
| univariate **clustering** | Normal-Gamma |       ~6 |     ~41 s |
| univariate **regression** | Student-t    |      ~13 |     ~51 s |
| univariate **regression** | Normal-Gamma |      ~14 |     ~58 s |

The pattern is not “one route always wins”:

- For **clustering**, the direct Student-t route mixes the partition far
  better (here about 7x the ESS per second). Augmenting with the latent
  `omega` couples the partition to a high-dimensional latent vector and
  slows mixing, even with `omega` slice-sampled (van Dyk & Meng 2001).
- For **regression**, the two are comparable. The conjugate Gibbs update
  of the `p` regression coefficients (which Normal-Gamma keeps and the
  direct t route gives up) is worth more here, roughly offsetting the
  augmentation penalty.

So: prefer **Student-t for clustering**; treat the two as
**interchangeable for regression**, choosing Normal-Gamma when you also
want the robustness weights `omega`. Always confirm on your own data
with the benchmark script rather than assuming.

For regression with heavy-tailed multivariate errors, note the
identifiability “pitfalls” with small `df` discussed by Fernández &
Steel (1999): small `df` can make some quantities weakly identified – a
property of the model, not a bug.
