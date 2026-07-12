# Most probable state path of a hidden-Markov mixture fit

Computes the Viterbi (maximum a posteriori joint) state sequence of a
`method = "hmm"` fit, by default at the posterior means of the state
parameters and transition matrix.

## Usage

``` r
viterbiPath(fit)
```

## Arguments

- fit:

  A `FitResult` from `nimixClust(..., method = "hmm")`.

## Value

Integer vector of length \\n\\: the decoded state per time point.

## Details

Note the difference from
[`binderPartition()`](https://madsyair.github.io/nimix/reference/binderPartition.md):
Viterbi gives the single jointly most probable *path* under the Markov
prior, while the Binder partition summarises marginal co-clustering
across all FFBS draws. They usually agree on well-separated regimes and
differ exactly where the state is genuinely uncertain.
