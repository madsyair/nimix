# Hidden-Markov mixture engine (regime switching in time)

Component labels follow a first-order Markov chain in the order of the
observations (which must therefore be time order), instead of being
independent: \\z_t \mid z\_{t-1} \sim \mathrm{Cat}(P\[z\_{t-1}, \])\\.
The state path is marginalised out of the likelihood by the forward
algorithm, so the MCMC samples only the continuous parameters;
allocation draws are then recovered exactly by forward-filter
backward-sampling (FFBS) per retained draw, which is what makes every
downstream tool (`relabel`, `psm`, `binderPartition`, plots) work
unchanged.

## Slots

- `transConc`:

  Positive scalar: symmetric Dirichlet concentration of the prior on
  each row of the transition matrix. Values above 1 favour
  persistence-agnostic rows; the default 1 is uniform on each row
  simplex.
