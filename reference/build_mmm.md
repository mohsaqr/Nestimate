# Fit a Mixed Markov Model

`build_mmm()` fits a mixture of first-order Markov chains to a set of
sequences by expectation-maximisation: it looks for `k` subgroups whose
members move between states in different ways, each with its own
transition matrix and initial-state distribution. Every sequence gets a
posterior probability of belonging to each subgroup and is assigned to
the most probable one. Covariates can shape the mixing proportions
during estimation or be related to the subgroups afterwards (see
`covariate_effect`).

EM finds a local maximum, so the fit is started `n_starts` times (one
start from k-means on the transition counts, the rest random), each is
screened with a short EM run, the most promising are continued for up to
`max_iter` iterations, and the one with the highest log-likelihood is
kept. [`print()`](https://rdrr.io/r/base/print.html) reports whether it
converged within `max_iter`; a fit that did not is not a result. The
information criteria (BIC, AIC, ICL) compare fits with different `k` –
[`compare_mmm()`](https://pak.dynasite.org/Nestimate/reference/compare_mmm.md)
fits a range of them – while average posterior probability (AvePP),
entropy and classification error describe how cleanly the sequences
separate within one fit. None of them picks `k` alone: also check that
the subgroups' transition matrices differ in ways that can be
interpreted, and that each subgroup holds enough sequences to estimate
its matrix.

## Usage

``` r
build_mmm(
  data,
  k = 2L,
  n_starts = 50L,
  max_iter = 200L,
  tol = 1e-06,
  smooth = 0.01,
  seed = NULL,
  covariates = NULL,
  covariate_effect = c("em", "posthoc"),
  estimator = c("auto", "firth", "multinom", "chisq")
)

# S3 method for class 'net_mmm'
print(x, digits = 3L, ...)

# S3 method for class 'net_mmm'
summary(object, ...)

# S3 method for class 'net_mmm'
plot(x, type = c("posterior", "covariates"), combined = TRUE, ...)

# S3 method for class 'net_mmm_clustering'
print(x, digits = 3L, ...)

# S3 method for class 'net_mmm_clustering'
plot(
  x,
  type = c("posterior", "covariates", "predictors"),
  combined = TRUE,
  ...
)
```

## Arguments

- data:

  The sequences: a data frame with one row per sequence and one column
  per position, a `netobject` or bare `cograph_network` built from
  sequences, or a `tna` model (whose stored states are decoded to their
  labels). Columns take the same roles as in
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md):
  an id or a covariate is not a position, and the columns set aside are
  named in a message and kept as metadata. Association networks (`cor`,
  `pcor`, `glasso`, `ising`, `mgm`) hold variables, not sequences, and
  raise an error of class `nestimate_bad_input`.

- k:

  The number of subgroups (mixture components), a whole number of at
  least 2. Default `2`.

- n_starts:

  How many starting configurations to try, a positive whole number. More
  starts make it likelier that the best local maximum is found. Default
  `50`.

- max_iter:

  The most EM iterations a start may run, a positive whole number.
  Default `200`.

- tol:

  The convergence tolerance: EM stops when the log-likelihood changes by
  less than this between iterations (with `covariate_effect = "em"`,
  also when its relative change does). A finite positive number; default
  `1e-6`.

- smooth:

  A pseudo-count added to every transition count before the matrices are
  estimated, so a transition unseen in a subgroup gets a small
  probability instead of zero. A finite non-negative number; default
  `0.01`.

- seed:

  An integer that makes the random starts reproducible. Default `NULL`.
  The session's own random stream is left as it was.

- covariates:

  Variables describing each sequence (such as a course or a grade),
  given as a column name, a character vector of names, a one-sided
  formula or a data frame, as in
  [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md).
  For a `netobject`, names are looked up in its metadata first, so a
  typical call is `build_mmm(net, k = 3, covariates = "session_label")`.
  How they enter the model is set by `covariate_effect`. Default `NULL`.

- covariate_effect:

  How `covariates` enter the model. `"em"` (default) makes the mixing
  proportions depend on them during EM, so they shape the subgroups
  found and add `(k - 1) * p` parameters for `p` covariate columns; rows
  with a missing covariate are dropped first, with a warning.
  `"posthoc"` fits a plain mixture on every sequence and relates the
  covariates to the subgroups only afterwards, with a multinomial logit,
  so neither their values nor their missingness change the subgroups.
  Ignored without `covariates`.

- estimator:

  How the after-fit multinomial logit is estimated; it never affects the
  EM. `"auto"` (default) uses Firth's penalised likelihood
  ([`brglm2::brmultinom()`](https://rdrr.io/pkg/brglm2/man/brmultinom.html),
  finite under separation) when a cell of the subgroup-by-covariate
  table has fewer than 5 sequences, and
  [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html)
  otherwise. `"firth"` and `"multinom"` force one of them (`"multinom"`
  warns about separation risk), and `"chisq"` runs descriptive tests
  instead of a logit.
  [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  has the details. Ignored without `covariates`.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: a
  `net_mmm` fit, or the `net_mmm_clustering` attribute of the
  `netobject_group` built from one.

- digits:

  For the [`print()`](https://rdrr.io/r/base/print.html) methods:
  decimal places of the statistics. Default `3`.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: not
  used; any extra argument raises an error.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `net_mmm` fit.

- type:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  methods: `"posterior"` (default), a histogram of each sequence's
  largest posterior probability, coloured by cluster; or `"covariates"`,
  a forest plot of the covariate effects, for a model fitted with
  `covariates`. On a `net_mmm_clustering`, `"predictors"` is an alias of
  `"covariates"`.

- combined:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  methods with `type = "covariates"`: `TRUE` (default) draws the forest
  panels as one faceted plot; `FALSE` returns a list of ggplots, one per
  cluster.

## Value

A `net_mmm` object. [`print()`](https://rdrr.io/r/base/print.html) shows
the fit statistics, classification quality and subgroup sizes;
[`summary()`](https://rdrr.io/r/base/summary.html) adds each subgroup's
transition matrix and returns one row per subgroup (mixing proportion,
sequences assigned, mean posterior, AvePP);
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
posterior probabilities or the covariate effects;
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
turns the fit into one network per subgroup (a `netobject_group`);
[`cluster_diagnostics()`](https://pak.dynasite.org/Nestimate/reference/cluster_diagnostics.md)
and
[`session_ids()`](https://pak.dynasite.org/Nestimate/reference/session_ids.md)
give the classification and each sequence's subgroup. The object's
fields:

- data:

  The full N-row sequence frame used for estimation.

- models:

  List of `netobject`s, one per component. Each component carries the
  rows assigned to that component in its `$data` slot, while its
  transition matrix is the EM-estimated component transition matrix.

- k:

  Number of components.

- mixing:

  Numeric vector of mixing proportions.

- posterior:

  N x k matrix of posterior probabilities.

- assignments:

  Integer vector of hard assignments (1..k).

- quality:

  List: `avepp` (per-class), `avepp_overall`, `entropy`,
  `relative_entropy`, `classification_error`, `class_entropy`.

- log_likelihood, BIC, AIC, ICL:

  Model fit statistics.

- n_params:

  Number of free parameters behind BIC/AIC/ICL. With
  `covariate_effect = "em"` it grows by `(k - 1) * p` for `p` covariate
  columns.

- iterations, converged:

  EM iterations used by the retained fit and whether it met `tol`.

- states:

  Character vector of state names.

- n_sequences:

  Number of sequences actually fitted (rows with missing covariates are
  dropped under `covariate_effect = "em"`).

- covariates:

  The post-hoc covariate analysis (list), or NULL when
  `covariates = NULL`.

- network_method, build_args, htna_partition:

  Provenance kept from `netobject` / HTNA input so per-cluster networks
  can be rebuilt the same way; NULL otherwise.

- metadata:

  The `netobject`'s per-sequence metadata, one row per fitted sequence
  in the row order of `posterior`, so
  [`session_ids`](https://pak.dynasite.org/Nestimate/reference/session_ids.md)
  can name each sequence; NULL for other input.

The [`print()`](https://rdrr.io/r/base/print.html) methods return their
input invisibly. For a model fitted with `covariates`,
[`summary()`](https://rdrr.io/r/base/summary.html) also prints the
covariate profiles and effects and returns, invisibly, a
`tidy_covariates` data frame of the effects (one row per subgroup other
than the reference and per predictor term, so a factor gives one row per
non-reference level; one row per covariate test with
`estimator = "chisq"`), with the per-subgroup table as its
`cluster_stats` attribute. The
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods print
the plot and return the `ggplot` invisibly; with `type = "covariates"`
and `combined = FALSE` they return, invisibly and without printing, a
list of `ggplot`s named by subgroup.

## Initial states

Each sequence's initial state is the value in its first position column,
read as it is: the function does not look ahead to the first observed
state, and it does not convert missing-value symbols as
[`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
does. The states are the distinct non-missing values of all position
columns, so a placeholder such as `"*"` or `"%"` for an empty cell
becomes a state, and a sequence starting with it starts in that state.
Recode such placeholders to `NA` first: a missing first position gives
the sequence no initial state, which leaves its subgroup to the
transitions alone. Or trim leading empty cells so every first position
holds an observed state.

A one-hot network is different: its units are windows of codes, not
positions, so its initial state is the one code active in the first
window, and none when that window holds several codes or none.

## Methods

[`print()`](https://rdrr.io/r/base/print.html) of a `net_mmm` shows a
header with the numbers of sequences, subgroups and states, the
log-likelihood and the information criteria; a quality line with the
overall AvePP, the normalised entropy (0 when every sequence is
classified with certainty, 1 when all posteriors are equal) and the
classification error (the expected share of sequences misclassified when
each is assigned to its most likely subgroup, the mean of one minus the
largest posterior; Vermunt & Magidson, 2005); a line when the EM did not
converge; and one row per subgroup with its size, mixing proportion and
AvePP. For a model fitted with `covariates`, a last line names them,
says whether they entered the EM (integrated) or were related post hoc,
and gives the number of predictors. The layout follows
[`print()`](https://rdrr.io/r/base/print.html) of a `net_clustering`
from
[`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md),
so dissimilarity-based and model-based clusterings read alike.

[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
on a fit, and
[`cluster_network()`](https://pak.dynasite.org/Nestimate/reference/cluster_network.md)
with `cluster_by = "mmm"`, return a `netobject_group` whose
`"clustering"` attribute is a `net_mmm_clustering`: the fit without its
subgroup networks. Its [`print()`](https://rdrr.io/r/base/print.html)
shows the same quality line, information criteria and subgroup table,
without the state count, the log-likelihood or the convergence line, and
its [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
same two views as
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a `net_mmm`.
The `"silhouette"`, `"mds"` and `"heatmap"` views of a
dissimilarity-based clustering need a distance matrix, which an MMM
clustering does not have, so asking for them raises an error.

## References

Helske, S., & Helske, J. (2019). Mixture hidden Markov models for
sequence data: The seqHMM package in R. *Journal of Statistical
Software*, 88(3), 1–32.
[doi:10.18637/jss.v088.i03](https://doi.org/10.18637/jss.v088.i03)

Vermunt, J. K., & Magidson, J. (2005). *Latent GOLD 4.0 User's Guide*.
Statistical Innovations.

## See also

[`compare_mmm`](https://pak.dynasite.org/Nestimate/reference/compare_mmm.md),
[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md)

## Examples

``` r
seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
                   V2 = sample(c("A","B","C"), 30, TRUE))
mmm <- build_mmm(seqs, k = 2, n_starts = 1, max_iter = 10, seed = 1)
mmm
#> Mixed Markov Model
#>   Sequences: 30  |  Clusters: 2  |  States: 3
#>   ICs: LL = -62.330  |  BIC = 182.480  |  AIC = 158.659  |  ICL = 184.652
#>   Quality: AvePP = 0.965  |  Entropy = 0.212  |  Class.Err = 3.5%
#> 
#>   Cluster  N           Mix%   AvePP
#>   1        24 (80.0%)  78.1%  0.966
#>   2        6 (20.0%)   21.9%  0.961
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:3], 30, TRUE), V2 = sample(LETTERS[1:3], 30, TRUE),
  V3 = sample(LETTERS[1:3], 30, TRUE), V4 = sample(LETTERS[1:3], 30, TRUE)
)
mmm <- build_mmm(seqs, k = 2, seed = 42)
print(mmm)
#> Mixed Markov Model
#>   Sequences: 30  |  Clusters: 2  |  States: 3
#>   ICs: LL = -120.823  |  BIC = 299.466  |  AIC = 275.646  |  ICL = 299.990
#>   Quality: AvePP = 0.992  |  Entropy = 0.035  |  Class.Err = 0.8%
#> 
#>   Cluster  N           Mix%   AvePP
#>   1        12 (40.0%)  40.7%  0.999
#>   2        18 (60.0%)  59.3%  0.987
summary(mmm)
#> Mixed Markov Model
#>   Sequences: 30  |  Clusters: 2  |  States: 3
#>   ICs: LL = -120.823  |  BIC = 299.466  |  AIC = 275.646  |  ICL = 299.990
#>   Quality: AvePP = 0.992  |  Entropy = 0.035  |  Class.Err = 0.8%
#> 
#>   Cluster  N           Mix%   AvePP
#>   1        12 (40.0%)  40.7%  0.999
#>   2        18 (60.0%)  59.3%  0.987
#> 
#> --- Cluster 1 (40.7%, n=12) ---
#>       A     B     C
#> A 0.003 0.365 0.632
#> B 0.002 0.001 0.997
#> C 0.176 0.309 0.515
#> 
#> --- Cluster 2 (59.3%, n=18) ---
#>       A     B     C
#> A 0.438 0.132 0.430
#> B 0.526 0.473 0.001
#> C 0.396 0.603 0.001
#> 
#>   component    prior n_assigned mean_posterior     avepp
#> 1         1 0.407233         12      0.9991433 0.9991433
#> 2         2 0.592767         18      0.9873740 0.9873740
# }
```
