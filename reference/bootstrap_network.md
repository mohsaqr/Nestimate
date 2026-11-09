# Bootstrap a Network Estimate

`bootstrap_network()` resamples the cases a network was estimated from,
re-estimates the network on every resample with the same method and
settings, and reports for each edge its bootstrap mean, standard
deviation, confidence interval and a p-value, together with the network
pruned to the edges that pass. It reads any network from
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
built-in method or registered estimator, as well as a `netobject_group`
or `mcml` (every network is bootstrapped) and a `wtna_mixed` (both of
its components). The cases are the sequences of a transition network
(for one-hot data, its actors or sessions) and the rows of the data
otherwise; with `actor`, whole actors are resampled (see *Nested data
and `actor`*).

The p-value comes from one of two tests. With `inference = "stability"`
it counts the replicates that fall outside a consistency range around
the original weight (by default 75% to 125% of it), so a significant
edge is one whose weight the resamples reproduce within that range. With
`inference = "threshold"` it counts the replicates whose absolute weight
falls below `edge_threshold`, so a significant edge stays above that
threshold. In both, the p-value is that count plus one, divided by the
number of replicates plus one, so it is never 0. Edges with a p-value
below `ci_level` are kept in the pruned network.

The p-value counts replicates, so it cannot fall below `1 / (iter + 1)`:
with too few iterations no edge can pass a small `ci_level`. The
stability test asks whether a weight is reproduced within a range, not
whether it differs from zero, and because the range is a proportion of
the weight, a weak edge has to be reproduced within a narrow band; read
the interval beside the test. The bootstrap describes how the estimate
varies over resamples of one sample; to compare groups, use
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md).

## Usage

``` r
bootstrap_network(
  x,
  iter = 1000L,
  ci_level = 0.05,
  inference = "stability",
  consistency_range = c(0.75, 1.25),
  edge_threshold = NULL,
  seed = NULL,
  boundary = c("inclusive", "strict"),
  ci_method = c("percentile", "basic"),
  actor = NULL
)

# S3 method for class 'net_bootstrap'
print(x, ...)

# S3 method for class 'net_bootstrap'
summary(object, ...)

# S3 method for class 'net_bootstrap_group'
print(x, ...)

# S3 method for class 'net_bootstrap_group'
summary(object, ...)

# S3 method for class 'wtna_boot_mixed'
print(x, ...)

# S3 method for class 'wtna_boot_mixed'
summary(object, ...)
```

## Arguments

- x:

  The network: a `netobject` that carries its data (any network from
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)),
  or a `cograph_network`, which is converted first. A `netobject_group`
  or `mcml` bootstraps every network it holds (for an `mcml`, the macro
  network and each cluster's), and a `wtna_mixed` bootstraps both of its
  components (see *Value*). A network without its data raises an error,
  and one that rests on a single sequence, or one-hot data with no actor
  or session, has nothing to resample and raises an error of class
  `nestimate_single_sequence`. For the
  [`print()`](https://rdrr.io/r/base/print.html) methods: the result of
  `bootstrap_network()`.

- iter:

  The number of bootstrap replicates, at least 2. Default `1000`.

- ci_level:

  The significance level, between 0 and 1: the confidence intervals
  cover `1 - ci_level`, and an edge is significant when its p-value is
  below `ci_level`. Default `0.05`.

- inference:

  The test behind the p-values. `"stability"` (default) counts the
  replicates outside the consistency range around the original weight;
  `"threshold"` counts the replicates whose absolute weight is below
  `edge_threshold`.

- consistency_range:

  For `inference = "stability"`: the two multipliers of the original
  weight that bound the consistency range. Default `c(0.75, 1.25)`. For
  a negative weight the bounds are swapped, so the range still runs from
  the lower to the upper value. Not used by `"threshold"`.

- edge_threshold:

  For `inference = "threshold"`: the absolute weight a replicate must
  reach. Default `NULL`, the 10th percentile of the non-zero absolute
  original weights. Not used by `"stability"`.

- seed:

  A number that makes the resamples reproducible. The session's own
  random stream is left as it was. For a `netobject_group` or `mcml`,
  each network is bootstrapped from the same seed; for a `wtna_mixed`,
  the seed is set once, before both components. Default `NULL`.

- boundary:

  For `inference = "stability"`: how a replicate that lands exactly on a
  bound of the consistency range counts. `"inclusive"` (default, as in
  tna) counts it as outside (\\\le\\ / \\\ge\\); `"strict"` counts only
  replicates strictly outside (\\\<\\ / \\\>\\).

- ci_method:

  The confidence interval for the edge weights. `"percentile"` (default)
  uses the empirical bootstrap quantiles (Efron). `"basic"` reflects
  those quantiles around the observed weight, \\(2\hat{\theta} -
  q\_{1-\alpha/2}, 2\hat{\theta} - q\_{\alpha/2})\\ (Davison & Hinkley
  1997, eq. 5.6), which corrects first-order bootstrap bias but can give
  bounds outside the natural weight range near its limits (for example
  below 0 for a transition probability close to 0).

- actor:

  The column identifying the actor each sequence belongs to (for example
  `"student_id"` for sessions nested in students, `"Group"` for students
  nested in teams), looked up in the network's metadata or its wide
  sequence data. With it, whole actors are resampled; see *Nested data
  and `actor`*. Default `NULL`: sequences are resampled one by one. It
  is available for sequence-based transition networks (`"relative"`,
  `"frequency"`, `"co_occurrence"`); a method outside the fast path (see
  *Details*) raises an error of class `nestimate_actor_unsupported`. A
  value that is not a single column name, or that names fewer than two
  actors, raises an error of class `nestimate_bad_actor`; a column that
  is absent or has missing values, one of class
  `nestimate_actor_missing`; and ids that do not give one actor per
  sequence, one of class `nestimate_actor_misaligned`.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods: ignored.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) methods:
  the result of `bootstrap_network()`.

## Value

A `net_bootstrap` object. [`print()`](https://rdrr.io/r/base/print.html)
lists the five strongest significant edges (by absolute bootstrap mean,
with their interval) and reports the iterations, the number of
significant edges, the interval level and the test, and with `actor` the
ICC and design effect;
[`summary()`](https://rdrr.io/r/base/summary.html) returns one row per
edge of the original network;
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
draws the bootstrapped network;
[`net_prune()`](https://pak.dynasite.org/Nestimate/reference/net_prune.md)
with `method = "bootstrap"` removes the non-significant edges from the
original network. The fields:

- original:

  The original `netobject`.

- mean, sd:

  The bootstrap mean and standard deviation of each weight, as matrices.

- p_values:

  The matrix of p-values.

- significant:

  The original weights where the p-value is below `ci_level`, and 0
  elsewhere, self-loops included.

- ci_lower, ci_upper:

  The matrices of interval bounds.

- cr_lower, cr_upper:

  The matrices of consistency-range bounds (used by `"stability"` only).

- summary:

  A data frame with one row per non-zero original edge (one row per
  unordered pair in an undirected network) and the columns `from`, `to`,
  `weight`, `mean`, `sd`, `p_value`, `sig`, `ci_lower` and `ci_upper`,
  plus `cr_lower` and `cr_upper` with `inference = "stability"`.

- model:

  The pruned `netobject`, whose weights are `significant`: its
  non-significant edges, self-loops included, are zero. (The network
  [`net_prune()`](https://pak.dynasite.org/Nestimate/reference/net_prune.md)
  returns keeps every self-loop.)

- method, params, iter, ci_level, inference, ci_method:

  The bootstrap settings.

- consistency_range, edge_threshold:

  The test settings.

- actor, n_actors:

  The `actor` column and its number of actors; `NULL` without `actor`.

- clustering:

  Only with `actor`: a one-row data frame of `n_sequences`, `n_actors`,
  `icc`, `icc_ci_lower`, `icc_ci_upper` and `deff_edges`.

- clustering_edges:

  Only with `actor`: one row per edge of `summary`, with `from`, `to`,
  `icc`, `sd_actor`, `sd_sequence` and `deff`.

A `netobject_group` or `mcml` gives a `net_bootstrap_group`, a named
list of `net_bootstrap` results: its
[`print()`](https://rdrr.io/r/base/print.html) reports each network's
significant and total edges and the significant edges all networks share
(when the networks have edges in common); its
[`summary()`](https://rdrr.io/r/base/summary.html) is described with the
methods below. A `wtna_mixed` gives a `wtna_boot_mixed` holding the
`transition` and `cooccurrence` results.

The [`print()`](https://rdrr.io/r/base/print.html) methods return their
input invisibly. The [`summary()`](https://rdrr.io/r/base/summary.html)
method of a `net_bootstrap` returns its per-edge table, the `summary`
field described above. For a `net_bootstrap_group` it stacks those
tables into one data frame, with a first column `group` naming the
network each row comes from; for a `wtna_boot_mixed` it returns a list
of the `transition` and `cooccurrence` tables.

## Details

The two tests and the inclusive boundary rule follow
[`tna::bootstrap()`](https://sonsoles.me/tna/reference/bootstrap.html).

Transition networks (`"relative"`, `"frequency"`, `"co_occurrence"`, and
`wtna` networks, which hold the same one-hot counts) take a fast path:
the count matrix of every sequence is computed once, and each replicate
sums the counts of the sequences drawn, then applies the network's
normalisation, scaling and threshold. With `params = list(concat = k)`,
the consecutive groups of up to `k` rows merged by the estimator are the
sequences resampled; joins inside each merged sequence are retained.
With `actor`, a merged sequence must belong entirely to one actor. Every
other method (the association methods, `"attention"`, `"ngram"`,
`"gap"`, `"reverse"` and registered estimators) is re-estimated in full
on resampled rows. A replicate whose estimation fails, or that loses a
node, is left out of the p-values without a warning, and the bootstrap
statistics are computed over the remaining replicates; if every
replicate fails, every p-value is 1.

A network built from an edge list (for example an `mcml` from
[`build_mcml()`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
with `from` and `to`) has no sequences, so each transition is resampled
as an independent case, with a warning that points to the alternatives.

## Nested data and `actor`

The bootstrap resamples sequences as independent units. When sequences
are nested in actors (sessions in students, students in teams), `actor`
names the column identifying the actor, and whole actors are resampled
with replacement, keeping all their sequences together: the cluster
bootstrap that resamples at the top level only (Davison & Hinkley, 1997,
section 3.8; Field & Welsh, 2007). The number of sequences per replicate
then varies with the actors drawn.

With `actor`, the result also reports the nesting effect. The ICC is the
proportion of the total variance that lies between actors (Shrout &
Fleiss, 1979); an ICC close to 0 indicates little evidence of a nesting
effect. It is computed as in
[`permutation`](https://pak.dynasite.org/Nestimate/reference/permutation.md).
The design effect is the ratio of the variance under the nested design
to the variance had the sequences been sampled independently (Kish,
1965): here, the variance of the edge weights over actor-level
replicates divided by their variance over sequence-level replicates
drawn in the same run, reported as the median over edges. `actor` is
available for transition networks (`"relative"`, `"frequency"`,
`"co_occurrence"`).

## References

Davison, A. C., & Hinkley, D. V. (1997). *Bootstrap Methods and Their
Application*. Cambridge University Press.

Field, C. A., & Welsh, A. H. (2007). Bootstrapping clustered data.
*Journal of the Royal Statistical Society: Series B*, 69(3), 369-390.

Kish, L. (1965). *Survey Sampling*. Wiley.

Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in
assessing rater reliability. *Psychological Bulletin*, 86(2), 420-428.

## See also

[`certainty`](https://pak.dynasite.org/Nestimate/reference/certainty.md)
for the closed-form Bayesian counterpart (same result layout, no
resampling);
[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md)

## Examples

``` r
net <- build_network(data.frame(V1 = c("A","B","C"), V2 = c("B","C","A")),
  method = "relative")
boot <- bootstrap_network(net, iter = 10)
# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
  V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
)
net <- build_network(seqs, method = "relative")
boot <- bootstrap_network(net, iter = 100)
print(boot)
#> Bootstrap Network  [Transition Network (relative) | directed]
#>   Iterations : 100  |  Nodes : 4
#>   Edges      : 0 significant / 16 total
#>   CI         : 95%  |  Inference: stability  |  CR [0.75, 1.25]
summary(boot)
#>    from to    weight      mean         sd   p_value   sig   ci_lower  ci_upper
#> 1     A  A 0.1304348 0.1274665 0.08074773 0.6831683 FALSE 0.00000000 0.3022826
#> 2     A  B 0.3043478 0.2997861 0.08278592 0.3663366 FALSE 0.14284091 0.4582168
#> 3     A  C 0.1304348 0.1418228 0.07482973 0.6732673 FALSE 0.00000000 0.3000000
#> 4     A  D 0.4347826 0.4309246 0.11390253 0.3267327 FALSE 0.24720280 0.6570000
#> 5     B  A 0.1724138 0.1737516 0.06474250 0.5346535 FALSE 0.07821429 0.3180000
#> 6     B  B 0.4137931 0.4205944 0.09336440 0.2970297 FALSE 0.23541667 0.5884848
#> 7     B  C 0.2068966 0.2050824 0.07346017 0.4752475 FALSE 0.06981818 0.3420833
#> 8     B  D 0.2068966 0.2005716 0.08216875 0.5049505 FALSE 0.05967023 0.3751667
#> 9     C  A 0.4000000 0.3857941 0.10757329 0.3861386 FALSE 0.17766798 0.6044118
#> 10    C  B 0.2000000 0.1884783 0.09391798 0.6534653 FALSE 0.04761905 0.3828755
#> 11    C  C 0.1500000 0.1571277 0.07338105 0.6039604 FALSE 0.00000000 0.2994885
#> 12    C  D 0.2500000 0.2686000 0.09868840 0.5544554 FALSE 0.10360963 0.4642308
#> 13    D  A 0.1666667 0.1651025 0.08394913 0.6336634 FALSE 0.02261905 0.3436275
#> 14    D  B 0.2222222 0.2362980 0.11192342 0.6930693 FALSE 0.00000000 0.4332589
#> 15    D  C 0.4444444 0.4507200 0.14203236 0.3861386 FALSE 0.18312325 0.7498039
#> 16    D  D 0.1666667 0.1478795 0.11004755 0.8118812 FALSE 0.00000000 0.3750000
#>      cr_lower  cr_upper
#> 1  0.09782609 0.1630435
#> 2  0.22826087 0.3804348
#> 3  0.09782609 0.1630435
#> 4  0.32608696 0.5434783
#> 5  0.12931034 0.2155172
#> 6  0.31034483 0.5172414
#> 7  0.15517241 0.2586207
#> 8  0.15517241 0.2586207
#> 9  0.30000000 0.5000000
#> 10 0.15000000 0.2500000
#> 11 0.11250000 0.1875000
#> 12 0.18750000 0.3125000
#> 13 0.12500000 0.2083333
#> 14 0.16666667 0.2777778
#> 15 0.33333333 0.5555556
#> 16 0.12500000 0.2083333

# Students nested in teams: resample whole teams
teams <- build_network(group_regulation_long, method = "relative",
                       actor = "Actor", action = "Action", time = "Time",
                       group = "Achiever")
bootstrap_network(teams, iter = 100, actor = "Group", seed = 1)
#>   Edge                   High      Low     
#>   --------------------------------------------
#>   synthesis→consensus   0.579     0.388   
#>   cohesion→consensus    0.536     0.451   
#>   adapt→consensus       0.519     0.458   
#>   consensus→plan        0.364     0.432   
#>   discuss→consensus     0.424     0.215   
#>   ... and 33 more shared significant edges
#> 
#> Grouped Bootstrap  [2 groups | 100 iterations | 95% CI]
#>   High                  41 sig / 76 total  |  Group: 100 actors, ICC -0.003, design effect 1.02
#>   Low                   44 sig / 75 total  |  Group: 100 actors, ICC -0.001, design effect 0.98
#>   Shared (all groups)   38 edges
# }
```
