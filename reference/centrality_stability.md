# Centrality Stability Coefficient (CS-coefficient)

`centrality_stability()` estimates how stable a network's centrality
indices are under case-dropping (Epskamp, Borsboom & Fried, 2018). For
each proportion in `drop_prop`, that share of the sequences (or, for
association networks, of the rows) is removed at random `iter` times,
the network is re-estimated on what is left, and its centralities are
correlated with those of the full network. The CS-coefficient of a
measure is the largest proportion that can be dropped while the
correlation stays at or above `threshold` in at least a `certainty`
share of the draws. It takes a `netobject` from
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
of any method, a bare `cograph_network`, or a `netobject_group` or
`mcml`, whose networks are assessed one by one.

The CS-coefficient says how well the nodes' centralities (their values,
or with `method = "spearman"` or `"kendall"` their ordering) survive the
loss of cases; Epskamp et al. (2018) suggest that it should not be below
0.25, and preferably above 0.5, before centrality differences are
interpreted. It does not say whether two nodes differ in centrality, nor
how precise single edges are: check the edges with
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md),
and for `glasso` networks the centrality differences with
[`boot_glasso()`](https://pak.dynasite.org/Nestimate/reference/boot_glasso.md).

## Usage

``` r
centrality_stability(
  x,
  measures = c("InStrength", "OutStrength", "Betweenness"),
  iter = 1000L,
  drop_prop = seq(0.1, 0.9, by = 0.1),
  threshold = 0.7,
  certainty = 0.95,
  method = "pearson",
  centrality_fn = NULL,
  loops = FALSE,
  normalize = FALSE,
  invert = TRUE,
  normalize_diffusion = TRUE,
  seed = NULL
)

# S3 method for class 'net_stability'
print(x, ...)

# S3 method for class 'net_stability_group'
print(x, ...)

# S3 method for class 'net_stability_group'
summary(object, ...)

# S3 method for class 'net_stability'
summary(object, ...)

# S3 method for class 'net_stability'
plot(x, ...)
```

## Arguments

- x:

  The network: a `netobject` from
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
  a bare `cograph_network` (converted to one), or a `netobject_group` or
  `mcml`, each of whose networks is assessed. Any other object, or a
  network built without its data, raises an error. For the
  [`print()`](https://rdrr.io/r/base/print.html) methods: a
  `net_stability` or `net_stability_group`; for the
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method: a
  `net_stability`.

- measures:

  The centrality measures to assess. The default is `"InStrength"`,
  `"OutStrength"` and `"Betweenness"`; `NULL` takes
  [`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)'s
  default (`"InStrength"`, `"Betweenness"`, `"Diffusion"`), and `"all"`
  every built-in measure: `"OutStrength"`, `"InStrength"`,
  `"ClosenessIn"`, `"ClosenessOut"`, `"Closeness"`, `"Betweenness"`,
  `"BetweennessRSP"`, `"Diffusion"` and `"Clustering"`
  ([`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
  defines them). `"InCloseness"` and `"OutCloseness"` are accepted for
  `"ClosenessIn"` and `"ClosenessOut"`. Other names need
  `centrality_fn`; without it they raise an error.

- iter:

  The number of random subsets drawn at each drop proportion, at
  least 2. Default `1000`.

- drop_prop:

  The proportions of cases to drop, each strictly between 0 and 1.
  Default `seq(0.1, 0.9, by = 0.1)`. A proportion that would leave fewer
  than 2 cases is skipped and counts as not reaching `threshold`. The
  CS-coefficient is the largest proportion that meets `certainty`, in
  whatever order they are given.

- threshold:

  The correlation with the full network's centralities that a draw must
  reach, between 0 and 1. Default `0.7`.

- certainty:

  The share of draws at a drop proportion that must reach `threshold`,
  between 0 and 1. Default `0.95`. It is a share of all draws: one with
  an undefined correlation counts as not reaching `threshold`.

- method:

  The correlation between the full and the subset centralities:
  `"pearson"` (default), `"spearman"` or `"kendall"`.

- centrality_fn:

  A function of your own that takes the weight matrix and returns a
  named list of centrality vectors, such as `list(my_metric = ...)`. It
  is called as `centrality_fn(mat)` and used only for the requested
  measures that are not built-in; list their names in `measures`.
  Default `NULL`.

- loops:

  Whether self-loops (staying in a state) count. With `FALSE` (default)
  they are left out of every measure; the stored matrix is not changed.

- normalize:

  Whether to range-normalise every measure to 0–1, as
  `tna::centralities(normalize = TRUE)` does. Default `FALSE`.

- invert:

  Whether a heavier edge is a shorter step in the shortest-path measures
  (distance = 1 / weight), as in `tna`. Default `TRUE`; with `FALSE` the
  weight is the distance. It affects betweenness and closeness only, not
  `"BetweennessRSP"`.

- normalize_diffusion:

  Whether to range-normalise `Diffusion` even when `normalize = FALSE`.
  Default `TRUE`.

- seed:

  An integer that makes the random subsets reproducible. Default `NULL`.
  The session's own random stream is left as it was; a group's networks
  are each started from the same seed.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: not
  used.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) methods: a
  `net_stability` or `net_stability_group`.

## Value

A `net_stability` object. [`print()`](https://rdrr.io/r/base/print.html)
shows the CS-coefficient of each measure;
[`summary()`](https://rdrr.io/r/base/summary.html) returns one row per
measure and drop proportion, with the mean and SD of the correlations
(over the draws where it is defined) and `prop_above`, the share of all
draws that reach `threshold`;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the mean
correlation against the proportion dropped, one line per measure with a
band of one SD and a dashed line at `threshold`, and returns the ggplot
invisibly. For a `netobject_group` or `mcml`, the result is a
`net_stability_group`, a named list of `net_stability` objects, one per
network, whose [`print()`](https://rdrr.io/r/base/print.html) shows the
CS-coefficients as a networks-by-measures table and whose
[`summary()`](https://rdrr.io/r/base/summary.html) stacks the
per-network summaries with a leading `group` column. The
[`print()`](https://rdrr.io/r/base/print.html) methods return their
input invisibly. The `net_stability`'s fields:

- cs:

  The CS-coefficients, a named numeric vector with one per retained
  measure.

- correlations:

  A named list, one per retained measure, of `iter` x
  `length(drop_prop)` matrices of correlations.

- measures:

  The measures assessed (see the rule for measures without variance
  under *Details*).

- drop_prop, threshold, certainty, iter, method:

  The settings used.

## Details

For transition methods, per-sequence count matrices are computed once
and each draw only re-sums them; the network's `scaling` and `threshold`
are applied to every re-estimated network as to the original. For the
other methods the estimator is re-run on each subset with the network's
`params` (after the multilevel decomposition, when the network was built
with a `level` and an id column and its estimator is undirected). The
built-in measures are computed from the weight matrix directly, without
igraph, as in
[`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md):
strength from row and column sums, betweenness and closeness from
Floyd-Warshall shortest paths. A network that rests on one sequence, or
one-hot data with no actor or session, has nothing to drop and raises an
error of class `nestimate_single_sequence`.

Every draw counts towards `certainty`: a draw whose correlation is
undefined (the measure is constant in it), and a draw whose network
could not be estimated, count as not reaching `threshold`. Failed
estimates are reported in a warning of class
`nestimate_resample_failed`.

Measures without variance are handled by two different rules, both
long-standing behaviour. When *some* requested measures do not vary
across the nodes of the original network (or are undefined on it), they
are dropped, and the result covers only the retained ones. When *every*
requested measure is so, a warning of class `nestimate_constant_measure`
is issued and *all* requested names are returned with a CS-coefficient
of 0 and correlations that are all `NA`. A measure counts as constant
when its values span no more than the square root of the machine epsilon
relative to their size, so `"OutStrength"` of a `"relative"` network
with `loops = TRUE`, which is 1 for every node up to rounding error, is
dropped rather than correlated on that error. The same rule makes a
measure constant in a draw an undefined correlation.

## References

Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
psychological networks and their accuracy: A tutorial paper. *Behavior
Research Methods* 50(1), 195-212.
[doi:10.3758/s13428-017-0862-1](https://doi.org/10.3758/s13428-017-0862-1)

## See also

[`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md),
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md),
[`network_reliability()`](https://pak.dynasite.org/Nestimate/reference/network_reliability.md),
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)

## Examples

``` r
seqs <- data.frame(
  T1 = c("plan", "code", "debug", "plan", "test", "code"),
  T2 = c("code", "debug", "code", "plan", "code", "test"),
  T3 = c("debug", "code", "plan", "code", "debug", "plan"),
  T4 = c("test", "plan", "test", "debug", "plan", "code")
)
net <- build_network(seqs, method = "relative")
cs <- centrality_stability(net, iter = 10, drop_prop = 0.3, seed = 1)
# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
  V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
)
net <- build_network(seqs, method = "relative")
cs <- centrality_stability(net, iter = 100, seed = 42,
  measures = c("InStrength", "OutStrength"))
print(cs)
#> Centrality Stability (100 iterations, threshold = 0.7)
#>   Drop proportions: 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9
#> 
#>   CS-coefficients:
#>     InStrength       0.00
#>     OutStrength      0.30
# }
```
