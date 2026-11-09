# Compare two or more networks

`compare_networks()` compares two or more networks cell by cell, node by
node and on global metrics, and returns one `net_network_comparison`
object. It compares every pair of networks, or every network against one
`reference`, and takes any mix of `netobject`s, groups of them, `tna`
models and square weight matrices, as long as they share the same nodes
in the same order. Each table carries a `pair` column, so one call to a
named verb returns it for every pair. The comparison is descriptive by
default; `test` adds permutation, Bayesian and bootstrap evidence to the
same tables.

Every difference is `diff = a - b`, network `a` minus network `b` of a
pair, between two sets of point estimates. The similarity metrics
(Pearson, cosine, Jaccard and the others in the global table) describe
how alike the two weight matrices are; neither they nor a large `diff`
say whether a difference is larger than sampling variation would
produce. For that, rerun with `test = "permutation"` (or `"bayes"` for
relative or frequency networks) and read the `sig` column, and use
`adjust` when many edges are read at once.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws one view
of the comparison per call.

## Usage

``` r
compare_networks(
  ...,
  reference = NULL,
  scaling = c("none", "minmax", "max", "rank", "zscore", "robust", "log", "log1p",
    "softmax", "quantile", "frobenius", "row"),
  measures = c("InStrength", "OutStrength", "Betweenness"),
  labels = NULL,
  test = "none",
  iter = 1000L,
  alpha = 0.05,
  adjust = "none",
  paired = FALSE,
  rope = NULL,
  seed = NULL,
  actor = NULL
)

# S3 method for class 'net_network_comparison'
print(x, digits = 2L, ...)

# S3 method for class 'net_network_comparison'
plot(
  x,
  type = c("networks", "difference", "edges", "nodes", "global", "heatmap", "scatter",
    "inference"),
  pair = NULL,
  combined = TRUE,
  top_n = 20L,
  measure = NULL,
  labels = TRUE,
  digits = 2L,
  what = NULL,
  ...
)
```

## Arguments

- ...:

  The networks, two or more, in any mix of `netobject`,
  `netobject_group` (its members are compared individually and keep
  their names), `cograph_network`, `psychnet`, `mcml`, `tna`,
  `group_tna`, square numeric matrices, or one unnamed list of these.
  Name the arguments to name the networks, as in
  `compare_networks(early = a, late = b)`; an unnamed network is named
  after its method, or `network_<i>`. Nested lists, non-square matrices
  and other classes raise an error of class
  `nestimate_compare_bad_input`. For the
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method:
  passed to
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
  in the `"networks"` and `"difference"` views (such as `layout`,
  `node_size` or `minimum`) and ignored by the other views; the
  [`print()`](https://rdrr.io/r/base/print.html) method ignores it.

- reference:

  The network every other one is compared against, by name or index; it
  is always `network_a` of its pairs, so `diff = reference - other`.
  Default `NULL`, which compares all pairs. A name or index that matches
  no network raises an error of class
  `nestimate_compare_reference_unknown`.

- scaling:

  A rescaling applied to every weight matrix before the comparison.
  `"none"` (default) keeps the weights as estimated; `"minmax"`,
  `"max"`, `"rank"`, `"zscore"`, `"robust"` (median and MAD), `"log"`,
  `"log1p"`, `"softmax"` and `"quantile"` (empirical percentile)
  transform all cells of a matrix together; `"frobenius"` divides by the
  Frobenius norm and `"row"` makes each row's absolute weights sum to 1.
  `"log"` needs strictly positive weights and `"log1p"` weights of at
  least -1, otherwise an error of class
  `nestimate_compare_scaling_domain` is raised. The tests are defined on
  the networks as estimated, so any other value than `"none"` with
  `test` raises an error of class `nestimate_compare_scaling_inference`.

- measures:

  The centrality measures compared in the node table. Default
  `c("InStrength", "OutStrength", "Betweenness")`. Any of
  `"OutStrength"`, `"InStrength"`, `"ClosenessIn"`, `"ClosenessOut"`,
  `"Closeness"`, `"Betweenness"`, `"BetweennessRSP"`, `"Diffusion"` and
  `"Clustering"`, or `"all"` for these nine (see
  [`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)).
  `NULL` or `character(0)` skips the node table. Unknown names are
  dropped with a warning of class `nestimate_compare_unknown_measure`.

- labels:

  Names for the networks, one per network after groups are expanded;
  they replace the argument names. Default `NULL`. A vector of the wrong
  length raises an error of class `nestimate_compare_labels_length`. For
  the [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method:
  whether to print the signed difference on the edge and node views and
  the p-values on the inference view. Default `TRUE`; the heatmap always
  prints its values.

- test:

  The inference to add, one or more of `"none"` (default),
  `"permutation"`, `"bayes"` and `"bootstrap"`. Each fills the columns
  it supports (see *Details*). The permutation and Bayesian tests need
  networks that carry their data and share one estimation method, and
  the Bayesian test needs `"relative"` or `"frequency"` networks;
  otherwise an error of class `nestimate_compare_test_unsupported` is
  raised.

- iter:

  The number of permutations, posterior draws or bootstrap replicates
  per pair, at least 2. Default `1000`.

- alpha:

  The significance level of the permutation and bootstrap tests, a
  number between 0 and 1; the Bayesian credible interval covers
  `1 - alpha`. Default `0.05`.

- adjust:

  The multiplicity correction for the permutation p-values, a method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html), applied
  within each pair and separately to the edge and node tables, as
  [`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
  does. Default `"none"`.

- paired:

  Whether the networks hold paired observations, so the permutation test
  swaps the two members of each pair instead of shuffling all labels;
  the networks then need the same number of sequences. It cannot be
  combined with `actor` (an error of class `nestimate_bad_actor`).
  Default `FALSE`.

- rope:

  The half-width of a region of practical equivalence on the difference
  scale, for the Bayesian test. It adds `bayes_p_rope` and a
  `bayes_decision` of `"different"`, `"equivalent"` or `"undecided"` to
  the edge table. Default `NULL` (none).

- seed:

  An integer that makes the tests reproducible. Pair `i` uses
  `seed + i`, so each pair's result does not depend on the others.
  Default `NULL`.

- actor:

  The column identifying the actor each sequence belongs to (such as
  `"student_id"` when sessions are nested in students, or `"Group"` for
  students nested in teams). It is passed to
  [`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md),
  which then shuffles whole actors; see its sections *Nested data and
  actor* and *ICC and design effect*. The ICC and design effects are
  added to the global table under the category `"Nesting"`. It needs
  `test = "permutation"`, otherwise an error of class
  `nestimate_compare_actor_needs_permutation` is raised. Default `NULL`.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: a
  `net_network_comparison` object.

- digits:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: the
  decimals shown in printed values. Default `2`.

- type:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: the view to draw, one per call. `"networks"` (default) draws
  each network once with
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html);
  `"difference"` draws the signed difference network of each pair;
  `"edges"` is a ranked dumbbell of the largest edge differences, and
  `"nodes"` the same for the centralities; `"global"` shows the 22
  comparison metrics as points, plus the inference rows with their
  p-values (and bootstrap intervals) when `test` was run; `"heatmap"`
  the signed difference matrix; `"scatter"` the weights of one network
  against the other's; and `"inference"` a forest of the edge
  differences on the difference scale, with the credible intervals when
  the Bayesian test ran and the p-value printed per edge. `"inference"`
  needs `test`, otherwise an error of class `nestimate_compare_no_test`
  is raised.

- pair:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: the comparisons to draw. Default `NULL`, all of them. Give the
  pair names as printed (`"A vs B"`), the two network names in either
  order (`c("A", "B")`), one network name (`"A"`, every pair it takes
  part in), or indices into the pair table (`1`, `c(1, 3)`). For
  `type = "networks"`, it selects the networks taking part in the chosen
  pairs. A selection that matches no comparison raises an error of class
  `nestimate_compare_unknown_pair`.

- combined:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: whether a view of several pairs is one figure. With `TRUE`
  (default), the `ggplot` views facet by pair and the `"networks"` and
  `"difference"` views draw all panels on one base-graphics page. With
  `FALSE`, the `ggplot` views return a named list of single-pair plots,
  for iterating over or saving, and the base-graphics views draw one
  panel per page.

- top_n:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: the number of edges shown per pair in the `"edges"` and
  `"inference"` views, largest absolute difference first. Default `20`.

- measure:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: the centrality measures shown in the `"nodes"` view. Default
  `NULL`, all measures of the comparison; a measure the object does not
  carry raises an error of class `nestimate_compare_unknown_measure`.

- what:

  Deprecated alias for `type`, kept so existing calls keep working.

## Value

A `net_network_comparison` object.
[`print()`](https://rdrr.io/r/base/print.html) shows one line per pair
(Pearson correlation, mean and largest absolute difference, the largest
change) and returns `x` invisibly;
[`summary()`](https://rdrr.io/r/base/summary.html) returns one row per
pair;
[`edge_differences()`](https://pak.dynasite.org/Nestimate/reference/comparison_tables.md),
[`node_differences()`](https://pak.dynasite.org/Nestimate/reference/comparison_tables.md),
[`global_differences()`](https://pak.dynasite.org/Nestimate/reference/comparison_tables.md)
and
[`network_metrics()`](https://pak.dynasite.org/Nestimate/reference/comparison_tables.md)
return the full tables (see
[comparison_tables](https://pak.dynasite.org/Nestimate/reference/comparison_tables.md));
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws one view
per call. The object's fields:

- `networks`: the input networks as `netobject`s, with scaled weights;

- `matrices`: the scaled weight matrices;

- `pairs`: one row per comparison, with `pair`, `network_a` and
  `network_b`;

- `edges`: one row per pair and cell: `pair`, `network_a`, `network_b`,
  `from`, `to`, `weight_a`, `weight_b`, `diff`, `abs_diff`, `rel_diff`,
  `ratio`, `log_ratio`, `rank_a`, `rank_b`, `rank_diff`,
  `percentile_diff`, `status` (`both`, `only_a`, `only_b` or `neither`),
  `higher` (the name of the higher network, or `equal`), and the
  inference columns;

- `nodes`: one row per pair, node and measure: `pair`, `network_a`,
  `network_b`, `node`, `measure`, `value_a`, `value_b`, `diff`,
  `abs_diff`, `rank_a`, `rank_b`, `higher`, and the inference columns;
  `NULL` without `measures`;

- `global`: one row per pair and metric (22 descriptive metrics in five
  categories, and the inference rows): `pair`, `network_a`, `network_b`,
  `category`, `metric`, `key`, `value`, and the inference columns;

- `network_metrics`: one row per network and structural metric:
  `network`, `metric`, `value`;

- `differences`: one `netdifference` object per pair (from
  [`subtract_networks()`](https://pak.dynasite.org/Nestimate/reference/subtract_networks.md));

- the settings `scaling`, `reference`, `measures`, `test`, `iter`,
  `alpha`, `adjust`, `paired`, `rope` and `actor`, with `directed` (one
  logical per network), `n_networks` and `n_pairs`.

For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method:
a `ggplot` for the `"edges"`, `"nodes"`, `"global"`, `"heatmap"`,
`"scatter"` and `"inference"` views, or, with `combined = FALSE`, a
named list of them, one per pair, returned invisibly. The `"networks"`
and `"difference"` views draw in base graphics (with
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
when cograph is installed, otherwise with a built-in circular drawer
that ignores `...`) and return `NULL` invisibly.

## Details

**Cells.** Every cell of the weight matrix is a row of the edge table,
the diagonal and the edges absent from one network (weight 0) included;
when both networks are undirected only the cells with `from <= to` are
kept. A difference smaller than `sqrt(.Machine$double.eps)` counts as no
difference (`higher = "equal"`).

**Guarding.** Ratios are never `Inf` or `NaN`, and no pseudo-counts are
added: `ratio` is `NA` when `weight_b` is 0, `rel_diff`
(`abs(diff) / (abs(a) + abs(b))`) is `NA` when both weights are 0, and
`log_ratio = log1p(a) - log1p(b)` is `NA` when either weight is
negative.

**Inference.** `"permutation"` runs
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
per pair and adds `perm_effect`, `perm_p` and `perm_sig` to the edge and
node tables, and two rows to the global table with permutation p-values:
`M`, the sum of the absolute differences over all cells, and `S`, the
largest absolute difference. With `actor`, whole actors are shuffled and
the global table also gains the rows `"ICC"`, `"Design effect (edges)"`
and `"Design effect (M)"`. `"bayes"` runs
[`bayes_compare()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
per pair (prior `0.5`, `iter` draws) and adds to the edge table
`bayes_diff` (posterior mean difference), `bayes_ci_lower`,
`bayes_ci_upper`, `bayes_pd` (probability of direction), `bayes_p`
(`2 * (1 - pd)`) and `bayes_sig`; with `rope`, also `bayes_p_rope` (a
normal approximation from the posterior mean and SD) and
`bayes_decision`. `"bootstrap"` runs the vertex bootstrap of
[`vertex_compare()`](https://pak.dynasite.org/Nestimate/reference/vertex_compare.md)
per pair and adds structural rows (density, mean weight, centralization
and, for directed networks, reciprocity) to the global table, with
`boot_observed_a`, `boot_observed_b`, `boot_se`, `boot_ci_lower`,
`boot_ci_upper`, `boot_z`, `boot_p` and `boot_sig`. It reads only the
weights, so plain matrices support it; the other two tests need
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
output with data, or `tna` models, which are rebuilt from their
sequences when these reproduce the weights. The `sig` and `evidence`
columns take the permutation result when it was run (edge and node
tables), and otherwise the Bayesian one (edge table only).

## Errors

The errors that `compare_networks()` and its verbs raise themselves have
the class `nestimate_compare_error` and one of these:

- `nestimate_compare_too_few`: fewer than two networks;

- `nestimate_compare_bad_input`: an input of an unsupported class, a
  non-square matrix or a nested list;

- `nestimate_compare_dim_mismatch`: different numbers of nodes;

- `nestimate_compare_node_mismatch`: node labels that differ in set or
  order;

- `nestimate_compare_na_weights`: missing weights;

- `nestimate_compare_reference_unknown`,
  `nestimate_compare_labels_length`, `nestimate_compare_scaling_domain`,
  `nestimate_compare_scaling_inference`,
  `nestimate_compare_test_unsupported` and
  `nestimate_compare_actor_needs_permutation`: see the arguments;

- `nestimate_compare_unknown_pair`: a `pair` that matches no comparison;

- `nestimate_compare_no_nodes`: the node view or table of a comparison
  without `measures`;

- `nestimate_compare_no_test`: `plot(type = "inference")` without
  `test`;

- `nestimate_compare_unknown_measure`: selecting a measure the object
  does not carry.

Plain argument-type checks raise unclassed errors, and errors from
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md),
[`bayes_compare()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
and
[`vertex_compare()`](https://pak.dynasite.org/Nestimate/reference/vertex_compare.md)
keep their own classes (for example `nestimate_bad_actor` for
`paired = TRUE` with `actor`).

An unknown name in `measures` gives a warning of class
`nestimate_compare_warning` and `nestimate_compare_unknown_measure`.

## Reading the figures

One colour contract holds in every view and every backend: `"#4A6FE3"`
marks `network_a` (the reference, when one is set) as the higher of the
two, `"#D33F6A"` marks `network_b`, and grey marks no difference; the
plotted quantity is always `diff = a - b`. Colour never carries the sign
alone: a solid line and a circular marker repeat "`a` higher", a dashed
line and a square marker repeat "`b` higher", and the printed value
carries its sign. When `test` was run, evidence is shown by opacity and
a starred (edge and node views) or annotated (inference view) label: a
non-significant difference is faded, never deleted.

## See also

[`compare_model()`](https://pak.dynasite.org/Nestimate/reference/compare_model.md)
(two-network predecessor),
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md),
[`bayes_compare()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md),
[`vertex_compare()`](https://pak.dynasite.org/Nestimate/reference/vertex_compare.md),
[`subtract_networks()`](https://pak.dynasite.org/Nestimate/reference/subtract_networks.md).

## Examples

``` r
# Regulation networks for the three courses, compared pairwise.
courses <- build_network(group_regulation_long, method = "relative",
                         actor = "Actor", action = "Action", time = "Time",
                         group = "Course")
cmp <- compare_networks(courses)
cmp
#> Network comparison (descriptive): 3 networks, 3 pairs, scaling = none
#>   networks: A (9 nodes, directed), B (9 nodes, directed), C (9 nodes, directed)
#> 
#>  pair   pearson mean|diff| max|diff| largest change                   
#>  A vs B 0.94    0.03       0.18      discuss -> consensus (A higher)  
#>  A vs C 0.91    0.04       0.22      synthesis -> consensus (A higher)
#>  B vs C 0.99    0.01       0.07      synthesis -> discuss (C higher)  
#> 
#> Tables: summary(x), edge_differences(x), node_differences(x), global_differences(x), network_metrics(x). Plot: plot(x, type = ...)
summary(cmp)
#>    pair network_a network_b n_cells n_differing share_higher_a share_higher_b
#>  A vs B         A         B      81          78           0.47           0.49
#>  A vs C         A         C      81          78           0.43           0.53
#>  B vs C         B         C      81          77           0.48           0.47
#>  mean_abs_diff max_abs_diff pearson spearman cosine jaccard
#>           0.03         0.18    0.94     0.93   0.97    0.78
#>           0.04         0.22    0.91     0.90   0.95    0.73
#>           0.01         0.07    0.99     0.98   0.99    0.88
#>                top_edge top_edge_higher
#>    discuss -> consensus               A
#>  synthesis -> consensus               A
#>    synthesis -> discuss               C
# `pair` takes the two network names, in either order.
edge_differences(cmp, pair = c("A", "B"))
#>    pair network_a network_b       from         to weight_a weight_b  diff
#>  A vs B         A         B      adapt      adapt        0        0     0
#>  A vs B         A         B      adapt   cohesion     0.26     0.28 -0.03
#>  A vs B         A         B      adapt  consensus     0.52     0.48  0.04
#>  A vs B         A         B      adapt coregulate        0     0.03 -0.03
#>  A vs B         A         B      adapt    discuss     0.03     0.05 -0.02
#>  A vs B         A         B      adapt    emotion     0.15     0.12  0.03
#>  A vs B         A         B      adapt    monitor     0.03     0.02  0.01
#>  A vs B         A         B      adapt       plan     0.02     0.01     0
#>  A vs B         A         B      adapt  synthesis        0        0     0
#>  A vs B         A         B   cohesion      adapt        0        0     0
#>  A vs B         A         B   cohesion   cohesion     0.05     0.01  0.03
#>  A vs B         A         B   cohesion  consensus     0.53     0.48  0.06
#>  A vs B         A         B   cohesion coregulate     0.08     0.16 -0.07
#>  A vs B         A         B   cohesion    discuss     0.04     0.08 -0.04
#>  A vs B         A         B   cohesion    emotion     0.12     0.09  0.03
#>  A vs B         A         B   cohesion    monitor     0.02     0.04 -0.02
#>  A vs B         A         B   cohesion       plan     0.15     0.14     0
#>  A vs B         A         B   cohesion  synthesis     0.01        0  0.01
#>  A vs B         A         B  consensus      adapt        0     0.01     0
#>  A vs B         A         B  consensus   cohesion     0.02     0.01  0.01
#>  A vs B         A         B  consensus  consensus     0.08     0.09 -0.01
#>  A vs B         A         B  consensus coregulate     0.17     0.21 -0.04
#>  A vs B         A         B  consensus    discuss     0.23     0.15  0.09
#>  A vs B         A         B  consensus    emotion     0.08     0.06  0.02
#>  A vs B         A         B  consensus    monitor     0.04     0.06 -0.02
#>  A vs B         A         B  consensus       plan     0.37     0.41 -0.05
#>  A vs B         A         B  consensus  synthesis     0.01     0.01     0
#>  A vs B         A         B coregulate      adapt     0.02     0.01  0.01
#>  A vs B         A         B coregulate   cohesion     0.04     0.03  0.01
#>  A vs B         A         B coregulate  consensus     0.11     0.15 -0.04
#>  A vs B         A         B coregulate coregulate     0.01     0.03 -0.02
#>  A vs B         A         B coregulate    discuss     0.24     0.28 -0.04
#>  A vs B         A         B coregulate    emotion     0.21     0.16  0.04
#>  A vs B         A         B coregulate    monitor     0.09     0.08  0.01
#>  A vs B         A         B coregulate       plan     0.26     0.24  0.03
#>  A vs B         A         B coregulate  synthesis     0.02     0.02     0
#>  A vs B         A         B    discuss      adapt     0.02     0.11 -0.09
#>  A vs B         A         B    discuss   cohesion     0.06     0.04  0.03
#>  A vs B         A         B    discuss  consensus     0.42     0.24  0.18
#>  A vs B         A         B    discuss coregulate     0.07     0.09 -0.02
#>  A vs B         A         B    discuss    discuss     0.17     0.21 -0.05
#>  A vs B         A         B    discuss    emotion     0.11     0.10  0.01
#>  A vs B         A         B    discuss    monitor     0.02     0.03 -0.01
#>  A vs B         A         B    discuss       plan     0.01     0.01     0
#>  A vs B         A         B    discuss  synthesis     0.11     0.17 -0.06
#>  A vs B         A         B    emotion      adapt        0        0     0
#>  A vs B         A         B    emotion   cohesion     0.33     0.33     0
#>  A vs B         A         B    emotion  consensus     0.34     0.29  0.05
#>  A vs B         A         B    emotion coregulate     0.02     0.04 -0.02
#>  A vs B         A         B    emotion    discuss     0.12     0.10  0.02
#>  A vs B         A         B    emotion    emotion     0.06     0.08 -0.02
#>  A vs B         A         B    emotion    monitor     0.03     0.04 -0.01
#>  A vs B         A         B    emotion       plan     0.09     0.11 -0.02
#>  A vs B         A         B    emotion  synthesis     0.01        0     0
#>  A vs B         A         B    monitor      adapt     0.01     0.01     0
#>  A vs B         A         B    monitor   cohesion     0.05     0.05     0
#>  A vs B         A         B    monitor  consensus     0.16     0.16     0
#>  A vs B         A         B    monitor coregulate     0.05     0.06 -0.01
#>  A vs B         A         B    monitor    discuss     0.37     0.38     0
#>  A vs B         A         B    monitor    emotion     0.09     0.09     0
#>  A vs B         A         B    monitor    monitor     0.02     0.02     0
#>  A vs B         A         B    monitor       plan     0.22     0.23     0
#>  A vs B         A         B    monitor  synthesis     0.02     0.01  0.01
#>  A vs B         A         B       plan      adapt        0        0     0
#>  A vs B         A         B       plan   cohesion     0.03     0.02  0.01
#>  A vs B         A         B       plan  consensus     0.29     0.28  0.02
#>  A vs B         A         B       plan coregulate     0.02     0.01  0.01
#>  A vs B         A         B       plan    discuss     0.06     0.07 -0.01
#>  A vs B         A         B       plan    emotion     0.18     0.13  0.05
#>  A vs B         A         B       plan    monitor     0.07     0.07     0
#>  A vs B         A         B       plan       plan     0.33     0.41 -0.08
#>  A vs B         A         B       plan  synthesis        0        0     0
#>  A vs B         A         B  synthesis      adapt     0.15     0.30 -0.15
#>  A vs B         A         B  synthesis   cohesion     0.03     0.04 -0.01
#>  A vs B         A         B  synthesis  consensus     0.58     0.41  0.17
#>  A vs B         A         B  synthesis coregulate     0.01     0.07 -0.05
#>  A vs B         A         B  synthesis    discuss     0.02     0.06 -0.04
#>  A vs B         A         B  synthesis    emotion     0.07     0.07 -0.01
#>  A vs B         A         B  synthesis    monitor        0     0.02 -0.02
#>  A vs B         A         B  synthesis       plan     0.14     0.04  0.10
#>  A vs B         A         B  synthesis  synthesis        0        0     0
#>  abs_diff rel_diff ratio log_ratio rank_a rank_b rank_diff percentile_diff
#>         0       NA    NA         0      3      3         0               0
#>      0.03     0.05  0.91     -0.02     70     72        -2           -0.02
#>      0.04     0.04  1.08      0.02     79     81        -2           -0.02
#>      0.03        1     0     -0.03      3     27       -24           -0.27
#>      0.02     0.27  0.58     -0.02  32.50     36     -3.50           -0.04
#>      0.03     0.12  1.27      0.03     60     56         4            0.05
#>      0.01     0.20  1.50      0.01  32.50     23      9.50            0.12
#>         0     0.11  1.25         0     18     17         1            0.01
#>         0       NA    NA         0      3      3         0               0
#>         0     0.43  2.54         0     10      9         1            0.01
#>      0.03     0.61  4.12      0.03     38     14        24            0.30
#>      0.06     0.06  1.12      0.04     80     80         0               0
#>      0.07     0.31  0.53     -0.07     49     61       -12           -0.15
#>      0.04     0.34  0.50     -0.04     37     45        -8           -0.10
#>      0.03     0.15  1.35      0.03     57     49         8            0.10
#>      0.02     0.38  0.45     -0.02  20.50     32    -11.50           -0.14
#>         0     0.01  1.02         0     59     58         1            0.01
#>      0.01        1    NA      0.01     12      3         9            0.09
#>         0     0.15  0.74         0      9     10        -1           -0.01
#>      0.01     0.35  2.08      0.01     24     13        11            0.14
#>      0.01     0.05  0.90     -0.01     47     48        -1           -0.01
#>      0.04     0.10  0.81     -0.03     64     65        -1           -0.01
#>      0.09     0.22  1.58      0.07     68     59         9            0.11
#>      0.02     0.16  1.37      0.02     48     38        10            0.12
#>      0.02     0.25  0.60     -0.02     35     39        -4           -0.05
#>      0.05     0.06  0.89     -0.03     76     79        -3           -0.04
#>         0     0.06  0.88         0     13     12         1            0.01
#>      0.01     0.28  1.76      0.01     25     16         9            0.11
#>      0.01     0.12  1.27      0.01     36     26        10            0.12
#>      0.04     0.15  0.74     -0.03     54     60        -6           -0.07
#>      0.02     0.51  0.32     -0.02     14     28       -14           -0.17
#>      0.04     0.07  0.86     -0.03     69     70        -1           -0.01
#>      0.04     0.12  1.26      0.04     66     63         3            0.04
#>      0.01     0.08  1.17      0.01     52     46         6            0.07
#>      0.03     0.06  1.12      0.02     71     68         3            0.04
#>         0     0.03  0.93         0  20.50     22     -1.50           -0.01
#>      0.09     0.68  0.19     -0.09     26     55       -29           -0.36
#>      0.03     0.28  1.78      0.03     43     29        14            0.17
#>      0.18     0.28  1.77      0.14     78     69         9            0.11
#>      0.02     0.13  0.77     -0.02     45     51        -6           -0.07
#>      0.05     0.12  0.78     -0.04     63     66        -3           -0.04
#>      0.01     0.07  1.15      0.01     55     53         2            0.02
#>      0.01     0.22  0.64     -0.01     19     25        -6           -0.07
#>         0     0.01  1.02         0     16     15         1            0.01
#>      0.06     0.20  0.66     -0.05     53     64       -11           -0.14
#>         0        1    NA         0      8      3         5            0.04
#>         0     0.01  0.99         0     73     75        -2           -0.02
#>      0.05     0.08  1.17      0.04     75     73         2            0.02
#>      0.02     0.32  0.51     -0.02     27     34        -7           -0.09
#>      0.02     0.10  1.23      0.02     56     52         4            0.05
#>      0.02     0.16  0.73     -0.02     42     47        -5           -0.06
#>      0.01     0.17  0.71     -0.01     31     33        -2           -0.02
#>      0.02     0.09  0.84     -0.02     51     54        -3           -0.04
#>         0     0.64  4.63         0     11      8         3            0.04
#>         0     0.24  1.65         0     17     11         6            0.07
#>         0     0.01  0.98         0     39     35         4            0.05
#>         0     0.01  1.02         0     62     62         0               0
#>      0.01     0.05  0.91         0     40     37         3            0.04
#>         0        0  0.99         0     77     76         1            0.01
#>         0     0.01  0.98         0     50     50         0               0
#>         0     0.07  1.15         0  22.50     21      1.50            0.02
#>         0     0.01  0.98         0     67     67         0               0
#>      0.01     0.19  1.48      0.01  22.50     18      4.50            0.06
#>         0     0.55  3.47         0      6      6         0               0
#>      0.01     0.15  1.36      0.01     34     24        10            0.12
#>      0.02     0.03  1.06      0.01     72     71         1            0.01
#>      0.01     0.22  1.56      0.01     29     19        10            0.12
#>      0.01     0.08  0.86     -0.01     41     42        -1           -0.01
#>      0.05     0.18  1.43      0.05     65     57         8            0.10
#>         0     0.01  0.99         0     46     44         2            0.02
#>      0.08     0.11  0.81     -0.06     74     78        -4           -0.05
#>         0     0.59  3.90         0      7      7         0               0
#>      0.15     0.33  0.50     -0.12     61     74       -13           -0.16
#>      0.01     0.15  0.73     -0.01     30  30.50     -0.50           -0.01
#>      0.17     0.18  1.43      0.12     81     77         4            0.05
#>      0.05     0.70  0.18     -0.05     15     41       -26           -0.32
#>      0.04     0.45  0.38     -0.04     28     40       -12           -0.15
#>      0.01     0.06  0.89     -0.01     44     43         1            0.01
#>      0.02        1     0     -0.02      3     20       -17           -0.19
#>      0.10     0.58  3.77      0.09     58  30.50     27.50            0.33
#>         0       NA    NA         0      3      3         0               0
#>   status higher
#>  neither  equal
#>     both      B
#>     both      A
#>   only_b      B
#>     both      B
#>     both      A
#>     both      A
#>     both      A
#>  neither  equal
#>     both      A
#>     both      A
#>     both      A
#>     both      B
#>     both      B
#>     both      A
#>     both      B
#>     both      A
#>   only_a      A
#>     both      B
#>     both      A
#>     both      B
#>     both      B
#>     both      A
#>     both      A
#>     both      B
#>     both      B
#>     both      B
#>     both      A
#>     both      A
#>     both      B
#>     both      B
#>     both      B
#>     both      A
#>     both      A
#>     both      A
#>     both      B
#>     both      B
#>     both      A
#>     both      A
#>     both      B
#>     both      B
#>     both      A
#>     both      B
#>     both      A
#>     both      B
#>   only_a      A
#>     both      B
#>     both      A
#>     both      B
#>     both      A
#>     both      B
#>     both      B
#>     both      B
#>     both      A
#>     both      A
#>     both      B
#>     both      A
#>     both      B
#>     both      B
#>     both      B
#>     both      A
#>     both      B
#>     both      A
#>     both      A
#>     both      A
#>     both      A
#>     both      A
#>     both      B
#>     both      A
#>     both      B
#>     both      B
#>     both      A
#>     both      B
#>     both      B
#>     both      A
#>     both      B
#>     both      B
#>     both      B
#>   only_b      B
#>     both      A
#>  neither  equal
global_differences(cmp, pair = c("A", "B"))
#>    pair network_a network_b             category               metric
#>  A vs B         A         B    Weight Deviations      Mean Abs. Diff.
#>  A vs B         A         B    Weight Deviations    Median Abs. Diff.
#>  A vs B         A         B    Weight Deviations            RMS Diff.
#>  A vs B         A         B    Weight Deviations       Max Abs. Diff.
#>  A vs B         A         B    Weight Deviations Rel. Mean Abs. Diff.
#>  A vs B         A         B    Weight Deviations             CV Ratio
#>  A vs B         A         B         Correlations              Pearson
#>  A vs B         A         B         Correlations             Spearman
#>  A vs B         A         B         Correlations              Kendall
#>  A vs B         A         B         Correlations             Distance
#>  A vs B         A         B      Dissimilarities            Euclidean
#>  A vs B         A         B      Dissimilarities            Manhattan
#>  A vs B         A         B      Dissimilarities             Canberra
#>  A vs B         A         B      Dissimilarities          Bray-Curtis
#>  A vs B         A         B      Dissimilarities            Frobenius
#>  A vs B         A         B         Similarities               Cosine
#>  A vs B         A         B         Similarities              Jaccard
#>  A vs B         A         B         Similarities                 Dice
#>  A vs B         A         B         Similarities              Overlap
#>  A vs B         A         B         Similarities                   RV
#>  A vs B         A         B Pattern Similarities       Rank Agreement
#>  A vs B         A         B Pattern Similarities       Sign Agreement
#>              key value
#>    mean_abs_diff  0.03
#>  median_abs_diff  0.02
#>         rms_diff  0.05
#>     max_abs_diff  0.18
#>     rel_mean_abs  0.25
#>         cv_ratio  1.09
#>          pearson  0.94
#>         spearman  0.93
#>          kendall  0.79
#>     distance_cor  0.94
#>        euclidean  0.41
#>        manhattan  2.25
#>         canberra 18.40
#>      bray_curtis  0.13
#>        frobenius  0.19
#>           cosine  0.97
#>          jaccard  0.78
#>             dice  0.87
#>          overlap  0.87
#>               rv  0.93
#>   rank_agreement  0.90
#>   sign_agreement  0.95
plot(cmp)                                    # each network once

plot(cmp, type = "difference")               # signed difference per pair

plot(cmp, type = "edges", pair = c("A", "B"))

# \donttest{
# High against low achievers, with a permutation test on every edge.
achievers <- build_network(group_regulation_long, method = "relative",
                           actor = "Actor", action = "Action",
                           time = "Time", group = "Achiever")
cmp_perm <- compare_networks(achievers, test = "permutation",
                             iter = 100, seed = 1)
summary(cmp_perm)
#>         pair network_a network_b n_cells n_differing share_higher_a
#>  High vs Low      High       Low      81          78           0.51
#>  share_higher_b mean_abs_diff max_abs_diff pearson spearman cosine jaccard
#>            0.46          0.03         0.21    0.92     0.92   0.95    0.75
#>              top_edge top_edge_higher n_sig_edges n_sig_nodes m_stat   m_p
#>  discuss -> consensus            High          43          16   2.61 0.010
#>  s_stat   s_p
#>    0.21 0.010
plot(cmp_perm, type = "inference", top_n = 10)

# }
```
