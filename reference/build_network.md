# Estimate a network

`build_network()` estimates a network from sequences, an event log,
one-hot codes or numeric variables, with the estimator named by
`method`, and returns a `netobject`: the object every other Nestimate
verb reads (bootstrap, permutation, centrality, clustering, plotting).
Each method's arguments are on its own page, named under `method`; they
are given here by name, as in
`build_network(df, method = "relative", actor = "id")`.

Column roles are declared, never read from names or values, and input
that fits none of the readings raises an error of class
`nestimate_bad_input`. The weights are point estimates: see how far they
move under resampling with
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md),
and compare groups with
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md).

## Usage

``` r
build_network(
  data,
  method,
  ...,
  group = NULL,
  scaling = NULL,
  threshold = 0,
  labels = NULL,
  inherit = FALSE
)

# S3 method for class 'netobject'
print(x, n = 10L, ...)

# S3 method for class 'netobject_group'
print(x, digits = 3L, ...)

# S3 method for class 'netobject_ml'
print(x, ...)

# S3 method for class 'netobject'
summary(object, ...)

# S3 method for class 'netobject_group'
summary(object, combined = TRUE, format = c("long", "wide"), ...)

# S3 method for class 'summary.netobject'
print(x, ...)

# S3 method for class 'summary.netobject_group'
x[i, j, ...]

# S3 method for class 'summary.netobject_group'
print(x, ...)
```

## Arguments

- data:

  The data the method reads (see its page), or a fitted clustering
  (`net_clustering`, `net_mmm`), to build one network per cluster with
  the clustering's method. Per-cluster `"latency"` networks are not
  available (a clustering keeps the sequences, not the event times):
  name another `method`, or the error has class
  `nestimate_latency_clusters`.

- method:

  The estimator. Each has its own page with its arguments:

  |  |  |  |
  |----|----|----|
  | *Sequences and event logs* |  |  |
  | `"relative"` | transition network analysis (TNA) | see [`build_tna()`](https://pak.dynasite.org/Nestimate/reference/build_tna.md) for details |
  | `"frequency"` | frequency transition network analysis (FTNA) | see [`build_ftna()`](https://pak.dynasite.org/Nestimate/reference/build_ftna.md) for details |
  | `"attention"` | attention transition network analysis (ATNA) | see [`build_atna()`](https://pak.dynasite.org/Nestimate/reference/build_atna.md) for details |
  | `"co_occurrence"` | co-occurrence network analysis (CNA) | see [`build_cna()`](https://pak.dynasite.org/Nestimate/reference/build_cna.md) for details |
  | `"ngram"` | n-gram transition network analysis | see [`build_ngram()`](https://pak.dynasite.org/Nestimate/reference/build_ngram.md) for details |
  | `"gap"` | gap transition network analysis | see [`build_gap()`](https://pak.dynasite.org/Nestimate/reference/build_gap.md) for details |
  | `"reverse"` | reverse (reply) transition network analysis | see [`build_reverse()`](https://pak.dynasite.org/Nestimate/reference/build_reverse.md) for details |
  | *Message logs* |  |  |
  | `"latency"` | latency network analysis | see [`build_latency()`](https://pak.dynasite.org/Nestimate/reference/build_latency.md) for details |
  | *One-hot codes* |  |  |
  | `"wtna"` | window transition network analysis (WTNA) | see [`build_wtna()`](https://pak.dynasite.org/Nestimate/reference/build_wtna.md) for details |
  | `"wtna_cooccurrence"` | window co-occurrence network analysis | see [`wtna()`](https://pak.dynasite.org/Nestimate/reference/wtna.md) for details |
  | *Numeric variables* |  |  |
  | `"cor"` | correlation network | see [`build_cor()`](https://pak.dynasite.org/Nestimate/reference/build_cor.md) for details |
  | `"pcor"` | partial correlation network | see [`build_pcor()`](https://pak.dynasite.org/Nestimate/reference/build_pcor.md) for details |
  | `"glasso"` | EBIC graphical lasso network | see [`build_glasso()`](https://pak.dynasite.org/Nestimate/reference/build_glasso.md) for details |
  | `"ising"` | Ising model | see [`build_ising()`](https://pak.dynasite.org/Nestimate/reference/build_ising.md) for details |
  | `"mgm"` | mixed graphical model (MGM) | see [`build_mgm()`](https://pak.dynasite.org/Nestimate/reference/build_mgm.md) for details |

  The other names
  [`list_estimators()`](https://pak.dynasite.org/Nestimate/reference/list_estimators.md)
  shows (`"tna"`, `"ftna"`, `"atna"`, `"cna"`, ...) also work, and an
  estimator added with
  [`register_estimator()`](https://pak.dynasite.org/Nestimate/reference/register_estimator.md)
  is called by its name. Required, except for a fitted clustering, whose
  own method is used.

- ...:

  The method's arguments, by name, as listed on its page, such as
  `actor`, `action`, `start` or `lambda`. For the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods: ignored.

- group:

  A column to split the data by: one network per value, returned as a
  `netobject_group`. Default `NULL`.

- scaling:

  Rescaling applied to the weights after estimation, in the order given:
  `"minmax"` rescales the weights to run from 0 to 1, over every cell
  for the transition methods (tna's rule, so with no zero cell the
  smallest weight becomes 0) and over the non-zero weights, into (0,
  1\], for the other methods, `"max"` divides by the largest absolute
  weight, `"rank"` replaces the non-zero weights by their ranks, and
  `"normalize"` makes each row's absolute weights sum to 1; an
  undirected network whose matrix it makes asymmetric becomes directed,
  so both of a pair's values are kept as edges. Zero weights stay zero
  under every option. Combine them as in `c("rank", "minmax")`. Default
  `NULL` (none).

- threshold:

  Weights whose absolute value is below this are set to zero. Default
  `0` (none).

- labels:

  New display labels for the nodes, applied after estimation: a
  two-column data frame (`name`, `label`), a named character vector
  `c(name = "label")` or a named list. It rewrites the node labels and
  the weight matrix's dimnames; nodes not listed keep their names.
  Default `NULL`.

- inherit:

  Only when `data` is a fitted clustering (`net_clustering` or
  `net_mmm`) of a network: whether the per-cluster networks also take
  that network's `threshold`, `scaling` and estimator settings. `FALSE`
  (default) builds them with its method and the arguments given here;
  `TRUE` uses its settings for any not given here, when the clusters are
  built with its method. The relative networks of a `net_mmm` are the
  fitted chains, which take only its `threshold` and `scaling`.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) methods: a
  `netobject`, `netobject_group` or `netobject_ml`, or the result of
  [`summary()`](https://rdrr.io/r/base/summary.html) on one.

- n:

  For [`print()`](https://rdrr.io/r/base/print.html) of a `netobject`:
  the number of edges listed, strongest first, self-loops included.
  Default `10`; `Inf` lists them all, as does
  [`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md)
  with `include_self = TRUE`;
  [`extract_transition_matrix()`](https://pak.dynasite.org/Nestimate/reference/extract_transition_matrix.md)
  the weight matrix and
  [`extract_initial_probs()`](https://pak.dynasite.org/Nestimate/reference/extract_initial_probs.md)
  the initial probabilities.

- digits:

  For [`print()`](https://rdrr.io/r/base/print.html) of a
  `netobject_group`: decimal places of the weight range. Default `3`.

- object:

  For [`summary()`](https://rdrr.io/r/base/summary.html): a `netobject`
  or `netobject_group`.

- combined:

  For [`summary()`](https://rdrr.io/r/base/summary.html) of a
  `netobject_group`: `TRUE` (default) joins the groups' metrics into one
  data frame; `FALSE` returns a named list of the per-group tables.

- format:

  For [`summary()`](https://rdrr.io/r/base/summary.html) of a
  `netobject_group` with `combined = TRUE`: `"long"` (default), one row
  per metric and group (`metric`, `group`, `value`), or `"wide"`, one
  row per metric with a column per group. Both print as a
  metric-by-group grid.

- i, j:

  For `[` on a grouped summary: rows and columns, as for a data frame. A
  subset is a plain data frame.

## Value

A `netobject` (class `c("netobject", "cograph_network")`). Print it for
an overview;
[`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md)
gives one row per edge,
[`extract_transition_matrix()`](https://pak.dynasite.org/Nestimate/reference/extract_transition_matrix.md)
and
[`extract_initial_probs()`](https://pak.dynasite.org/Nestimate/reference/extract_initial_probs.md)
the weights and initial-state probabilities, and
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
draws it. [`summary()`](https://rdrr.io/r/base/summary.html) returns the
descriptive metrics, one row per metric (a `netobject_group` gives one
row per metric and group, with a `group` column). The
[`print()`](https://rdrr.io/r/base/print.html) methods return their
input invisibly. The `netobject`'s fields:

- data:

  The state columns of the cleaned input data, as a data frame.

- metadata:

  Data frame of the non-state columns of the cleaned input (and, for
  long input, the per-sequence metadata), or NULL.

- weights:

  The estimated network weight matrix.

- nodes:

  Data frame with columns `id`, `label`, `name`, `x`, `y`.

- edges:

  Data frame of non-zero edges with integer `from`/`to` (node IDs) and
  numeric `weight`.

- directed:

  Logical. Whether the network is directed.

- method:

  The resolved method name.

- params:

  The params list used (for reproducibility).

- scaling:

  The scaling applied (or NULL).

- threshold:

  The threshold applied.

- n_nodes:

  Number of nodes.

- n_edges:

  Number of non-zero edges.

- level:

  Decomposition level used, or NULL; for `mgm` built without a
  decomposition, the category count of each variable.

- build_args:

  The resolved column/format arguments (`actor`, `action`, `time`,
  `session`, `order`, `codes`, `format`, `window_size`, `mode`) used for
  this build.

- meta:

  List with `source`, `layout`, and `tna` metadata (cograph-compatible).

- node_groups:

  Node groupings data frame, or NULL.

- predictability:

  Named numeric vector of R-squared predictability values per node, for
  `cor`, `pcor` and `glasso` when `predictability = TRUE`; NULL
  otherwise.

The estimator's own results are kept beside these fields, depending on
the method: for example `frequency_matrix` and `initial` for `relative`,
`frequency`, `ngram`, `gap` and `reverse` (only `initial` for
`attention`, neither for `co_occurrence`), `cor_matrix` and the sample
size `n` for `cor`, `pcor` and `glasso`, `precision_matrix` for `pcor`
and `glasso`, and `lambda_selected` for `glasso`. A `latency` network
keeps `initial`, `n_replies` and `median_latency` (matrices laid out as
the weights), `latency_cap_seconds` (the resolved global cap, also in
`params`), `latency_unit_caps` (the per-unit caps of a non-global
`latency_scope`, named by unit, also in `params`) and `latency_time`,
the event times in seconds (and, for the `"actor"` and `"session"`
scopes, each sequence's unit), one row per row of `data`.

With `level = "both"`, the result is a `netobject_ml` holding the
between-person and within-person networks (`between` and `within`, each
a `netobject`) and the `method`. `level = "between"` or `"within"`
returns a single `netobject` estimated on the decomposed data.

With `group`, or when `data` is a `net_clustering` or `net_mmm` object,
the result is a `netobject_group`: a named list of `netobject`s, one per
group, that records the grouping column in its `group_col` attribute.

## Methods

Printing a `netobject` shows the kind of network, the range of its
weights (with the numbers of positive and negative edges for `cor`,
`pcor`, `glasso` and `ising`, the mean weight otherwise), the weight
matrix under the node labels, the initial probabilities and the
predictability values when present, the sample size for `cor`, `pcor`,
`glasso` and `ising`, the tuning values of `glasso` and `ising`, and any
scaling or threshold. Printing a `netobject_group` gives a header naming
its source (a clustering from `build_network()` on a `net_clustering` or
`net_mmm` fit, or
[`cluster_network()`](https://pak.dynasite.org/Nestimate/reference/cluster_network.md),
or a plain split by the grouping column) and one row per group with its
node and edge counts and weight range; when the groups are clusters, it
also shows the number and share of sequences in each, in the layout of
the
[`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
and
[`build_mmm()`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
prints. A `netobject_ml` prints the node, edge and sample counts of its
between- and within-person networks.

[`summary()`](https://rdrr.io/r/base/summary.html) gives the descriptive
metrics
[`tna::summary.tna()`](https://sonsoles.me/tna/reference/summary.tna.html)
reports: node and edge counts, density, mean shortest-path distance, the
mean and SD of out- and in-strength, the mean and SD of out-degree, out-
and in-degree centralisation (Freeman) and reciprocity. A self-loop
counts as an edge, and an undirected edge is counted once. Density
divides the edge count, self-loops included, by the n(n - 1)
off-diagonal pairs (n(n - 1) / 2 when undirected), as
[`tna::summary.tna()`](https://sonsoles.me/tna/reference/summary.tna.html)
and `igraph::edge_density(loops = FALSE)` do, and is capped at 1.
Strengths sum absolute weights. Distances take each edge's absolute
weight as its length, not its inverse, and pairs with no path between
them are left out of the mean. Degree centralisation ignores self-loops,
and reciprocity is the share of off-diagonal edges whose reverse is also
present (1 for an undirected network). Mean distance and reciprocity are
`NaN` when no pair qualifies. Printing a summary shows the counts as
integers and the other metrics to four significant digits.

## References

Saqr, M., López-Pernas, S., Törmänen, T., Kaliisa, R., Misiejuk, K., &
Tikka, S. (2025). Transition network analysis: A novel framework for
modeling, visualizing, and identifying the temporal patterns of learners
and learning processes. In *Proceedings of the 15th International
Learning Analytics and Knowledge Conference (LAK25)*. ACM.
[doi:10.1145/3706468.3706513](https://doi.org/10.1145/3706468.3706513)

## See also

[`list_estimators()`](https://pak.dynasite.org/Nestimate/reference/list_estimators.md)
and
[`register_estimator()`](https://pak.dynasite.org/Nestimate/reference/register_estimator.md)
for the estimator registry;
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md),
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
and
[`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
for what to do with the network.

## Examples

``` r
build_network(trajectories, method = "relative")
#> Transition Network (relative probabilities) [directed]
#>   Weights: [0.035, 0.698]  |  mean: 0.333
#> 
#>   Strongest edges (9 of 9):
#>    from       to         weight
#>    Active     Active     0.698 
#>    Average    Average    0.610 
#>    Disengaged Disengaged 0.483 
#>    Disengaged Average    0.397 
#>    Active     Average    0.267 
#>    Average    Active     0.204 
#>    Average    Disengaged 0.186 
#>    Disengaged Active     0.120 
#>    Active     Disengaged 0.035  
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix;
#>   extract_initial_probs() the initial probabilities.
build_network(group_regulation_long, method = "relative", actor = "Actor",
              action = "Action", time = "Time", group = "Course")
#> Group Networks (3 groups, group_col: Course)
#> 
#>   Group  Nodes  Edges  Weights
#>   A      9      76     [0.002, 0.581]
#>   B      9      76     [0.000, 0.479]
#>   C      9      74     [0.001, 0.440]
build_network(trajectories, method = "attention", lambda = 2)
#> Attention Network (decay-weighted transitions) [directed]
#>   Weights: [44.056, 608.921]  |  mean: 263.021
#> 
#>   Strongest edges (9 of 9):
#>    from       to         weight 
#>    Average    Average    608.921
#>    Active     Active     605.510
#>    Active     Average    255.056
#>    Average    Active     227.679
#>    Disengaged Disengaged 200.905
#>    Average    Disengaged 188.805
#>    Disengaged Average    179.156
#>    Disengaged Active      57.105
#>    Active     Disengaged  44.056 
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix;
#>   extract_initial_probs() the initial probabilities.
build_network(srl_strategies, method = "glasso", gamma = 0.25)
#> Partial Correlation Network (EBICglasso) [undirected]
#>   Sample size: 250
#>   Weights: [-0.144, 0.431]  |  +17 / -12 edges
#> 
#>   Strongest edges (10 of 29):
#>    from         to           weight
#>    Time_Mgmt    Effort_Reg   0.431 
#>    Monitoring   Evaluating   0.373 
#>    Elaboration  Organization 0.332 
#>    Planning     Time_Mgmt    0.293 
#>    Evaluating   Organization 0.287 
#>    Planning     Monitoring   0.285 
#>    Evaluating   Effort_Reg   0.222 
#>    Planning     Evaluating   0.217 
#>    Organization Rehearsal    0.216 
#>    Elaboration  Rehearsal    0.198  
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix.
#> 
#>   Predictability (R²):
#>   Evaluating    0.343  ████████████████████████████████████████
#>   Monitoring    0.309  ████████████████████████████████████
#>   Organization  0.283  █████████████████████████████████
#>   Time_Mgmt     0.278  █████████████████████████████████
#>   Planning      0.277  ████████████████████████████████
#>   Effort_Reg    0.247  █████████████████████████████
#>   Elaboration   0.228  ███████████████████████████
#>   Rehearsal     0.157  ██████████████████
#>   Help_Seeking  0.072  ████████
#> 
#>   Gamma: 0.25  |  Lambda: 0.0312
```
