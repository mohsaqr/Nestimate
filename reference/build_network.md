# Build a Network

`build_network()` estimates a network from sequences, an event log,
one-hot codes or numeric variables, and returns a `netobject`: the
object every other Nestimate verb reads (bootstrap, permutation,
centrality, clustering, plotting). One call covers two kinds of network.
Transition networks (`"relative"`, `"frequency"`, `"co_occurrence"` and
their variants) are directed or undirected graphs of states, where an
edge says how often or how likely one state follows another. Association
networks (`"cor"`, `"pcor"`, `"glasso"`, `"ising"`, `"mgm"`) are
undirected graphs of variables, where an edge is a correlation or a
(regularised) partial correlation. The *Estimators* section lists every
method and its aliases; estimators added with
[`register_estimator()`](https://pak.dynasite.org/Nestimate/reference/register_estimator.md)
are called the same way.

Column roles are declared, never read from names or values. An event log
needs `action`, and usually `actor`, `time` or `session`; a wide table
needs position columns named as a numbered run (`T1`, `T2`, ...) or
listed in `state_cols`; one-hot data needs `codes` or
`format = "onehot"`. Input that fits none of these raises an error of
class `nestimate_bad_input` instead of being read one way or another.

The weights are point estimates. Before reading single edges, see how
far they move under resampling with
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md),
and compare groups with
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
rather than by eye. An event log that yields a single sequence is
reported with a message of class `nestimate_single_sequence`: that
network describes one sequence and cannot be resampled.

## Usage

``` r
build_network(
  data,
  method,
  actor = NULL,
  action = NULL,
  time = NULL,
  session = NULL,
  order = NULL,
  codes = NULL,
  group = NULL,
  format = "auto",
  window_size = 3L,
  mode = c("non-overlapping", "overlapping"),
  scaling = NULL,
  threshold = 0,
  level = NULL,
  time_threshold = 900,
  timezone = "UTC",
  predictability = TRUE,
  state_cols = NULL,
  metadata_cols = NULL,
  start = FALSE,
  end = FALSE,
  params = list(),
  labels = NULL,
  ...
)

# S3 method for class 'netobject'
print(x, ...)

# S3 method for class 'netobject_group'
print(x, digits = 3L, ...)

# S3 method for class 'netobject_ml'
print(x, ...)

# S3 method for class 'netobject'
summary(object, ...)

# S3 method for class 'netobject_group'
summary(object, combined = TRUE, ...)

# S3 method for class 'summary.netobject'
print(x, ...)

# S3 method for class 'summary.netobject_group'
print(x, ...)
```

## Arguments

- data:

  The data to estimate from: a data frame of sequences (one row per
  sequence, one column per position), an event log (one row per event,
  with `action` set), one-hot code columns, numeric variables (one row
  per observation) for the association methods, or a square symmetric
  correlation or covariance matrix. A fitted `net_clustering` or
  `net_mmm` object also works: its per-cluster networks are built, with
  the `scaling`, `threshold`, `start`, `end` and `params` given here,
  and returned as a `netobject_group`. The relative networks of a
  `net_mmm` are the fitted chains: `scaling` and `threshold` are applied
  to them, and any other argument raises an error of class
  `nestimate_bad_input`.

- method:

  The estimator: a built-in method or alias (see the *Estimators*
  section), or a name added with
  [`register_estimator()`](https://pak.dynasite.org/Nestimate/reference/register_estimator.md);
  [`list_estimators()`](https://pak.dynasite.org/Nestimate/reference/list_estimators.md)
  lists them. Required, except when `data` is a `net_clustering` or
  `net_mmm` object, whose own method is used.

- actor:

  The column naming who each row belongs to (a person, a team). In an
  event log, each actor's events form a sequence, split further by
  `session` or `time_threshold`. With `level`, it is the id the
  multilevel decomposition groups by. Default `NULL`.

- action:

  The column holding each event's state; setting it declares the data an
  event log. Column names give no role: `action`, `actor`, `time` and
  `session` are read only when set. A column called `Action` with
  nothing declared raises an error of class `nestimate_bad_input`, and
  an undeclared time or session column beside `action` is named in a
  message of class `nestimate_column_role`. Default `NULL`.

- time:

  The timestamp column of an event log. Events are ordered by it, and a
  gap longer than `time_threshold` starts a new sequence. Default
  `NULL`.

- session:

  The column identifying sessions in an event log; each session is its
  own sequence, also without `actor`. Default `NULL`.

- order:

  The column that orders events sharing a timestamp; without it, those
  events keep their row order. Default `NULL`.

- codes:

  The one-hot code columns (one per state, 1 when the state is active);
  setting them declares one-hot data. Default `NULL`.

- group:

  A column to split the data by: one network is estimated per value and
  the result is a `netobject_group`. Default `NULL`.

- format:

  How to read `data`. `"auto"` (default) reads an event log when
  `action` is set, one-hot data when `codes` is set, and wide positions
  otherwise; it never infers one-hot data from the values. `"wide"`,
  `"long"` and `"onehot"` force a reading. The association methods and
  the `wtna` methods read 0/1 columns as they always do.

- window_size:

  For one-hot data: how many consecutive rows form one window (default
  `3`). Windows on sequence data are set through `params` instead.

- mode:

  For one-hot data: `"non-overlapping"` (default) or `"overlapping"`
  windows;
  [`wtna()`](https://pak.dynasite.org/Nestimate/reference/wtna.md)
  describes how overlapping windows are counted. Wide and long sequence
  data ignore it.

- scaling:

  Rescaling applied to the weights after estimation, in the order given:
  `"minmax"` maps the non-zero weights into (0, 1\], `"max"` divides by
  the largest absolute weight, `"rank"` replaces the non-zero weights by
  their ranks, and `"normalize"` makes each row's absolute weights sum
  to 1. Zero weights stay zero under every option. Combine them as in
  `c("rank", "minmax")`. Default `NULL` (none).

- threshold:

  Weights whose absolute value is below this are set to zero. Default
  `0` (none). For `method = "mgm"`, `"LW"` or `"none"` instead selects
  mgm's coefficient-thresholding rule (a `threshold` in `params` takes
  precedence); any other non-numeric value raises an error.

- level:

  Multilevel decomposition for the undirected association methods
  (`cor`, `pcor`, `glasso`, `mgm`): `"between"` estimates the network of
  person means, `"within"` the network of deviations from them, and
  `"both"` returns the two as a `netobject_ml`. It needs an id column,
  given as `actor` or as `params$id`; a directed method raises an error.
  Default `NULL` (no decomposition). For `method = "mgm"`, a numeric
  vector instead gives one category count per modelled variable (`1` for
  continuous ones), after group, actor and metadata columns are set
  aside; it is checked and kept for compatibility, and the fit reads the
  categories from the data.

- time_threshold:

  The longest gap, in seconds, between two events of the same sequence
  in an event log; a longer gap starts a new sequence. `FALSE` turns the
  split off, so each actor (or actor-session) is one sequence. Default
  `900`.

- timezone:

  The Olson time zone for timestamps without an offset (passed to
  [`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md));
  timestamps with an offset, such as `...Z` or `+02:00`, are converted
  from it. Default `"UTC"`.

- predictability:

  Whether to compute each node's predictability (the R-squared of the
  node on its neighbours) for the undirected association methods `cor`,
  `pcor` and `glasso`. It is stored on the network, and
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
  draws it as a ring around each node. Default `TRUE`.

- state_cols:

  The columns of a wide table that are sequence positions; every other
  column of the cleaned input goes to the network's metadata. Declared
  columns are positions even when they hold only 0 and 1; one-hot data
  is declared with `codes` or `format = "onehot"`. With nothing
  declared, the frame must be one column or hold position columns named
  as a numbered run (a shared prefix and a number, `T1`, `T2`, ..., or
  `1`, `2`, ...); numeric columns beside the run and character or factor
  columns outside it are then covariates. Any other frame on a sequence
  method, and columns of only 0 and 1, raise an error of class
  `nestimate_bad_input`: nothing says whether they are positions or
  one-hot codes. Default `NULL`.

- metadata_cols:

  Columns to keep as the network's metadata (covariates), never as
  positions; they cannot overlap `state_cols`. A column left out of a
  numbered run of positions is named in a message of class
  `nestimate_column_role`, and listing it here or in `state_cols`
  settles it. Default `NULL`.

- start:

  A start state put before every sequence, so the network shows where
  sequences begin (an edge from it to each sequence's first state).
  `FALSE` (default) adds none, `TRUE` adds one labelled `"Start"`, and a
  single string adds one with that label. Only for the transition
  methods `relative`, `frequency`, `co_occurrence`, `attention`,
  `ngram`, `gap` and `reverse`; any other method, `wtna` included,
  raises an error.

- end:

  An end state put in the cell after each sequence's last observed
  state, so the network shows where sequences stop. It has no outgoing
  edges and no self-loop, unlike
  [`mark_terminal_state()`](https://pak.dynasite.org/Nestimate/reference/mark_terminal_state.md),
  which fills all trailing missing cells with an absorbing state.
  `FALSE` (default) adds none, `TRUE` adds one labelled `"End"`, and a
  single string adds one with that label. Same methods as `start`.

- params:

  Settings passed to the estimator, as a named list, such as
  `list(gamma = 0.5)` for `glasso` or `list(n_gram = 3)` for `ngram`.
  The list is stored on the network, so
  [`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  and the other resampling verbs re-estimate with exactly the same
  settings. The transition methods take tna's sequence options
  `weighted` and `concat`, and the low-level `begin_state` / `end_state`
  that `start` and `end` set. `window_size` (with `window_span`, or
  `windowed = TRUE` for a size of 2) counts `relative`, `frequency` and
  `co_occurrence` between consecutive non-overlapping blocks of that
  many positions, as
  [`tna::tna()`](https://sonsoles.me/tna/reference/build_model.html)
  does; other transition methods raise an error of class
  `nestimate_bad_input`, and `mode = "overlapping"` a warning of class
  `nestimate_ignored_mode`. Column roles given here (`action`, `id`,
  `id_col`, `actor`, `time`, `session`, `order`, `cols`, `codes`,
  `group`) must name existing columns, and must agree with the same role
  given directly.

- labels:

  New display labels for the nodes, applied after estimation: a
  two-column data frame (`name`, `label`), a named character vector
  `c(name = "label")` or a named list. It rewrites the node labels and
  the weight matrix's dimnames; nodes not listed keep their names.

- ...:

  Further arguments for the estimator. For event logs,
  [`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md)'s
  time options `unix_time_unit`, `is_unix_time` and `custom_format` go
  to
  [`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md)
  instead (e.g. `unix_time_unit = "milliseconds"`). For the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods: ignored.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) methods: a
  `netobject`, `netobject_group` or `netobject_ml`, or the result of
  [`summary()`](https://rdrr.io/r/base/summary.html) on one.

- digits:

  For [`print()`](https://rdrr.io/r/base/print.html) of a
  `netobject_group`: decimal places of the weight range. Default `3`.

- object:

  For [`summary()`](https://rdrr.io/r/base/summary.html): a `netobject`
  or `netobject_group`.

- combined:

  For [`summary()`](https://rdrr.io/r/base/summary.html) of a
  `netobject_group`: `TRUE` (default) joins the groups' metrics into one
  data frame with a column per group; `FALSE` returns a named list of
  the per-group tables.

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
column per group). The [`print()`](https://rdrr.io/r/base/print.html)
methods return their input invisibly. The `netobject`'s fields:

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

  Decomposition level used, or NULL; for `mgm`, the category count of
  each variable.

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
and `glasso`, and `lambda_selected` for `glasso`.

With `level = "both"`, the result is a `netobject_ml` holding the
between-person and within-person networks (`between` and `within`, each
a `netobject`) and the `method`. `level = "between"` or `"within"`
returns a single `netobject` estimated on the decomposed data.

With `group`, or when `data` is a `net_clustering` or `net_mmm` object,
the result is a `netobject_group`: a named list of `netobject`s, one per
group, that records the grouping column in its `group_col` attribute.

## Details

`build_network()` resolves the method alias, checks the declared column
arguments before the format is resolved, and takes the estimator from
the registry. For an association method with `level`, it first
decomposes the data into person means (`"between"`) or deviations from
them (`"within"`). It then calls the estimator with `data` and `params`,
applies `scaling` and then `threshold` to the weight matrix, and
extracts the edges.

In an event log, `action` alone makes all rows one sequence in row
order; `session`, or a `time` gap above `time_threshold`, splits them
without an `actor`. When the log yields a single sequence, a message of
class `nestimate_single_sequence` says so, because such a network cannot
be validated by bootstrap or the other confirmatory tests.

For `method = "mgm"`, missing continuous values are filled with the
observed column mean and missing categorical values with the modal
level. This keeps the rows and is psychnets' documented
`na_method = "pairwise"` policy; it is not pairwise deletion.

For a network with non-negative weights, the
[`summary()`](https://rdrr.io/r/base/summary.html) values equal those of
[`tna::summary.tna()`](https://sonsoles.me/tna/reference/summary.tna.html)
from tna 1.5.4 on.

## Estimators

- `"relative"`:

  Transition probabilities (row-normalised counts); directed. Aliases:
  `"tna"`, `"transition"`.

- `"frequency"`:

  Transition counts; directed. Aliases: `"ftna"`, `"counts"`.

- `"co_occurrence"`:

  Co-occurrence within sequences; undirected. Aliases: `"cna"`,
  `"wcna"`, `"co-occurrence"`.

- `"attention"`:

  Decay-weighted transitions; directed. Alias: `"atna"`.

- `"ngram"`:

  Adjacent pairs counted once per n-gram window; directed. Aliases:
  `"n-gram"`, `"n_gram"`.

- `"gap"`:

  Pairs up to `max_gap + 1` positions apart, weighted `1 / distance`;
  directed.

- `"reverse"`:

  Reply network, the transpose of `"frequency"`; directed.

- `"wtna"`:

  Window-based transitions on one-hot data; directed. Alias:
  `"wtna_transition"`.

- `"wtna_cooccurrence"`:

  Window-based co-occurrence on one-hot data; undirected.

- `"cor"`:

  Pearson correlations; undirected. Aliases: `"corr"`, `"correlation"`.

- `"pcor"`:

  Partial correlations, unregularised; undirected. Alias: `"partial"`.

- `"glasso"`:

  EBIC graphical lasso partial correlations; undirected. Aliases:
  `"ebicglasso"`, `"regularized"`.

- `"ising"`:

  Ising model (L1-penalised logistic regression); undirected. Alias:
  `"isingfit"`.

- `"mgm"`:

  Mixed graphical model; undirected. Aliases: `"mixed"`,
  `"mixed_graphical"`.

`"ngram"`, `"gap"` and `"reverse"` are
[`tna::build_model()`](https://sonsoles.me/tna/reference/build_model.html)'s
types `"n-gram"`, `"gap"` and `"reverse"`. Their settings go in
`params`: `list(n_gram = 2)`, `list(max_gap = 1)` and
`list(weighted = FALSE)` are the defaults. All three return raw weights;
add `scaling = "normalize"` for row probabilities.

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
counts as an edge, an undirected edge is counted once, and density is
capped at 1. Strengths sum absolute weights. Distances take each edge's
absolute weight as its length, not its inverse, and pairs with no path
between them are left out of the mean. Degree centralisation ignores
self-loops, and reciprocity is the share of off-diagonal edges whose
reverse is also present (1 for an undirected network). Mean distance and
reciprocity are `NaN` when no pair qualifies. Printing a summary shows
the counts as integers and the other metrics to four significant digits.

## References

Saqr, M., López-Pernas, S., Törmänen, T., Kaliisa, R., Misiejuk, K., &
Tikka, S. (2025). Transition network analysis: A novel framework for
modeling, visualizing, and identifying the temporal patterns of learners
and learning processes. In *Proceedings of the 15th International
Learning Analytics and Knowledge Conference (LAK25)*. ACM.
[doi:10.1145/3706468.3706513](https://doi.org/10.1145/3706468.3706513)

Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse
covariance estimation with the graphical lasso. *Biostatistics*, 9(3),
432–441.
[doi:10.1093/biostatistics/kxm045](https://doi.org/10.1093/biostatistics/kxm045)

Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria
for Gaussian graphical models. *Advances in Neural Information
Processing Systems*, 23, 604–612.

van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F., Boschloo,
L., Schoevers, R. A., & Waldorp, L. J. (2014). A new method for
constructing networks from binary data. *Scientific Reports*, 4, 5918.
[doi:10.1038/srep05918](https://doi.org/10.1038/srep05918)

Haslbeck, J. M. B., & Waldorp, L. J. (2020). mgm: Estimating
time-varying mixed graphical models in high-dimensional data. *Journal
of Statistical Software*, 93(8), 1–46.
[doi:10.18637/jss.v093.i08](https://doi.org/10.18637/jss.v093.i08)

## See also

[`register_estimator`](https://pak.dynasite.org/Nestimate/reference/register_estimator.md),
[`list_estimators`](https://pak.dynasite.org/Nestimate/reference/list_estimators.md),
[`bootstrap_network`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","C","A"), V2 = c("B","C","A","B"))
net <- build_network(seqs, method = "relative")
net
#> Transition Network (relative probabilities) [directed]
#>   Weights: [1.000, 1.000]  |  mean: 1.000
#> 
#>   Weight matrix:
#>     A B C
#>   A 0 1 0
#>   B 0 0 1
#>   C 1 0 0 
#> 
#>   Initial probabilities:
#>   A             0.500  ████████████████████████████████████████
#>   B             0.250  ████████████████████
#>   C             0.250  ████████████████████
# \donttest{
# Transition network (relative probabilities)
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
  V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
)
net <- build_network(seqs, method = "relative")
print(net)
#> Transition Network (relative probabilities) [directed]
#>   Weights: [0.143, 0.429]  |  mean: 0.250
#> 
#>   Weight matrix:
#>         A     B     C     D
#>   A 0.310 0.310 0.172 0.207
#>   B 0.353 0.176 0.176 0.294
#>   C 0.261 0.174 0.217 0.348
#>   D 0.429 0.238 0.143 0.190 
#> 
#>   Initial probabilities:
#>   C             0.367  ████████████████████████████████████████
#>   A             0.300  █████████████████████████████████
#>   D             0.233  █████████████████████████
#>   B             0.100  ███████████

# Association network (glasso)
freq_data <- convert_sequence_format(seqs, format = "frequency")
net_glasso <- build_network(freq_data, method = "glasso",
                             params = list(gamma = 0.5, nlambda = 50))

# With scaling
net_scaled <- build_network(seqs, method = "relative",
                             scaling = c("rank", "minmax"))
# }
```
