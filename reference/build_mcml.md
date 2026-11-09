# Build MCML from Raw Transition Data

`build_mcml()` builds a Multi-Cluster Multi-Level (MCML) model from
transitions: it recodes every state to its cluster and counts the
transitions that occur, which gives the Markov chain over clusters (the
macro network) and one network per cluster of the transitions among its
own states. It reads sequences (one row per sequence), event logs (with
`action`), edge lists (with `from` and `to`), `tna` models and
`netobject`s that carry their data, and it re-partitions an existing
`mcml` (see `combine` and `expand`). A weight matrix has no transitions
to count, so it is passed to
[`cluster_summary()`](https://pak.dynasite.org/Nestimate/reference/cluster_summary.md),
which aggregates the weights instead.

With the defaults (`type = "tna"`, `method = "sum"`), a macro edge is
the probability that the next state lies in cluster j given that the
current one lies in cluster i; its diagonal is the probability of
staying in a cluster, whether in the same state or a different one. Both
levels depend on the partition, which is the analyst's: compare a few
partitions with `combine` and `expand` rather than reading one. The
weights are point estimates;
[`as_tna()`](https://pak.dynasite.org/Nestimate/reference/as_tna.md)
turns the levels into networks that carry their sequences, so
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
can show how far each edge moves under resampling.

## Usage

``` r
build_mcml(
  x,
  clusters = NULL,
  method = c("sum", "mean", "median", "max", "min", "density", "geomean"),
  type = c("tna", "frequency", "cooccurrence", "raw"),
  directed = TRUE,
  compute_within = TRUE,
  actor = NULL,
  action = NULL,
  time = NULL,
  order = NULL,
  session = NULL,
  time_threshold = 900,
  exclude = NULL,
  trim = NULL,
  end = FALSE,
  end_by = NULL,
  labels = NULL,
  combine = NULL,
  expand = NULL,
  from = NULL,
  to = NULL,
  weight = NULL,
  cols = NULL
)

# S3 method for class 'mcml_layer'
print(x, ...)

# S3 method for class 'mcml'
print(x, ...)

# S3 method for class 'mcml'
summary(object, ...)
```

## Arguments

- x:

  The data. Accepted forms:

  data frame with `from` and `to` set

  :   An edge list: one row per transition, with an optional `weight`
      column.

  data frame with `action` set

  :   An event log, turned into sequences by
      [`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md).

  data frame of positions

  :   Sequences, one row per sequence and one column per time step:
      either a numbered run of columns (`T1`, `T2`, ...) and nothing
      else, or the columns named in `cols`. Each pair of adjacent
      columns (t, t + 1) whose two states are both observed is a
      transition; a missing cell breaks the sequence.

  `tna` model

  :   Its stored sequences are decoded and counted; a model without data
      is passed to
      [`cluster_summary()`](https://pak.dynasite.org/Nestimate/reference/cluster_summary.md).

  `netobject` or bare `cograph_network`

  :   A network of a counting method (`relative`, `frequency`,
      `co_occurrence` and their aliases) is re-counted from the data it
      carries. A network of another method (such as an association or
      attention network), a one-hot network, or one without data is
      passed to
      [`cluster_summary()`](https://pak.dynasite.org/Nestimate/reference/cluster_summary.md),
      which aggregates its weights.

  `mcml`

  :   Returned unchanged, unless `clusters`, `combine` or `expand`
      changes its partition: it is then re-estimated from the sequences
      it carries (see `combine`).

  square numeric matrix

  :   Passed to
      [`cluster_summary()`](https://pak.dynasite.org/Nestimate/reference/cluster_summary.md).

  non-square or character matrix

  :   Read as sequences.

  A data frame that declares no role and is not a numbered run of
  positions, or whose columns hold only 0 and 1, raises an error of
  class `nestimate_bad_input`. For the
  [`print()`](https://rdrr.io/r/base/print.html) methods: an `mcml`
  object or one of its layers (an `mcml_layer`).

- clusters:

  The cluster of each state. Accepted forms:

  named list

  :   Cluster names mapped to their states, such as
      `list(A = c("N1", "N2"), B = c("N3", "N4"))`.

  data frame

  :   State names in the first column, cluster names in the second, such
      as
      `data.frame(node = c("N1", "N2", "N3"), group = c("A", "A", "B"))`.

  membership vector

  :   One cluster per state, such as `c("A", "A", "B", "B")`, matched to
      the states of the data in alphabetical order (a named vector is
      matched by name).

  column name

  :   For an edge list only: a column holding each row's cluster, given
      to both its endpoints. It suits edge lists whose edges never cross
      clusters; a node that gets two clusters raises an error, so pass a
      list or a data frame for any other edge list.

  `NULL`

  :   (default) Read from a `netobject`'s node table (a column
      `clusters`, `cluster`, `groups` or `group`) or from its
      `node_groups` (aligned to the node table by label). The values are
      then taken in the order of the node table, which must list the
      states in the order
      [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
      gives them (alphabetical for sequence data). Other input needs an
      explicit assignment and raises an error without one.

  Every state must belong to exactly one cluster, the `end` state
  included. On the sequence and edge-list routes the clusters of the
  result come in alphabetical order.

- method:

  How the transitions between (or within) two clusters are combined into
  one weight, with
  [`net_aggregate_weights()`](https://pak.dynasite.org/Nestimate/reference/net_aggregate_weights.md):
  `"sum"` (default), `"mean"`, `"median"`, `"max"`, `"min"`, `"density"`
  or `"geomean"`. On sequences and event logs every transition weighs 1,
  so only `"sum"` keeps the counts: `"mean"`, `"median"`, `"max"`,
  `"min"` and `"geomean"` give 1 for every pair of clusters with a
  transition, and `"density"` divides each count by `n_i * n_j`. The
  other methods suit weighted edge lists and weight matrices, where each
  row is already a measurement rather than one observation. On sequences
  and event logs the within-cluster networks likewise get 1 for each
  observed pair of states under every method but `"sum"` (`"density"`
  too, which is not divided by `n_i * n_j` there); on an edge list they
  hold the aggregated weights.

- type:

  How the aggregated counts are turned into weights:

  `"tna"`

  :   (default) Each row is normalised to sum to 1, giving first-order
      Markov transition probabilities.

  `"raw"`

  :   The counts as they are.

  `"frequency"`

  :   The same counts as `"raw"`, under the name frequency networks use.

  `"cooccurrence"`

  :   The mean of the count matrix and its transpose (undirected
      co-occurrence); the result is recorded as undirected.

  `"semi_markov"` raises an error: the package has no semi-Markov
  (holding-time) construction. For input passed to
  [`cluster_summary()`](https://pak.dynasite.org/Nestimate/reference/cluster_summary.md)
  (a weight matrix, a `tna` model or `netobject` without data, or a
  `netobject` of a method that is not a count, such as `"attention"` or
  `"glasso"`) `type` has no effect; setting it gives a warning for a
  matrix, or for a `tna` model or `netobject` without data.

- directed:

  Whether transitions are directed. Default `TRUE`. `FALSE` replaces the
  counts by the mean of the count matrix and its transpose before `type`
  is applied.

- compute_within:

  Whether to build the within-cluster networks. Default `TRUE`. An
  `mcml` built with `FALSE` keeps no per-cluster sequences, so it cannot
  be re-partitioned (see `combine`).

- actor, action, time, order, session, time_threshold:

  For an event log: the columns and the session gap (in seconds, default
  `900`) passed to
  [`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md),
  which forms the sequences that are then counted, as
  `prepare(...) |> build_network() |> build_mcml()` would. Setting
  `action` declares the data an event log.

- exclude:

  States to drop before anything is counted, such as a technical marker.
  In an event log the matching events are removed before sequences are
  formed, so a dropped state never takes a position; in wide sequences
  the matching cells become missing. Default `NULL`. Only for wide
  sequences and event logs: an edge list, a `netobject` or a `tna` model
  raises an error. The same holds for `trim` and `end`. On an existing
  `mcml` they raise an error together with a re-partition (`clusters`,
  `combine` or `expand`); without one the `mcml` is returned unchanged
  and they have no effect (see `combine`).

- trim:

  How much of each sequence to keep. `NULL` (default) keeps every time
  step. A fraction in (0, 1) keeps the columns up to that quantile of
  the sequence lengths (the last observed position of each row), so
  `trim = 0.95` cuts only the longest 5% of sequences; a value of 1 or
  more keeps that many leading columns (`trim = 10` keeps the first 10).
  It changes the data, so it only ever cuts the end of a sequence; the
  fraction in
  [`sequence_plot()`](https://pak.dynasite.org/Nestimate/reference/sequence_plot.md)'s
  `trim` works differently because it changes only the drawing.

- end:

  A terminal state put in the cell after each sequence's last observed
  state. `FALSE` (default) adds none, `TRUE` adds one labelled `"End"`,
  and a single string adds one with that label; assign it to a cluster
  in `clusters`. It is added after `trim`, so trimming never removes it.

- end_by:

  Columns grouping the terminal marker at a coarser unit than the
  sequence: only the sequence holding each group's last event is marked.
  For example,
  `session = "AttemptID", end = "Conclude", end_by = "SkillID"` closes
  each skill once, not each attempt. It needs `end` and an event log,
  and raises an error otherwise. Default `NULL` (every sequence is
  marked).

- labels:

  New display labels for the states of the within-cluster networks,
  applied after estimation: a two-column data frame (`name`, `label`), a
  named character vector `c(name = "label")` or a named list. States not
  listed keep their names. The new labels also replace the state names
  in the cluster membership and the edge table; the macro network keeps
  its cluster names.

- combine, expand:

  Change the partition. `combine` merges clusters into one: a character
  vector merges one group, and a list merges several, labelled by its
  names (default label `"A + B"`); the merged cluster lists its states
  in cluster-name order. `expand` then splits the named clusters (or
  every cluster, with `"all"` or `TRUE`) into one cluster per state,
  named by the state. On new input they apply to `clusters` (in any of
  its forms, auto-detection included) before estimation, so
  `build_mcml(data, clusters = cl, combine = c("A", "B"))` equals a
  build with `A` and `B` merged in `cl`; this works for every input,
  matrices included. On an existing `mcml` they re-estimate the model
  (macro network, within-cluster networks and stored sequences) from the
  sequences it carries, with its original `type`, `method` and
  `directed` unless these are passed; a new `clusters` does the same
  under that partition, and an unchanged partition returns an equal
  object. Any session split, `exclude`, `trim` or `end` applied at build
  time is already in those sequences, so passing one of them (or
  `actor`, `action`, `time`, `order`, `session`, `time_threshold`,
  `end_by`, `labels`, `from`, `to`, `weight`, `cols`) with a
  re-partition raises an error, as does passing `clusters` together with
  `combine` or `expand`. An `mcml` built from a matrix, from an edge
  list (it keeps only the within-cluster edges) or with
  `compute_within = FALSE` raises an error of class
  `nestimate_mcml_no_sequences`.
  [`sequence_plot()`](https://pak.dynasite.org/Nestimate/reference/sequence_plot.md)
  takes the same two arguments as display options: its `combine` draws
  what it draws for `build_mcml(x, combine = )`, while its `expand`
  opens clusters in the Summary panel only and keeps one panel per
  cluster.

- from, to, weight:

  The columns of an edge list: its source, its target and an optional
  non-negative weight (each transition counts 1 without it). Names give
  no role: a frame is read as an edge list only when `from` and `to` are
  set, and setting them with `action` or `cols` raises an error of class
  `nestimate_bad_input`.

- cols:

  For a wide data frame: the columns holding the positions, needed when
  the frame is anything but a numbered run (`T1`, `T2`, ...). Default
  `NULL`.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods: not used;
  any extra argument raises an error.

- object:

  For [`summary()`](https://rdrr.io/r/base/summary.html): an `mcml`
  object.

## Value

An `mcml` object with the layout described under
[`cluster_summary()`](https://pak.dynasite.org/Nestimate/reference/cluster_summary.md).
Print it for the type, the method, the node and cluster counts, the
clusters, the number of transitions between and within clusters and the
macro weights; printing one layer (an `mcml_layer`) shows its node and
edge counts, weight range, weight matrix, initial probabilities and the
size of the data it carries.
[`summary()`](https://rdrr.io/r/base/summary.html) returns a data frame
with one row per cluster, with columns `cluster`, `size`, `within_total`
(the total absolute weight of the within-cluster network, `NA` when it
was not built), `between_out` and `between_in` (the absolute macro
weight leaving and entering the cluster, diagonal excluded); for an
undirected macro network, `between_out` is the total incident weight and
`between_in` is `NA`.
[`as_tna()`](https://pak.dynasite.org/Nestimate/reference/as_tna.md)
turns the levels into a `netobject_group`,
[`as_htna()`](https://pak.dynasite.org/Nestimate/reference/as_htna.md)
rebuilds the node-level network with its clusters,
[`macro_network()`](https://pak.dynasite.org/Nestimate/reference/macro_network.md)
opens chosen clusters into their states, and
[`sequence_plot()`](https://pak.dynasite.org/Nestimate/reference/sequence_plot.md)
draws the sequences by cluster. The
[`print()`](https://rdrr.io/r/base/print.html) methods return their
input invisibly.

On the sequence and edge-list routes, the object differs from the matrix
route's in these fields:

- edges:

  A data frame with one row per observed node-level transition and
  columns `from`, `to`, `weight`, `cluster_from`, `cluster_to` and
  `type` (`"within"` or `"between"`).

- macro:

  Its `inits` are the shares of the clusters among the first-column
  states (missing cells left out), or the normalised column sums for an
  edge list; its `data` holds the sequences recoded to cluster labels,
  or one `from`/`to` row per transition for an edge list.

- clusters:

  Each layer counts only the transitions between states of its own
  cluster. Its `inits` are the shares of its states among the
  first-column states that belong to it (normalised column sums for an
  edge list), and its `data` the sequences with the states of other
  clusters set missing (for an edge list, its within-cluster
  transitions).

- meta:

  `source` is `"transitions"`, `type` records the `type` used, and
  `directed` whether the stored weights are directed (`FALSE` for
  `type = "cooccurrence"`).

Input passed to
[`cluster_summary()`](https://pak.dynasite.org/Nestimate/reference/cluster_summary.md)
gives its result: `source` `"matrix"` in `meta`, and no edge table.

## See also

[`cluster_summary`](https://pak.dynasite.org/Nestimate/reference/cluster_summary.md)
for matrix-based aggregation,
[`as_tna`](https://pak.dynasite.org/Nestimate/reference/as_tna.md) to
promote the layers to netobjects,
[`macro_network`](https://pak.dynasite.org/Nestimate/reference/macro_network.md)
for the cluster-level network with one cluster expanded

## Examples

``` r
# Edge list with clusters: its columns are declared
edges <- data.frame(
  source = c("A", "A", "B", "C", "C", "D"),
  target = c("B", "C", "A", "D", "D", "A"),
  count  = c(1, 2, 1, 3, 1, 2)
)
clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))
build_mcml(edges, clusters, from = "source", to = "target", weight = "count")
#> MCML Network
#> ============
#> Type: tna  | Method: sum 
#> Nodes: 4  | Clusters: 2 
#> Transitions: 6 
#>   Macro: 2  | Per-cluster: 4 
#> 
#> Clusters:
#>   G1 (2): A, B
#>   G2 (2): C, D
#> 
#> Macro (cluster-level) weights:
#>        G1     G2
#> G1 0.5000 0.5000
#> G2 0.3333 0.6667

# Sequence data with clusters
seqs <- data.frame(
  T1 = c("A", "C", "B"),
  T2 = c("B", "D", "A"),
  T3 = c("C", "C", "D"),
  T4 = c("D", "A", "C")
)
cs <- build_mcml(seqs, clusters, type = "raw")
cs
#> MCML Network
#> ============
#> Type: raw  | Method: sum 
#> Nodes: 4  | Clusters: 2 
#> Transitions: 9 
#>   Macro: 3  | Per-cluster: 6 
#> 
#> Clusters:
#>   G1 (2): A, B
#>   G2 (2): C, D
#> 
#> Macro (cluster-level) weights:
#>    G1 G2
#> G1  2  2
#> G2  1  4
summary(cs)
#>   cluster size within_total between_out between_in
#> 1      G1    2            2           2          1
#> 2      G2    2            4           1          2

# Change the partition while building ...
three <- build_mcml(seqs, list(G1 = "A", G2 = "B", G3 = c("C", "D")))
build_mcml(seqs, list(G1 = "A", G2 = "B", G3 = c("C", "D")),
           combine = c("G1", "G2"))
#> MCML Network
#> ============
#> Type: tna  | Method: sum 
#> Nodes: 4  | Clusters: 2 
#> Transitions: 9 
#>   Macro: 3  | Per-cluster: 6 
#> 
#> Clusters:
#>   G1 + G2 (2): A, B
#>   G3 (2): C, D
#> 
#> Macro (cluster-level) weights:
#>         G1 + G2  G3
#> G1 + G2     0.5 0.5
#> G3          0.2 0.8

# ... or re-partition an existing mcml: the model is re-estimated
build_mcml(three, combine = c("G1", "G2"))           # G1 + G2 as one cluster
#> MCML Network
#> ============
#> Type: tna  | Method: sum 
#> Nodes: 4  | Clusters: 2 
#> Transitions: 9 
#>   Macro: 3  | Per-cluster: 6 
#> 
#> Clusters:
#>   G1 + G2 (2): A, B
#>   G3 (2): C, D
#> 
#> Macro (cluster-level) weights:
#>         G1 + G2  G3
#> G1 + G2     0.5 0.5
#> G3          0.2 0.8
build_mcml(three, combine = list(AB = c("G1", "G2"))) # named merge
#> MCML Network
#> ============
#> Type: tna  | Method: sum 
#> Nodes: 4  | Clusters: 2 
#> Transitions: 9 
#>   Macro: 3  | Per-cluster: 6 
#> 
#> Clusters:
#>   AB (2): A, B
#>   G3 (2): C, D
#> 
#> Macro (cluster-level) weights:
#>     AB  G3
#> AB 0.5 0.5
#> G3 0.2 0.8
build_mcml(three, expand = "G3")                      # C and D as clusters
#> MCML Network
#> ============
#> Type: tna  | Method: sum 
#> Nodes: 4  | Clusters: 4 
#> Transitions: 9 
#>   Macro: 9  | Per-cluster: 0 
#> 
#> Clusters:
#>   C (1): C
#>   D (1): D
#>   G1 (1): A
#>   G2 (1): B
#> 
#> Macro (cluster-level) weights:
#>      C      D     G1  G2
#> C  0.0 0.6667 0.3333 0.0
#> D  1.0 0.0000 0.0000 0.0
#> G1 0.0 0.5000 0.0000 0.5
#> G2 0.5 0.0000 0.5000 0.0
```
