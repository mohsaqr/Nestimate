# Build MCML from Raw Transition Data

Builds a Multi-Cluster Multi-Level (MCML) model from raw transition data
(edge lists or sequences) by recoding node labels to cluster labels and
counting actual transitions. Unlike
[`cluster_summary`](https://saqr.me/Nestimate/reference/cluster_summary.md)
which aggregates a pre-computed weight matrix, this function works from
the original transition data to produce the TRUE Markov chain over
cluster states.

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
  labels = NULL
)
```

## Arguments

- x:

  Input data. Accepts multiple formats:

  data.frame with from/to columns

  :   Edge list. Columns named from/source/src/v1/node1/i and
      to/target/tgt/v2/node2/j are auto-detected. Optional weight column
      (weight/w/value/strength).

  data.frame without from/to columns

  :   Sequence data. Each row is a sequence, columns are time steps.
      Consecutive pairs (t, t+1) become transitions.

  tna object

  :   If `x$data` is non-NULL, uses sequence path on the raw data.
      Otherwise falls back to
      [`cluster_summary`](https://saqr.me/Nestimate/reference/cluster_summary.md).

  netobject

  :   If `x$data` is non-NULL, detects edge list vs sequence data.
      Otherwise falls back to
      [`cluster_summary`](https://saqr.me/Nestimate/reference/cluster_summary.md).

  mcml

  :   An `mcml` object is returned unchanged.

  square numeric matrix

  :   Falls back to
      [`cluster_summary`](https://saqr.me/Nestimate/reference/cluster_summary.md).

  non-square or character matrix

  :   Treated as sequence data.

- clusters:

  Cluster/group assignments. Accepts:

  named list

  :   Direct mapping. List names = cluster names, values = character
      vectors of node labels. Example:
      `list(A = c("N1","N2"), B = c("N3","N4"))`

  data.frame

  :   A data frame where the first column contains node names and the
      second column contains group/cluster names. Example:
      `data.frame(node = c("N1","N2","N3"), group = c("A","A","B"))`

  membership vector

  :   Character or numeric vector. Node names are extracted from the
      data. Example: `c("A","A","B","B")`

  column name string

  :   For edge list data.frames, the name of a column containing cluster
      labels. The mapping is built from unique (node, group) pairs in
      both from and to columns. **Limitation:** this mode assigns the
      row's group label to *both* endpoints, so it only makes sense for
      edge lists where source and target nodes always share the same
      group (within-group edges only). For general edge lists where a
      single node may be source in some rows and target in others, or
      where source and target belong to different groups, pass an
      explicit named list (`list(G1 = c("N1","N2"), ...)`) or a
      two-column data frame `data.frame(node, group)` instead.

  NULL

  :   Auto-detect from a `netobject`'s node table (a
      `clusters`/`cluster`/`groups`/`group` column) or from its
      `node_groups`. Only `netobject` input can auto-detect; data frames
      and matrices require an explicit assignment.

- method:

  Aggregation method for combining edge weights: "sum", "mean",
  "median", "max", "min", "density", "geomean". Default "sum". For raw
  sequence/event-log inputs the function is counting observed
  transitions, so `"sum"` is the only interpretation that preserves the
  count semantics – the other methods are useful when aggregating
  weighted edge lists or pre-existing weight matrices, where each row
  already represents a measurement rather than a single observation.

- type:

  Post-processing of the aggregated count matrix. One of:

  "tna"

  :   (default) Row-normalize so each row sums to 1 (first-order Markov
      transition probabilities).

  "raw"

  :   Return the un-normalized count matrix.

  "frequency"

  :   Explicit alias of `"raw"` – identical raw count construction (kept
      as a synonym for callers using frequency-network terminology).

  "cooccurrence"

  :   Symmetrize the matrix (undirected co-occurrence).

  `"semi_markov"` is *not* accepted: the package does not implement a
  semi-Markov / holding-time construction, so passing it errors rather
  than silently aliasing `"tna"`.

- directed:

  Logical. If `TRUE` (default), treat transitions as directed. If
  `FALSE`, symmetrize sequence- and edge-derived weights before
  returning raw/frequency weights or before row-normalizing transition
  probabilities.

- compute_within:

  Logical. Compute within-cluster matrices? Default TRUE.

- actor, action, time, order, session, time_threshold:

  Long-format event-log shortcut. When `action` is supplied on a
  data.frame input, the data is passed through
  [`prepare()`](https://saqr.me/Nestimate/reference/prepare.md) to
  derive a wide sequence, which is then routed to the existing sequence
  path. Behaves identically to
  `prepare(...) |> build_network() |> build_mcml()`.

- exclude:

  Optional character vector of state labels to drop before the network
  is built (e.g. a technical-void marker). On long-format input the
  matching events are removed before sequences are formed, so a dropped
  state never occupies a sequence position; on wide sequence input the
  matching cells are set to `NA`.

- trim:

  Optional truncation of each sequence. `NULL` (default) keeps every
  time point. A fraction in `(0, 1)` keeps the columns covering that
  quantile of sequence lengths (`trim = 0.95` keeps the shortest 95\\
  (`trim = 10` keeps the first 10 time points). Same semantics as
  [`sequence_plot`](https://saqr.me/Nestimate/reference/sequence_plot.md)'s
  `trim`.

- end:

  Terminal state appended after each sequence's last observed state.
  `FALSE` (default) adds none, `TRUE` adds one labelled `"End"`, and a
  string supplies the label. Applied *after* `trim`, so trimming can
  never remove the marker.

- end_by:

  Optional column name(s) grouping the terminal marker at a coarser unit
  than the sequence. `NULL` (default) marks every sequence. When
  supplied, only the *last* sequence of each group is marked – e.g.
  `session = "AttemptID", end = "Conclude", end_by = "SkillID"` closes
  each skill once, not each attempt. Requires long-format input, since
  the grouping column lives there.

- labels:

  Optional name -\> label remap applied to within-cluster nodes (the
  macro layer is left untouched because its labels are cluster names).
  Accepts a 2-column data.frame `(name, label)`, a named character
  vector `c(name = "label")`, or a named list. Unmapped names pass
  through unchanged.

## Value

An `mcml` object with the same layout as the return value of
[`cluster_summary`](https://saqr.me/Nestimate/reference/cluster_summary.md)
(`macro`, `clusters`, `cluster_members`, `edges`, `meta`). On the
sequence and edge-list paths `meta$source` is `"transitions"`,
`meta$type` records the `type` post-processing, and `edges` is a tidy
data frame with one row per observed node-level transition and columns
`from`, `to`, `weight`, `cluster_from`, `cluster_to`, `type`
(`"within"`/`"between"`). Matrix input falls through to
[`cluster_summary`](https://saqr.me/Nestimate/reference/cluster_summary.md),
so `meta$source` is `"matrix"` and `edges` is `NULL`. Works with
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`as_tna`](https://saqr.me/Nestimate/reference/as_tna.md),
[`as_htna`](https://saqr.me/Nestimate/reference/as_htna.md) and
[`macro_network`](https://saqr.me/Nestimate/reference/macro_network.md).

## See also

[`cluster_summary`](https://saqr.me/Nestimate/reference/cluster_summary.md)
for matrix-based aggregation,
[`as_tna`](https://saqr.me/Nestimate/reference/as_tna.md) to promote the
layers to netobjects,
[`macro_network`](https://saqr.me/Nestimate/reference/macro_network.md)
for the cluster-level network with one cluster expanded

## Examples

``` r
# Edge list with clusters
edges <- data.frame(
  from = c("A", "A", "B", "C", "C", "D"),
  to   = c("B", "C", "A", "D", "D", "A"),
  weight = c(1, 2, 1, 3, 1, 2)
)
clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))
build_mcml(edges, clusters)
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
```
