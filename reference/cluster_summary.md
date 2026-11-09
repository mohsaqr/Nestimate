# Cluster Summary Statistics

`cluster_summary()` aggregates a node-level weight matrix to clusters of
nodes and returns an `mcml` object with two levels: a macro network
whose entry (i, j) combines every edge from a node of cluster i to a
node of cluster j (the diagonal combines the edges inside a cluster),
and one within-cluster network per cluster, the submatrix of its nodes.
It takes a weight matrix, a `netobject`, a bare `cograph_network` or a
`tna` model, and the clusters as a membership vector, a node-to-cluster
data frame or a named list. This is the matrix route of Multi-Cluster
Multi-Level (MCML) analysis;
[`build_mcml()`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
calls it for matrix input.

The macro weights are arithmetic aggregates (see `method`) of weights
that already exist. They are not row-normalised and are not the
transition probabilities of a process over clusters: the transitions
between clusters can be counted only from the sequences, which
[`build_mcml()`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
does (its default `type = "tna"` gives transition probabilities between
clusters).
[`as_tna()`](https://pak.dynasite.org/Nestimate/reference/as_tna.md)
keeps the macro weights as they are and does not normalise them.

## Usage

``` r
cluster_summary(
  x,
  clusters = NULL,
  method = c("sum", "mean", "median", "max", "min", "density", "geomean"),
  directed = TRUE,
  compute_within = TRUE
)
```

## Arguments

- x:

  The node-level network: a numeric square matrix whose row and column
  names (identical, in the same order) label the nodes, with finite
  weights such as transition counts, co-occurrence frequencies or
  probabilities; a `netobject` from
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  or a bare `cograph_network`, whose weight matrix is used; or a `tna`
  model, whose weights are used. An `mcml` object is returned unchanged.
  A matrix without names labels its nodes `"1"`, `"2"`, ....

- clusters:

  The cluster of each node, required (`NULL`, the default, raises an
  error;
  [`build_mcml()`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
  is the verb that reads clusters from a `netobject`'s node table).
  Accepted forms:

  membership vector

  :   One value per node, in the order of the matrix rows, such as
      `c(1, 1, 2, 2, 3, 3)` or `c("A", "A", "B")`. A named vector is
      matched to the nodes by name instead. Numeric values give clusters
      in sorted order, character values in order of first appearance.

  data frame

  :   Node names in the first column, cluster names in the second, such
      as
      `data.frame(node = c("A", "B", "C"), group = c("G1", "G1", "G2"))`;
      the clusters come in alphabetical order.

  named list

  :   Cluster names mapped to their nodes, such as
      `list(Alpha = c("A", "B"), Beta = c("C", "D"))`; the clusters keep
      the list's order.

  Every node must belong to exactly one cluster: an unknown node, a node
  listed twice, a node left out or an empty cluster raises an error.

- method:

  How the edges between (or within) two clusters are combined into one
  weight, with
  [`net_aggregate_weights()`](https://pak.dynasite.org/Nestimate/reference/net_aggregate_weights.md);
  zero edges are left out first.

  `"sum"`

  :   (default) The total weight, which keeps the total flow. Suits
      counts such as transition frequencies.

  `"mean"`

  :   The mean weight, so a cluster of 5 nodes does not get 5 times the
      weight of a cluster of 1. Suits a matrix of probabilities (rows
      summing to 1), where a sum would grow with cluster size.

  `"median"`

  :   The median weight, less affected by a few extreme weights.

  `"max"`

  :   The strongest connection.

  `"min"`

  :   The weakest connection.

  `"density"`

  :   The total weight divided by the number of possible edges between
      the two clusters (`n_i * n_j`). It equals `"mean"` when the block
      has no zero edges and is smaller than it when it has some.

  `"geomean"`

  :   The geometric mean of the positive weights, for multiplicative
      processes.

- directed:

  Whether the network is directed, so that A-\>B and B-\>A are separate
  edges. Default `TRUE`. `FALSE` replaces the matrix by the mean of it
  and its transpose before aggregating, which makes every layer
  symmetric.

- compute_within:

  Whether to build the within-cluster networks, one n_i x n_i matrix per
  cluster. Default `TRUE`; `FALSE` skips them when only the macro
  network is needed, and
  [`summary()`](https://rdrr.io/r/base/summary.html) then reports their
  totals as `NA`.

## Value

An `mcml` object. Print it for the aggregation method, the node and
cluster counts, the clusters and the macro weights;
[`summary()`](https://rdrr.io/r/base/summary.html) returns a data frame
with one row per cluster with its size, total absolute within-cluster
weight and between-cluster weight sent and received (described under
[`build_mcml()`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md));
[`as_tna()`](https://pak.dynasite.org/Nestimate/reference/as_tna.md)
turns the layers into a `netobject_group`;
[`simplicial_features()`](https://pak.dynasite.org/Nestimate/reference/simplicial_features.md)
summarises the topology of each within-cluster network.
[`macro_network()`](https://pak.dynasite.org/Nestimate/reference/macro_network.md),
which opens chosen clusters into their states, needs an `mcml` built
from sequences by
[`build_mcml()`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md);
on this result it raises an error of class `nestimate_no_expand_source`.
The object's fields:

- macro:

  The cluster-level layer (class `mcml_layer`): `weights`, the k x k
  matrix of aggregated weights, diagonal included and rows not
  normalised; `inits`, the column sums of that matrix normalised to sum
  to 1 (equal shares when every weight is zero); `labels`, the cluster
  names; and `data`, `NULL` on this route.

- clusters:

  A named list of one `mcml_layer` per cluster: the n_i x n_i submatrix
  of its nodes, self-loops included, and its normalised column sums as
  `inits` (1 for a single-node cluster). `NULL` when
  `compute_within = FALSE`.

- cluster_members:

  A named list of each cluster's node labels, such as
  `list(A = c("n1", "n2"), B = c("n3", "n4", "n5"))`.

- edges:

  `NULL`: a matrix carries no node-level transitions. The sequence and
  edge-list routes of
  [`build_mcml()`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
  fill it.

- meta:

  A list with `type` (always `"aggregate"` here), `method`, `directed`,
  `n_nodes`, `n_clusters`, `cluster_sizes` (a named integer vector) and
  `source` (`"matrix"`).

## Details

### Workflow

A typical MCML analysis:


    # 1. Create the node-level network
    net <- build_network(data, method = "relative")

    # 2. Aggregate its edges to cluster level (arithmetic aggregation)
    cs <- cluster_summary(net, clusters = group_assignments, method = "sum")

    # 3. Read the result
    print(cs)      # macro weights
    summary(cs)    # one row per cluster

    # 4. Promote the layers to netobjects for downstream verbs
    nets <- as_tna(cs)

### The macro matrix

The macro matrix has clusters as both rows and columns. An off-diagonal
entry (i, j) is the aggregated weight from cluster i to cluster j; a
diagonal entry (i, i) aggregates the edges inside cluster i.

Rows are not normalised: each entry is the aggregate `method` produces.
Row-normalising an arithmetic aggregate here (the old `type = "tna"`
combined with `method = "min"`, `"mean"` and so on) gave rows that sum
to 1 but are not a probability distribution over any process, which is
why `type` was removed from the matrix route. The sequence and edge-list
routes of
[`build_mcml()`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
keep `type`: there the aggregate is always a count, and `type` chooses
between well-defined network constructions.

### Choosing method

|  |  |  |
|----|----|----|
| **Input data** | **Suggested method** | **Reason** |
| Edge counts | `"sum"` | Keeps the total flow between clusters |
| Transition matrix | `"mean"` | Avoids cluster-size bias |
| Correlation matrix | `"mean"` | Average correlations |
| Dense weighted | `"max"` / `"median"` | Less sensitive to single weights |

## See also

[`build_mcml`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
to build an mcml from raw transitions instead of a weight matrix,
[`as_tna`](https://pak.dynasite.org/Nestimate/reference/as_tna.md) to
promote the layers to netobjects,
[`macro_network`](https://pak.dynasite.org/Nestimate/reference/macro_network.md)
for the cluster-level network with one cluster expanded back into its
member states

## Examples

``` r
# -----------------------------------------------------
# Basic usage with matrix and cluster vector
# -----------------------------------------------------
set.seed(1)
mat <- matrix(runif(100), 10, 10)
rownames(mat) <- colnames(mat) <- LETTERS[1:10]

cs <- cluster_summary(mat, c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3))
cs            # cluster-level (macro) weights
#> MCML Network
#> ============
#> Type: aggregate  | Method: sum 
#> Nodes: 10  | Clusters: 3 
#> 
#> Clusters:
#>   1 (3): A, B, C
#>   2 (3): D, E, F
#>   3 (4): G, H, I, J
#> 
#> Macro (cluster-level) weights:
#>        1      2      3
#> 1 4.0786 5.6031 5.6788
#> 2 4.4388 3.9691 6.6812
#> 3 6.7692 5.8665 8.6995
summary(cs)   # one row per cluster
#>   cluster size within_total between_out between_in
#> 1       1    3     4.078562    11.28182   11.20801
#> 2       2    3     3.969112    11.12004   11.46957
#> 3       3    4     8.699460    12.63571   12.35999

# -----------------------------------------------------
# Named list clusters (more readable)
# -----------------------------------------------------
clusters <- list(
  Alpha = c("A", "B", "C"),
  Beta = c("D", "E", "F"),
  Gamma = c("G", "H", "I", "J")
)
cluster_summary(mat, clusters)
#> MCML Network
#> ============
#> Type: aggregate  | Method: sum 
#> Nodes: 10  | Clusters: 3 
#> 
#> Clusters:
#>   Alpha (3): A, B, C
#>   Beta (3): D, E, F
#>   Gamma (4): G, H, I, J
#> 
#> Macro (cluster-level) weights:
#>        Alpha   Beta  Gamma
#> Alpha 4.0786 5.6031 5.6788
#> Beta  4.4388 3.9691 6.6812
#> Gamma 6.7692 5.8665 8.6995

# -----------------------------------------------------
# A netobject as input: its weight matrix is aggregated
# -----------------------------------------------------
seqs <- data.frame(
  T1 = c("A", "C", "B", "D"), T2 = c("B", "D", "A", "C"),
  T3 = c("C", "A", "D", "B")
)
net <- build_network(seqs, method = "relative")
cluster_summary(net, list(G1 = c("A", "B"), G2 = c("C", "D")))
#> MCML Network
#> ============
#> Type: aggregate  | Method: sum 
#> Nodes: 4  | Clusters: 2 
#> 
#> Clusters:
#>   G1 (2): A, B
#>   G2 (2): C, D
#> 
#> Macro (cluster-level) weights:
#>    G1 G2
#> G1  1  1
#> G2  1  1

# -----------------------------------------------------
# Different aggregation methods
# -----------------------------------------------------
summary(cluster_summary(mat, clusters, method = "sum"))   # total flow
#>   cluster size within_total between_out between_in
#> 1   Alpha    3     4.078562    11.28182   11.20801
#> 2    Beta    3     3.969112    11.12004   11.46957
#> 3   Gamma    4     8.699460    12.63571   12.35999
summary(cluster_summary(mat, clusters, method = "mean"))  # average
#>   cluster size within_total between_out between_in
#> 1   Alpha    3     4.078562    1.095792   1.057301
#> 2    Beta    3     3.969112    1.049971   1.111438
#> 3   Gamma    4     8.699460    1.052976   1.029999
summary(cluster_summary(mat, clusters, method = "max"))   # strongest
#>   cluster size within_total between_out between_in
#> 1   Alpha    3     4.078562    1.774085   1.900114
#> 2    Beta    3     3.969112    1.800406   1.655449
#> 3   Gamma    4     8.699460    1.786146   1.805074

# -----------------------------------------------------
# Skip within-cluster computation for speed
# -----------------------------------------------------
cluster_summary(mat, clusters, compute_within = FALSE)
#> MCML Network
#> ============
#> Type: aggregate  | Method: sum 
#> Nodes: 10  | Clusters: 3 
#> 
#> Clusters:
#>   Alpha (3): A, B, C
#>   Beta (3): D, E, F
#>   Gamma (4): G, H, I, J
#> 
#> Macro (cluster-level) weights:
#>        Alpha   Beta  Gamma
#> Alpha 4.0786 5.6031 5.6788
#> Beta  4.4388 3.9691 6.6812
#> Gamma 6.7692 5.8665 8.6995

# -----------------------------------------------------
# Promote the layers to netobjects
# (as_tna() stores the aggregated weights as they are)
# -----------------------------------------------------
as_tna(cluster_summary(mat, clusters, method = "sum"))
#> Group Networks (4 groups)
#> 
#>   Group  Nodes  Edges  Weights
#>   macro  3      9      [3.969, 8.699]
#>   Alpha  3      9      [0.177, 0.935]
#>   Beta   3      9      [0.071, 0.827]
#>   Gamma  4      16     [0.084, 0.961]
```
