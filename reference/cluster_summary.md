# Cluster Summary Statistics

Aggregates node-level network weights to cluster-level summaries.
Computes both between-cluster transitions (how clusters connect to each
other) and within-cluster transitions (how nodes connect within each
cluster).

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

  Network input. Accepts multiple formats:

  matrix

  :   Numeric adjacency/weight matrix. Row and column names are used as
      node labels. Values represent edge weights (e.g., transition
      counts, co-occurrence frequencies, or probabilities).

  netobject

  :   A network object built by
      [`build_network`](https://saqr.me/Nestimate/reference/build_network.md).
      Its weight matrix `x$weights` is aggregated; a bare
      `cograph_network` is coerced to a netobject first.

  tna

  :   A tna object from the tna package. Extracts `x$weights`.

  mcml

  :   An `mcml` object is returned unchanged.

- clusters:

  Cluster/group assignments for nodes. Accepts multiple formats:

  NULL

  :   (default) Not usable here: `clusters` is required and `NULL`
      raises an error. Auto-detection of a cluster column in a
      `netobject`'s node table happens in
      [`build_mcml`](https://saqr.me/Nestimate/reference/build_mcml.md),
      which then calls this function with the detected assignment.

  vector

  :   Cluster membership for each node, in the same order as the matrix
      rows/columns. Can be numeric (1, 2, 3) or character ("A", "B").
      Cluster names will be derived from unique values. Example:
      `c(1, 1, 2, 2, 3, 3)` assigns first two nodes to cluster 1.

  data.frame

  :   A data frame where the first column contains node names and the
      second column contains group/cluster names. Example:
      `data.frame(node = c("A", "B", "C"), group = c("G1", "G1", "G2"))`

  named list

  :   Explicit mapping of cluster names to node labels. List names
      become cluster names, values are character vectors of node labels
      that must match matrix row/column names. Example:
      `list(Alpha = c("A", "B"), Beta = c("C", "D"))`

- method:

  Aggregation method for combining edge weights within/between clusters.
  Controls how multiple node-to-node edges are summarized:

  "sum"

  :   (default) Sum of all edge weights. Best for count data (e.g.,
      transition frequencies). Preserves total flow.

  "mean"

  :   Average edge weight. Best when cluster sizes differ and you want
      to control for size. Note: when input is already a transition
      matrix (rows sum to 1), "mean" avoids size bias. Example: cluster
      with 5 nodes won't have 5x the weight of cluster with 1 node.

  "median"

  :   Median edge weight. Robust to outliers.

  "max"

  :   Maximum edge weight. Captures strongest connection.

  "min"

  :   Minimum edge weight. Captures weakest connection.

  "density"

  :   Sum of (non-zero) edge weights divided by the number of possible
      edges between the two clusters (`n_i * n_j`). Normalizes by
      cluster size combinations. Because zero/`NA` edges are stripped
      before aggregation, this equals `"mean"` exactly when the
      cluster-pair block is fully dense (no zero edges), and is strictly
      smaller than `"mean"` when zero edges are present (it divides by
      the larger possible-edge count).

  "geomean"

  :   Geometric mean of positive weights. Useful for multiplicative
      processes.

- directed:

  Logical. If `TRUE` (default), treat network as directed. A-\>B and
  B-\>A are separate edges. If `FALSE`, edges are undirected and the
  matrix is symmetrized before processing.

- compute_within:

  Logical. If `TRUE` (default), compute within-cluster transition
  matrices for each cluster. Each cluster gets its own n_i x n_i matrix
  showing internal node-to-node transitions. Set to `FALSE` to skip this
  computation for better performance when only between-cluster summary
  is needed.

## Value

An `mcml` object (S3 class): a list with

- macro:

  An `mcml_layer` holding the cluster-level network: `$weights`, the k x
  k matrix whose entry (i, j) is the aggregation (per `method`) of all
  edges from nodes in cluster i to nodes in cluster j, with the diagonal
  holding the within-cluster edges – pure arithmetic, no row
  normalization; `$inits`, the length-k column sums of that matrix
  normalized to sum to 1; `$labels`, the cluster names; and `$data`,
  `NULL` on this path.

- clusters:

  Named list with one `mcml_layer` per cluster. Its `$weights` is the
  n_i x n_i submatrix of the nodes in that cluster and its `$inits` the
  normalized column sums of that submatrix. `NULL` when
  `compute_within = FALSE`.

- cluster_members:

  Named list mapping cluster names to their member node labels, e.g.
  `list(A = c("n1", "n2"), B = c("n3", "n4", "n5"))`.

- edges:

  `NULL` on this path – a matrix carries no node-level transitions. The
  sequence and edge-list paths of
  [`build_mcml`](https://saqr.me/Nestimate/reference/build_mcml.md) fill
  in a tidy edge table here.

- meta:

  List with `type` (always `"aggregate"` here), `method`, `directed`,
  `n_nodes`, `n_clusters`, `cluster_sizes` (named integer vector) and
  `source` (`"matrix"`).

## Details

This is the core function for Multi-Cluster Multi-Level (MCML) analysis.
Use [`as_tna()`](https://saqr.me/Nestimate/reference/as_tna.md) to
convert results to tna objects for further analysis with the tna
package.

### Workflow

Typical MCML analysis workflow:


    # 1. Create the node-level network
    net <- build_network(data, method = "relative")

    # 2. Aggregate its edges to cluster level (arithmetic aggregation)
    cs <- cluster_summary(net, clusters = group_assignments, method = "sum")

    # 3. Read the result
    print(cs)      # macro weights
    summary(cs)    # one row per cluster

    # 4. Promote the layers to netobjects for downstream verbs
    nets <- as_tna(cs)

### Between-Cluster Matrix Structure

The `macro$weights` matrix has clusters as both rows and columns:

- Off-diagonal (row i, col j): Aggregated weight from cluster i to
  cluster j

- Diagonal (row i, col i): Within-cluster total (aggregation of internal
  edges)

Rows are NOT normalized. Entries are elementwise aggregates produced by
`method`. If the caller wants probabilities, they should normalize
downstream (e.g. via
[`as_tna()`](https://saqr.me/Nestimate/reference/as_tna.md)). Mixing an
arithmetic aggregation with row-normalization here (the old
`type = "tna"` combined with `method = "min"` / `"mean"` etc.) produces
numbers that sum to 1 per row but are not a probability distribution
over any process; that silently-wrong combination is why `type` was
removed from the matrix path. The sequence and edgelist paths of
[`build_mcml()`](https://saqr.me/Nestimate/reference/build_mcml.md) keep
`type`, where the aggregation is always counts and the post-processing
chooses between well-defined network constructions.

### Choosing method

|  |  |  |
|----|----|----|
| **Input data** | **Recommended method** | **Reason** |
| Edge counts | `"sum"` | Preserves total flow between clusters |
| Transition matrix | `"mean"` | Avoids cluster size bias |
| Correlation matrix | `"mean"` | Average correlations |
| Dense weighted | `"max"` / `"median"` | Robust summary |

## See also

[`build_mcml`](https://saqr.me/Nestimate/reference/build_mcml.md) to
build an mcml from raw transitions instead of a weight matrix,
[`as_tna`](https://saqr.me/Nestimate/reference/as_tna.md) to promote the
layers to netobjects,
[`macro_network`](https://saqr.me/Nestimate/reference/macro_network.md)
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
#> 1       1    3     2.984822    11.28182   11.20801
#> 2       2    3     3.153709    11.12004   11.46957
#> 3       3    4     6.980503    12.63571   12.35999

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
#> 1   Alpha    3     2.984822    11.28182   11.20801
#> 2    Beta    3     3.153709    11.12004   11.46957
#> 3   Gamma    4     6.980503    12.63571   12.35999
summary(cluster_summary(mat, clusters, method = "mean"))  # average
#>   cluster size within_total between_out between_in
#> 1   Alpha    3     2.984822    1.095792   1.057301
#> 2    Beta    3     3.153709    1.049971   1.111438
#> 3   Gamma    4     6.980503    1.052976   1.029999
summary(cluster_summary(mat, clusters, method = "max"))   # strongest
#>   cluster size within_total between_out between_in
#> 1   Alpha    3     2.984822    1.774085   1.900114
#> 2    Beta    3     3.153709    1.800406   1.655449
#> 3   Gamma    4     6.980503    1.786146   1.805074

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
