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
  type = c("tna", "cooccurrence", "semi_markov", "raw"),
  directed = TRUE,
  compute_within = TRUE
)

csum(
  x,
  clusters = NULL,
  method = c("sum", "mean", "median", "max", "min", "density", "geomean"),
  type = c("tna", "cooccurrence", "semi_markov", "raw"),
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

  :   A cograph network object. The function extracts the weight matrix
      from `x$weights` or converts via `to_matrix()`. Clusters can be
      auto-detected from node attributes.

  tna

  :   A tna object from the tna package. Extracts `x$weights`.

  cluster_summary

  :   If already a cluster_summary, returns unchanged.

- clusters:

  Cluster/group assignments for nodes. Accepts multiple formats:

  NULL

  :   (default) Auto-detect from netobject. Looks for columns named
      'clusters', 'cluster', 'groups', or 'group' in `x$nodes`. Throws
      an error if no cluster column is found. This option only works
      when `x` is a netobject.

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

  :   Sum divided by number of possible edges. Normalizes by cluster
      size combinations.

  "geomean"

  :   Geometric mean of positive weights. Useful for multiplicative
      processes.

- type:

  Post-processing applied to aggregated weights. Determines the
  interpretation of the resulting matrices:

  "tna"

  :   (default) Row-normalize so each row sums to 1. Creates transition
      probabilities suitable for Markov chain analysis. Interpretation:
      "Given I'm in cluster A, what's the probability of transitioning
      to cluster B?" Required for use with tna package functions.
      Diagonal represents within-cluster transition probability.

  "raw"

  :   No normalization. Returns aggregated counts/weights as-is. Use for
      frequency analysis or when you need raw counts. Compatible with
      igraph's contract + simplify output.

  "cooccurrence"

  :   Symmetrize the matrix: (A + t(A)) / 2. For undirected
      co-occurrence analysis.

  "semi_markov"

  :   Row-normalize with duration weighting. For semi-Markov process
      analysis.

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

A `cluster_summary` object (S3 class) containing:

- between:

  List with two elements:

  weights

  :   k x k matrix of cluster-to-cluster weights, where k is the number
      of clusters. Row i, column j contains the aggregated weight from
      cluster i to cluster j. Diagonal contains within-cluster totals.
      Processing depends on `type`.

  inits

  :   Numeric vector of length k. Initial state distribution across
      clusters, computed from column sums of the original matrix.
      Represents the proportion of incoming edges to each cluster.

- within:

  Named list with one element per cluster. Each element contains:

  weights

  :   n_i x n_i matrix for nodes within that cluster. Shows internal
      transitions between nodes in the same cluster.

  inits

  :   Initial distribution within the cluster.

  NULL if `compute_within = FALSE`.

- clusters:

  Named list mapping cluster names to their member node labels. Example:
  `list(A = c("n1", "n2"), B = c("n3", "n4", "n5"))`

- meta:

  List of metadata:

  type

  :   The `type` argument used ("tna", "raw", etc.)

  method

  :   The `method` argument used ("sum", "mean", etc.)

  directed

  :   Logical, whether network was treated as directed

  n_nodes

  :   Total number of nodes in original network

  n_clusters

  :   Number of clusters

  cluster_sizes

  :   Named vector of cluster sizes

See `cluster_summary`.

## Details

This is the core function for Multi-Cluster Multi-Level (MCML) analysis.
Use
[`as_tna()`](https://mohsaqr.github.io/Nestimate/reference/as_tna.md) to
convert results to tna objects for further analysis with the tna
package.

### Workflow

Typical MCML analysis workflow:

    # 1. Create network
    net <- build_network(data, method = "relative")
    net$nodes$clusters <- group_assignments

    # 2. Compute cluster summary
    cs <- cluster_summary(net, type = "tna")

    # 3. Convert to tna models
    tna_models <- as_tna(cs)

    # 4. Analyze/visualize
    plot(tna_models$macro)
    tna::centralities(tna_models$macro)

### Between-Cluster Matrix Structure

The `macro$weights` matrix has clusters as both rows and columns:

- Off-diagonal (row i, col j): Aggregated weight from cluster i to
  cluster j

- Diagonal (row i, col i): Within-cluster total (sum of internal edges
  in cluster i)

When `type = "tna"`, rows sum to 1 and diagonal values represent
"retention rate" - the probability of staying within the same cluster.

### Choosing method and type

|                    |                           |                                                   |
|--------------------|---------------------------|---------------------------------------------------|
| **Input data**     | **Recommended**           | **Reason**                                        |
| Edge counts        | method="sum", type="tna"  | Preserves total flow, normalizes to probabilities |
| Transition matrix  | method="mean", type="tna" | Avoids cluster size bias                          |
| Frequencies        | method="sum", type="raw"  | Keep raw counts for analysis                      |
| Correlation matrix | method="mean", type="raw" | Average correlations                              |

## See also

[`as_tna()`](https://mohsaqr.github.io/Nestimate/reference/as_tna.md) to
convert results to tna objects,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) for two-layer
visualization, [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
for flat cluster visualization

## Examples

``` r
# -----------------------------------------------------
# Basic usage with matrix and cluster vector
# -----------------------------------------------------
mat <- matrix(runif(100), 10, 10)
rownames(mat) <- colnames(mat) <- LETTERS[1:10]

clusters <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3)
cs <- cluster_summary(mat, clusters)

# Access results
cs$macro$weights    # 3x3 cluster transition matrix
#>           1         2         3
#> 1 0.3096851 0.3024384 0.3878765
#> 2 0.3022261 0.2590855 0.4386884
#> 3 0.2565129 0.2628054 0.4806818
cs$macro$inits      # Initial distribution
#>         1         2         3 
#> 0.2868688 0.2746885 0.4384427 
cs$clusters$`1`$weights # Within-cluster 1 transitions
#>            A         B         C
#> A 0.35884842 0.6027075 0.0384441
#> B 0.07374815 0.4482965 0.4779554
#> C 0.31270223 0.3622318 0.3250660
cs$meta               # Metadata
#> $type
#> [1] "tna"
#> 
#> $method
#> [1] "sum"
#> 
#> $directed
#> [1] TRUE
#> 
#> $n_nodes
#> [1] 10
#> 
#> $n_clusters
#> [1] 3
#> 
#> $cluster_sizes
#> 1 2 3 
#> 3 3 4 
#> 

# -----------------------------------------------------
# Named list clusters (more readable)
# -----------------------------------------------------
clusters <- list(
  Alpha = c("A", "B", "C"),
  Beta = c("D", "E", "F"),
  Gamma = c("G", "H", "I", "J")
)
cs <- cluster_summary(mat, clusters, type = "tna")
cs$macro$weights    # Rows/cols named Alpha, Beta, Gamma
#>           Alpha      Beta     Gamma
#> Alpha 0.3096851 0.3024384 0.3878765
#> Beta  0.3022261 0.2590855 0.4386884
#> Gamma 0.2565129 0.2628054 0.4806818
cs$clusters$Alpha       # Within Alpha cluster
#> $weights
#>            A         B         C
#> A 0.35884842 0.6027075 0.0384441
#> B 0.07374815 0.4482965 0.4779554
#> C 0.31270223 0.3622318 0.3250660
#> 
#> $inits
#>         A         B         C 
#> 0.2358790 0.4452836 0.3188374 
#> 
#> $labels
#> [1] "A" "B" "C"
#> 
#> $data
#> NULL
#> 

# -----------------------------------------------------
# Auto-detect clusters from netobject
# -----------------------------------------------------
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:10], 30, TRUE), V2 = sample(LETTERS[1:10], 30, TRUE),
  V3 = sample(LETTERS[1:10], 30, TRUE)
)
net <- build_network(seqs, method = "relative")
cs2 <- cluster_summary(net, c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3))
# }

# -----------------------------------------------------
# Different aggregation methods
# -----------------------------------------------------
cs_sum <- cluster_summary(mat, clusters, method = "sum")   # Total flow
cs_mean <- cluster_summary(mat, clusters, method = "mean") # Average
cs_max <- cluster_summary(mat, clusters, method = "max")   # Strongest

# -----------------------------------------------------
# Raw counts vs TNA probabilities
# -----------------------------------------------------
cs_raw <- cluster_summary(mat, clusters, type = "raw")
cs_tna <- cluster_summary(mat, clusters, type = "tna")

rowSums(cs_raw$macro$weights)  # Various sums
#>    Alpha     Beta    Gamma 
#> 16.59946 14.45172 19.78782 
rowSums(cs_tna$macro$weights)  # All equal to 1
#> Alpha  Beta Gamma 
#>     1     1     1 

# -----------------------------------------------------
# Skip within-cluster computation for speed
# -----------------------------------------------------
cs_fast <- cluster_summary(mat, clusters, compute_within = FALSE)
cs_fast$clusters  # NULL
#> NULL

# -----------------------------------------------------
# Convert to tna objects for tna package
# -----------------------------------------------------
cs <- cluster_summary(mat, clusters, type = "tna")
tna_models <- as_tna(cs)
# tna_models$macro      # tna object
# tna_models$clusters$Alpha # tna object
# \donttest{
mat <- matrix(runif(16), 4, 4)
rownames(mat) <- colnames(mat) <- LETTERS[1:4]
csum(mat, c(1, 1, 2, 2))
#> MCML Network
#> ============
#> Type: tna  | Method: sum 
#> Nodes: 4  | Clusters: 2 
#> 
#> Clusters:
#>   1 (2): A, B
#>   2 (2): C, D
#> 
#> Macro (cluster-level) weights:
#>        1      2
#> 1 0.5383 0.4617
#> 2 0.4832 0.5168
# }
```
