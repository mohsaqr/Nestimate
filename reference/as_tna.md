# Convert cluster_summary to tna Objects

Converts a `cluster_summary` object to proper tna objects that can be
used with all functions from the tna package. Creates both a
between-cluster tna model (cluster-level transitions) and within-cluster
tna models (internal transitions within each cluster).

## Usage

``` r
as_tna(x)

# S3 method for class 'mcml'
as_tna(x)

# Default S3 method
as_tna(x)
```

## Arguments

- x:

  A `cluster_summary` object created by
  [`cluster_summary`](https://saqr.me/Nestimate/reference/cluster_summary.md).
  The aggregated weights are passed to
  [`tna::tna()`](http://sonsoles.me/tna/reference/build_model.md), which
  row-normalises them as needed.

## Value

A `cluster_tna` object (S3 class) containing:

- between:

  A tna object representing cluster-level transitions. Contains
  `$weights` (k x k transition matrix), `$inits` (initial distribution),
  and `$labels` (cluster names). Use this for analyzing how
  learners/entities move between high-level groups or phases.

- within:

  Named list of tna objects, one per cluster. Each tna object represents
  internal transitions within that cluster. Contains `$weights` (n_i x
  n_i matrix), `$inits` (initial distribution), and `$labels` (node
  labels). Clusters with single nodes or zero-row nodes are excluded
  (tna requires positive row sums).

A `netobject_group` with data preserved from each sub-network.

A `tna` object constructed from the input.

## Details

This is the final step in the MCML workflow, enabling full integration
with the tna package for centrality analysis, bootstrap validation,
permutation tests, and visualization.

### Requirements

The tna package must be installed. If not available, the function throws
an error with installation instructions.

### Workflow


    # Full MCML workflow
    net <- build_network(data, method = "relative")
    net$nodes$clusters <- group_assignments
    cs <- cluster_summary(net)
    tna_models <- as_tna(cs)

    # Now use tna package functions
    plot(tna_models$macro)
    tna::centralities(tna_models$macro)
    tna::bootstrap(tna_models$macro, iter = 1000)

    # Analyze within-cluster patterns
    plot(tna_models$clusters$ClusterA)
    tna::centralities(tna_models$clusters$ClusterA)

### Excluded Clusters

A within-cluster network is dropped only when row-normalisation would
fail. Specifically, when the recorded `net_method` is `"relative"`
(row-stochastic transitions) and any node in the cluster has zero
outgoing weight, that cluster is excluded from `$clusters` and a
[`warning()`](https://rdrr.io/r/base/warning.html) is emitted listing
the dropped cluster names. For `net_method = "frequency"` (raw counts),
a zero-row node is a legitimate sink and the cluster is retained. The
macro / between-cluster network always includes every cluster regardless
of the per-cluster drop decisions.

If a cluster you expect to see is missing from the returned `$clusters`,
check the warning output and consider building with `type = "raw"`
(which carries through to a frequency-method netobject and skips the
drop) or inspect `rowSums(x$clusters[[cl]]$weights)`.

## See also

[`cluster_summary`](https://saqr.me/Nestimate/reference/cluster_summary.md)
to create the input object,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) for
visualization without conversion,
[`tna::tna`](http://sonsoles.me/tna/reference/build_model.md) for the
underlying tna constructor

## Examples

``` r
mat <- matrix(runif(36), 6, 6)
rownames(mat) <- colnames(mat) <- LETTERS[1:6]
clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"), G3 = c("E", "F"))
cs <- cluster_summary(mat, clusters)
tna_models <- as_tna(cs)
tna_models
#> Group Networks (4 groups)
#> 
#>   Group  Nodes  Edges  Weights
#>   macro  3      9      [0.907, 3.076]
#>   G1     2      4      [0.018, 0.777]
#>   G2     2      4      [0.056, 0.736]
#>   G3     2      4      [0.272, 0.905]
tna_models$macro$weights
#>           G1       G2       G3
#> G1 0.9065801 2.585793 2.645868
#> G2 2.5352275 1.534746 2.120124
#> G3 2.2492125 3.076226 2.460135
```
