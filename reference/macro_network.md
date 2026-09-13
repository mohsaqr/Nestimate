# Cluster-Level Network, With One Cluster Expanded

Returns the macro (cluster-level) network of an `mcml`, optionally with
one or more clusters expanded back into their member states. Every other
cluster stays collapsed to a single node, so the result is a network at
mixed resolution: the cluster of interest in detail, its context in
summary.

## Usage

``` r
macro_network(x, expand = NULL, method = "relative", ...)
```

## Arguments

- x:

  An `mcml` built from sequence data. A matrix-derived `mcml` carries no
  node-level data and cannot be expanded.

- expand:

  Names of clusters to expand into their member states. `NULL` (default)
  collapses every cluster, reproducing the macro layer; `"all"` or
  `TRUE` expands every cluster, so each state is its own node.

- method:

  Estimator passed to
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md).
  Default `"relative"` (row-normalised transitions).

- ...:

  Further arguments passed to
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md).

## Value

A `netobject` (also a `cograph_network`) whose nodes are the collapsed
clusters plus the member states of any expanded cluster, with weights
re-counted from the sequence data by
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md).
`$node_groups` is a two-column data frame (`node`, `group`) mapping
every node to its cluster, and the same labels are a factor in
`$nodes$groups`, so the result plots grouped; an expanded cluster's
states each map to that cluster, a collapsed cluster maps to itself.
`$expanded` records the cluster names that were expanded (`NULL` when
none were).

## See also

[`build_mcml`](https://saqr.me/Nestimate/reference/build_mcml.md),
[`as_tna`](https://saqr.me/Nestimate/reference/as_tna.md)

## Examples

``` r
seqs <- data.frame(
  t1 = c("A", "C", "A", "B"), t2 = c("B", "D", "C", "A"),
  t3 = c("C", "A", "D", "C"), stringsAsFactors = FALSE
)
mc <- build_mcml(seqs, clusters = list(G1 = c("A", "B"), G2 = c("C", "D")))
macro_network(mc)                    # every cluster collapsed
#> Transition Network (relative probabilities) [directed]
#>   Weights: [0.333, 0.667]  |  mean: 0.500
#> 
#>   Weight matrix:
#>         G1    G2
#>   G1 0.400 0.600
#>   G2 0.333 0.667 
#> 
#>   Initial probabilities:
#>   G1            0.750  ████████████████████████████████████████
#>   G2            0.250  █████████████
macro_network(mc, expand = "G2")     # G2 shown as C and D
#> Transition Network (relative probabilities) [directed]
#>   Weights: [0.400, 1.000]  |  mean: 0.750
#> 
#>   Weight matrix:
#>        C D  G1
#>   C  0.0 1 0.0
#>   D  0.0 0 1.0
#>   G1 0.6 0 0.4 
#> 
#>   Initial probabilities:
#>   G1            0.750  ████████████████████████████████████████
#>   C             0.250  █████████████
#>   D             0.000  
```
