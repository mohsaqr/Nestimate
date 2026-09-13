# Promote the Layers of an mcml to Networks

Converts an `mcml` object into a `netobject_group`: one netobject for
the cluster-level (macro) layer, and one per cluster for the
within-cluster layers. The stored weights are carried over as they are –
nothing is re-normalised here, so the aggregation chosen when the `mcml`
was built is what the networks hold.

## Usage

``` r
as_tna(x, ...)

# S3 method for class 'mcml'
as_tna(x, expand = NULL, ...)

# Default S3 method
as_tna(x, ...)
```

## Arguments

- x:

  An `mcml` object created by
  [`cluster_summary`](https://saqr.me/Nestimate/reference/cluster_summary.md)
  or [`build_mcml`](https://saqr.me/Nestimate/reference/build_mcml.md).

- ...:

  Passed to methods.

- expand:

  For the `mcml` method, names of clusters whose member states replace
  the collapsed cluster node in the `macro` layer (see
  [`macro_network`](https://saqr.me/Nestimate/reference/macro_network.md)).
  `NULL` (default) keeps the macro fully collapsed. The per-cluster
  layers are unaffected.

## Value

A `netobject_group`: a named list whose first element is `macro` (the k
x k cluster-level network) followed by one element per cluster, each a
`netobject`/`cograph_network` carrying `$weights`, `$inits`, `$nodes`,
`$edges` and the recorded `$method` (`"relative"` for an mcml whose
weights are already row-normalised, `"frequency"` otherwise).

The `mcml` method returns that `netobject_group`, each layer keeping the
data the corresponding `mcml` layer carried. With `expand`, its `macro`
element is the mixed-resolution network of
[`macro_network`](https://saqr.me/Nestimate/reference/macro_network.md)
rather than the fully collapsed one.

The default method returns the input unchanged when it already inherits
from `tna`, and otherwise raises an error.

## Details

This is the step that lets an MCML result flow into the verbs that take
a group of networks (printing, network-metric summaries, rendering with
cograph).

### Workflow


    # Full MCML workflow
    net <- build_network(data, method = "relative")
    cs   <- cluster_summary(net, clusters = group_assignments)
    nets <- as_tna(cs)

    # Every layer is an ordinary netobject
    print(nets)      # one line per layer
    summary(nets)    # network metrics per layer

### Zero-out-degree (sink) nodes

Every cluster is returned, regardless of its row sums. A node with zero
outgoing weight is a legitimate sink (a terminal state); its row in the
wrapped network is left all-zero. This holds for both
`net_method = "relative"` and `"frequency"` – the stored weights are
never re-normalised, so a sink row needs no special handling. Inspect
`rowSums(x$clusters[[cl]]$weights)` to find sink nodes.

## See also

[`cluster_summary`](https://saqr.me/Nestimate/reference/cluster_summary.md)
and [`build_mcml`](https://saqr.me/Nestimate/reference/build_mcml.md) to
create the input object,
[`macro_network`](https://saqr.me/Nestimate/reference/macro_network.md)
for a macro layer with one cluster expanded,
[`as_networks`](https://saqr.me/Nestimate/reference/as_networks.md) for
the psychometric-network counterpart

## Examples

``` r
set.seed(1)
mat <- matrix(runif(36), 6, 6)
rownames(mat) <- colnames(mat) <- LETTERS[1:6]
clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"), G3 = c("E", "F"))
cs <- cluster_summary(mat, clusters)
nets <- as_tna(cs)
nets
#> Group Networks (4 groups)
#> 
#>   Group  Nodes  Edges  Weights
#>   macro  3      9      [1.076, 2.706]
#>   G1     2      4      [0.266, 0.945]
#>   G2     2      4      [0.212, 0.935]
#>   G3     2      4      [0.340, 0.870]
summary(nets)
#> Network metrics by group:
#>                       metric  macro     G1      G2     G3
#>                   Node Count      3      2       2      2
#>                   Edge Count      9      4       4      4
#>              Network Density      1      1       1      1
#>                Mean Distance  1.863 0.6584  0.7162 0.5839
#>            Mean Out-Strength  6.181  1.122   1.207  1.353
#>              SD Out-Strength 0.8432 0.6844 0.08534 0.2021
#>             Mean In-Strength  6.181  1.122   1.207  1.353
#>               SD In-Strength 0.5072 0.1253  0.7034 0.4867
#>              Mean Out-Degree      3      2       2      2
#>                SD Out-Degree      0      0       0      0
#>  Centralization (Out-Degree)      0      0       0      0
#>   Centralization (In-Degree)      0      0       0      0
#>                  Reciprocity      1      1       1      1
```
