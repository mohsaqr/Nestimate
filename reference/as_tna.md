# Promote the Layers of an mcml to Networks

`as_tna()` turns an `mcml` object into a `netobject_group`: one network
for the cluster-level (macro) layer, followed by one per cluster for the
within-cluster layers. The stored weights are carried over as they are –
nothing is re-normalised – so the networks hold the aggregation and
`type` chosen when the `mcml` was built. Despite its name it returns
Nestimate networks, not `tna` models; a `tna` model passed to it is
returned unchanged.

This is the step that lets an MCML result reach the verbs that take a
group of networks: printing, the network metrics of
[`summary()`](https://rdrr.io/r/base/summary.html),
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
(for an `mcml` built from sequences, whose layers keep them; the layers
of a
[`cluster_summary()`](https://pak.dynasite.org/Nestimate/reference/cluster_summary.md)
result carry none) and drawing with cograph.

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

  An `mcml` object from
  [`cluster_summary()`](https://pak.dynasite.org/Nestimate/reference/cluster_summary.md)
  or
  [`build_mcml()`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md).
  The default method also takes a `tna` model.

- ...:

  Not used.

- expand:

  For the `mcml` method: clusters whose member states replace the
  collapsed cluster node in the macro network, as in
  [`macro_network()`](https://pak.dynasite.org/Nestimate/reference/macro_network.md).
  Default `NULL` (the macro network stays fully collapsed). The
  within-cluster networks are not affected. It needs an `mcml` that
  [`build_mcml()`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
  built from a wide sequence frame, a sequence matrix or an event log;
  one from
  [`cluster_summary()`](https://pak.dynasite.org/Nestimate/reference/cluster_summary.md),
  an edge list, a `netobject` or a `tna` model raises an error of class
  `nestimate_no_expand_source`.

## Value

A `netobject_group`: a named list whose first element, `macro`, is the k
x k cluster-level network, followed by one network per cluster. Print it
for one row per layer (nodes, edges, weight range), call
[`summary()`](https://rdrr.io/r/base/summary.html) for the network
metrics of each layer, and draw a layer with
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html).
Each element is a `netobject` (class
`c("netobject", "cograph_network")`) with the layer's weights, initial
probabilities, node and edge tables and data, and records
`method = "relative"` when the `mcml` holds row-normalised weights
(`type = "tna"`), `method = "co_occurrence"` for symmetrised counts
(`type = "cooccurrence"`), and `method = "frequency"` for counts (`type`
`"raw"` or `"frequency"`) and for the aggregates of
[`cluster_summary()`](https://pak.dynasite.org/Nestimate/reference/cluster_summary.md).
With `expand`, the `macro` element is the mixed-resolution network of
[`macro_network()`](https://pak.dynasite.org/Nestimate/reference/macro_network.md)
instead of the fully collapsed one. The default method returns a `tna`
model unchanged and raises an error for any other class.

## Details

### Workflow


    # Full MCML workflow
    net <- build_network(data, method = "relative")
    cs   <- cluster_summary(net, clusters = group_assignments)
    nets <- as_tna(cs)

    # Every layer is an ordinary netobject
    print(nets)      # one line per layer
    summary(nets)    # network metrics per layer

### Sink nodes

Every cluster is returned, whatever its row sums. A node with no
outgoing weight is a legitimate sink (a terminal state), and its row
stays all zero in the network, for `method = "relative"` and
`"frequency"` alike: the stored weights are never re-normalised, so a
sink row needs no special handling.

## See also

[`cluster_summary`](https://pak.dynasite.org/Nestimate/reference/cluster_summary.md)
and
[`build_mcml`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
to create the input object,
[`macro_network`](https://pak.dynasite.org/Nestimate/reference/macro_network.md)
for a macro layer with one cluster expanded,
[`as_networks`](https://pak.dynasite.org/Nestimate/reference/as_networks.md)
for the psychometric-network counterpart

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
#>              SD Out-Strength 0.5072 0.1253  0.7034 0.4867
#>             Mean In-Strength  6.181  1.122   1.207  1.353
#>               SD In-Strength 0.8432 0.6844 0.08534 0.2021
#>              Mean Out-Degree      3      2       2      2
#>                SD Out-Degree      0      0       0      0
#>  Centralization (Out-Degree)      0      0       0      0
#>   Centralization (In-Degree)      0      0       0      0
#>                  Reciprocity      1      1       1      1
```
