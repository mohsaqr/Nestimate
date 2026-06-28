# Build a grouped node-level network (htna) from data and a clustering

Builds the full node-level network from the original data and attaches a
cluster grouping, producing a single `netobject` in which every actor is
a node and cluster membership labels the actors. This is the node-level
counterpart of
[`build_mcml`](https://saqr.me/Nestimate/reference/build_mcml.md): where
`build_mcml` collapses the network to a cluster-level (macro) summary,
`as_htna` keeps every node and every transition — including the
between-cluster transitions an mcml only retains in aggregate.

## Usage

``` r
as_htna(x, clusters = NULL, method = "relative", ...)

# S3 method for class 'mcml'
as_htna(x, clusters = NULL, method = "relative", data = NULL, ...)

# Default S3 method
as_htna(x, clusters = NULL, method = "relative", ...)
```

## Arguments

- x:

  Data accepted by
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md)
  (sequence data frame, edgelist, transition matrix, `netobject`, or
  `tna`); or an `mcml` object, in which case the original `data` must
  also be supplied and the mcml provides the cluster membership.

- clusters:

  Cluster assignment: a named list of node-name vectors, a per-node
  membership vector, or a two-column data frame. When `NULL` and `x`
  carries node groups (or is an `mcml`), those are used.

- method:

  Estimator passed to
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md).
  Default `"relative"` (row-normalized transitions).

- ...:

  Further arguments forwarded to
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md)
  (e.g. `actor`, `action`, `time` for long-format data).

- data:

  For the `mcml` method, the original data the mcml was built from
  (sequence/edgelist/etc.). Optional when the mcml was built from
  sequence/edgelist data:
  [`build_mcml()`](https://saqr.me/Nestimate/reference/build_mcml.md)
  stashes that source (and the `actor`/`action`/`time` roles), so
  `as_htna(mcml)` works on its own. Required for an mcml built from a
  matrix/aggregate, which retains no node-level data.

## Value

A `netobject` (`cograph_network`) over all nodes, with a `cluster`
column on `$nodes`, `$node_groups` populated, and the membership stored
in the `"cluster_members"` attribute.

## Details

**Why this rebuilds from data.** An `mcml` stores cluster-level data
(the macro sequences are recoded to cluster labels, and the per-cluster
data is filtered to within-cluster nodes), so it does not retain a
faithful node-level transition network. The only faithful source of
node-level between-cluster transitions is the original data. `as_htna()`
therefore rebuilds from data via
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md);
an `mcml` can supply the cluster membership, but the data must be
provided.

The result is a genuine `netobject`, so it supports inference
([`bootstrap_network`](https://saqr.me/Nestimate/reference/bootstrap_network.md),
centrality,
[`permutation`](https://saqr.me/Nestimate/reference/permutation.md)) and
plots directly as a grouped network with cograph:
`cograph::plot_htna(as_htna(data, clusters))`.

## See also

[`build_mcml`](https://saqr.me/Nestimate/reference/build_mcml.md),
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md);
plot with
[`cograph::plot_htna()`](https://sonsoles.me/cograph/reference/plot_htna.html).

## Examples

``` r
seqs <- data.frame(
  t1 = c("A", "C", "E", "B"), t2 = c("B", "D", "F", "A"),
  t3 = c("C", "A", "E", "D"), stringsAsFactors = FALSE
)
clusters <- list(C1 = c("A", "B"), C2 = c("C", "D"), C3 = c("E", "F"))
net <- as_htna(seqs, clusters)
net$nodes$cluster
#> [1] "C1" "C1" "C2" "C2" "C3" "C3"
if (FALSE) { # \dontrun{
cograph::plot_htna(net)
} # }
```
