# Edge Betweenness Network

`net_edge_betweenness()` returns a copy of a network in which each
edge's weight is replaced by its betweenness: the number of shortest
paths between ordered pairs of distinct nodes that pass through that
edge, with a fractional share when several shortest paths tie (Girvan &
Newman 2002). It takes a single network or a `netobject_group`, which is
scored one network at a time.

With `invert = TRUE` (the default) a heavier edge is a shorter step, so
the shortest path between two states of a transition network is the most
probable route, not the one with the fewest hops. Only positive weights
count as edges, and a self-loop never lies on a shortest path, so its
betweenness is 0. A directed network gives an asymmetric betweenness
matrix; an undirected one a symmetric matrix, each cell holding the
betweenness of that edge. The result keeps the source network's
directedness (directed when it records none).

An edge with high betweenness carries many of the shortest routes
between states; Girvan & Newman (2002) read such edges as the bridges
between groups of states. To compare the edge betweenness of two groups,
pass two results to
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
or
[`bayes_compare()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md),
which resample the source networks behind them; both need two results
with the same source method and `invert`, and
[`bayes_compare()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
needs `"relative"` source networks.

## Usage

``` r
net_edge_betweenness(x, invert = TRUE, ...)

# S3 method for class 'netobject'
net_edge_betweenness(x, invert = TRUE, ...)

# S3 method for class 'netobject_group'
net_edge_betweenness(x, invert = TRUE, ...)

# Default S3 method
net_edge_betweenness(x, invert = TRUE, ...)

# S3 method for class 'net_edge_betweenness'
plot(x, style = c("bar", "forest", "delta"), top_n = NULL, labels = TRUE, ...)
```

## Arguments

- x:

  A `netobject` or `netobject_group`; any other input raises an error.
  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: a `net_edge_betweenness` network.

- invert:

  Whether to turn weights into distances by `1 / weight` before finding
  shortest paths, which suits probability and frequency networks. With
  `FALSE`, the weight itself is the distance; use it when the weights
  already are distances. Default `TRUE`.

- ...:

  Not used.

- style:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: the kind of chart. `"bar"` (default) draws one horizontal bar
  per edge; `"forest"` a stem from zero to a point, with a dashed line
  at the mean betweenness of the edges shown; `"delta"` each edge's
  deviation from that mean as a bar pointing right (above the mean,
  blue) or left (below, red).

- top_n:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: the number of highest edges to keep. Default `NULL` (every
  edge with non-zero betweenness).

- labels:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: whether to print each edge's value (its deviation, for
  `"delta"`) beside it. Default `TRUE`.

## Value

For a `netobject`: a network of class
`c("net_edge_betweenness", "netobject", "cograph_network")` with method
`"edge_betweenness"`, whose weights are the edge betweenness scores. It
prints like any network;
[`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md)
gives one row per edge,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) ranks the edges
and
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
draws the network. The source network's data, metadata, parameters,
weights and initial probabilities are kept on it, which is what lets
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
and
[`bayes_compare()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
test differences in edge betweenness between two such networks. For a
`netobject_group`: a `netobject_group` of these networks, one per group.

For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method:
a `ggplot` object.

## Details

Shortest paths come from all-pairs Floyd-Warshall on the distance matrix
(1 / weight, or the weight with `invert = FALSE`), counting tied paths
with a relative tolerance so that floating-point rounding does not break
a tie. The values agree with
[`tna::betweenness_network()`](https://sonsoles.me/tna/reference/betweenness_network.html)
for transition networks, and on an undirected network each cell equals
[`igraph::edge_betweenness()`](https://r.igraph.org/reference/betweenness.html)
for that edge; the function is named differently to avoid masking
either.

## Methods

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) ranks the edges
with non-zero betweenness, highest first, one row per edge
(`from -> to`). It needs no cograph;
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
draws the same object as a node-link diagram. A network whose edges all
have zero betweenness raises an error.

## References

Girvan, M. & Newman, M. E. J. (2002). Community structure in social and
biological networks. *Proceedings of the National Academy of Sciences*,
99(12), 7821–7826.
[doi:10.1073/pnas.122653799](https://doi.org/10.1073/pnas.122653799)

## Examples

``` r
seqs <- data.frame(
  V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
  V3 = c("C","A","C","B"))
net <- build_network(seqs, method = "relative")
eb  <- net_edge_betweenness(net)
extract_edges(eb)
#>   from to weight
#> 1    C  A      3
#> 2    A  B      3
#> 3    B  C      3
```
