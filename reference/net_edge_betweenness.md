# Edge Betweenness Network

Builds a network in which each edge's weight is replaced by its
betweenness: the number of shortest paths between all node pairs that
traverse that edge (fractional when shortest paths tie). This is the
Nestimate counterpart of
[`tna::betweenness_network()`](https://sonsoles.me/tna/reference/betweenness_network.html)
and produces identical values for transition networks; the name differs
to avoid a clash with
[`tna::betweenness_network()`](https://sonsoles.me/tna/reference/betweenness_network.html)
and
[`igraph::edge_betweenness()`](https://r.igraph.org/reference/betweenness.html).

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

  A `netobject` or `netobject_group`. For the
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method: an
  object of class `net_edge_betweenness`.

- invert:

  Logical. Invert weights to distances by `1/w` before computing
  shortest paths? Default `TRUE` (correct for probability and frequency
  networks).

- ...:

  Additional arguments (ignored). In `plot.net_edge_betweenness()`:
  Additional arguments (ignored).

- style:

  Plot style. `"bar"` (default) draws one horizontal bar per edge;
  `"forest"` draws a forest/lollipop chart (a stem from zero to a point)
  with a dashed reference line at the mean betweenness; `"delta"` draws
  each edge's deviation from the mean edge betweenness as a diverging
  bar (above the mean in blue, below in red).

- top_n:

  Integer or `NULL`. Keep only the `top_n` highest edges. Default `NULL`
  (all edges with non-zero betweenness).

- labels:

  Logical. Print the betweenness value beside each edge. Default `TRUE`.

## Value

For a `netobject`: a new network of class
`c("net_edge_betweenness", "netobject", "cograph_network")` whose
`$weights` are the edge-betweenness scores, with
`method = "edge_betweenness"`. Call
[`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md)
on it for a tidy per-edge table, or
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) to render it.
The object preserves source-network metadata so
[`permutation`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
can test edge-betweenness differences by permuting the source networks.
For a `netobject_group`: a `netobject_group` of such networks, one per
group.

In `plot.net_edge_betweenness()`: A `ggplot` object.

## Details

For a probability/transition network the edge weights are transition
probabilities, so they are inverted to distances (`invert = TRUE`)
before path-finding: the geodesic between two states is then the most
probable route rather than the one with the fewest hops. Pass
`invert = FALSE` when the weights already represent distances.

Directedness is taken from the network itself. A directed network yields
an asymmetric betweenness matrix; an undirected (symmetric) network
yields a symmetric one.

## Methods

- `plot.net_edge_betweenness()`: Draws the edges of a
  `net_edge_betweenness` network ranked by their betweenness, as a
  horizontal bar chart. This is the tidy, cograph-free companion to the
  node-link diagram: render the diagram with `cograph::splot(eb)` and
  the ranking with `plot(eb)`.

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
#> 4    B  A      0
#> 5    C  B      0
#> 6    A  C      0
```
