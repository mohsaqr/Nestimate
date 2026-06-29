# Plot edge-betweenness scores

Draws the edges of a
[`net_edge_betweenness`](https://saqr.me/Nestimate/reference/net_edge_betweenness.md)
network ranked by their betweenness, as a horizontal bar chart. This is
the tidy, cograph-free companion to the node-link diagram: render the
diagram with `cograph::splot(eb)` and the ranking with `plot(eb)`.

## Usage

``` r
# S3 method for class 'net_edge_betweenness'
plot(x, style = c("bar", "forest"), top_n = NULL, labels = TRUE, ...)
```

## Arguments

- x:

  A `net_edge_betweenness` network from
  [`net_edge_betweenness`](https://saqr.me/Nestimate/reference/net_edge_betweenness.md).

- style:

  Plot style. `"bar"` (default) draws one horizontal bar per edge;
  `"forest"` draws a forest/lollipop chart (a stem from zero to a point)
  with a dashed reference line at the mean betweenness.

- top_n:

  Integer or `NULL`. Keep only the `top_n` highest edges. Default `NULL`
  (all edges with non-zero betweenness).

- labels:

  Logical. Print the betweenness value beside each edge. Default `TRUE`.

- ...:

  Additional arguments (ignored).

## Value

A `ggplot` object.

## Examples

``` r
seqs <- data.frame(
  V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
  V3 = c("C","A","C","B"))
eb <- net_edge_betweenness(build_network(seqs, method = "relative"))
plot(eb)

plot(eb, style = "forest")
```
