# Transductive label spreading on a hypergraph

Semi-supervised classification of hypergraph nodes by the regularization
framework of Zhou et al. (2006): given labels for a subset of nodes, the
scores `F = (1 - xi) * (I - xi * S)^{-1} Y` spread the labels over the
hypergraph, where `S = I - L` is the normalized similarity operator of
the chosen Laplacian and `Y` is the label indicator matrix. Each node is
assigned the class with the highest score. This is the non-neural
ancestor of hypergraph-attention text classifiers: with documents as
hyperedges over words (or vice versa) it classifies unlabeled nodes from
a handful of labeled ones.

## Usage

``` r
hypergraph_transduction(
  hg,
  labels,
  xi = 0.99,
  type = c("zhou", "random_walk"),
  edge_weights = NULL
)
```

## Arguments

- hg:

  A connected `net_hypergraph`.

- labels:

  Node labels. Either a named vector (names = node names, values = class
  labels) covering a subset of nodes, or a full-length vector aligned
  with `hg$nodes` with `NA` for unlabeled nodes. At least two distinct
  classes must be labeled.

- xi:

  Numeric in `(0, 1)`. Spreading coefficient (default `0.99`); larger
  values weight the hypergraph structure more relative to the initial
  labels.

- type, edge_weights:

  Passed to
  [`hypergraph_laplacian()`](https://saqr.me/Nestimate/reference/hypergraph_laplacian.md).

## Value

An object of class `net_hypergraph_transduction`: a list with
`$predictions` (data.frame, one row per node: `node`, `label` (given,
`NA` if unlabeled), `predicted`, `score` (winning class score), `margin`
(winning minus runner-up score)), `$classes`, `$scores` (node x class
score matrix), `$xi`, `$type`, `$n_labeled`, `$n_nodes` and `$params`
(the `edge_weights` used). Has `print`, `summary`, `plot` and
`as.data.frame` methods; `as.data.frame(x, what = "scores")` returns the
tidy long score table.

## References

Zhou, D., Huang, J., & Scholkopf, B. (2006). Learning with hypergraphs:
Clustering, classification, and embedding. *NeurIPS 19*.

## Examples

``` r
events <- data.frame(
  person = c("a", "b", "c", "a", "b", "c", "d", "e", "f",
             "d", "e", "f", "c", "d"),
  meeting = c("m1", "m1", "m1", "m2", "m2", "m2", "m3", "m3", "m3",
              "m4", "m4", "m4", "m5", "m5")
)
hg <- bipartite_groups(events, player = "person", group = "meeting")
tr <- hypergraph_transduction(hg, labels = c(a = "x", d = "y"))
tr
#> Hypergraph transductive label spreading (zhou Laplacian, xi = 0.99)
#>   Nodes: 6 (2 labeled) | Classes: x, y
#>   Predicted: x = 1, y = 5
as.data.frame(tr)
#>   node label predicted     score      margin
#> 1    a     x         x 0.1654808 0.003987167
#> 2    b  <NA>         y 0.1614936 0.006012833
#> 3    c  <NA>         y 0.2037821 0.019834793
#> 4    d     y         y 0.2321154 0.070621782
#> 5    e  <NA>         y 0.1839473 0.055966489
#> 6    f  <NA>         y 0.1839473 0.055966489
```
