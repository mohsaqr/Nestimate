# Transductive label spreading on a hypergraph

`hypergraph_transduction()` classifies the unlabelled nodes of a
hypergraph from a few labelled ones, by the label-spreading framework of
Zhou et al. (2006). The labels spread over the hypergraph as the scores
\$\$F = (1 - \xi) (I - \xi S)^{-1} Y,\$\$ where \\S = I - L\\ is the
similarity operator of the chosen Laplacian (see
[`hypergraph_laplacian()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_laplacian.md))
and \\Y\\ marks each labelled node's class, and each node is assigned
the class with the highest score. This is the non-neural ancestor of
hypergraph-attention text classifiers: with documents as hyperedges over
words (or the reverse), it classifies unlabelled nodes from a handful of
labelled ones.

The `margin` of a node, its winning minus its runner-up score, says how
clearly it was assigned; assignments with a small margin are the ones to
distrust. The given labels are not held fixed, so a labelled node can be
predicted into another class, and a class's scores grow with its number
of labelled nodes. As `xi` approaches 1, every node's scores become
proportional to the same vector, so in the limit every node is assigned
the same class. The default `0.99` is the value of Zhou et al. (2006)
and suits large hypergraphs; on a small one it over-smooths. In the
example below it assigns five of the six nodes to one class, where
`xi = 0.9` splits them three and three, so on a small hypergraph pass a
smaller `xi` and compare.

## Usage

``` r
hypergraph_transduction(
  hg,
  labels,
  xi = 0.99,
  type = c("zhou", "random_walk"),
  edge_weights = NULL
)

# S3 method for class 'net_hypergraph_transduction'
print(x, ...)

# S3 method for class 'net_hypergraph_transduction'
summary(object, ...)

label_predictions(x)

label_scores(x)

# S3 method for class 'net_hypergraph_transduction'
plot(x, ...)
```

## Arguments

- hg:

  The hypergraph: a connected `net_hypergraph`, as for
  [`hypergraph_laplacian()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_laplacian.md).

- labels:

  The known class labels, in one of two forms: a vector named by node,
  covering some of the nodes, or an unnamed vector with one entry per
  node, in the hypergraph's node order, and `NA` for the unlabelled
  ones. Required. It must hold at least two distinct classes; a name
  that is not a node, or an unnamed vector of the wrong length, raises
  an error.

- xi:

  The spreading coefficient, a number strictly between 0 and 1. Larger
  values give the hypergraph structure more weight relative to the given
  labels. Default `0.99`.

- type, edge_weights:

  The Laplacian and its hyperedge weights, passed to
  [`hypergraph_laplacian()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_laplacian.md).
  Defaults `"zhou"` and `NULL`.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods and
  for `label_predictions()` and `label_scores()`: a
  `net_hypergraph_transduction` object.

- ...:

  For the methods: ignored.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `net_hypergraph_transduction` object.

## Value

A `net_hypergraph_transduction` object.
[`print()`](https://rdrr.io/r/base/print.html) reports the Laplacian,
`xi`, the classes and how many nodes each was predicted for;
[`summary()`](https://rdrr.io/r/base/summary.html) returns one row per
class; [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws
the node-by-class score heatmap. The object's fields:

- predictions:

  A data frame with one row per node: `node`, `label` (the given label,
  `NA` for an unlabelled node), `predicted`, `score` (the winning
  class's score) and `margin` (the winning minus the runner-up score).

- classes:

  The classes, sorted.

- scores:

  The node-by-class score matrix.

- xi, type:

  The spreading coefficient and the Laplacian type.

- n_labeled, n_nodes:

  The numbers of labelled and of all nodes.

- params:

  The `edge_weights` used.

The [`print()`](https://rdrr.io/r/base/print.html) method returns its
input invisibly. The [`summary()`](https://rdrr.io/r/base/summary.html)
method returns a data frame with one row per class: `class`,
`n_labeled`, `n_predicted` and `mean_margin` (the mean margin of the
nodes predicted into the class, `NA` when there are none).
`label_predictions()` returns a data frame with one row per node:
`node`, `label` (`NA` for an unlabelled node), `predicted`, `score` and
`margin`. `label_scores()` returns a data frame with one row per node
and class: `node`, `class` and `score`. The
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method returns
the score heatmap as a `ggplot`, visibly, so `plot(x)` draws it: one row
per node, grouped by predicted class and ordered by margin within it,
one column per class, the tiles shaded and labelled by score. Labelled
nodes have a black border on their given class, and a dot marks each
node's winning class, so seeds, scores and assignments can be compared
in one panel.

## References

Zhou, D., Huang, J., & Schölkopf, B. (2006). Learning with hypergraphs:
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
tr <- hypergraph_transduction(hg, labels = c(a = "x", d = "y"), xi = 0.9)
tr
#> Hypergraph transductive label spreading (zhou Laplacian, xi = 0.9)
#>   Nodes: 6 (2 labeled) | Classes: x, y
#>   Predicted: x = 3, y = 3
label_predictions(tr)
#>   node label predicted     score     margin
#> 1    a     x         x 0.3062500 0.21439413
#> 2    b  <NA>         x 0.2062500 0.11439413
#> 3    c  <NA>         x 0.2143304 0.06433035
#> 4    d     y         y 0.3500000 0.25814413
#> 5    e  <NA>         y 0.2143304 0.15808035
#> 6    f  <NA>         y 0.2143304 0.15808035
label_scores(tr)
#>    node class      score
#> 1     a     x 0.30625000
#> 2     b     x 0.20625000
#> 3     c     x 0.21433035
#> 4     d     x 0.09185587
#> 5     e     x 0.05625000
#> 6     f     x 0.05625000
#> 7     a     y 0.09185587
#> 8     b     y 0.09185587
#> 9     c     y 0.15000000
#> 10    d     y 0.35000000
#> 11    e     y 0.21433035
#> 12    f     y 0.21433035
```
