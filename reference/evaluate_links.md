# Evaluate Link Predictions Against Known Edges

`evaluate_links()` measures how well each index of a
[`predict_links()`](https://pak.dynasite.org/Nestimate/reference/predict_links.md)
result recovers a set of known edges. For each index it returns the area
under the ROC curve (AUC), the average precision, and the precision
among the top `k` predictions.

The AUC is the probability that a known edge scores higher than a pair
that is not one (ties count half), so 0.5 is chance level. Average
precision is the mean of the precisions at the rank of each known edge,
and precision\\@\\k the share of known edges among the `k`
highest-scoring pairs. Only the pairs listed in the predictions are
evaluated: with `exclude_existing = TRUE`, the default of
[`predict_links()`](https://pak.dynasite.org/Nestimate/reference/predict_links.md),
the network's current edges are not scored, so the known edges should be
ones held out of the network (a later period, a removed subset), and
`top_n` limits the evaluation to the pairs it kept. The AUC is `NA` when
the evaluated pairs include no known edge, or nothing but known edges.

## Usage

``` r
evaluate_links(pred, true_edges, k = c(5L, 10L, 20L))
```

## Arguments

- pred:

  A
  [`predict_links()`](https://pak.dynasite.org/Nestimate/reference/predict_links.md)
  result.

- true_edges:

  The known edges: a data frame with columns `from` and `to` (every row
  is an edge, whatever its other columns, such as a weight, hold), or a
  matrix whose non-zero cells are edges. Its row names name both its
  rows and its columns (column names are not read), so its columns must
  be in the same order as its rows; without row names, the nodes of
  `pred` are used. For a directed prediction, pairs are matched by
  `from` and `to`, direction included; for an undirected one, in either
  order, so `(C, A)` matches the listed pair `(A, C)`.

- k:

  The cut-offs for precision\\@\\k. Default `c(5, 10, 20)`. A `k` larger
  than an index's number of predictions gives `NA`.

## Value

A data frame with one row per index: `method`, `auc`,
`average_precision` (`0` when no evaluated pair is a known edge), and
one `precision_at_<k>` column per value of `k`. A prediction with no
listed pairs (every possible link already exists) gives `NULL`.

## Examples

``` r
seqs <- data.frame(
  V1 = c("A", "B", "C", "D", "A", "C", "E", "B"),
  V2 = c("B", "C", "D", "E", "C", "E", "A", "D"),
  V3 = c("C", "D", "E", "A", "D", "A", "B", "E")
)
net <- build_network(seqs, method = "relative")
pred <- predict_links(net, exclude_existing = FALSE)

# Evaluate against the network's own edges as the known truth
evaluate_links(pred, extract_edges(net, threshold = 0.001))
#>                    method       auc average_precision precision_at_5
#> 1        common_neighbors 0.5833333         0.4618602            0.4
#> 2     resource_allocation 0.5833333         0.4618602            0.4
#> 3             adamic_adar 0.5833333         0.4618602            0.4
#> 4                 jaccard 0.3333333         0.3242384            0.0
#> 5 preferential_attachment 0.6875000         0.7002999            0.8
#> 6                    katz 0.7916667         0.6196293            0.4
#>   precision_at_10 precision_at_20
#> 1             0.4             0.4
#> 2             0.4             0.4
#> 3             0.4             0.4
#> 4             0.3             0.4
#> 5             0.5             0.4
#> 6             0.7             0.4
```
