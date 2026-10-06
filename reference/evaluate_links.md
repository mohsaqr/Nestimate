# Evaluate Link Predictions Against Known Edges

Computes AUC-ROC, precision\\@\\k, and average precision for link
predictions against a set of known true edges.

## Usage

``` r
evaluate_links(pred, true_edges, k = c(5L, 10L, 20L))
```

## Arguments

- pred:

  A `net_link_prediction` object.

- true_edges:

  A data frame with columns `from` and `to`, or a binary matrix where 1
  indicates a true edge.

- k:

  Integer vector. Values of k for precision\\@\\k. Default:
  `c(5, 10, 20)`.

## Value

A data frame with columns: method, auc, average_precision, and one
precision_at_k column per k value.

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
#> 4                 jaccard 0.3333333         0.3240105            0.0
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
