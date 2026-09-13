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
set.seed(42)
seqs <- data.frame(
  V1 = sample(LETTERS[1:5], 50, TRUE),
  V2 = sample(LETTERS[1:5], 50, TRUE),
  V3 = sample(LETTERS[1:5], 50, TRUE)
)
net <- build_network(seqs, method = "relative")
pred <- predict_links(net, exclude_existing = FALSE)

# Evaluate against the network's own edges as the known truth
evaluate_links(pred, extract_edges(net, threshold = 0.001))
#>                    method       auc average_precision precision_at_5
#> 1        common_neighbors 0.3421053         0.9536904              1
#> 2     resource_allocation 0.3421053         0.9536904              1
#> 3             adamic_adar 0.3421053         0.9536904              1
#> 4                 jaccard 0.6578947         0.9817801              1
#> 5 preferential_attachment 1.0000000         1.0000000              1
#> 6                    katz 1.0000000         1.0000000              1
#>   precision_at_10 precision_at_20
#> 1             0.9            0.95
#> 2             0.9            0.95
#> 3             0.9            0.95
#> 4             1.0            0.95
#> 5             1.0            0.95
#> 6             1.0            0.95
```
