# Print Method for net_link_prediction

Print Method for net_link_prediction

## Usage

``` r
# S3 method for class 'net_link_prediction'
print(x, ...)
```

## Arguments

- x:

  A `net_link_prediction` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
seqs <- data.frame(
  V1 = c("A", "B", "C", "D", "A", "C", "E", "B"),
  V2 = c("B", "C", "D", "E", "C", "E", "A", "D"),
  V3 = c("C", "D", "E", "A", "D", "A", "B", "E")
)
net <- build_network(seqs, method = "relative")
pred <- predict_links(net)
print(pred)
#> Link Prediction  [directed | weighted | 5 nodes | 8 existing edges]
#>   Methods: common_neighbors, resource_allocation, adamic_adar, jaccard, preferential_attachment, katz
#> 
#>   Top predicted links (consensus across 6 methods):
#>     1. C -> B  (avg rank: 4.5, agreed: 6/6)
#>     2. B -> A  (avg rank: 4.8, agreed: 6/6)
#>     3. C -> A  (avg rank: 5.0, agreed: 6/6)
#>     4. D -> C  (avg rank: 5.5, agreed: 6/6)
#>     5. A -> D  (avg rank: 5.8, agreed: 6/6)
#>     6. D -> A  (avg rank: 6.2, agreed: 6/6)
#>     7. E -> D  (avg rank: 6.5, agreed: 6/6)
#>     8. E -> C  (avg rank: 7.5, agreed: 6/6)
#>     9. B -> E  (avg rank: 7.7, agreed: 6/6)
#>     10. A -> E  (avg rank: 8.0, agreed: 6/6)
#>     ... and 2 more predictions
```
