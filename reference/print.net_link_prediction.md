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
  V1 = sample(LETTERS[1:4], 30, TRUE),
  V2 = sample(LETTERS[1:4], 30, TRUE),
  V3 = sample(LETTERS[1:4], 30, TRUE)
)
net <- build_network(seqs, method = "relative")
pred <- predict_links(net)
print(pred)
#> Link Prediction  [directed | weighted | 4 nodes | 12 existing edges]
#>   Methods: common_neighbors, resource_allocation, adamic_adar, jaccard, preferential_attachment, katz
```
