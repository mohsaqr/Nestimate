# Summary Method for net_link_prediction

Summary Method for net_link_prediction

## Usage

``` r
# S3 method for class 'net_link_prediction'
summary(object, ...)
```

## Arguments

- object:

  A `net_link_prediction` object.

- ...:

  Additional arguments (ignored).

## Value

A data frame, one row per method, with columns `method`,
`n_predictions`, `score_mean`, `score_sd`, `score_max` and `score_min`.
A method with no predictions (every possible link already exists) has
`n_predictions = 0` and `NA` scores.

## Examples

``` r
seqs <- data.frame(
  V1 = c("A", "B", "C", "D", "A", "C", "E", "B"),
  V2 = c("B", "C", "D", "E", "C", "E", "A", "D"),
  V3 = c("C", "D", "E", "A", "D", "A", "B", "E")
)
net <- build_network(seqs, method = "relative")
pred <- predict_links(net)
summary(pred)
#>                    method n_predictions score_mean   score_sd score_max
#> 1        common_neighbors            12 0.11284722 0.18530385 0.4722222
#> 2     resource_allocation            12 0.03477045 0.06055055 0.1574074
#> 3             adamic_adar            12 0.09626854 0.16577482 0.4298352
#> 4                 jaccard            12 0.14555751 0.11315182 0.3251534
#> 5 preferential_attachment            12 2.25000000 1.13818037 4.0000000
#> 6                    katz            12 1.71583160 0.85510259 3.3856655
#>    score_min
#> 1 0.00000000
#> 2 0.00000000
#> 3 0.00000000
#> 4 0.01408451
#> 5 1.00000000
#> 6 0.64215374
```
