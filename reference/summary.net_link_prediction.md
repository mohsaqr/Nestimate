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

A data frame with per-method summary statistics, invisibly.

## Examples

``` r
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 30, TRUE),
  V2 = sample(LETTERS[1:4], 30, TRUE),
  V3 = sample(LETTERS[1:4], 30, TRUE)
)
net <- build_network(seqs, method = "relative")
pred <- predict_links(net)
summary(pred)
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to min; returning Inf
#>                    method n_predictions score_mean score_sd score_max score_min
#> 1        common_neighbors             0        NaN       NA      -Inf       Inf
#> 2     resource_allocation             0        NaN       NA      -Inf       Inf
#> 3             adamic_adar             0        NaN       NA      -Inf       Inf
#> 4                 jaccard             0        NaN       NA      -Inf       Inf
#> 5 preferential_attachment             0        NaN       NA      -Inf       Inf
#> 6                    katz             0        NaN       NA      -Inf       Inf
```
