# Summary Method for net_reliability

Summary Method for net_reliability

## Usage

``` r
# S3 method for class 'net_reliability'
summary(object, ...)
```

## Arguments

- object:

  A `net_reliability` object.

- ...:

  Ignored.

## Value

A tidy data frame with columns `model`, `metric`, `mean`, `sd`
summarising the split-half iterations.

## Examples

``` r
net <- build_network(data.frame(V1 = c("A","B","C","A"),
  V2 = c("B","C","A","B")), method = "relative")
rel <- network_reliability(net, iter = 10)
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
  V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
)
net <- build_network(seqs, method = "relative")
rel <- network_reliability(net, iter = 100, seed = 42)
print(rel)
#> Split-Half Reliability (100 iterations, split = 50%)
#>   Mean Abs. Diff.     mean = 0.1498  sd = 0.0326
#>   Median Abs. Diff.   mean = 0.1315  sd = 0.0380
#>   Pearson             mean = -0.3774  sd = 0.1721
#>   Max Abs. Diff.      mean = 0.3721  sd = 0.0810
# }
```
