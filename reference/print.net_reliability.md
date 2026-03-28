# Print Method for net_reliability

Print Method for net_reliability

## Usage

``` r
# S3 method for class 'net_reliability'
print(x, ...)
```

## Arguments

- x:

  A `net_reliability` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(c("A","B","C"), 30, TRUE),
  V2 = sample(c("A","B","C"), 30, TRUE),
  V3 = sample(c("A","B","C"), 30, TRUE)
)
net <- build_network(seqs, method = "relative")
rel <- reliability(net, iter = 20, seed = 1)
print(rel)
#> Split-Half Reliability (20 iterations, split = 50%)
#>   Mean Abs. Dev.      mean = 0.1710  sd = 0.0498
#>   Median Abs. Dev.    mean = 0.1503  sd = 0.0679
#>   Correlation         mean = 0.0523  sd = 0.2687
#>   Max Abs. Dev.       mean = 0.3745  sd = 0.1078
# }
```
