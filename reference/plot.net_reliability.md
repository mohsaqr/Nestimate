# Plot Method for net_reliability

Density plots of split-half metrics faceted by metric type. Multi-model
comparisons show overlaid densities colored by model.

## Usage

``` r
# S3 method for class 'net_reliability'
plot(x, ...)
```

## Arguments

- x:

  A `net_reliability` object.

- ...:

  Additional arguments (ignored).

## Value

A `ggplot` object (invisibly).

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
plot(rel)

# }
```
