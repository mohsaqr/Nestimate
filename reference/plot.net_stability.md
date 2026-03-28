# Plot Method for net_stability

Plots mean correlation vs drop proportion for each centrality measure.
The CS-coefficient is marked where the curve crosses the threshold.

## Usage

``` r
# S3 method for class 'net_stability'
plot(x, ...)
```

## Arguments

- x:

  A `net_stability` object.

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
stab <- centrality_stability(net, measures = c("InStrength","OutStrength"),
                              iter = 10)
plot(stab)

# }
```
