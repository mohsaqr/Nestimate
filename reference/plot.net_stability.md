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
seqs <- data.frame(
  T1 = c("plan", "code", "debug", "plan", "test", "code"),
  T2 = c("code", "debug", "code", "plan", "code", "test"),
  T3 = c("debug", "code", "plan", "code", "debug", "plan"),
  T4 = c("test", "plan", "test", "debug", "plan", "code")
)
net <- build_network(seqs, method = "relative")
cs <- centrality_stability(net, iter = 10,
  drop_prop = c(0.1, 0.3, 0.5), seed = 1)
plot(cs)

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
