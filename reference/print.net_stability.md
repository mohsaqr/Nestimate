# Print Method for net_stability

Print Method for net_stability

## Usage

``` r
# S3 method for class 'net_stability'
print(x, ...)
```

## Arguments

- x:

  A `net_stability` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
seqs <- data.frame(
  T1 = c("plan", "code", "debug", "plan", "test", "code"),
  T2 = c("code", "debug", "code", "plan", "code", "test"),
  T3 = c("debug", "code", "plan", "code", "debug", "plan"),
  T4 = c("test", "plan", "test", "debug", "plan", "code")
)
net <- build_network(seqs, method = "relative")
cs <- centrality_stability(net, iter = 10, drop_prop = 0.3, seed = 1)
print(cs)
#> Centrality Stability (10 iterations, threshold = 0.7)
#>   Drop proportions: 0.3
#> 
#>   CS-coefficients:
#>     InStrength       0.30
#>     OutStrength      0.30
#>     Betweenness      0.00
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
print(stab)
#> Centrality Stability (10 iterations, threshold = 0.7)
#>   Drop proportions: 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9
#> 
#>   CS-coefficients:
#>     InStrength       0.10
#>     OutStrength      0.10
# }
```
