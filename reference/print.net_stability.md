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
net <- build_network(data.frame(V1 = c("A","B","C","A"),
  V2 = c("B","C","A","B")), method = "relative")
cs <- centrality_stability(net, iter = 10, drop_prop = 0.3)
#> Warning: All centrality measures have zero variance. No stability can be assessed.
print(cs)
#> Centrality Stability (10 iterations, threshold = 0.7)
#>   Drop proportions: 0.3
#> 
#>   CS-coefficients:
#>     InStrength       0.00
#>     OutStrength      0.00
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
