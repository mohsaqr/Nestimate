# Re-apply Network Pruning

Re-applies a previously computed pruning that was undone by
[`net_deprune`](https://saqr.me/Nestimate/reference/net_deprune.md),
without recomputation.

## Usage

``` r
net_reprune(x, ...)

# S3 method for class 'netobject'
net_reprune(x, ...)

# S3 method for class 'netobject_group'
net_reprune(x, ...)

# Default S3 method
net_reprune(x, ...)
```

## Arguments

- x:

  A depruned `netobject` or `netobject_group`.

- ...:

  Ignored.

## Value

The network (or group) with pruned weights re-applied and its pruning
marked active.

## See also

[`net_prune`](https://saqr.me/Nestimate/reference/net_prune.md),
[`net_deprune`](https://saqr.me/Nestimate/reference/net_deprune.md)

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"),
                   V3 = c("C","A","B"))
net    <- build_network(seqs, method = "relative")
pruned <- net_prune(net, threshold = 0.2)
undone <- net_deprune(pruned)
net_reprune(undone)
#> Transition Network (relative probabilities) [directed]
#>   Weights: [1.000, 1.000]  |  mean: 1.000
#> 
#>   Weight matrix:
#>     A B C
#>   A 0 1 0
#>   B 0 0 1
#>   C 1 0 0 
#> 
#>   Initial probabilities:
#>   A             0.333  ████████████████████████████████████████
#>   B             0.333  ████████████████████████████████████████
#>   C             0.333  ████████████████████████████████████████
```
