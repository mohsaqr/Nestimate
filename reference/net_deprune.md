# Undo Network Pruning

Restores the original (pre-pruning) weights of a network pruned by
[`net_prune`](https://saqr.me/Nestimate/reference/net_prune.md), without
recomputation. The pruning record is kept, so
[`net_reprune`](https://saqr.me/Nestimate/reference/net_reprune.md) can
re-apply it.

## Usage

``` r
net_deprune(x, ...)

# S3 method for class 'netobject'
net_deprune(x, ...)

# S3 method for class 'netobject_group'
net_deprune(x, ...)

# Default S3 method
net_deprune(x, ...)
```

## Arguments

- x:

  A pruned `netobject` or `netobject_group`.

- ...:

  Ignored.

## Value

The network (or group) with original weights restored and its pruning
marked inactive.

## See also

[`net_prune`](https://saqr.me/Nestimate/reference/net_prune.md),
[`net_reprune`](https://saqr.me/Nestimate/reference/net_reprune.md)

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"),
                   V3 = c("C","A","B"))
net    <- build_network(seqs, method = "relative")
pruned <- net_prune(net, threshold = 0.2)
net_deprune(pruned)
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
