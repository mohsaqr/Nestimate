# Undo Network Pruning

`net_deprune()` restores the weights a network had before
[`net_prune()`](https://pak.dynasite.org/Nestimate/reference/net_prune.md),
without recomputing anything, and marks its pruning inactive. The
pruning record stays on the network, so
[`net_reprune()`](https://pak.dynasite.org/Nestimate/reference/net_reprune.md)
can re-apply it and
[`net_pruning_details()`](https://pak.dynasite.org/Nestimate/reference/net_pruning_details.md)
still lists what the pruning removed. A `netobject_group` is restored
network by network.

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

  A pruned `netobject` or `netobject_group`. A network that was never
  pruned, or whose pruning is already inactive, raises an error, and so
  does another class.

- ...:

  Ignored.

## Value

The input, of the same class, with the original weights, its edge table
and edge count updated to them, and its pruning marked inactive.

## See also

[`net_prune`](https://pak.dynasite.org/Nestimate/reference/net_prune.md),
[`net_reprune`](https://pak.dynasite.org/Nestimate/reference/net_reprune.md)

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
