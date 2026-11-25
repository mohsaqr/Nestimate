# Re-apply the pruning of a network

`net_reprune()` re-applies a pruning that
[`net_deprune()`](https://pak.dynasite.org/Nestimate/reference/net_deprune.md)
undid, without recomputing it: the pruned weights stored by
[`net_prune()`](https://pak.dynasite.org/Nestimate/reference/net_prune.md)
come back and the pruning is marked active again. A `netobject_group` is
handled network by network.

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

  A depruned `netobject` or `netobject_group`. A network that was never
  pruned, or whose pruning is already active, raises an error, and so
  does another class.

- ...:

  Ignored.

## Value

The input, of the same class, with the pruned weights, its edge table
and edge count updated to them, and its pruning marked active.

## See also

[`net_prune`](https://pak.dynasite.org/Nestimate/reference/net_prune.md),
[`net_deprune`](https://pak.dynasite.org/Nestimate/reference/net_deprune.md)

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
#>   Strongest edges (3 of 3):
#>    from to weight
#>    A    B  1     
#>    B    C  1     
#>    C    A  1      
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix;
#>   extract_initial_probs() the initial probabilities.
```
