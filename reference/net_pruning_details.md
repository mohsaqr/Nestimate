# Report Network Pruning Details

`net_pruning_details()` lists the edges that
[`net_prune()`](https://pak.dynasite.org/Nestimate/reference/net_prune.md)
removed from a network, one row per edge, and prints the method with its
setting and the numbers of edges removed and retained above the table.
It reads the pruning record, so it works on a network whose pruning was
undone by
[`net_deprune()`](https://pak.dynasite.org/Nestimate/reference/net_deprune.md)
as well.

## Usage

``` r
net_pruning_details(x, ...)

# S3 method for class 'netobject'
net_pruning_details(x, ...)

# S3 method for class 'netobject_group'
net_pruning_details(x, ...)

# Default S3 method
net_pruning_details(x, ...)

# S3 method for class 'net_pruning_details'
print(x, ...)
```

## Arguments

- x:

  A pruned `netobject` or `netobject_group`; a network that was never
  pruned, or another class, raises an error. For the
  [`print()`](https://rdrr.io/r/base/print.html) method: a
  `net_pruning_details`.

- ...:

  Ignored.

## Value

For a `netobject`: a `net_pruning_details` data frame of the removed
edges, with columns `from`, `to` and `weight` (the weight before
pruning). The method, its setting, the cut-off and the counts are
attributes; [`print()`](https://rdrr.io/r/base/print.html) shows the
method, its setting and the counts. The counts are edges: non-zero cells
of either sign (self-loops included), one per pair on an undirected
network. The retained count is the non-zero edges left after pruning.
For a `netobject_group`: a named list of such data frames, one per
network. The [`print()`](https://rdrr.io/r/base/print.html) method
returns its input invisibly.

## See also

[`net_prune`](https://pak.dynasite.org/Nestimate/reference/net_prune.md)

## Examples

``` r
seqs <- data.frame(
  V1 = c("A","B","A","C","B"), V2 = c("B","C","B","A","C"),
  V3 = c("C","A","C","B","A"))
net <- build_network(seqs, method = "relative")
net_pruning_details(net_prune(net, threshold = 0.2))
#> Pruning details
#>   Method:  user-specified threshold (0.2)
#>   Removed: 0 edges
#>   Retained: 3 edges
#> 
#> (no edges removed)
```
