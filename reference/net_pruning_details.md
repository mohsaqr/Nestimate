# Report Network Pruning Details

Returns the edges removed by
[`net_prune`](https://pak.dynasite.org/Nestimate/reference/net_prune.md)
as a tidy one-row-per-edge data frame, with the method, cut-off, and
retained/removed counts attached as attributes and shown by its print
method.

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

  A pruned `netobject` or `netobject_group`. For the
  `net_pruning_details()` and
  [`print()`](https://rdrr.io/r/base/print.html) methods: an object of
  class `netobject`, `netobject_group`, `default` or
  `net_pruning_details`.

- ...:

  Ignored. For the S3 methods: further arguments passed to or from other
  methods.

## Value

For a `netobject`: a `net_pruning_details` data frame (columns `from`,
`to`, `weight`) of removed edges. For a `netobject_group`: a named list
of such data frames.

In `print.net_pruning_details()`: `x`, invisibly.

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
