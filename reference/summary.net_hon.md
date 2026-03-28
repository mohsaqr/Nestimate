# Summary Method for net_hon

Summary Method for net_hon

## Usage

``` r
# S3 method for class 'net_hon'
summary(object, ...)
```

## Arguments

- object:

  A `net_hon` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
# \donttest{
seqs <- data.frame(
  V1 = c("A","B","C","A","B"),
  V2 = c("B","C","A","B","C"),
  V3 = c("C","A","B","C","A")
)
hon <- build_hon(seqs, max_order = 2L)
summary(hon)
#> Higher-Order Network (HON) Summary
#>   Nodes: 3 | Edges: 3 | Trajectories: 5
#>   First-order states: A, B, C
#>   Max order observed: 1 (requested: 2)
#>   Min frequency: 1
#>   Node order distribution:
#>     Order 1: 3 nodes
# }
```
