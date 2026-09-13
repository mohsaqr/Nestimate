# Subtract one network from another

Returns `x - y` as a `netdifference` object: the element-wise difference
of the two weight matrices. Works on any pair of networks; for an
edge-betweenness difference, subtract two
[`net_edge_betweenness`](https://saqr.me/Nestimate/reference/net_edge_betweenness.md)
results. Draw the signed difference network with `cograph::splot(d)` or
`cograph::plot_difference(d)`; cograph handles the colouring and node
palette.

## Usage

``` r
subtract_networks(x, y)

# S3 method for class 'netdifference'
print(x, max_print = 12L, ...)
```

## Arguments

- x, y:

  A `netobject`, `cograph_network`, or numeric square matrix. Both must
  share the same nodes in the same order.

- max_print:

  Integer. Rows to show in
  [`print()`](https://rdrr.io/r/base/print.html). Default `12`.

- ...:

  Ignored.

## Value

A `netdifference` object: a `netobject` whose `$weights` and
`$difference_matrix` are `x - y`, carrying the source matrices `$x` and
`$y`.

## Examples

``` r
early <- data.frame(
  V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
  V3 = c("C","A","C","B"))
late <- data.frame(
  V1 = c("B","A","C","B"), V2 = c("C","B","A","C"),
  V3 = c("A","C","B","A"))
a <- build_network(early, method = "relative")
b <- build_network(late, method = "relative")
subtract_networks(a, b)
#> Network difference (x - y): 3 nodes, 0 differing edges
#> Plot: cograph::splot(d) or cograph::plot_difference(d)
#> 
#>  from to x y difference
#>     C  A 1 1          0
#>     A  B 1 1          0
#>     B  C 1 1          0
# edge-betweenness difference:
subtract_networks(net_edge_betweenness(a), net_edge_betweenness(b))
#> Network difference (x - y): 3 nodes, 0 differing edges
#> Plot: cograph::splot(d) or cograph::plot_difference(d)
#> 
#>  from to x y difference
#>     C  A 3 3          0
#>     A  B 3 3          0
#>     B  C 3 3          0
```
