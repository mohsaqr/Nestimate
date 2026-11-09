# Subtract one network from another

`subtract_networks()` subtracts one network's weights from another's,
cell by cell, and returns the difference `x - y` as a `netdifference`
object: a network whose edge weights are the signed differences. It
takes any pair of networks on the same nodes (`netobject`s,
`cograph_network`s or square numeric matrices); for an edge-betweenness
difference, subtract two
[`net_edge_betweenness`](https://pak.dynasite.org/Nestimate/reference/net_edge_betweenness.md)
results. `cograph::splot(d)` or `cograph::plot_difference(d)` draws the
signed difference network, with cograph's colouring and node palette.

The difference is taken between two point estimates and carries no
uncertainty, so a large difference on a rarely observed edge can be
noise. To see which edges differ beyond sampling variation, compare the
two networks with
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
or
[`bayes_compare()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md);
[`as_netdifference()`](https://pak.dynasite.org/Nestimate/reference/as_netdifference.md)
turns a
[`bayes_compare()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
result into the same kind of object, keeping, by default, only the
credible differences.

## Usage

``` r
subtract_networks(x, y)

# S3 method for class 'netdifference'
print(x, max_print = 12L, ...)
```

## Arguments

- x, y:

  The two networks: `netobject`s, `cograph_network`s or square numeric
  matrices. They must have the same number of nodes, with the same
  labels in the same order, or an error is raised; any other kind of
  object raises an error too. A matrix without dimnames is labelled
  `V1`, `V2`, ..., so it matches only a network with those labels. The
  difference is directed when `x` is a directed `netobject`; when `x` is
  not a `netobject`, it is directed when `x`'s weights are not
  symmetric. For the [`print()`](https://rdrr.io/r/base/print.html)
  method: `x` is the `netdifference` object.

- max_print:

  For the [`print()`](https://rdrr.io/r/base/print.html) method: the
  most edges to list; the number left out is reported. Default `12`.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html) method:
  ignored.

## Value

A `netdifference` object, which is also a `netobject` (method
`"difference"`). [`print()`](https://rdrr.io/r/base/print.html) lists
the edges present in either network, largest absolute difference first,
with both weights and their difference (for an undirected difference,
each pair once);
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
draws it. Beside the usual `netobject` fields, its `weights` and
`difference_matrix` both hold `x - y`, and `x` and `y` hold the two
weight matrices that were subtracted.

The [`print()`](https://rdrr.io/r/base/print.html) method returns its
input invisibly.

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
