# Print Method for Group Network Object

Compact summary of a `netobject_group`. Header surfaces the source (a
clustering attached by
[`cluster_network`](https://mohsaqr.github.io/Nestimate/reference/cluster_network.md)
or
[`cluster_mmm`](https://mohsaqr.github.io/Nestimate/reference/cluster_mmm.md),
or a plain split by `group_col`). The per-group table carries node and
edge counts, weight range, and – when a clustering attribute is present
– N and percentage of sequences per cluster (matching the layout used by
[`print.net_clustering`](https://mohsaqr.github.io/Nestimate/reference/print.net_clustering.md)
and
[`print.net_mmm`](https://mohsaqr.github.io/Nestimate/reference/print.net_mmm.md)).

## Usage

``` r
# S3 method for class 'netobject_group'
print(x, digits = 3L, ...)
```

## Arguments

- x:

  A `netobject_group`.

- digits:

  Integer. Decimal places for the weight summary. Default `3`.
  Non-breaking: `print(x)` keeps the same shape as before, with the
  addition of a weight-range column.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","A","B"), V2 = c("B","A","B","A"),
                   grp = c("X","X","Y","Y"))
nets <- build_network(seqs, method = "relative", group = "grp")
print(nets)
#> Group Networks (2 groups, group_col: grp)
#> 
#>   Group  Nodes  Edges  Weights
#>   X      2      2      [1.000, 1.000]
#>   Y      2      2      [1.000, 1.000]
# \donttest{
seqs <- data.frame(
  V1 = c("A","B","A","C","B","A"),
  V2 = c("B","C","B","A","C","B"),
  V3 = c("C","A","C","B","A","C"),
  grp = c("X","X","X","Y","Y","Y")
)
nets <- build_network(seqs, method = "relative", group = "grp")
print(nets)
#> Group Networks (2 groups, group_col: grp)
#> 
#>   Group  Nodes  Edges  Weights
#>   X      3      3      [1.000, 1.000]
#>   Y      3      3      [1.000, 1.000]
# }
```
