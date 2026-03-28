# Print Method for Group Network Object

Print Method for Group Network Object

## Usage

``` r
# S3 method for class 'netobject_group'
print(x, ...)
```

## Arguments

- x:

  A `netobject_group`.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
# \donttest{
seqs <- data.frame(
  V1 = c("A","B","A","C","B","A"),
  V2 = c("B","C","B","A","C","B"),
  V3 = c("C","A","C","B","A","C"),
  grp = c("X","X","X","Y","Y","Y")
)
nets <- build_network(seqs, method = "relative", group = "grp")
print(nets)
#> Group Networks (2 groups)
#>   X: 3 nodes, 3 edges
#>   Y: 3 nodes, 3 edges
# }
```
