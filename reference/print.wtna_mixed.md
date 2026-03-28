# Print Method for wtna_mixed

Print Method for wtna_mixed

## Usage

``` r
# S3 method for class 'wtna_mixed'
print(x, ...)
```

## Arguments

- x:

  A `wtna_mixed` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
# \donttest{
oh <- data.frame(
  A = c(1,0,1,0,1,0,1,0),
  B = c(0,1,0,1,0,1,0,1),
  C = c(1,1,0,0,1,1,0,0)
)
mixed <- wtna(oh, method = "both")
print(mixed)
#> Mixed Window TNA (transition + co-occurrence)
#> -- Transition (directed) --
#>   Nodes: 3  |  Edges: 6
#> -- Co-occurrence (undirected) --
#>   Nodes: 3  |  Edges: 2
# }
```
