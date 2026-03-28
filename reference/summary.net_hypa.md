# Summary Method for net_hypa

Summary Method for net_hypa

## Usage

``` r
# S3 method for class 'net_hypa'
summary(object, ...)
```

## Arguments

- object:

  A `net_hypa` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
# \donttest{
seqs <- data.frame(
  V1 = c("A","B","C","A","B","C","A","B","C","A"),
  V2 = c("B","C","A","B","C","A","B","C","A","B"),
  V3 = c("C","A","B","C","A","B","C","A","B","C"),
  V4 = c("A","B","C","A","B","C","A","B","C","A")
)
hypa <- build_hypa(seqs, k = 2L)
summary(hypa)
#> HYPA Summary
#> 
#>   Order: 2 | Nodes: 3 | Edges: 3
#>   Alpha: 0.05
#> 
#>   No anomalous paths detected.
# }
```
