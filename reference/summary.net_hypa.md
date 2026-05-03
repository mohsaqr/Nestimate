# Summary Method for net_hypa

Summary Method for net_hypa

## Usage

``` r
# S3 method for class 'net_hypa'
summary(object, n = 10L, type = c("all", "over", "under"), ...)
```

## Arguments

- object:

  A `net_hypa` object.

- n:

  Integer. Maximum number of paths to display per category (default:
  10).

- type:

  Character. Which anomalies to show: `"all"` (default), `"over"`, or
  `"under"`.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
seqs <- list(c("A","B","C"), c("B","C","A"), c("A","C","B"), c("A","B","C"))
hyp <- build_hypa(seqs, k = 2)
summary(hyp)
#> HYPA Summary
#> 
#>   Order: 2 | Nodes: 5 | Edges: 3
#>   Alpha: 0.05 | p_adjust: BH
#>   Anomalous: 0 (over: 0, under: 0)
#> 
#>   No anomalous paths detected.
#> [1] path      observed  expected  ratio     p_value   direction
#> <0 rows> (or 0-length row.names)

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
#>   Alpha: 0.05 | p_adjust: BH
#>   Anomalous: 0 (over: 0, under: 0)
#> 
#>   No anomalous paths detected.
#> [1] path      observed  expected  ratio     p_value   direction
#> <0 rows> (or 0-length row.names)
summary(hypa, type = "over", n = 5)
#> HYPA Summary
#> 
#>   Order: 2 | Nodes: 3 | Edges: 3
#>   Alpha: 0.05 | p_adjust: BH
#>   Anomalous: 0 (over: 0, under: 0)
#> 
#>   No anomalous paths detected.
#> [1] path      observed  expected  ratio     p_value   direction
#> <0 rows> (or 0-length row.names)
# }
```
