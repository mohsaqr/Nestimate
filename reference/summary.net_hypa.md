# Summary Method for net_hypa

Summary Method for net_hypa

## Usage

``` r
# S3 method for class 'net_hypa'
summary(
  object,
  n = 10L,
  type = c("all", "over", "under"),
  order_by = c("sig", "freq", "frequency", "ratio", "path"),
  ...
)
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

- order_by:

  Character. Ranking used within each anomaly direction: `"sig"`
  (default) ranks by the active tail probability, `"freq"` (or its alias
  `"frequency"`) by observed count, `"ratio"` by observed/expected
  ratio, and `"path"` alphabetically.

- ...:

  Additional arguments (ignored).

## Value

A data frame of the reported anomalies, at most `n` rows per direction,
with columns `order` (the De Bruijn order the path was found at),
`path`, `observed`, `expected`, `ratio`, `p_tail` (the raw tail
probability in the reported direction: `p_over` for over-represented
paths, `p_under` for under-represented ones) and `direction`
(`"over"`/`"under"`). Returned visibly; the summary text and the top-`n`
tables are printed as a side effect. When no anomalies were detected, a
zero-row data frame with the same columns except `order` is returned.

## Examples

``` r
seqs <- list(c("A","B","C"), c("B","C","A"), c("A","C","B"), c("A","B","C"))
hyp <- build_hypa(seqs, order = 2)
summary(hyp)
#> HYPA Summary
#> 
#>   Order(s): 2 | Edges: 3
#>   Alpha: 0.05 | p_adjust: BH
#>   Anomalous: 0 (over: 0, under: 0) | order_by: sig
#> 
#>   No anomalous paths detected.
#> [1] path      observed  expected  ratio     p_tail    direction
#> <0 rows> (or 0-length row.names)

# \donttest{
seqs <- data.frame(
  V1 = c("A","B","C","A","B","C","A","B","C","A"),
  V2 = c("B","C","A","B","C","A","B","C","A","B"),
  V3 = c("C","A","B","C","A","B","C","A","B","C"),
  V4 = c("A","B","C","A","B","C","A","B","C","A")
)
hypa <- build_hypa(seqs, order = 2L)
summary(hypa)
#> HYPA Summary
#> 
#>   Order(s): 2 | Edges: 3
#>   Alpha: 0.05 | p_adjust: BH
#>   Anomalous: 0 (over: 0, under: 0) | order_by: sig
#> 
#>   No anomalous paths detected.
#> [1] path      observed  expected  ratio     p_tail    direction
#> <0 rows> (or 0-length row.names)
summary(hypa, type = "over", n = 5)
#> HYPA Summary
#> 
#>   Order(s): 2 | Edges: 3
#>   Alpha: 0.05 | p_adjust: BH
#>   Anomalous: 0 (over: 0, under: 0) | order_by: sig
#> 
#>   No anomalous paths detected.
#> [1] path      observed  expected  ratio     p_tail    direction
#> <0 rows> (or 0-length row.names)
# }
```
