# Compute Centrality Measures for a Network

Computes centrality measures from a `netobject`, `netobject_group`,
`mcml`, or `cograph_network`. The built-in measures match
[`tna::centralities()`](http://sonsoles.me/tna/reference/centralities.md)
without importing `tna` or `igraph`. The only intentional default
difference is that `Diffusion` is range-normalized by default.

## Usage

``` r
net_centrality(
  x,
  measures = NULL,
  loops = FALSE,
  normalize = FALSE,
  invert = TRUE,
  normalize_diffusion = TRUE,
  centrality_fn = NULL,
  ...
)
```

## Arguments

- x:

  A `netobject`, `netobject_group`, or `cograph_network`.

- measures:

  Character vector. Centrality measures to compute. Defaults to
  `c("InStrength", "Betweenness", "Diffusion")`. Pass `"all"` for every
  built-in measure: `"OutStrength"`, `"InStrength"`, `"ClosenessIn"`,
  `"ClosenessOut"`, `"Closeness"`, `"Betweenness"`, `"BetweennessRSP"`,
  `"Diffusion"`, and `"Clustering"`. The legacy aliases `"InCloseness"`
  and `"OutCloseness"` are also accepted.

- loops:

  Logical. Include self-loops (diagonal) in computation? Default:
  `FALSE`.

- normalize:

  Logical. Range-normalize all requested measures using the same
  transformation as `tna::centralities(normalize = TRUE)`. Default:
  `FALSE`.

- invert:

  Logical. Invert weights for shortest-path measures? Default: `TRUE`,
  matching `tna`.

- normalize_diffusion:

  Logical. Range-normalize `Diffusion` even when `normalize = FALSE`.
  Default: `TRUE`.

- centrality_fn:

  Optional function. Custom centrality function that takes a weight
  matrix and returns a named list of centrality vectors.

- ...:

  Additional arguments (ignored).

## Value

For a `netobject`: a `net_centrality` data frame with node names as
rows, a `state` column, and one column per centrality measure. For a
`netobject_group`: a `net_centrality_group` list of such data frames.

## Examples

``` r
seqs <- data.frame(
  V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
  V3 = c("C","A","C","B"))
net <- build_network(seqs, method = "relative")
net_centrality(net)
#> centralities computed excluding loops (diagonal). Pass `loops = TRUE` to include self-transitions.
#>   state InStrength Betweenness Diffusion
#> A     A          1           1       NaN
#> B     B          1           1       NaN
#> C     C          1           1       NaN
```
