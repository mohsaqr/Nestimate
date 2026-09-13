# Compute Centrality Measures for a Network

Computes centrality measures from a `netobject`, `netobject_group`,
`mcml`, or `cograph_network`. The built-in measures match
[`tna::centralities()`](http://sonsoles.me/tna/reference/centralities.md)
without importing `tna` or `igraph`: strength is taken from the weight
matrix directly, and the path-based measures (betweenness, closeness)
come from all-pairs shortest paths computed in-package by
Floyd-Warshall. The only intentional default difference from `tna` is
that `Diffusion` is range-normalized by default.

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

  A `netobject`, `netobject_group`, `mcml`, or `cograph_network`.

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

For a `netobject` or `cograph_network`: a `net_centrality` data frame,
one row per node, with a `state` column and one further column per
requested measure (node names are also the row names). For a
`netobject_group` or an `mcml`: a `net_centrality_group` list of such
data frames, one per group.

## References

Freeman, L. C. (1978). Centrality in social networks: conceptual
clarification. *Social Networks*, 1(3), 215–239. (betweenness,
closeness)

Opsahl, T., Agneessens, F. & Skvoretz, J. (2010). Node centrality in
weighted networks: generalizing degree and shortest paths. *Social
Networks*, 32(3), 245–251. (weighted strength and geodesics)

Kivimaki, I., Lebichot, B., Saramaki, J. & Saerens, M. (2016). Two
betweenness centrality measures based on randomized shortest paths.
*Scientific Reports*, 6, 19668. (`BetweennessRSP`)

Banerjee, A., Chandrasekhar, A. G., Duflo, E. & Jackson, M. O. (2013).
The diffusion of microfinance. *Science*, 341(6144), 1236498.
(`Diffusion`)

Onnela, J.-P., Saramaki, J., Kertesz, J. & Kaski, K. (2005). Intensity
and coherence of motifs in weighted complex networks. *Physical Review
E*, 71, 065103. (`Clustering`)

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
