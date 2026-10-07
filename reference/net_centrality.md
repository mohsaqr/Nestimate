# Compute Centrality Measures for a Network

Computes centrality measures from a `netobject`, `netobject_group`,
`mcml`, or `cograph_network`. The built-in measures match
[`tna::centralities()`](https://sonsoles.me/tna/reference/centralities.html)
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

# S3 method for class 'net_centrality'
plot(
  x,
  reorder = TRUE,
  ncol = 3L,
  type = c("bar", "line", "heatmap"),
  scales = c("free_x", "fixed"),
  profile_scale = c("measure", "none"),
  labels = TRUE,
  drop_zero = FALSE,
  ...
)

# S3 method for class 'net_centrality_group'
plot(
  x,
  reorder = TRUE,
  ncol = 3L,
  type = c("bar", "line", "delta"),
  scales = c("free_x", "fixed"),
  palette = "Set2",
  profile_scale = c("measure", "none"),
  labels = FALSE,
  drop_zero = FALSE,
  ...
)
```

## Arguments

- x:

  A `netobject`, `netobject_group`, `mcml`, or `cograph_network`. For
  the [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method:
  an object of class `net_centrality` or `net_centrality_group`.

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

  Additional arguments (ignored). In `plot.net_centrality()` and
  `plot.net_centrality_group()`: Additional arguments ignored.

- reorder:

  In `plot.net_centrality()`: Logical. Reorder states within each
  centrality panel by centrality value. Default: `TRUE`. In
  `plot.net_centrality_group()`: Logical. Reorder states by their mean
  value within each centrality panel. Default: `TRUE`.

- ncol:

  Integer. Number of facet columns. Default: `3`.

- type:

  In `plot.net_centrality()`: Plot type. `"bar"` shows one faceted
  horizontal bar chart per measure; `"line"` shows state profiles as
  lines across measures (the value `"profile"` is still accepted as an
  alias); `"heatmap"` shows a states-by-measures tile grid, each measure
  scaled to 0–1 for cross-measure comparability with the raw value
  printed in the tile. Default: `"bar"`. In
  `plot.net_centrality_group()`: Plot type. `"bar"` shows grouped bars
  within each measure; `"line"` facets by state and draws one line per
  group across centrality measures (`"profile"` is accepted as an
  alias); `"delta"` draws a diverging bar of group differences. With two
  groups it is the per-state difference (second group minus first); with
  three or more groups it is each group's deviation from the per-state
  group mean, so the largest gaps stand out either way. Default:
  `"bar"`.

- scales:

  Facet scale mode. `"free_x"` (default) uses free centrality axes;
  `"fixed"` keeps a common centrality axis.

- profile_scale:

  Scaling used by `type = "line"`. `"measure"` (default) rescales each
  centrality measure to 0–1 before drawing cross-measure profiles;
  `"none"` uses raw values.

- labels:

  In `plot.net_centrality()`: Logical. Add compact value labels.
  Default: `TRUE`. In `plot.net_centrality_group()`: Logical. Add
  compact value labels. Default: `FALSE`.

- drop_zero:

  In `plot.net_centrality()`: Logical. Drop measures whose values are
  all (near) zero so empty panels do not waste space. Default: `FALSE`
  (every requested measure is shown). In `plot.net_centrality_group()`:
  Logical. Drop measures whose values are all (near) zero so empty
  panels do not waste space. Default: `FALSE`.

- palette:

  Brewer palette for groups. Default: `"Set2"`.

## Value

For a `netobject` or `cograph_network`: a `net_centrality` data frame,
one row per node, with a `state` column and one further column per
requested measure (node names are also the row names). For a
`netobject_group` or an `mcml`: a `net_centrality_group` list of such
data frames, one per group.

In `plot.net_centrality()` and `plot.net_centrality_group()`: A `ggplot`
object.

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
