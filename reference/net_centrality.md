# Compute Centrality Measures for a Network

`net_centrality()` scores how central each node is in a network, for one
or more measures, and returns one row per node. Strength sums a node's
edge weights (incoming or outgoing); closeness and betweenness read the
shortest paths between nodes, where a heavier edge is a shorter step;
diffusion counts the weighted walks that start at a node; and clustering
measures how strongly a node's neighbours are tied to each other. It
takes a single network or a group of them (a `netobject_group`, the
clusters of an `mcml`), scoring each separately.

Centralities describe the estimated network, not the population it was
estimated from. Rankings built on a few weak edges change from sample to
sample: check them with
[`centrality_stability()`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md)
before reporting which node is most central.

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

  The network: a `netobject`, `netobject_group`, `mcml` or
  `cograph_network`. For the
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: the
  result of `net_centrality()`.

- measures:

  The measures to compute. The default is `"InStrength"`,
  `"Betweenness"` and `"Diffusion"`; `"all"` computes every built-in
  one: `"OutStrength"` and `"InStrength"` (sum of outgoing or incoming
  weights), `"ClosenessOut"`, `"ClosenessIn"` and `"Closeness"` (inverse
  of the summed shortest-path distances from, to, or ignoring the
  direction of edges), `"Betweenness"` (how many shortest paths pass
  through the node), `"BetweennessRSP"` (betweenness over randomised
  shortest paths, which also counts near-shortest routes), `"Diffusion"`
  (the summed weights of walks of length 1 to n from the node) and
  `"Clustering"` (weighted clustering coefficient). `"InCloseness"` and
  `"OutCloseness"` are accepted for `"ClosenessIn"` and
  `"ClosenessOut"`. With `centrality_fn`, name its measures here to have
  them returned.

- loops:

  Whether self-loops (staying in a state) count. With `FALSE` (default)
  they are left out of every measure, and a message says so.

- normalize:

  Whether to range-normalise every measure to 0–1, as
  `tna::centralities(normalize = TRUE)` does. Default `FALSE`.

- invert:

  Whether a heavier edge is a shorter step in the shortest-path measures
  (distance = 1 / weight), as in `tna`. With `FALSE`, the weight itself
  is the distance. Default `TRUE`.

- normalize_diffusion:

  Whether to range-normalise `Diffusion` even when `normalize = FALSE`.
  Default `TRUE`.

- centrality_fn:

  A function of your own that takes the weight matrix and returns a
  named list of centrality vectors, one per measure; list their names in
  `measures`. Default `NULL`.

- ...:

  Not used.

- reorder:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  methods: whether to sort the nodes in each measure's panel by their
  value (for a group, by their mean value across networks). Default
  `TRUE`.

- ncol:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  methods: the number of panel columns. Default `3`.

- type:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  methods: the kind of plot. `"bar"` (default) draws one panel of bars
  per measure, grouped bars for several networks. `"line"` draws each
  node's profile across the measures (for several networks, one panel
  per node and one line per network); `"profile"` is accepted for it. A
  single network also takes `"heatmap"`, a nodes-by-measures grid with
  each measure scaled to 0–1 for the colour and its raw value printed in
  the tile. Several networks also take `"delta"`, bars of the
  differences between networks: the second minus the first with two,
  each network's deviation from the per-node mean with more.

- scales:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  methods: `"free_x"` (default) gives each measure its own value axis,
  `"fixed"` one shared axis.

- profile_scale:

  For `type = "line"`: `"measure"` (default) rescales each measure to
  0–1 so the profiles compare across measures; `"none"` draws raw
  values.

- labels:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  methods: whether to print the values on the plot. Default `TRUE` for
  one network, `FALSE` for several.

- drop_zero:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  methods: whether to leave out measures that are zero for every node.
  Default `FALSE`.

- palette:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of
  several networks: the ColorBrewer palette for the networks. Default
  `"Set2"`. Ignored by `type = "delta"` with two networks, whose bars
  are coloured by which network is higher.

## Value

For a `netobject` or `cograph_network`: a `net_centrality` data frame
with one row per node, a `state` column and one column per requested
measure (the node names are also the row names). For a `netobject_group`
or an `mcml`: a `net_centrality_group`, a list of such data frames, one
per network. [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
draws the measures per node.

For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
methods: a `ggplot` object.

## Details

Strength is read from the weight matrix directly, and closeness and
betweenness from all-pairs shortest paths (Floyd-Warshall), without
igraph; the values agree with
[`tna::centralities()`](https://sonsoles.me/tna/reference/centralities.html).
Diffusion is range-normalised by default, which `tna` does not do.

Range normalisation of a measure that is the same for every node is
undefined and gives `NaN`, with a warning of class
`nestimate_constant_measure`. Diffusion on a `"relative"` network with
`loops = TRUE` is such a case: each row sums to 1, so every node has the
same diffusion; use `normalize_diffusion = FALSE` there.

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
#> Warning: Diffusion is the same for every node, so it cannot be range-normalised and is NaN. Use `normalize_diffusion = FALSE` for the raw values.
#>   state InStrength Betweenness Diffusion
#> A     A          1           1       NaN
#> B     B          1           1       NaN
#> C     C          1           1       NaN
```
