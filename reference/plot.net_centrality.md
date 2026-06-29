# Plot centrality measures

Plot centrality measures

## Usage

``` r
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
```

## Arguments

- x:

  A `net_centrality` object returned by
  [`net_centrality`](https://saqr.me/Nestimate/reference/net_centrality.md).

- reorder:

  Logical. Reorder states within each centrality panel by centrality
  value. Default: `TRUE`.

- ncol:

  Integer. Number of facet columns.

- type:

  Plot type. `"bar"` shows one faceted horizontal bar chart per measure;
  `"line"` shows state profiles as lines across measures (the value
  `"profile"` is still accepted as an alias); `"heatmap"` shows a
  states-by-measures tile grid, each measure scaled to 0–1 for
  cross-measure comparability with the raw value printed in the tile.

- scales:

  Facet scale mode. `"free_x"` uses free centrality axes; `"fixed"`
  keeps a common centrality axis.

- profile_scale:

  Scaling used by `type = "line"`. `"measure"` rescales each centrality
  measure to 0–1 before drawing cross-measure profiles; `"none"` uses
  raw values.

- labels:

  Logical. Add compact value labels.

- drop_zero:

  Logical. Drop measures whose values are all (near) zero so empty
  panels do not waste space. Default: `FALSE` (every requested measure
  is shown).

- ...:

  Additional arguments ignored.

## Value

A `ggplot` object.
