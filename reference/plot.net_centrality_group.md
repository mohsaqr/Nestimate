# Plot grouped centrality measures

Plot grouped centrality measures

## Usage

``` r
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

  A `net_centrality_group` object returned by
  [`net_centrality`](https://saqr.me/Nestimate/reference/net_centrality.md).

- reorder:

  Logical. Reorder states by their mean value within each centrality
  panel. Default: `TRUE`.

- ncol:

  Integer. Number of facet columns.

- type:

  Plot type. `"bar"` shows grouped bars within each measure; `"line"`
  facets by state and draws one line per group across centrality
  measures (`"profile"` is accepted as an alias); `"delta"` draws a
  diverging bar of group differences. With two groups it is the
  per-state difference (second group minus first); with three or more
  groups it is each group's deviation from the per-state group mean, so
  the largest gaps stand out either way.

- scales:

  Facet scale mode. `"free_x"` uses free centrality axes; `"fixed"`
  keeps a common centrality axis.

- palette:

  Brewer palette for groups.

- profile_scale:

  Scaling used by `type = "line"`. `"measure"` rescales each centrality
  measure to 0–1 before drawing cross-measure profiles; `"none"` uses
  raw values.

- labels:

  Logical. Add compact value labels.

- drop_zero:

  Logical. Drop measures whose values are all (near) zero so empty
  panels do not waste space. Default: `FALSE`.

- ...:

  Additional arguments ignored.

## Value

A `ggplot` object.
