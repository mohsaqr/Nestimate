# Plot State Frequency Distributions

Visualise state (node) frequency distributions across groups for any
Nestimate object that carries sequence data: a single `netobject`, a
`netobject_group`, an `mcml` model, or an `htna` network.

## Usage

``` r
plot_state_frequencies(
  x,
  style = c("marimekko", "bars", "mosaic"),
  metric = c("prop", "freq"),
  label = c("prop", "freq", "both", "state", "all", "none"),
  legend = c("bottom", "right", "top", "left", "none", "per_facet"),
  legend_dir = c("auto", "horizontal", "vertical"),
  legend_frame = c("none", "border"),
  sort_states = c("frequency", "alpha", "none"),
  colors = NULL,
  label_size = 3,
  include_macro = FALSE,
  combine = TRUE,
  ncol = NULL,
  node_groups = NULL,
  ...
)

# S3 method for class 'netobject'
plot_state_frequencies(
  x,
  style = c("marimekko", "bars", "mosaic"),
  metric = c("prop", "freq"),
  label = c("prop", "freq", "both", "state", "all", "none"),
  legend = c("bottom", "right", "top", "left", "none", "per_facet"),
  legend_dir = c("auto", "horizontal", "vertical"),
  legend_frame = c("none", "border"),
  sort_states = c("frequency", "alpha", "none"),
  colors = NULL,
  label_size = 3,
  include_macro = FALSE,
  combine = TRUE,
  ncol = NULL,
  node_groups = NULL,
  ...
)

# S3 method for class 'htna'
plot_state_frequencies(
  x,
  style = c("marimekko", "bars", "mosaic"),
  metric = c("prop", "freq"),
  label = c("prop", "freq", "both", "state", "all", "none"),
  legend = c("per_facet", "bottom", "right", "top", "left", "none"),
  legend_dir = c("auto", "horizontal", "vertical"),
  legend_frame = c("border", "none"),
  sort_states = c("frequency", "alpha", "none"),
  colors = NULL,
  label_size = 3,
  include_macro = FALSE,
  combine = TRUE,
  ncol = NULL,
  node_groups = NULL,
  ...
)

# S3 method for class 'mcml'
plot_state_frequencies(
  x,
  style = c("marimekko", "bars", "mosaic"),
  metric = c("prop", "freq"),
  label = c("prop", "freq", "both", "state", "all", "none"),
  legend = c("per_facet", "bottom", "right", "top", "left", "none"),
  legend_dir = c("auto", "horizontal", "vertical"),
  legend_frame = c("none", "border"),
  sort_states = c("frequency", "alpha", "none"),
  colors = NULL,
  label_size = 3,
  include_macro = FALSE,
  combine = TRUE,
  ncol = NULL,
  node_groups = NULL,
  ...
)

# S3 method for class 'netobject_group'
plot_state_frequencies(
  x,
  style = c("marimekko", "bars", "mosaic"),
  metric = c("prop", "freq"),
  label = c("prop", "freq", "both", "state", "all", "none"),
  legend = c("bottom", "right", "top", "left", "none", "per_facet"),
  legend_dir = c("auto", "horizontal", "vertical"),
  legend_frame = c("none", "border"),
  sort_states = c("frequency", "alpha", "none"),
  colors = NULL,
  label_size = 3,
  include_macro = FALSE,
  combine = TRUE,
  ncol = NULL,
  node_groups = NULL,
  ...
)

# Default S3 method
plot_state_frequencies(
  x,
  style = c("marimekko", "bars", "mosaic"),
  metric = c("prop", "freq"),
  label = c("prop", "freq", "both", "state", "all", "none"),
  legend = c("bottom", "right", "top", "left", "none"),
  sort_states = c("frequency", "alpha", "none"),
  colors = NULL,
  label_size = 3,
  include_macro = FALSE,
  ...
)
```

## Arguments

- x:

  A `netobject`, `netobject_group`, `mcml`, or `htna` object.

- style:

  Either `"marimekko"` (default) or `"bars"`.

- metric:

  For `style = "bars"`: which value the bar length encodes – `"prop"`
  (default) or `"freq"`. Treemap and hierarchical-marimekko areas always
  encode proportion within group.

- label:

  Inline tile / bar annotation. All formats render on a single line.

  - `"prop"` (default) – proportion only, e.g. `"66%"`

  - `"freq"` – count only, e.g. `"1,234"`

  - `"both"` – count + proportion, e.g. `"1,234 (66%)"`

  - `"state"` – state name only, e.g. `"Average"`

  - `"all"` – state + proportion, e.g. `"Average (66%)"`

  - `"none"` – no inline labels

- legend:

  Legend position. `"bottom"` (default), `"top"`, `"right"`, `"left"`,
  `"none"` (hide; pair with `label = "state"` or `"all"` so state names
  show on tiles), or `"per_facet"` – each group renders as its own panel
  with an isolated legend showing only the states present in that panel.
  For `htna` this gives each actor (AI, Human) its own legend; for
  `mcml` each cluster gets its own. Requires the gridExtra package and
  returns a `gtable`.

- legend_dir:

  Legend internal layout: `"auto"` (default – horizontal for top/bottom,
  vertical for left/right), or force `"horizontal"` or `"vertical"`
  regardless of position.

- legend_frame:

  `"none"` (default) for an unframed legend, or `"border"` to draw a
  thin grey rectangle around the legend ("legend enclosed in a square").

- sort_states:

  One of `"frequency"` (default – most frequent first), `"alpha"`, or
  `"none"`.

- colors:

  Optional character vector overriding the default Okabe-Ito state
  palette. Length must be at least the number of unique states.

- label_size:

  Numeric size of inline labels.

- include_macro:

  For `mcml` only: prepend a `"macro"` reference column showing
  aggregate state frequencies across all clusters. Default `FALSE`.

- combine:

  For `legend = "per_facet"` only. When `TRUE` (default), per-panel
  ggplots are arranged into a single gtable via gridExtra. When `FALSE`,
  returns a list of ggplots that are printed one-per-figure by knitr (so
  each panel uses the full chunk `fig.width` / `fig.height`).

- ncol:

  For `legend = "per_facet"` with `combine = TRUE`: number of columns in
  the grid arrangement. `NULL` (default) picks 1, 2, or 3 columns based
  on the number of panels.

- ...:

  Reserved for future use.

## Value

A `ggplot` object.

## Details

The marimekko layout is dispatched per class:

- For `mcml`, where states partition cleanly into clusters, the chart is
  a hierarchical 2D marimekko: cluster columns of width proportional to
  cluster total, segments stacked vertically with heights proportional
  to within-cluster state proportions.

- For all other classes (`netobject`, `netobject_group`, `htna`), each
  group is rendered as its own panel containing a squarified treemap:
  each state becomes a rectangular tile whose AREA is exactly
  proportional to the state's share within that group. Single-panel when
  no groups exist; faceted when groups are present.

The bar style produces horizontal bars (state on the y-axis), faceted by
group when groups exist. All variants use the Okabe-Ito palette.

## Examples

``` r
# \donttest{
if (requireNamespace("ggplot2", quietly = TRUE)) {
  data(group_regulation_long, package = "Nestimate")
  nw <- build_network(group_regulation_long,
                      method = "relative", format = "long",
                      actor = "Actor", action = "Action",
                      order = "Time", group = "Course")
  plot_state_frequencies(nw)
  plot_state_frequencies(nw, style = "bars")
}

# }
```
