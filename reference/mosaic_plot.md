# Draw a mosaic plot of transition or co-occurrence counts

`mosaic_plot()` draws a mosaic display of a contingency table and fills
each tile with its standardised residual under independence: blue tiles
hold more cases than independence predicts, red tiles fewer. It reads
the transition or co-occurrence counts of a `netobject` or `htna`
network, the state counts of each group of a `netobject_group`, the
cluster-level or per-cluster counts of an `mcml`, or a `table` or count
`matrix` given directly. A data frame in long form is drawn as a plain
marimekko of shares instead, coloured by the stacked variable (or the
`fill` column), without residuals.

In a network's mosaic each column is a target state, its width
proportional to the transitions into it, and its segments are the source
states, stacked in proportion to their share of those transitions. For a
`netobject_group` the columns are the groups and the segments the
states, counted from each group's sequences, so the mosaic shows how the
state mix differs between groups. For a `table` or `matrix`, its rows
are the columns of the mosaic and its columns the stacked segments. An
`htna` network has its states ordered by actor, with each state's label
set in bold in its actor's colour.

The residuals describe where the counts depart from independence; they
are not tests of single cells. The default residuals treat every
transition or state as an independent draw, but those of one sequence
(or one student, or one team) are not. With `nested = TRUE` the
residuals allow for that nesting; `report = TRUE` adds the overall test
as the subtitle and prints it; see *Nested data*. To test edges or group
differences, use
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md);
for a chi-square or Fisher test with Cramér's V on two variables of a
data frame, use
[`mosaic_analysis()`](https://pak.dynasite.org/Nestimate/reference/mosaic_analysis.md).

## Usage

``` r
mosaic_plot(x, ...)

# Default S3 method
mosaic_plot(x, ...)

# S3 method for class 'netobject'
mosaic_plot(
  x,
  xlab = NULL,
  ylab = NULL,
  range = NULL,
  top_angle = NULL,
  left_angle = NULL,
  residuals = c("permutation", "asymptotic"),
  n_perm = 500L,
  seed = NULL,
  values = FALSE,
  style = c("classic", "flat"),
  ...,
  nested = FALSE,
  nested_method = c("rao_scott", "donner_klar", "cmh", "gee", "glmm"),
  actor = NULL,
  report = FALSE
)

# S3 method for class 'htna'
mosaic_plot(
  x,
  xlab = NULL,
  ylab = NULL,
  range = NULL,
  top_angle = NULL,
  left_angle = NULL,
  residuals = c("permutation", "asymptotic"),
  n_perm = 500L,
  seed = NULL,
  values = FALSE,
  style = c("classic", "flat"),
  ...,
  nested = FALSE,
  nested_method = c("rao_scott", "donner_klar", "cmh", "gee", "glmm"),
  actor = NULL,
  report = FALSE
)

# S3 method for class 'mcml'
mosaic_plot(
  x,
  level = c("macro", "clusters"),
  xlab = NULL,
  ylab = NULL,
  range = NULL,
  top_angle = NULL,
  left_angle = NULL,
  residuals = c("permutation", "asymptotic"),
  n_perm = 500L,
  seed = NULL,
  ncol = 2L,
  values = FALSE,
  style = c("classic", "flat"),
  ...
)

# S3 method for class 'netobject_group'
mosaic_plot(
  x,
  xlab = NULL,
  ylab = NULL,
  range = NULL,
  top_angle = NULL,
  left_angle = NULL,
  residuals = c("permutation", "asymptotic"),
  n_perm = 500L,
  seed = NULL,
  ncol = 2L,
  values = FALSE,
  style = c("classic", "flat"),
  ...,
  nested = FALSE,
  nested_method = c("rao_scott", "donner_klar", "cmh", "gee", "glmm"),
  actor = NULL,
  report = FALSE
)

# S3 method for class 'table'
mosaic_plot(
  x,
  xlab = "Row",
  ylab = "Column",
  range = NULL,
  top_angle = NULL,
  left_angle = NULL,
  residuals = c("permutation", "asymptotic"),
  n_perm = 500L,
  seed = NULL,
  values = FALSE,
  style = c("classic", "flat"),
  ...,
  nested = FALSE,
  nested_method = c("rao_scott", "donner_klar", "cmh", "gee", "glmm"),
  report = FALSE
)

# S3 method for class 'matrix'
mosaic_plot(x, ...)

# S3 method for class 'array'
mosaic_plot(x, ...)

# S3 method for class 'data.frame'
mosaic_plot(
  x,
  x_var,
  y_var,
  weight,
  fill = "y",
  colors = NULL,
  show_labels = TRUE,
  label_size = 3.5,
  xlab = NULL,
  ylab = NULL,
  ...
)
```

## Arguments

- x:

  The object to draw: a `netobject`, `netobject_group`, `mcml` or `htna`
  network, a contingency `table` (two-way, or three-way with the
  clusters as the third dimension), a numeric count `matrix`, or a data
  frame in long form (see `x_var`, `y_var` and `weight`). A `tna` model
  or a `group_tna` (package tna) is first rebuilt as a network from its
  own sequences (see
  [`as_netobject()`](https://pak.dynasite.org/Nestimate/reference/as_netobject.md)).

- ...:

  With `style = "flat"`: styling arguments for the flat mosaic of a
  network, table or matrix; an unknown name raises an error that lists
  the valid ones. Raises `nestimate_bad_input` when the data-frame
  method, or the `mcml` method with `style = "classic"`, gets an
  argument it does not take; the other methods ignore such an argument
  with a warning.

- xlab, ylab:

  Axis titles. The network and `mcml` methods draw none by default
  (`NULL`); the `table` and `matrix` methods default to `"Row"` and
  `"Column"`; a data frame defaults to the `x_var` name and
  `"Proportion"`. Not drawn with `style = "flat"`.

- range:

  The lower and upper limits of the residual colour scale, two
  increasing numbers; residuals beyond them take the end colours.
  Default `NULL`, which fits the symmetric range `c(-M, M)`, with `M`
  the largest absolute residual (at least 1), so no tile is clipped.
  `c(-4, 4)` matches tna's display.

- top_angle, left_angle:

  Rotation, in degrees, of the column labels on top and the segment
  labels on the left. Default `NULL`: 90 when an axis has more than
  three levels, 0 otherwise. Not used with `style = "flat"` or for the
  per-cluster `mcml` panels.

- residuals:

  How the standardised residuals are computed. `"permutation"` (default)
  shuffles one variable's labels against the other `n_perm` times, which
  keeps both margins, and gives each cell `(observed - mean) / sd` over
  the shuffles; a cell that never varies gets 0. It does not rely on
  large expected counts, so it suits sparse tables. `"asymptotic"` gives
  the adjusted residuals of
  [`stats::chisq.test()`](https://rdrr.io/r/stats/chisq.test.html),
  `(O - E) / sqrt(E (1 - p_row) (1 - p_col))`, as vcd and tna do. Used
  when `nested = FALSE`; given with `nested = TRUE` it is ignored with a
  warning of class `nestimate_ignored_argument`.

- n_perm:

  The number of shuffles for `residuals = "permutation"`, at least 2.
  Default `500`; use 1000 or more for stable values in the tails.

- seed:

  A seed for the shuffles, so the residuals repeat exactly. The
  session's own random stream is left as it was. Default `NULL`. Ignored
  with `residuals = "asymptotic"`.

- values:

  Whether to print each tile's residual (one decimal) inside it. Default
  `FALSE`. In the classic style the text is white on strongly coloured
  tiles (absolute residual above 1.5) and dark grey otherwise; with
  `style = "flat"` it is white on tiles whose absolute residual exceeds
  0.55 of the largest, and tiles too small for a label stay blank.

- style:

  `"classic"` (default) draws black-bordered tiles with axis labels.
  `"flat"` draws the mosaic of
  [`mosaic_analysis()`](https://pak.dynasite.org/Nestimate/reference/mosaic_analysis.md):
  white gutters between tiles, labels beside them, no axes, and a
  ColorBrewer red-blue residual ramp; `values = TRUE` then prints the
  residuals in the tiles, and its styling arguments (`tile_label`,
  `col_label_side`, `row_label_side`, `legend_position`, `legend_size`,
  `label_size`, ...) go in `...`. An `mcml` with `level = "clusters"`
  and more than one panel is drawn classic, with a message; an `htna`
  drawn flat loses its actor-coloured labels, with a warning.

- nested:

  Whether the counts are nested in actors (transitions in sequences,
  sequences in students or teams) and the residuals should allow for it.
  Default `FALSE`: every count is an independent draw, as before. `TRUE`
  draws the residuals of `nested_method`; see *Nested data*.

- nested_method:

  With `nested = TRUE`, the method: `"rao_scott"` (default),
  `"donner_klar"`, `"cmh"`, `"gee"` or `"glmm"`; see *Nested data*. A
  method that cannot run on the design raises an error of class
  `nestimate_bad_input`.

- actor:

  With `nested = TRUE`, on a network or `netobject_group`: `NULL`
  (default) makes every sequence its own actor; the name of a column of
  the network's metadata (such as `"Actor"` or a team column, as `actor`
  in
  [`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
  and
  [`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md))
  pools the sequences that share its value; several names give nested
  levels, outermost first, which `"glmm"` models level by level and the
  other methods reduce to the outermost. A name that is not a column
  raises an error of class `nestimate_actor_missing`. `nested_method`
  and `actor` given with `nested = FALSE` are ignored with a warning of
  class `nestimate_ignored_argument`.

- report:

  Whether to report the overall test: `TRUE` adds it as the plot
  subtitle and prints it as a one-row data frame with the columns of
  [`mosaic_test()`](https://pak.dynasite.org/Nestimate/reference/mosaic_test.md)
  but `note`. The test is `nested_method` with `nested = TRUE`,
  Pearson's chi-square otherwise. Default `FALSE`. Not available for an
  `mcml`.

- level:

  For an `mcml`: `"macro"` (default) draws one mosaic of the
  cluster-by-cluster counts; `"clusters"` draws one mosaic per cluster,
  as panels of one plot that share the residual scale.

- ncol:

  For an `mcml` with `level = "clusters"`: the number of panel columns.
  Default `2`. Accepted by the `netobject_group` method, whose single
  mosaic does not use it.

- x_var, y_var, weight:

  For a data frame: the column whose levels form the mosaic's columns
  (each as wide as its total `weight`), the column whose levels are
  stacked within each column (in proportion to their share), and the
  column holding each row's weight, such as a count. Levels keep their
  factor order, or sort alphabetically. `x_var`, `y_var`, `weight`,
  `fill`, `colors` and `show_labels` given with any other input raise an
  error of class `nestimate_bad_input`, as does `label_size` outside
  `style = "flat"`.

- fill:

  For a data frame: `"y"` (default) colours the segments by `y_var`; the
  name of another column maps that column to the fill instead (the first
  value per cell).

- colors:

  For a data frame with `fill = "y"`: fill colours, either an unnamed
  vector in level order (at least one per level) or a named lookup such
  as `c(plan = "#0072B2")` that overrides only the levels named. Default
  `NULL`, the Okabe-Ito palette, recycled.

- show_labels, label_size:

  For a data frame: whether to print each segment's percentage inside it
  (default `TRUE`; segments under 4% stay unlabelled), and the label
  size (default `3.5`). With `style = "flat"`, `label_size` also sizes
  the flat mosaic's tile labels.

## Value

A `ggplot` object, one rectangle per cell (empty cells are left out of
the flat and data-frame mosaics), filled by the standardised residual
(by category for a data frame). Print it to draw it, or add ggplot2
layers and themes. An `mcml` with `level = "clusters"` gives one faceted
`ggplot` with a panel per cluster and one shared fill scale. An object
of an unsupported class raises an error.

## Details

A mosaic needs counts, never probabilities. A network is drawn from the
transition counts stored with it at estimation; failing that from its
weights, when its method counts (such as `"frequency"`) and they are not
scaled, thresholded or weighted; and failing that from first-order
transitions recounted from its sequence data. A network that records no
method is drawn from its weights if they are whole numbers. A
`netobject_group` is counted from each group's sequences, or from the
column sums of its weights when they are counts. A single network with
none of these, such as an association network (`cor`, `pcor`, `glasso`,
`ising`, `mgm`), raises an error.

A co-occurrence network (`method = "co_occurrence"`) is undirected: it
stores each unordered pair of positions in both mirrored cells and the
same-state pairs once on the diagonal, so its matrix is not a table of
row-by-column counts. The mosaic draws its table of pairs instead, each
pair counted once and split equally between its two mirrored cells (so a
cell can hold a half; Newman, 2003), and the residuals compare it with
random pairing of the pair ends, which keeps how often each state ends a
pair: `"permutation"` re-pairs the ends at random `n_perm` times, and
`"asymptotic"` gives each unordered pair's multinomial residual under
that null, as for genotypes under Hardy-Weinberg proportions (Guo &
Thompson, 1992). The pairs of a sequence share its positions, so they
are not independent counts, and `report = TRUE` or `nested = TRUE` on a
co-occurrence network raises an error of class
`nestimate_mosaic_undirected`. So does any mosaic of an undirected
network whose diagonal holds each state's own occurrences rather than
pairs: a windowed or one-hot co-occurrence network,
`"wtna_cooccurrence"` and
[`cooccurrence()`](https://pak.dynasite.org/Nestimate/reference/cooccurrence.md).

A `netobject_group` must hold sequence networks: the group method counts
whatever its networks' data hold, so a group of association networks is
drawn without an error, each distinct numeric value counted as a state.
An `mcml` is counted the same way as a network, layer by layer; a
cluster with a single state has no table to draw: when it holds
self-transitions it is left out of the per-cluster figure with a
warning, and when it holds no transition at all the call stops with an
error. When a row or column of the table is empty, its cells have no
asymptotic residual and are drawn at 0; the other cells keep theirs. The
geometry and the asymptotic residuals follow
[`tna::plot_mosaic()`](https://sonsoles.me/tna/reference/plot_mosaic.html);
the function has its own name so the two packages can be attached
together.

## Nested data

A network's counts are summed over sequences. With `nested = TRUE` the
counts are taken sequence by sequence (as the network counted them),
summed within the actors named by `actor` (each sequence on its own by
default), and `nested_method` reads the table of each actor. With
`report = TRUE` the subtitle names the test, its statistic, degrees of
freedom and p-value, the number of actors, and the design effect or
intracluster correlation where the method gives one. For a network, the
source state is treated as the predictor and the target as the response;
for a `netobject_group`, the group is the predictor and the state the
response. The `nested_method` choices:

- `"rao_scott"` (default):

  Rao-Scott second-order corrected chi-square, an F statistic on
  fractional degrees of freedom (as `survey::svychisq()`); any cluster
  sizes. Residual: the standardised residual divided by the square root
  of the cell's design effect. With one observation per actor it equals
  the default asymptotic residual, so it corrects only as much as the
  data are clustered.

- `"donner_klar"`:

  Pearson's chi-square divided by the cluster-size inflation factor
  `1 + (m - 1) * rho`. Needs one variable that is constant within every
  cluster, as the group of a `netobject_group`; a single network raises
  an error.

- `"cmh"`:

  Generalised Cochran-Mantel-Haenszel test with each cluster as a
  stratum; both variables must vary within clusters, so a
  `netobject_group` raises an error.

- `"gee"`:

  Wald chi-square of a marginal model (identity link, independence
  working correlation) with the cluster-robust covariance; refused when
  that covariance is rank-deficient, and liberal unless the actors are
  many relative to the degrees of freedom. Residual: conditional minus
  pooled share, over its cluster-robust standard error.

- `"glmm"`:

  Wald chi-square on the interaction of a log-linear Poisson mixed model
  with cluster random effects at every level of `actor`; needs lme4.
  Residual: observed minus expected count under independence, over its
  cluster-robust standard error.

The methods are described in full under *Nested data* in
[`mosaic_analysis()`](https://pak.dynasite.org/Nestimate/reference/mosaic_analysis.md).
`nested = TRUE` needs per-sequence counts: a network built with
`method = "relative"` or `"frequency"` that keeps its sequence data, or
a `netobject_group` whose groups keep theirs. Other networks, and
per-sequence counts that do not add up to the drawn table, raise an
error of class `nestimate_bad_input`; a co-occurrence network raises one
of class `nestimate_mosaic_undirected` (see Details). An `mcml` takes no
`nested`. A `table` or `matrix` (or a plain array) takes it as a
three-way table, rows x columns x actors (e.g.
`table(var1, var2, actor)`); the mosaic draws the table summed over
actors. A three-way table without `nested = TRUE`, or a two-way table
with it, raises an error. Residuals a method leaves undefined are drawn
at 0. For the tests as a table, all methods side by side, use
[`mosaic_test()`](https://pak.dynasite.org/Nestimate/reference/mosaic_test.md).

## References

Friendly, M. (1994). Mosaic displays for multi-way contingency tables.
*Journal of the American Statistical Association*, 89(425), 190–200.
[doi:10.1080/01621459.1994.10476460](https://doi.org/10.1080/01621459.1994.10476460)

Guo, S. W., & Thompson, E. A. (1992). Performing the exact test of
Hardy-Weinberg proportion for multiple alleles. *Biometrics*, 48(2),
361–372. [doi:10.2307/2532296](https://doi.org/10.2307/2532296)

Newman, M. E. J. (2003). Mixing patterns in networks. *Physical Review
E*, 67(2), 026126.
[doi:10.1103/PhysRevE.67.026126](https://doi.org/10.1103/PhysRevE.67.026126)

Rao, J. N. K., & Scott, A. J. (1984). On chi-squared tests for multiway
contingency tables with cell proportions estimated from survey data.
*The Annals of Statistics*, 12(1), 46–60.
[doi:10.1214/aos/1176346391](https://doi.org/10.1214/aos/1176346391)

Donner, A., & Klar, N. (1994). Methods for comparing event rates in
intervention studies when the unit of allocation is a cluster. *American
Journal of Epidemiology*, 140(3), 279–289.
[doi:10.1093/oxfordjournals.aje.a117247](https://doi.org/10.1093/oxfordjournals.aje.a117247)

Landis, J. R., Heyman, E. R., & Koch, G. G. (1978). Average partial
association in three-way contingency tables: A review and discussion of
alternative tests. *International Statistical Review*, 46(3), 237–254.
[doi:10.2307/1402373](https://doi.org/10.2307/1402373)

Liang, K.-Y., & Zeger, S. L. (1986). Longitudinal data analysis using
generalized linear models. *Biometrika*, 73(1), 13–22.
[doi:10.1093/biomet/73.1.13](https://doi.org/10.1093/biomet/73.1.13)

## Examples

``` r
data(group_regulation_long, package = "Nestimate")
net <- build_network(group_regulation_long, method = "frequency",
                     format = "long", actor = "Actor", action = "Action",
                     order = "Time")
mosaic_plot(net, seed = 1)


# Transitions are nested in students' sequences, students in teams:
mosaic_plot(net, nested = TRUE, actor = "Group", report = TRUE)
#>     method      test statistic       df      df2 p_value     n clusters
#>  rao_scott Rao-Scott  227.7646 39.41762 7844.106       0 25533      200
#>  design_effect icc
#>       1.058882  NA


# A tidy data frame: columns by group, segments by state
df <- data.frame(
  group = rep(c("A", "B", "C"), each = 3),
  state = rep(c("s1", "s2", "s3"), 3),
  count = c(10, 5, 2,  7, 8, 3,  4, 6, 12)
)
mosaic_plot(df, x_var = "group", y_var = "state", weight = "count")
```
