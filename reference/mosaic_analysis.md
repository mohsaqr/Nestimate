# Test and draw the association of two categorical variables

`mosaic_analysis()` tests the association between two categorical
columns of a data frame and draws it as a mosaic. It builds their
contingency table, drops sparse categories (fewer than `min_count`
rows), runs a Pearson chi-square test (or Fisher's exact test), computes
Cramer's V with an effect-size label, and draws a flat ggplot2 mosaic in
which tile area is the cell count and the fill is the cell's
standardised residual. The per-cell results come as a data frame, one
row per cell.

The p-value tests the null hypothesis that the two variables are
independent; with large samples almost any departure is significant, so
read Cramer's V for the size of the association. The label follows
Cohen's (1988) thresholds for the table's degrees of freedom,
`min(rows, columns) - 1`: 0.10 / 0.30 / 0.50 for small / medium / large
at 1, 0.07 / 0.21 / 0.35 at 2, 0.06 / 0.17 / 0.29 at 3, and 0.05 / 0.15
/ 0.25 above (Cohen's w divided by the square root of the degrees of
freedom); a V below the small threshold is labelled `"negligible"`. The
residuals are
[`chisq.test()`](https://rdrr.io/r/stats/chisq.test.html)'s standardised
residuals, roughly standard normal under independence, so one beyond
about +-2 marks a cell that holds more (or fewer) rows than independence
predicts; with many cells, some pass that line by chance.

## Usage

``` r
mosaic_analysis(
  data,
  var1,
  var2,
  min_count = 10L,
  test = c("chisq", "fisher"),
  percentage_base = c("total", "row", "column"),
  tile_label = c("count", "percent", "residual", "category", "none"),
  title = "",
  ...,
  nested = FALSE,
  nested_method = c("rao_scott", "donner_klar", "cmh", "gee", "glmm"),
  actor = NULL,
  report = FALSE
)

# S3 method for class 'mosaic_analysis'
plot(x, ...)

# S3 method for class 'mosaic_analysis'
print(x, ...)

# S3 method for class 'mosaic_analysis'
summary(object, ...)
```

## Arguments

- data:

  A data frame holding the two variables. Rows missing either of them
  (or an `actor` column) are left out; a `var1` or `var2` that is not a
  column of `data`, or no complete rows, raises an error.

- var1:

  The name of the first variable, drawn as the mosaic's columns (column
  width is its marginal share).

- var2:

  The name of the second variable, stacked within each column.

- min_count:

  The smallest number of rows a category must have to be kept.
  Categories of either variable with fewer rows in the full table are
  dropped, in one pass, before the test. Default `10`. If fewer than two
  categories of either variable remain, an error is raised.

- test:

  The test reported: `"chisq"` (default) for Pearson's chi-square test,
  or `"fisher"` for Fisher's exact test: exact for a 2 x 2 table, and
  otherwise with a p-value simulated from 10,000 tables (so it changes
  with the random seed and cannot fall below 1/10,001). Cramer's V and
  the residuals come from Pearson's chi-square (V from the statistic
  without Yates' correction). Used when `nested = FALSE`; given with
  `nested = TRUE` it is ignored with a warning of class
  `nestimate_ignored_argument`.

- percentage_base:

  What the `pct` column and the `"percent"` tile label are a percentage
  of: `"total"` (default) the whole table, `"row"` each level of `var1`,
  `"column"` each level of `var2`.

- tile_label:

  What to print inside each tile: `"count"` (default), `"percent"`,
  `"residual"`, `"category"` (the `var2` level) or `"none"`.

- title:

  The plot title. Default `""` (none).

- ...:

  Further styling for the flat mosaic:

  - labels: `col_label_side`, `row_label_side`, `col_label_angle`,
    `row_label_angle`, `label_size`, `min_label_h` and `min_label_w`;

  - legend: `show_legend`, `legend_position`, `legend_size` and
    `legend_title`;

  - fill: `range` and `palette` (default the ColorBrewer RdBu ramp);

  - `pct_base`, which overrides `percentage_base` in the plot only.

  Column labels turn vertical when there are more than 6 columns;
  `col_label_angle` sets an angle. Any other name raises an error. For
  the [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method:
  the same styling and `title`, overriding what the fit stored; for the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods: ignored.

- nested:

  Whether the rows are nested in actors and the test should allow for
  it. Default `FALSE`: every row is independent, as before. `TRUE`
  reports `nested_method` and needs `actor`; see *Nested data*.

- nested_method:

  With `nested = TRUE`, the test: `"rao_scott"` (default),
  `"donner_klar"`, `"cmh"`, `"gee"` or `"glmm"`. A test that cannot run
  on the design raises `nestimate_bad_input`; `"glmm"` needs the lme4
  package.

- actor:

  With `nested = TRUE` (required then): the name of the column whose
  values identify the actors, or several names for nested levels,
  outermost first (e.g. `c("Group", "Actor")`). `nested_method` and
  `actor` given with `nested = FALSE` are ignored with a warning of
  class `nestimate_ignored_argument`.

- report:

  Whether to report the test: `TRUE` draws it as the plot subtitle and
  prints it as a one-row data frame with the columns of
  [`mosaic_test()`](https://pak.dynasite.org/Nestimate/reference/mosaic_test.md)
  but `note`. Default `FALSE`. The test is always in the `stats` table
  and in [`print()`](https://rdrr.io/r/base/print.html).

- x:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`print()`](https://rdrr.io/r/base/print.html) methods: a
  `mosaic_analysis`.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `mosaic_analysis`.

## Value

A `mosaic_analysis` object.
[`print()`](https://rdrr.io/r/base/print.html) shows the sample size,
the test, Cramer's V with its label and the dropped categories, and
returns the object invisibly;
[`summary()`](https://rdrr.io/r/base/summary.html) returns the per-cell
data frame (the `counts` table below) with the one-row test summary
attached as the attribute `"stats"`;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
mosaic on the active device, with new styling if given, and returns the
`ggplot` invisibly. The object is a list with:

- plot:

  The flat mosaic, a `ggplot`. Building the object does not draw it.

- counts:

  A data frame with one row per (`var1`, `var2`) cell: the two variables
  (named after them), `observed`, `expected`, `residual` (the
  standardised residual, or the nested test's residual, `NA` where it is
  undefined and drawn at 0) and `pct` (on `percentage_base`); `expected`
  and `residual` are rounded to 3 decimals, `pct` to 2.

- stats:

  A one-row data frame: `test`, `statistic`, `df`, `df2`, `p_value`,
  `cramers_v`, `effect_size`, `n`, `clusters`, `design_effect` and
  `icc`, rounded (`statistic` and `cramers_v` to 3 decimals, `p_value`
  to 4 significant digits). `statistic` and `df` are `NA` for Fisher's
  test. `df2` (the denominator degrees of freedom of the Rao-Scott F),
  `clusters`, `design_effect` and `icc` are `NA` where the test has
  none.

- test:

  The `htest` object of the test reported; `NULL` for a nested test.

- actor:

  The `actor` column names, or `NULL`.

- cramers_v, effect_size:

  Cramer's V, unrounded, and its label. V is computed from Pearson's
  chi-square without Yates' continuity correction, as Cramer defined it,
  so a perfect 2 x 2 association gives 1; the 2 x 2 test reported keeps
  the correction.

- table:

  The filtered contingency table.

- removed:

  The dropped categories, as a list with `var1` and `var2`.

- n_original, n_filtered:

  The number of complete rows before and after the sparse categories are
  dropped.

- vars:

  The two variable names, as `c(var1 = , var2 = )`.

- plot_parts, plot_args:

  The residual matrix, the table, the test line drawn as the subtitle
  with `report = TRUE`, and the styling, kept so
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) can redraw
  without testing again.

## Nested data

Pearson's test treats every row as an independent observation. When rows
are nested in actors (events within students, students within classes),
rows of one actor resemble each other and the test is too liberal. Set
`nested = TRUE` and name the actor column(s) in `actor` (as in
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)).
The test is then `nested_method`, one of five that allow for the
nesting; each gives a test and its own residual per cell, which fills
the mosaic.

- `"rao_scott"` (default):

  Rao-Scott second-order corrected chi-square, reported as an F
  statistic on fractional degrees of freedom (the default of
  `survey::svychisq()`, which it reproduces for an unweighted cluster
  sample). Works for any cluster sizes and when both variables vary
  within clusters. `design_effect` is the mean generalised design
  effect; 1 means the clustering costs nothing. Residual: the
  standardised residual divided by the square root of the cell's design
  effect (cluster-robust over simple-random-sampling variance of
  `p_ij - p_i. p_.j`); with one row per cluster it is
  [`chisq.test()`](https://rdrr.io/r/stats/chisq.test.html)'s, and it is
  `NA` where the cell's design effect is zero. With empty cells the
  design effects are those of the estimable contrasts. The test is
  liberal with few actors when one variable is constant within actors:
  in simulations of 3 x 4 tables its type I error at alpha = 0.05 was
  0.22 with 10 actors, 0.10 with 30 and 0.06 with 100.

- `"donner_klar"`:

  Pearson's chi-square divided, group by group, by the inflation factor
  `1 + (m - 1) * rho` (its weighted form for unequal cluster sizes),
  where `rho` is the ANOVA intracluster correlation of `var2`, pooled
  within the levels of the other variable. For designs where one
  variable is constant within each cluster, as in a cluster-randomised
  trial; otherwise an error is raised. `icc` is `rho`; a negative
  estimate is reported and used as 0. Residual: the standardised
  residual divided by the square root of its group's factor.

- `"cmh"`:

  Generalised Cochran-Mantel-Haenszel test with each cluster as a
  stratum: is there an association within clusters? Suits a few clusters
  with enough rows each; both variables must vary within clusters, and
  clusters with fewer than two rows are left out (`clusters` counts
  those that inform the test). Residual: summed observed minus expected
  over the square root of the summed hypergeometric variances. This
  tests conditional association, whereas Rao-Scott and GEE test marginal
  association; the null hypotheses can differ. Stratification does not
  correct arbitrary serial dependence among observations within an
  actor.

- `"gee"`:

  Marginal model for `var2` given `var1` with an identity link and
  independence working correlation; Wald chi-square with the
  cluster-robust (sandwich) covariance for equal distributions of `var2`
  across `var1`, on `(rows - 1) (columns - 1)` degrees of freedom. When
  the covariance has lower rank (fewer actors than degrees of freedom,
  or cells that never occur) the test would depend on the order of rows
  and columns, and an error of class `nestimate_bad_input` is raised.
  The Wald test is liberal unless the actors are many relative to the
  degrees of freedom: in simulations of 3 x 4 tables its type I error at
  alpha = 0.05 was 0.28 with 30 actors and 0.11 with 100. Residual: the
  share of `var2` within the level of `var1` minus the pooled share,
  over its cluster-robust standard error.

- `"glmm"`:

  Log-linear Poisson mixed model of the cluster-by-cell counts with
  random intercepts for cluster-by-`var1` and cluster-by-`var2` at every
  level of `actor`, so several columns give nested random effects
  (outermost first); Wald chi-square on the interaction. Fitted by
  `lme4::glmer(nAGQ = 0)`, a faster, less exact approximation than the
  default Laplace fit; needs lme4, and a failed fit raises an error of
  class `nestimate_bad_input`. With an empty cell the Wald test loses
  power (its interaction term runs to minus infinity). Residual:
  observed minus expected count under independence, over its
  cluster-robust standard error (the GLMM supplies the test, the
  residual is model-free). With several actor levels, this standard
  error clusters at the outermost level; it is undefined with fewer than
  two outer actors. With about 10 actors this standard error is small
  and 7 to 8% of null cells pass 1.96.

With several `actor` columns, the four other tests use the outermost
(the first), which carries the dependence of every level below it.
Cramer's V is the pooled-table effect size whatever the test: the
clustering changes the uncertainty, not the size of the association.

## Methods

- `plot.mosaic_analysis()` redraws the flat mosaic from the stored
  contingency table and residuals, so the styling can change without the
  test being run again. Any flat-mosaic styling argument (`tile_label`,
  `pct_base`, `col_label_side`, `legend_size`, ...) can be overridden
  through `...`.

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Lawrence Erlbaum.

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

Agresti, A. (2013). *Categorical data analysis* (3rd ed.), chapter 13.
Wiley.

## See also

[`mosaic_plot`](https://pak.dynasite.org/Nestimate/reference/mosaic_plot.md)
for the network/table mosaic (which also accepts `style = "flat"`);
[`mosaic_test`](https://pak.dynasite.org/Nestimate/reference/mosaic_test.md)
for every test, clustered or not, side by side in one table.

## Examples

``` r
data(group_regulation_long, package = "Nestimate")
res <- mosaic_analysis(group_regulation_long, "Course", "Action",
                       min_count = 20)
res
#> Mosaic analysis: Course x Action
#>   N = 27533 (filtered from 27533); table 3 x 9
#>   Chi-square: X2 = 233.560, df = 16, p < 2.2e-16
#>   Cramer's V = 0.065 (negligible)
head(summary(res))
#>   Course   Action observed expected residual  pct
#> 1      A    adapt      140  249.303   -9.430 0.51
#> 2      B    adapt      261  193.688    6.059 0.95
#> 3      C    adapt      153  111.009    4.502 0.56
#> 4      A cohesion      923  827.560    4.631 3.35
#> 5      B cohesion      590  642.945   -2.680 2.14
#> 6      C cohesion      326  368.495   -2.563 1.18
# \donttest{
plot(res, tile_label = "percent")


# Actions are nested in students: Rao-Scott test by student
nested <- mosaic_analysis(group_regulation_long, "Course", "Action",
                          nested = TRUE, actor = "Actor")
nested
#> Mosaic analysis: Course x Action
#>   N = 27533 (filtered from 27533); table 3 x 9
#>   Rao-Scott: F(13.45, 26887.5) = 14.99, p < 1e-16; 2000 actors, design effect 0.97
#>   Actor: Actor
#>   Cramer's V = 0.065 (negligible)
# }
```
