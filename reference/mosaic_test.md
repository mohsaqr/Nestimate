# Test the association of a mosaic with and without clustering

`mosaic_test()` runs the tests behind
[`mosaic_plot()`](https://pak.dynasite.org/Nestimate/reference/mosaic_plot.md)
and
[`mosaic_analysis()`](https://pak.dynasite.org/Nestimate/reference/mosaic_analysis.md)
and returns their results as a data frame. Give it the same object you
would draw: a network, a `netobject_group`, a three-way table, or a data
frame with two variables. By default it runs Pearson's chi-square, which
treats every count as independent, next to the five tests that allow for
counts nested in actors (sequences, students, teams), so the table shows
what the nesting changes.

The tests are those of *Nested data* in
[`mosaic_analysis()`](https://pak.dynasite.org/Nestimate/reference/mosaic_analysis.md):
`"rao_scott"` (Rao-Scott corrected F), `"donner_klar"` (Pearson divided
by the cluster-size inflation factor), `"cmh"` (generalised
Cochran-Mantel-Haenszel with each actor as a stratum), `"gee"` (Wald
test of a marginal model with cluster-robust covariance) and `"glmm"`
(Wald test of a log-linear mixed model, needs lme4). `"chisq"` is
Pearson's test on the pooled table. As in the mosaics, a network's
source state and a group are the predictors.

Not every test fits every design: Donner-Klar needs one variable that is
constant within actors (the group of a `netobject_group`), and CMH needs
both to vary within actors (the states of a network). When several
methods are asked for, a test that cannot run on the design gives a row
of `NA` with the reason in `note`; asked for alone, it raises the error.

## Usage

``` r
mosaic_test(x, ...)

# Default S3 method
mosaic_test(x, ...)

# S3 method for class 'netobject'
mosaic_test(
  x,
  method = c("chisq", "rao_scott", "donner_klar", "cmh", "gee", "glmm"),
  actor = NULL,
  level = c("test", "cells"),
  ...
)

# S3 method for class 'htna'
mosaic_test(
  x,
  method = c("chisq", "rao_scott", "donner_klar", "cmh", "gee", "glmm"),
  actor = NULL,
  level = c("test", "cells"),
  ...
)

# S3 method for class 'netobject_group'
mosaic_test(
  x,
  method = c("chisq", "rao_scott", "donner_klar", "cmh", "gee", "glmm"),
  actor = NULL,
  level = c("test", "cells"),
  ...
)

# S3 method for class 'array'
mosaic_test(x, ...)

# S3 method for class 'table'
mosaic_test(
  x,
  method = c("chisq", "rao_scott", "donner_klar", "cmh", "gee", "glmm"),
  level = c("test", "cells"),
  ...
)

# S3 method for class 'data.frame'
mosaic_test(
  x,
  var1,
  var2,
  actor,
  method = c("chisq", "rao_scott", "donner_klar", "cmh", "gee", "glmm"),
  level = c("test", "cells"),
  ...
)
```

## Arguments

- x:

  A `netobject` or `htna` network (relative or frequency, with its
  sequences), a `netobject_group`, a three-way `table` or array (rows x
  columns x actors), or a data frame (with `var1`, `var2` and `actor`).
  A co-occurrence network raises an error of class
  `nestimate_mosaic_undirected`: it counts each pair of states in both
  mirrored cells and the pairs of a sequence share its positions, so its
  counts are not independent and the tests would not hold their level
  (under independence Pearson's test rejected 200 of 200 simulated
  tables at 0.05).
  [`mosaic_plot()`](https://pak.dynasite.org/Nestimate/reference/mosaic_plot.md)
  still draws it. A `tna` model or a `group_tna` (package tna) is first
  rebuilt as a network from its own sequences (see
  [`as_netobject()`](https://pak.dynasite.org/Nestimate/reference/as_netobject.md)).
  `actor` is read from the metadata the model keeps
  (`attr(model$data, "meta")`, kept when the model is built from
  [`tna::prepare_data()`](https://sonsoles.me/tna/reference/prepare_data.html));
  a model without it raises `nestimate_actor_missing`.

- ...:

  Not used; a name given here raises an error of class
  `nestimate_bad_input`.

- method:

  The tests to run, any of `"chisq"`, `"rao_scott"`, `"donner_klar"`,
  `"cmh"`, `"gee"` and `"glmm"`. Default: all six.

- actor:

  The actor column(s), as in
  [`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md):
  for a network or `netobject_group`, `NULL` (default) makes every
  sequence its own actor, and the name of a metadata column (such as a
  team) pools the sequences that share its value; for a data frame, the
  column(s) that identify the actors (required). Several names give
  nested levels, outermost first: `"glmm"` models each level, the other
  tests use the outermost. Not used for a three-way table, whose third
  dimension is the actor.

- level:

  `"test"` (default) for one row per test, `"cells"` for one row per
  cell per test with the cell's residual.

- var1, var2:

  For a data frame: the names of the two variables; `var1` is the
  predictor.

## Value

A data frame. With `level = "test"`, one row per test: `method`, `test`
(its name), `statistic` (F for Rao-Scott, chi-square otherwise), `df`,
`df2` (the Rao-Scott denominator degrees of freedom), `p_value`, `n`
(counts in the table), `clusters` (actors the test used),
`design_effect` (Rao-Scott's mean generalised design effect, or
Donner-Klar's mean inflation factor), `icc` (Donner-Klar's intracluster
correlation) and `note` (why a test did not run, the caution on the GEE
Wald test, else `NA`). With `level = "cells"`, one row per cell and
test: `method`, the two variables (`from`/`to` for a network,
`group`/`state` for a group, the variable names for a data frame,
`row`/`column` for a table), `observed` and `residual` (the residual the
mosaic draws for that method; `NA` where undefined). An unsupported
input raises an error of class `nestimate_bad_input`.

## References

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

## See also

[`mosaic_plot()`](https://pak.dynasite.org/Nestimate/reference/mosaic_plot.md),
[`mosaic_analysis()`](https://pak.dynasite.org/Nestimate/reference/mosaic_analysis.md)

## Examples

``` r
data(group_regulation_long, package = "Nestimate")
net <- build_network(group_regulation_long, method = "frequency",
                     format = "long", actor = "Actor", action = "Action",
                     order = "Time")
mosaic_test(net, method = c("chisq", "rao_scott"), actor = "Group")
#>      method      test  statistic       df      df2 p_value     n clusters
#> 1     chisq   Pearson 15435.2605 64.00000       NA       0 25533       NA
#> 2 rao_scott Rao-Scott   227.7646 39.41762 7844.106       0 25533      200
#>   design_effect icc note
#> 1            NA  NA <NA>
#> 2      1.058882  NA <NA>
# \donttest{
# Several tests, sequences nested in teams (the default runs them all,
# the mixed model included)
mosaic_test(net, method = c("rao_scott", "donner_klar", "cmh", "gee"),
            actor = "Group")
#>        method                    test  statistic       df      df2 p_value
#> 1   rao_scott               Rao-Scott   227.7646 39.41762 7844.106       0
#> 2 donner_klar             Donner-Klar         NA       NA       NA      NA
#> 3         cmh Cochran-Mantel-Haenszel 15167.3307 64.00000       NA       0
#> 4         gee                GEE Wald         NA       NA       NA      NA
#>       n clusters design_effect icc
#> 1 25533      200      1.058882  NA
#> 2 25533       NA            NA  NA
#> 3 25533      200            NA  NA
#> 4 25533       NA            NA  NA
#>                                                                                          note
#> 1                                                                                        <NA>
#> 2                        Donner-Klar needs one variable that is constant within every cluster
#> 3                                                                                        <NA>
#> 4 GEE: the cluster-robust covariance has rank 63, below the 64 degrees of freedom of the test

# Cell residuals of one test
head(mosaic_test(net, method = "rao_scott", actor = "Group",
                 level = "cells"))
#>      method  from         to observed   residual
#> 1 rao_scott adapt      adapt        0 -1.9966253
#> 2 rao_scott adapt   cohesion      139 17.4827797
#> 3 rao_scott adapt  consensus      243 11.5347947
#> 4 rao_scott adapt coregulate       11 -4.6447682
#> 5 rao_scott adapt    discuss       30 -6.2361428
#> 6 rao_scott adapt    emotion       61  0.8545643
# }
```
