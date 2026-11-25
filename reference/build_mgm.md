# Estimate a mixed graphical model (MGM)

`build_network(method = "mgm")`, or for short `build_mgm()`, estimates a
mixed graphical model (MGM) of continuous and categorical variables. The
edge between two variables is the strength of their interaction, fitted
by penalised node-wise regression with the penalty chosen by EBIC
(Haslbeck & Waldorp, 2020).

Factor and character columns are categorical, and an integer column with
at most 10 values is categorical with a warning. Missing continuous
values are filled with the observed column mean and missing categorical
values with the modal level; this keeps the rows and is psychnets'
documented `na_method = "pairwise"` policy, not pairwise deletion.

## Usage

``` r
build_mgm(
  data,
  actor = NULL,
  level = NULL,
  lambdaGam = 0.25,
  ruleReg = c("AND", "OR"),
  metadata_cols = NULL,
  group = NULL,
  scaling = NULL,
  threshold = 0,
  labels = NULL,
  params = list(),
  ...
)
```

## Arguments

- data:

  The variables: a data frame with one row per observation and one
  column per variable, continuous or categorical.

- actor:

  The id column naming whose observation each row is (a person). With
  `level`, the decomposition groups by it; without, it is set aside, not
  modelled. Default `NULL`.

- level:

  With `actor`, a multilevel decomposition (`"between"`, `"within"` or
  `"both"`), as for
  [`build_cor()`](https://pak.dynasite.org/Nestimate/reference/build_cor.md).
  A numeric vector instead gives one category count per modelled
  variable (`1` for continuous ones); it is checked and kept for
  compatibility, and the fit reads the categories from the data. Default
  `NULL`.

- lambdaGam:

  The EBIC hyperparameter. Default `0.25`.

- ruleReg:

  `"AND"` (default) keeps an edge when both regressions select it,
  `"OR"` when either does.

- metadata_cols:

  Columns to keep as the network's metadata (covariates), never read as
  variables. Default `NULL`.

- group:

  A column to split the data by: one network per value, returned as a
  `netobject_group`. Default `NULL`.

- scaling:

  Rescaling applied to the weights after estimation, in the order given:
  `"minmax"` rescales the weights to run from 0 to 1, over every cell
  for the transition methods (tna's rule, so with no zero cell the
  smallest weight becomes 0) and over the non-zero weights, into (0,
  1\], for the other methods, `"max"` divides by the largest absolute
  weight, `"rank"` replaces the non-zero weights by their ranks, and
  `"normalize"` makes each row's absolute weights sum to 1; an
  undirected network whose matrix it makes asymmetric becomes directed,
  so both of a pair's values are kept as edges. Zero weights stay zero
  under every option. Combine them as in `c("rank", "minmax")`. Default
  `NULL` (none).

- threshold:

  Weights whose absolute value is below this are set to zero (default
  `0`); or `"LW"` or `"none"`, mgm's rule for thresholding the
  coefficients (mgm's default is `"LW"`).

- labels:

  New display labels for the nodes, applied after estimation: a
  two-column data frame (`name`, `label`), a named character vector
  `c(name = "label")` or a named list. It rewrites the node labels and
  the weight matrix's dimnames; nodes not listed keep their names.
  Default `NULL`.

- params:

  Settings for the estimator, as a named list. The network stores them.
  Default [`list()`](https://rdrr.io/r/base/list.html).

- ...:

  Further arguments for
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md).
  `method` cannot be given (error of class `nestimate_bad_input`).

## Value

A `netobject` (a `netobject_group` with `group`), as
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
returns it.

## References

Haslbeck, J. M. B., & Waldorp, L. J. (2020). mgm: Estimating
time-varying mixed graphical models in high-dimensional data. *Journal
of Statistical Software*, 93(8), 1–46.
[doi:10.18637/jss.v093.i08](https://doi.org/10.18637/jss.v093.i08)

## See also

[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
for every method;
[`build_ising()`](https://pak.dynasite.org/Nestimate/reference/build_ising.md)
for binary variables,
[`build_glasso()`](https://pak.dynasite.org/Nestimate/reference/build_glasso.md)
for continuous ones.

## Examples

``` r
if (requireNamespace("glmnet", quietly = TRUE)) {
  set.seed(1)
  mixed <- data.frame(x = rnorm(100), y = rnorm(100),
                      g = factor(sample(c("a", "b"), 100, TRUE)))
  build_network(mixed, method = "mgm")
  build_mgm(mixed)
}
#> Network (method: mgm) [undirected]
#> 
#>   No edges: every weight is zero.
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix.
```
