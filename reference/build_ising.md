# Estimate an Ising model

`build_network(method = "ising")`, or for short `build_ising()`,
estimates an Ising model of binary variables. The edge between two
variables is their interaction, fitted by L1-penalised logistic
regression of each variable on the others, with the penalty chosen by
EBIC, and the two estimates of each edge combined into one undirected
weight (van Borkulo et al., 2014). It needs the glmnet package.

An edge is a conditional association between two variables given all the
others; an absent edge means the penalty removed it.

## Usage

``` r
build_ising(
  data,
  actor = NULL,
  ...,
  gamma = 0.25,
  nlambda = 100L,
  rule = c("AND", "OR"),
  metadata_cols = NULL,
  group = NULL,
  scaling = NULL,
  threshold = 0,
  labels = NULL,
  params = list()
)
```

## Arguments

- data:

  The binary variables: a data frame with one row per observation and
  one column per variable, holding only 0 and 1. Non-numeric columns and
  a column named `rid` are set aside; columns with non-syntactic names,
  all-missing and constant columns, and rows with a missing value are
  dropped, each with a message; fewer than 2 variables or 3 complete
  rows, or any value other than 0 and 1, raise an error.

- actor:

  An id column naming whose observation each row is; it is set aside,
  not modelled. Default `NULL`.

- ...:

  Further arguments for
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md).
  `method` cannot be given (error of class `nestimate_bad_input`).

- gamma:

  The EBIC hyperparameter. Default `0.25`.

- nlambda:

  The number of penalties tried. Default `100`.

- rule:

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

  Weights whose absolute value is below this are set to zero. Default
  `0` (none).

- labels:

  New display labels for the nodes, applied after estimation: a
  two-column data frame (`name`, `label`), a named character vector
  `c(name = "label")` or a named list. It rewrites the node labels and
  the weight matrix's dimnames; nodes not listed keep their names.
  Default `NULL`.

- params:

  Settings for the estimator, as a named list. The network stores them,
  so the resampling verbs re-estimate with the same settings. Default
  [`list()`](https://rdrr.io/r/base/list.html).

## Value

A `netobject` (a `netobject_group` with `group`), as
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
returns it.

## References

van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F., Boschloo,
L., Schoevers, R. A., & Waldorp, L. J. (2014). A new method for
constructing networks from binary data. *Scientific Reports*, 4, 5918.
[doi:10.1038/srep05918](https://doi.org/10.1038/srep05918)

## See also

[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
for every method;
[`build_mgm()`](https://pak.dynasite.org/Nestimate/reference/build_mgm.md)
for mixed continuous and categorical variables.

## Examples

``` r
if (requireNamespace("glmnet", quietly = TRUE)) {
  set.seed(1)
  binary <- data.frame(matrix(rbinom(200, 1, 0.5), ncol = 5))
  build_network(binary, method = "ising")
  build_ising(binary, rule = "OR")
}
#> Ising Model Network [undirected]
#>   Sample size: 40
#> 
#>   No edges: every weight is zero.
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix.
#> 
#>   Gamma: 0.25  |  Rule: OR
#>   Thresholds: [-0.201, 0.405]
```
