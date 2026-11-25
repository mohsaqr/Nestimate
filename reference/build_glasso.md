# Estimate an EBIC graphical lasso network

`build_network(method = "glasso")`, or for short `build_glasso()`,
estimates an EBIC graphical lasso network, a Gaussian graphical model,
from numeric variables. The edge between two variables is their partial
correlation given all the others, shrunk by the graphical lasso, with
the penalty chosen by EBIC (Friedman et al., 2008; Foygel & Drton,
2010).

The penalty is chosen along a path of `nlambda` values; an L1 penalty
sets small partial correlations to exactly zero. An absent edge means
the penalty removed it at the selected lambda, not that the two
variables are shown to be independent. Check the stability of the edges
and centralities with
[`boot_glasso()`](https://pak.dynasite.org/Nestimate/reference/boot_glasso.md).

## Usage

``` r
build_glasso(
  data,
  actor = NULL,
  level = NULL,
  ...,
  gamma = 0.5,
  nlambda = 100L,
  cor_method = c("pearson", "spearman", "kendall"),
  predictability = TRUE,
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

  The variables: a data frame with one row per observation and one
  column per variable, or a square symmetric correlation or covariance
  matrix, which needs the sample size in `params` (as in
  `params = list(n = 250)`). Non-numeric columns are set aside; a frame
  of two or more character or factor columns and fewer than two numeric
  ones is read as sequences and turned into per-row state counts first.
  A column named `rid` is set aside; columns with non-syntactic names,
  all-missing and constant columns, and rows with a missing value are
  dropped, each with a message; fewer than 2 variables or 3 complete
  rows raise an error.

- actor:

  The id column naming whose observation each row is (a person). With
  `level`, the decomposition groups by it; without, it is set aside, not
  modelled. Default `NULL`.

- level:

  Multilevel decomposition, with `actor`: `"between"` estimates the
  network of person means, `"within"` the network of deviations from
  them, and `"both"` returns the two as a `netobject_ml`. Default `NULL`
  (no decomposition).

- ...:

  Further arguments for
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md).
  `method` cannot be given (error of class `nestimate_bad_input`).

- gamma:

  The EBIC hyperparameter: larger values give sparser networks. Default
  `0.5`.

- nlambda:

  The number of penalties tried. Default `100`.

- cor_method:

  `"pearson"` (default), `"spearman"` or `"kendall"`, passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html).

- predictability:

  Whether to compute each node's predictability (the R-squared of the
  node on its neighbours). It is stored on the network, and
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
  draws it as a ring around each node. Default `TRUE`.

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

  Settings for the estimator, as a named list: `lambda.min.ratio`
  (default `0.01`), `penalize.diagonal` (default `FALSE`), `refit`
  (default `FALSE`; `TRUE` re-estimates the selected edges without the
  penalty), and `n`, the sample size when `data` is a correlation or
  covariance matrix. The network stores them, so the resampling verbs
  re-estimate with the same settings. Default
  [`list()`](https://rdrr.io/r/base/list.html).

## Value

A `netobject` (a `netobject_group` with `group`, a `netobject_ml` with
`level = "both"`), as
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
returns it.

## References

Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse
covariance estimation with the graphical lasso. *Biostatistics*, 9(3),
432–441.
[doi:10.1093/biostatistics/kxm045](https://doi.org/10.1093/biostatistics/kxm045)

Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria
for Gaussian graphical models. *Advances in Neural Information
Processing Systems*, 23, 604–612.

## See also

[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
for every method;
[`build_pcor()`](https://pak.dynasite.org/Nestimate/reference/build_pcor.md)
for the unregularised network;
[`boot_glasso()`](https://pak.dynasite.org/Nestimate/reference/boot_glasso.md)
and [`nct()`](https://pak.dynasite.org/Nestimate/reference/nct.md) for
inference.

## Examples

``` r
build_network(srl_strategies, method = "glasso")
#> Partial Correlation Network (EBICglasso) [undirected]
#>   Sample size: 250
#>   Weights: [0.020, 0.325]  |  +13 / -0 edges
#> 
#>   Strongest edges (10 of 13):
#>    from         to           weight
#>    Time_Mgmt    Effort_Reg   0.325 
#>    Monitoring   Evaluating   0.300 
#>    Elaboration  Organization 0.262 
#>    Planning     Monitoring   0.229 
#>    Elaboration  Rehearsal    0.158 
#>    Evaluating   Organization 0.154 
#>    Organization Rehearsal    0.152 
#>    Planning     Time_Mgmt    0.130 
#>    Planning     Evaluating   0.130 
#>    Evaluating   Effort_Reg   0.093  
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix.
#> 
#>   Predictability (R²):
#>   Monitoring    0.173  ████████████████████████████████████████
#>   Evaluating    0.173  ████████████████████████████████████████
#>   Organization  0.137  ████████████████████████████████
#>   Time_Mgmt     0.133  ███████████████████████████████
#>   Effort_Reg    0.122  ████████████████████████████
#>   Planning      0.119  ████████████████████████████
#>   Elaboration   0.115  ██████████████████████████
#>   Rehearsal     0.067  ███████████████
#>   Help_Seeking  0.006  █
#> 
#>   Gamma: 0.50  |  Lambda: 0.1319
build_glasso(srl_strategies)
#> Partial Correlation Network (EBICglasso) [undirected]
#>   Sample size: 250
#>   Weights: [0.020, 0.325]  |  +13 / -0 edges
#> 
#>   Strongest edges (10 of 13):
#>    from         to           weight
#>    Time_Mgmt    Effort_Reg   0.325 
#>    Monitoring   Evaluating   0.300 
#>    Elaboration  Organization 0.262 
#>    Planning     Monitoring   0.229 
#>    Elaboration  Rehearsal    0.158 
#>    Evaluating   Organization 0.154 
#>    Organization Rehearsal    0.152 
#>    Planning     Time_Mgmt    0.130 
#>    Planning     Evaluating   0.130 
#>    Evaluating   Effort_Reg   0.093  
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix.
#> 
#>   Predictability (R²):
#>   Monitoring    0.173  ████████████████████████████████████████
#>   Evaluating    0.173  ████████████████████████████████████████
#>   Organization  0.137  ████████████████████████████████
#>   Time_Mgmt     0.133  ███████████████████████████████
#>   Effort_Reg    0.122  ████████████████████████████
#>   Planning      0.119  ████████████████████████████
#>   Elaboration   0.115  ██████████████████████████
#>   Rehearsal     0.067  ███████████████
#>   Help_Seeking  0.006  █
#> 
#>   Gamma: 0.50  |  Lambda: 0.1319

build_glasso(srl_strategies, gamma = 0.75)
#> Partial Correlation Network (EBICglasso) [undirected]
#>   Sample size: 250
#>   Weights: [0.020, 0.325]  |  +13 / -0 edges
#> 
#>   Strongest edges (10 of 13):
#>    from         to           weight
#>    Time_Mgmt    Effort_Reg   0.325 
#>    Monitoring   Evaluating   0.300 
#>    Elaboration  Organization 0.262 
#>    Planning     Monitoring   0.229 
#>    Elaboration  Rehearsal    0.158 
#>    Evaluating   Organization 0.154 
#>    Organization Rehearsal    0.152 
#>    Planning     Time_Mgmt    0.130 
#>    Planning     Evaluating   0.130 
#>    Evaluating   Effort_Reg   0.093  
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix.
#> 
#>   Predictability (R²):
#>   Monitoring    0.173  ████████████████████████████████████████
#>   Evaluating    0.173  ████████████████████████████████████████
#>   Organization  0.137  ████████████████████████████████
#>   Time_Mgmt     0.133  ███████████████████████████████
#>   Effort_Reg    0.122  ████████████████████████████
#>   Planning      0.119  ████████████████████████████
#>   Elaboration   0.115  ██████████████████████████
#>   Rehearsal     0.067  ███████████████
#>   Help_Seeking  0.006  █
#> 
#>   Gamma: 0.75  |  Lambda: 0.1319
```
