# Estimate a correlation network

`build_network(method = "cor")`, or for short `build_cor()`, estimates a
correlation network of numeric variables. The edge between two variables
is their correlation (Epskamp & Fried, 2018).

A correlation edge is a marginal association: it can come from a shared
link with a third variable. To separate direct from indirect links, use
[`build_pcor()`](https://pak.dynasite.org/Nestimate/reference/build_pcor.md)
or
[`build_glasso()`](https://pak.dynasite.org/Nestimate/reference/build_glasso.md).
Check the stability of the edges with
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md).

## Usage

``` r
build_cor(
  data,
  actor = NULL,
  level = NULL,
  ...,
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

  Settings for the estimator, as a named list, such as `n`, the sample
  size when `data` is a correlation or covariance matrix. The network
  stores them, so
  [`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  re-estimates with the same settings. Default
  [`list()`](https://rdrr.io/r/base/list.html).

## Value

A `netobject` (a `netobject_group` with `group`, a `netobject_ml` with
`level = "both"`), as
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
returns it: print it for an overview, read the edges with
[`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md),
and draw it with
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html).

## References

Epskamp, S., & Fried, E. I. (2018). A tutorial on regularized partial
correlation networks. *Psychological Methods*, 23(4), 617–634.
[doi:10.1037/met0000167](https://doi.org/10.1037/met0000167)

Saqr, M., Beck, E., & López-Pernas, S. (2024). Psychological networks: A
modern approach to analysis of learning and complex learning processes.
In M. Saqr & S. López-Pernas (Eds.), *Learning Analytics Methods and
Tutorials: A Practical Guide Using R*. Springer.
[doi:10.1007/978-3-031-54464-4_19](https://doi.org/10.1007/978-3-031-54464-4_19)

## See also

[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
for every method;
[`build_pcor()`](https://pak.dynasite.org/Nestimate/reference/build_pcor.md),
[`build_glasso()`](https://pak.dynasite.org/Nestimate/reference/build_glasso.md)
and
[`build_ising()`](https://pak.dynasite.org/Nestimate/reference/build_ising.md)
for the other association networks;
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
and [`nct()`](https://pak.dynasite.org/Nestimate/reference/nct.md) for
inference.

## Examples

``` r
build_network(srl_strategies, method = "cor")
#> Correlation Network [undirected]
#>   Sample size: 250
#>   Weights: [-0.130, 0.485]  |  +26 / -10 edges
#> 
#>   Strongest edges (10 of 36):
#>    from         to           weight
#>    Monitoring   Evaluating   0.485 
#>    Time_Mgmt    Effort_Reg   0.467 
#>    Elaboration  Organization 0.432 
#>    Planning     Monitoring   0.423 
#>    Planning     Evaluating   0.358 
#>    Elaboration  Rehearsal    0.341 
#>    Organization Rehearsal    0.339 
#>    Evaluating   Organization 0.313 
#>    Planning     Time_Mgmt    0.285 
#>    Evaluating   Effort_Reg   0.250  
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix.
#> 
#>   Predictability (R²):
#>   Evaluating    0.456  ████████████████████████████████████████
#>   Monitoring    0.393  ███████████████████████████████████
#>   Organization  0.387  ██████████████████████████████████
#>   Time_Mgmt     0.377  █████████████████████████████████
#>   Planning      0.372  █████████████████████████████████
#>   Effort_Reg    0.336  █████████████████████████████
#>   Elaboration   0.293  ██████████████████████████
#>   Rehearsal     0.212  ███████████████████
#>   Help_Seeking  0.110  ██████████
build_cor(srl_strategies)
#> Correlation Network [undirected]
#>   Sample size: 250
#>   Weights: [-0.130, 0.485]  |  +26 / -10 edges
#> 
#>   Strongest edges (10 of 36):
#>    from         to           weight
#>    Monitoring   Evaluating   0.485 
#>    Time_Mgmt    Effort_Reg   0.467 
#>    Elaboration  Organization 0.432 
#>    Planning     Monitoring   0.423 
#>    Planning     Evaluating   0.358 
#>    Elaboration  Rehearsal    0.341 
#>    Organization Rehearsal    0.339 
#>    Evaluating   Organization 0.313 
#>    Planning     Time_Mgmt    0.285 
#>    Evaluating   Effort_Reg   0.250  
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix.
#> 
#>   Predictability (R²):
#>   Evaluating    0.456  ████████████████████████████████████████
#>   Monitoring    0.393  ███████████████████████████████████
#>   Organization  0.387  ██████████████████████████████████
#>   Time_Mgmt     0.377  █████████████████████████████████
#>   Planning      0.372  █████████████████████████████████
#>   Effort_Reg    0.336  █████████████████████████████
#>   Elaboration   0.293  ██████████████████████████
#>   Rehearsal     0.212  ███████████████████
#>   Help_Seeking  0.110  ██████████

build_cor(srl_strategies, cor_method = "spearman", threshold = 0.2)
#> Correlation Network [undirected]
#>   Sample size: 250
#>   Weights: [0.209, 0.476]  |  +11 / -0 edges
#> 
#>   Strongest edges (10 of 11):
#>    from         to           weight
#>    Monitoring   Evaluating   0.476 
#>    Time_Mgmt    Effort_Reg   0.455 
#>    Elaboration  Organization 0.432 
#>    Planning     Monitoring   0.420 
#>    Organization Rehearsal    0.353 
#>    Elaboration  Rehearsal    0.351 
#>    Planning     Evaluating   0.340 
#>    Evaluating   Organization 0.322 
#>    Planning     Time_Mgmt    0.278 
#>    Evaluating   Effort_Reg   0.271  
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix.
#> 
#>   Predictability (R²):
#>   Evaluating    0.411  ████████████████████████████████████████
#>   Time_Mgmt     0.312  ██████████████████████████████
#>   Organization  0.308  ██████████████████████████████
#>   Monitoring    0.302  █████████████████████████████
#>   Planning      0.269  ██████████████████████████
#>   Effort_Reg    0.260  █████████████████████████
#>   Elaboration   0.232  ███████████████████████
#>   Rehearsal     0.173  █████████████████
#>   Help_Seeking  0.044  ████
#>   Threshold: 0.2
```
