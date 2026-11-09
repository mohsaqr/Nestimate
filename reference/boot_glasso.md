# Bootstrap for Regularised Partial Correlation Networks

`boot_glasso()` assesses how accurately an EBIC graphical lasso network
is estimated and how stable its centralities are, by resampling the
observations (Epskamp, Borsboom & Fried, 2018). One call runs a
nonparametric bootstrap, which gives percentile intervals and inclusion
probabilities for every edge, intervals for each node's centralities and
predictability, and bootstrap tests of the difference between every pair
of edges and every pair of nodes; and a case-dropping bootstrap, which
gives each centrality measure's correlation-stability (CS) coefficient.
It takes the raw data or a network built with `method = "glasso"`.

An interval shows how precisely an edge weight is estimated: a wide one
means the weight changes a good deal between samples. Epskamp et al.
(2018) caution that, because the lasso shrinks estimates towards zero,
these intervals are not tests of whether an edge differs from zero; the
thresholded network (edges whose interval excludes zero) and the
"significant" count in [`print()`](https://rdrr.io/r/base/print.html)
read them that way, and the inclusion probability (how often an edge is
non-zero across resamples) is a separate summary. The CS-coefficient is
the largest proportion of cases that can be dropped while more than 95%
of the subsamples keep a correlation above 0.7 with the original
centralities; Epskamp et al. (2018) suggest it should not be below 0.25
and preferably above 0.5 before the centrality order is interpreted, and
[`print()`](https://rdrr.io/r/base/print.html) labels the values
accordingly. The difference tests are not corrected for the number of
pairs compared.

## Usage

``` r
boot_glasso(
  x,
  iter = 1000L,
  cs_iter = 500L,
  cs_drop = seq(0.1, 0.9, by = 0.1),
  alpha = 0.05,
  gamma = 0.5,
  nlambda = 100L,
  centrality = c("strength", "expected_influence", "betweenness", "closeness"),
  centrality_fn = NULL,
  cor_method = "pearson",
  ncores = 1L,
  seed = NULL
)

# S3 method for class 'boot_glasso'
print(x, ...)

# S3 method for class 'boot_glasso'
summary(object, type = "edges", ...)

# S3 method for class 'boot_glasso'
plot(x, type = "edges", measure = NULL, ...)
```

## Arguments

- x:

  The data: a data frame or numeric matrix with one row per observation
  and one column per variable, cleaned as in
  [`build_glasso()`](https://pak.dynasite.org/Nestimate/reference/build_glasso.md):
  non-numeric columns and a column named `rid` are set aside; columns
  with non-syntactic names, all-missing and constant columns, and rows
  with a missing value are dropped, each with a message; fewer than 2
  variables or 3 complete rows, or no more rows than variables, raise an
  error. A square symmetric matrix with matching row and column names is
  read as a correlation or covariance matrix and refused (class
  `nestimate_bad_input`): it has no observations to resample. Or a
  `netobject` (or `cograph_network`) built with `method = "glasso"`,
  whose data are resampled; any `gamma`, `nlambda` or `cor_method` set
  in its `params` replaces the arguments of the same name, and its
  `lambda.min.ratio`, `penalize.diagonal` and `refit` are used too, so
  `original_pcor` is the network's weights. Any other network raises an
  error. A correlation or covariance matrix is not accepted: a matrix is
  always read as data, its rows as observations. For the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: a
  `boot_glasso`.

- iter:

  The number of nonparametric bootstrap resamples, at least 2. Default
  `1000`.

- cs_iter:

  The total number of case-dropping iterations, at least 1. Default
  `500`. As in bootnet, each iteration draws one proportion from
  `cs_drop` at random, so the iterations are spread across the
  proportions (about `cs_iter / length(cs_drop)` each) rather than
  repeated `cs_iter` times at each one.

- cs_drop:

  The proportions of cases to drop, each in (0, 1). Default
  `seq(0.1, 0.9, by = 0.1)`.

- alpha:

  The significance level: the intervals run from the `alpha / 2` to the
  `1 - alpha / 2` quantile, and
  [`print()`](https://rdrr.io/r/base/print.html) and the `"edge_diff"` /
  `"centrality_diff"` plots mark the difference tests' p-values below it
  as differing. Default `0.05`.

- gamma:

  The EBIC hyperparameter; larger values give sparser networks. Default
  `0.5`. Replaced by the value stored in a `netobject`.

- nlambda:

  The number of lambda values in the regularisation path, at least 2.
  Default `100`. Replaced by the value stored in a `netobject`.

- centrality:

  The centrality measures to compute. Default all four built-in ones:
  `"strength"` (sum of absolute edge weights), `"expected_influence"`
  (sum of signed edge weights), `"betweenness"` and `"closeness"` (on
  shortest paths with distance 1 / \|weight\|). Without `centrality_fn`,
  names may be abbreviated. Any other name needs `centrality_fn`, which
  then computes it; the built-in measures are always computed
  internally.

- centrality_fn:

  A function for measures beyond the four built-in ones: it is called as
  `centrality_fn(mat)` on the partial correlation matrix and returns a
  named list of centrality vectors, such as `list(my_metric = ...)`.
  Default `NULL`; a measure in `centrality` that is neither built in nor
  returned by it raises an error.

- cor_method:

  The correlation method: `"pearson"` (default), `"spearman"` or
  `"kendall"`. Replaced by the value stored in a `netobject`.

- ncores:

  The number of cores for the resampling, used through
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html),
  which forks the R process; on Windows it accepts only 1, and any
  larger value raises an error. Default `1` (sequential).

- seed:

  A seed for reproducible resampling; the session's random state is
  restored when the call returns. Default `NULL` uses the current random
  stream. With `ncores > 1` the resamples are reproducible only under
  `RNGkind("L'Ecuyer-CMRG")`.

- ...:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method with `type = "edge_diff"` or `"centrality_diff"`: `order`,
  `"sample"` (default) to sort the edges by their absolute weight or the
  nodes by their centrality, or `"id"` to sort them alphabetically; with
  `type = "edge_diff"` also `pos_color` and `neg_color`, the fills for a
  cell where the edge on the x-axis has a larger or a smaller mean
  bootstrap weight than the edge on the y-axis (defaults `"#C0392B"` and
  `"#2C6E8A"`). Ignored otherwise, and by the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `boot_glasso`.

- type:

  What to return or draw. For the
  [`summary()`](https://rdrr.io/r/base/summary.html) method: `"edges"`
  (default), the edge intervals sorted by decreasing absolute weight;
  `"centrality"`, the centrality intervals; `"cs"`, the case-dropping
  table; `"predictability"`, the predictability intervals; or `"all"`,
  the four together. For the
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method:
  `"edges"` (default), each edge's weight and interval; `"stability"`,
  the mean correlation with the original centralities against the
  proportion dropped, with a dashed line at 0.7; `"edge_diff"`, a
  heatmap of the mean bootstrap difference between every pair of edges,
  faint where the pair does not differ at `alpha` and starred where it
  does, with the edge weights on the diagonal; `"centrality_diff"`, a
  heatmap of which node pairs differ on `measure`, with the node values
  on the diagonal; or `"inclusion"`, each edge's inclusion probability,
  with a dashed line at 0.5.

- measure:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method with `type = "centrality_diff"`: the centrality measure to
  compare, one of those computed (names may be abbreviated). Default
  `NULL`, the first measure computed.

## Value

A `boot_glasso` object. [`print()`](https://rdrr.io/r/base/print.html)
reports the data size, `alpha` and `gamma`, the number of edges whose
interval excludes zero, the mean inclusion probability, the failed
resamples when there are any, each measure's CS-coefficient with its
label, the number of edge pairs that differ, and the timing;
[`summary()`](https://rdrr.io/r/base/summary.html) returns the table
chosen by `type` (edges, centralities, case-dropping results or
predictability; described at the end of this section);
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws one
figure chosen by `type`. The object's fields:

- original_pcor:

  The partial correlation matrix of the original data.

- original_precision:

  The precision matrix of the original data.

- original_centrality:

  A named list with one vector of node values per centrality measure.

- original_predictability:

  A named vector of each node's R-squared.

- edge_ci:

  A data frame with one row per edge: `edge`, `weight`, `ci_lower`,
  `ci_upper` and `inclusion`.

- edge_inclusion:

  A named vector of edge inclusion probabilities.

- thresholded_pcor:

  The partial correlation matrix with every edge whose interval includes
  zero set to zero.

- centrality_ci:

  A named list with one data frame per centrality measure: `node`,
  `value`, `ci_lower`, `ci_upper`.

- cs_coefficient:

  A named vector with the CS-coefficient of each centrality measure.

- cs_data:

  A data frame of the case-dropping results, one row per drop proportion
  and measure: `drop_prop`, `measure`, `mean_cor`, `prop_above` (the
  share of that proportion's iterations correlating above 0.7) and
  `n_samples` (the iterations that drew that proportion).

- edge_diff_p:

  A symmetric matrix of the edge-difference p-values; `NULL` when the
  network has more than 500 edges.

- centrality_diff_p:

  A named list with one symmetric matrix of node-difference p-values per
  centrality measure.

- predictability_ci:

  A data frame with one row per node: `node`, `r2`, `ci_lower`,
  `ci_upper`.

- boot_edges:

  An `iter` x edges matrix of the bootstrap edge weights.

- boot_centrality:

  A named list with one `iter` x nodes matrix of bootstrap values per
  centrality measure.

- boot_predictability:

  An `iter` x nodes matrix of bootstrap R-squared values.

- nodes:

  The node names.

- n:

  The sample size.

- p:

  The number of variables.

- iter, cs_iter, cs_drop, alpha, gamma, nlambda:

  The settings used.

- centrality_measures:

  The centrality measures computed.

- cor_method:

  The correlation method used.

- lambda_path:

  The lambda path used for every fit.

- lambda_selected:

  The lambda selected on the original data.

- n_failed, cs_n_failed:

  The numbers of bootstrap and case-dropping resamples whose estimation
  failed and that were left out.

- timing:

  A named vector of elapsed seconds: `total`, `phase1` (the original
  fit), `bootstrap`, `case_drop` and `statistics`.

[`print()`](https://rdrr.io/r/base/print.html) returns its input
invisibly. [`summary()`](https://rdrr.io/r/base/summary.html) returns a
data frame: for `type = "edges"` the `edge_ci` table (`edge`, `weight`,
`ci_lower`, `ci_upper`, `inclusion`) sorted by decreasing absolute
weight, for `"cs"` the `cs_data` table and for `"predictability"` the
`predictability_ci` table; `"centrality"` returns a named list with one
data frame per measure (`node`, `value`, `ci_lower`, `ci_upper`), and
`"all"` a named list of the four (`edges`, `centrality`, `cs`,
`predictability`).
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns a
`ggplot` object, which is drawn when the call is made at the top level.
The edge-difference plot raises an error when the test was skipped (more
than 500 edges).

## Details

The network is estimated as
[`build_glasso()`](https://pak.dynasite.org/Nestimate/reference/build_glasso.md)
estimates it: with a `netobject`, with its own `params`
(`lambda.min.ratio`, `penalize.diagonal`, `refit`); with data, with that
function's defaults (smallest lambda 0.01 times the largest, no penalty
on the diagonal, the selected penalty re-solved exactly, no
unregularised refit). The lambda path is set once from the original data
and reused for every resample, which is fitted with the same settings;
the penalty is chosen on it by EBIC with `gamma`. Strength is the sum of
a node's absolute edge weights and expected influence the sum of its
signed weights; betweenness and closeness use shortest paths with
distance 1 / \|weight\|, computed without igraph. Predictability is
\\1 - 1 / \Theta\_{ii}\\, from the diagonal of the precision matrix,
bounded to 0–1.

The case-dropping bootstrap follows bootnet: each of the `cs_iter`
iterations draws one proportion from `cs_drop` at random, keeps that
share of the rows (at least `p + 1` of them) without replacement, and
correlates the subsample's centralities with the original ones; a
centrality that is constant across nodes counts as a correlation of 0.
The difference tests give, for each pair, the two-sided bootstrap
p-value `2 * min(P(d >= 0), P(d <= 0))` of the difference `d` between
the two edges or nodes, capped at 1, with ties counted in both tails.
The edge-difference test is skipped for networks with more than 500
edges. A resample whose estimation fails (for example, a variable
constant in the resample) is kept as missing and left out of the
intervals and the case-dropping results; the failures are counted in
`n_failed` and `cs_n_failed`, shown by
[`print()`](https://rdrr.io/r/base/print.html), and reported in a
warning of class `nestimate_resample_failed`. The function was written
as a faster, single-call alternative to running bootnet several times.

## Methods

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws one
  figure of the bootstrap results per call, chosen with `type`.

## References

Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
psychological networks and their accuracy: A tutorial paper. *Behavior
Research Methods* 50(1), 195-212.
[doi:10.3758/s13428-017-0862-1](https://doi.org/10.3758/s13428-017-0862-1)
(source of the CS-coefficient and of the case-dropping, edge-difference
and centrality-difference procedures reproduced here.)

## See also

[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
[`bootstrap_network`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)

## Examples

``` r
set.seed(1)
dat <- as.data.frame(matrix(rnorm(60), ncol = 3))
net <- build_network(dat, method = "glasso")
bg <- boot_glasso(net, iter = 10, cs_iter = 5, centrality = "strength")
# \donttest{
set.seed(42)
mat <- matrix(rnorm(60), ncol = 4)
colnames(mat) <- LETTERS[1:4]
net <- build_network(as.data.frame(mat), method = "glasso")
# iter = 20 keeps the example fast; a real analysis uses 1000 or more.
boot <- boot_glasso(net, iter = 20, cs_iter = 10, seed = 42,
  centrality = c("strength", "expected_influence"))
print(boot)
#> GLASSO Bootstrap (20 iterations, 10 case-drop iterations)
#>   Data: 15 x 4  |  Alpha: 0.05  |  Gamma: 0.50
#>   Edges: 0/6 significant (CI excludes zero)
#>   Mean inclusion probability: 0.34
#> 
#>   Centrality Stability (CS-coefficient):
#>     strength:              0.00 [Unstable]
#>     expected_influence:    0.00 [Unstable]
#> 
#>   Edge differences: 0/15 pairs significantly different
#>   Timing: 1.0s (bootstrap: 0.6s, case-drop: 0.3s)
summary(boot, type = "edges")
#>     edge weight    ci_lower  ci_upper inclusion
#> 1 A -- B      0 -0.11158199 0.2237159      0.35
#> 2 A -- C      0 -0.35185567 0.0000000      0.45
#> 3 B -- C      0 -0.30946900 0.2065199      0.25
#> 4 A -- D      0  0.00000000 0.3570903      0.50
#> 5 B -- D      0 -0.15048304 0.3105727      0.35
#> 6 C -- D      0 -0.03724162 0.1787669      0.15
# }
```
