# Compare two networks descriptively

`compare_model()` compares two networks, or two weight matrices, edge by
edge and returns a battery of descriptive agreement measures: weight
deviations (mean, median, root-mean-square and maximum absolute
difference, relative mean absolute difference, coefficient-of-variation
ratio), four correlations (Pearson, Spearman, Kendall, distance
correlation), five dissimilarities (Euclidean, Manhattan, Canberra,
Bray-Curtis, Frobenius), five similarities (cosine, Jaccard, Dice,
overlap, RV coefficient) and two pattern agreements (rank and sign
agreement), with the two networks' metrics side by side and, on request,
their centrality differences. It takes two `netobject`s,
`cograph_network`s or square numeric matrices, or two members of a
`netobject_group`. Both weight matrices are scaled with the same
`scaling` before comparison, which decides how weights from different
estimators are put on a common footing.

The measures describe how alike the two estimated networks are; they say
nothing about whether a difference is larger than sampling variation
would produce. For that, test the difference with
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
(two transition networks) or with
[`nct()`](https://pak.dynasite.org/Nestimate/reference/nct.md) on the
two groups' observations (psychological networks).
[`compare_networks()`](https://pak.dynasite.org/Nestimate/reference/compare_networks.md)
compares any number of networks with tidy tables and a
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) that draws one
view per call.

## Usage

``` r
compare_model(x, ...)

# Default S3 method
compare_model(x, ...)

# S3 method for class 'netobject'
compare_model(
  x,
  y,
  scaling = "none",
  measures = character(0),
  network = TRUE,
  ...
)

# S3 method for class 'cograph_network'
compare_model(
  x,
  y,
  scaling = "none",
  measures = character(0),
  network = TRUE,
  ...
)

# S3 method for class 'matrix'
compare_model(
  x,
  y,
  scaling = "none",
  measures = character(0),
  network = TRUE,
  ...
)

# S3 method for class 'netobject_group'
compare_model(
  x,
  i = 1L,
  j = 2L,
  scaling = "none",
  measures = character(0),
  network = TRUE,
  ...
)

# S3 method for class 'net_comparison'
print(x, ...)

# S3 method for class 'net_comparison'
plot(
  x,
  type = c("scatter", "heatmap", "diff_hist", "weight_dist", "all"),
  combined = TRUE,
  ...
)
```

## Arguments

- x:

  The first network: a `netobject`, a `cograph_network` or a numeric
  square matrix, or a `netobject_group` whose members `i` and `j` are
  compared. For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: a
  `net_comparison` object.

- ...:

  Not used.

- y:

  The second network, in the same forms (not used for a
  `netobject_group`). The two must hold the same nodes and have no
  missing weights. Named inputs are aligned by node name, so the node
  order may differ; different node sets, or missing, empty or duplicated
  names, raise an error. A matrix without names is matched to the other
  input by position and takes its node names; when neither input has
  names the nodes are labelled `V1`, `V2`, ....

- scaling:

  The scaling applied to each weight matrix before comparison. Default
  `"none"`. One of:

  `"none"`

  :   The weights as they are.

  `"minmax"`

  :   \\(w - \min) / (\max - \min)\\, which maps them to \\\[0, 1\]\\; a
      constant matrix is left as it is.

  `"max"`

  :   \\w / \max\|w\|\\, dividing by the largest absolute weight, so
      signs are kept; a matrix of zeros is left as it is.

  `"rank"`

  :   Average ranks, mapped to \\\[0, 1\]\\ by min-max; an ordinal
      scaling. A constant matrix keeps its (equal) ranks.

  `"zscore"`

  :   \\(w - \bar w) / s_w\\, the standard score; a constant matrix is
      only centred.

  `"robust"`

  :   \\(w - \mathrm{med}(w)) / \mathrm{mad}(w)\\, a standard score from
      the median and the median absolute deviation, less affected by
      outliers; with a MAD of 0 the weights are only centred.

  `"log"`

  :   \\\log(w)\\; needs \\w \> 0\\: a weight at or below 0 raises an
      error of class `nestimate_bad_input` (use `"log1p"` for networks
      with zeros).

  `"log1p"`

  :   \\\log(1 + w)\\; defined for \\w \> -1\\: a weight of -1 or below
      raises an error of class `nestimate_bad_input`.

  `"softmax"`

  :   A numerically stable softmax over all the cells.

  `"quantile"`

  :   The empirical cumulative distribution of the cells.

  `"frobenius"`

  :   The matrix divided by its Frobenius norm \\\\W\\\_F = \sqrt{\sum
      w\_{ij}^2}\\, a matrix-level normalisation.

  `"row"`

  :   Each row divided by the sum of its absolute values, so the rows
      sum to 1 for non-negative weights; rows summing to zero are left
      unchanged.

  Scalings that produce negative weights (`"zscore"`, `"robust"`) work
  with `network = TRUE`: the side-by-side distances are computed on the
  absolute weights.

- measures:

  Centrality measures to compare. Default `character(0)` (none). Any
  built-in measure of
  [`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
  can be named: `"OutStrength"`, `"InStrength"`, `"ClosenessIn"`,
  `"ClosenessOut"`, `"Closeness"`, `"Betweenness"`, `"BetweennessRSP"`,
  `"Diffusion"` and `"Clustering"` (`"InCloseness"` and `"OutCloseness"`
  are accepted for `"ClosenessIn"` and `"ClosenessOut"`). Unknown names
  are dropped with a warning. The centralities are computed on the
  scaled matrices as directed networks, without self-loops.

- network:

  Whether to add the two networks' metrics side by side (the metrics
  [`summary()`](https://rdrr.io/r/base/summary.html) reports for a
  `netobject`). Default `TRUE`.

- i, j:

  For a `netobject_group`: the two members to compare, by index or name.
  Defaults `1` and `2`; the two must differ, and a member that does not
  exist raises an error.

- type:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): the
  view. `"scatter"` (default) plots each cell's weight in `y` against
  its weight in `x`, with the line of perfect agreement (dashed), an OLS
  fit with its 95% band and the Pearson, Spearman and Kendall
  correlations. `"heatmap"` draws the difference matrix `x - y` on the
  diverging palette: blue where `x` is larger, red where `y` is larger,
  white where they agree. `"diff_hist"` is a histogram of the absolute
  differences `|x - y|`, with a density curve, a rug of the cells and
  lines at the mean and median. `"weight_dist"` overlays the
  distributions of `|x|` and `|y|`. `"all"` draws these four views in a
  2 x 2 grid.

- combined:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html) with
  `type = "all"`: `TRUE` (default) arranges the four views into one 2 x
  2 gtable, which needs the gridExtra package; `FALSE` returns them as a
  named list of ggplots, to print, save or lay out separately. Ignored
  for the other views.

## Value

A `net_comparison` object. Print it for the agreement measures, the
side-by-side network metrics and the centrality correlations;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
edge-weight scatter or another view (see `type`). The object is a list
with:

- matrices:

  The two scaled weight matrices, `x` and `y`, aligned by node.

- difference_matrix:

  The scaled `x` minus the scaled `y`.

- edge_metrics:

  A data frame with one row per cell (`source`, `target`, self-loops
  included): the two weights and their raw, absolute, squared, relative
  and rank differences, percentile difference, ratio
  (`similarity_strength_index`), relative change (`difference_index`),
  log ratio of 1 + weight (`logarithmic_ratio`, `NA` where a weight is
  -1 or below), standardised weights and their ratio. A ratio is `NA`
  where its denominator is zero, and the standardised weights are `NA`
  when a matrix is constant.

- summary_metrics:

  A data frame of the 22 agreement measures, with columns `category`,
  `metric` and `value`.

- network_metrics:

  When `network = TRUE`: a data frame with columns `metric`, `x` and
  `y`, one row per network metric.

- centrality_differences, centrality_correlations:

  When `measures` is not empty (empty data frames when none of its names
  is valid): one row per node and measure (`state`, `centrality`, `x`,
  `y`, `difference`), and one row per measure with the correlation of
  the two networks' values.

[`print()`](https://rdrr.io/r/base/print.html) returns its input
invisibly. [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
returns a ggplot; with `type = "all"`, it draws the 2 x 2 grid and
returns it invisibly as a gtable, or with `combined = FALSE` returns the
four ggplots invisibly as a named list.

## Details

The measures, scalings and edge metrics mirror
[`tna::compare()`](https://sonsoles.me/tna/reference/compare.html) (tna
1.5.4 and later) numerically. The distance correlation is that of
Szekely, Rizzo and Bakirov (2007), `NA` when one network is constant,
and the RV coefficient that of Robert and Escoufier (1976). A network's
own `directed` flag decides how its network metrics are computed; a
plain matrix has none, so an asymmetric matrix is read as directed and a
symmetric one as undirected.

## References

Szekely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and
testing dependence by correlation of distances. *The Annals of
Statistics*, 35(6), 2769–2794.
[doi:10.1214/009053607000000505](https://doi.org/10.1214/009053607000000505)

Robert, P., & Escoufier, Y. (1976). A unifying tool for linear
multivariate statistical methods: The RV-coefficient. *Journal of the
Royal Statistical Society, Series C (Applied Statistics)*, 25(3),
257–265. [doi:10.2307/2347233](https://doi.org/10.2307/2347233)

## Examples

``` r
nets <- build_network(group_regulation_long, method = "relative",
                      actor = "Actor", action = "Action", time = "Time",
                      group = "Achiever")
compare_model(nets)
#> Network comparison
#> ==================
#> Summary metrics:
#>              category               metric   value
#>     Weight Deviations      Mean Abs. Diff. 0.03225
#>     Weight Deviations    Median Abs. Diff. 0.01813
#>     Weight Deviations            RMS Diff. 0.05217
#>     Weight Deviations       Max Abs. Diff.  0.2103
#>     Weight Deviations Rel. Mean Abs. Diff.  0.2902
#>     Weight Deviations             CV Ratio   1.103
#>          Correlations              Pearson  0.9211
#>          Correlations             Spearman  0.9153
#>          Correlations              Kendall  0.7672
#>          Correlations             Distance  0.9156
#>       Dissimilarities            Euclidean  0.4695
#>       Dissimilarities            Manhattan   2.612
#>       Dissimilarities             Canberra   19.76
#>       Dissimilarities          Bray-Curtis  0.1451
#>       Dissimilarities            Frobenius  0.2213
#>          Similarities               Cosine   0.954
#>          Similarities              Jaccard  0.7466
#>          Similarities                 Dice  0.8549
#>          Similarities              Overlap  0.8549
#>          Similarities                   RV  0.8983
#>  Pattern Similarities       Rank Agreement  0.8194
#>  Pattern Similarities       Sign Agreement  0.9383
#> 
#> Network metrics (x vs y):
#>                       metric       x       y
#>                   Node Count       9       9
#>                   Edge Count      76      75
#>              Network Density       1       1
#>                Mean Distance 0.04229 0.05596
#>            Mean Out-Strength       1       1
#>              SD Out-Strength       0       0
#>             Mean In-Strength       1       1
#>               SD In-Strength  0.9141  0.7186
#>              Mean Out-Degree   8.444   8.333
#>                SD Out-Degree    1.13   0.866
#>  Centralization (Out-Degree) 0.04688  0.0625
#>   Centralization (In-Degree) 0.04688  0.0625
#>                  Reciprocity  0.9565  0.9412
```
