# Compare two networks descriptively

Computes a battery of descriptive comparison metrics between two
networks or two weight matrices: weight deviations (mean / median / RMS
/ max absolute difference, relative mean absolute difference,
coefficient-of- variation ratio), four correlation measures (Pearson,
Spearman, Kendall, distance correlation), five dissimilarity measures
(Euclidean, Manhattan, Canberra, Bray-Curtis, Frobenius), five
similarity measures (Cosine, Jaccard, Dice, Overlap, RV), pattern
agreements, and side-by-side network metrics. Optionally adds centrality
differences and centrality correlations.

## Usage

``` r
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

  A `netobject`, `cograph_network`, or numeric square matrix, or a
  `netobject_group` whose members `i` and `j` are compared. For the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: an
  object of class `net_comparison`.

- ...:

  Ignored. For the S3 methods: further arguments passed to or from other
  methods.

- y:

  A `netobject`, `cograph_network`, or numeric square matrix.

- scaling:

  Scaling applied to both weight matrices before comparison. One of:

  `"none"`

  :   Identity (default).

  `"minmax"`

  :   \\(w - \min) / (\max - \min)\\; maps to \\\[0, 1\]\\.

  `"max"`

  :   \\w / \max(\|w\|)\\; preserves sign.

  `"rank"`

  :   Min-max of average ranks; ordinal scaling.

  `"zscore"`

  :   \\(w - \bar w) / s_w\\; standard score.

  `"robust"`

  :   \\(w - \mathrm{med}(w)) / \mathrm{mad}(w)\\; Huber-style robust
      z-score, resists outliers.

  `"log"`

  :   \\\log(w)\\; requires \\w \> 0\\.

  `"log1p"`

  :   \\\log(1 + w)\\; admits \\w \ge 0\\.

  `"softmax"`

  :   Numerically stable softmax over the flattened vector.

  `"quantile"`

  :   Empirical CDF of the flattened vector.

  `"frobenius"`

  :   Divide the matrix by its Frobenius norm \\\\W\\\_F = \sqrt{\sum
      w\_{ij}^2}\\; matrix-level normalisation.

  `"row"`

  :   Row-stochastic normalisation (each row's absolute values sum to
      1). Only meaningful for non-negative matrices; rows summing to
      zero are left unchanged.

  Scalings that produce negative weights (`zscore`, `robust`) are
  compatible with `network = TRUE` because the side-by-side metrics use
  Nestimate's base-R Floyd-Warshall, which handles negative weights.

- measures:

  Character vector of centrality measures to compare. Empty by default
  (no centrality block). Any built-in measure is valid: `"OutStrength"`,
  `"InStrength"`, `"ClosenessIn"`, `"ClosenessOut"`, `"Closeness"`,
  `"Betweenness"`, `"BetweennessRSP"`, `"Diffusion"`, `"Clustering"`,
  `"InCloseness"`, `"OutCloseness"`. Unknown names are ignored with a
  warning.

- network:

  Logical. Include side-by-side network metrics from
  [`summary()`](https://rdrr.io/r/base/summary.html)? Default `TRUE`.

- i, j:

  For a `netobject_group`: index or name of the two member networks to
  compare. Defaults `1L` and `2L`.

- type:

  Character. One of `"scatter"` (default - edge-weight scatter with OLS
  fit and correlation overlay), `"heatmap"` (n by n grid of x - y
  differences using the diverging palette), `"diff_hist"` (histogram of
  \|x - y\| absolute differences with rug + density), `"weight_dist"`
  (overlaid distributions of \|x\| and \|y\| edge weights), or `"all"`
  (2 by 2 grid of all four panels; requires the gridExtra package).

- combined:

  When `type = "all"` and `combined = TRUE` (default), the four panels
  are stitched into a 2x2 gtable. When `FALSE`, returns a named list of
  the four ggplots so each can be printed, saved, or re-laid-out
  independently. Ignored for other `type` values.

## Value

A `net_comparison` object: a named list with `matrices`,
`difference_matrix`, `edge_metrics`, `summary_metrics`, optionally
`network_metrics`, `centrality_differences`, `centrality_correlations`.

In `plot.net_comparison()`: A `ggplot` object; for `type = "all"` with
`combined = TRUE` a `gtable` arranged 2 by 2; for `type = "all"` with
`combined = FALSE` a named list of four ggplots.

## Details

Mirrors
[`tna::compare()`](https://sonsoles.me/tna/reference/compare.html)
numerically. Inputs are converted to weight matrices and scaled before
comparison; the choice of scaling determines how weights from different
estimators are placed on a common footing.

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
#>          Correlations             Distance  0.8382
#>       Dissimilarities            Euclidean  0.4695
#>       Dissimilarities            Manhattan   2.612
#>       Dissimilarities             Canberra   14.76
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
#>              SD Out-Strength  0.9141  0.7186
#>             Mean In-Strength       1       1
#>               SD In-Strength       0       0
#>              Mean Out-Degree   8.444   8.333
#>                SD Out-Degree    1.13   0.866
#>  Centralization (Out-Degree) 0.04688  0.0625
#>   Centralization (In-Degree) 0.04688  0.0625
#>                  Reciprocity  0.9565  0.9412
```
