# Compare Subsequence Patterns Between Groups

Extracts all k-gram patterns (subsequences of length k) from sequences
in each group, computes standardized residuals against the independence
model, and optionally runs a permutation or chi-square test of group
differences.

## Usage

``` r
sequence_compare(
  x,
  group = NULL,
  sub = 3:5,
  min_freq = 5L,
  test = c("permutation", "chisq", "none"),
  iter = 1000L,
  adjust = "fdr"
)

# S3 method for class 'net_sequence_comparison'
print(x, ...)

# S3 method for class 'net_sequence_comparison'
summary(object, ...)

# S3 method for class 'net_sequence_comparison'
plot(
  x,
  top_n = 10L,
  style = c("auto", "pyramid", "heatmap"),
  sort = c("statistic", "frequency"),
  alpha = 0.05,
  show_residuals = FALSE,
  ...
)
```

## Arguments

- x:

  A `netobject_group` (from grouped `build_network`), a `netobject`
  (requires `group`), or a wide-format `data.frame` (requires `group`).
  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: an
  object of class `net_sequence_comparison`.

- group:

  Character or vector. Column name or vector of group labels. Not needed
  for `netobject_group`.

- sub:

  Integer vector. Pattern lengths to analyze. Default: `3:5`.

- min_freq:

  Integer. Minimum frequency in each group for a pattern to be included:
  a pattern is kept only when its count reaches this threshold in
  *every* group. Default: 5.

- test:

  Character. Inference method: one of `"permutation"` (default),
  `"chisq"`, or `"none"`. See Details.

- iter:

  Integer. Permutation iterations. Only used when
  `test = "permutation"`. Default: 1000.

- adjust:

  Character. P-value correction method (see
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html)). Default:
  `"fdr"`.

- ...:

  For the S3 methods: further arguments passed to or from other methods.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: an
  object of class `net_sequence_comparison`.

- top_n:

  Integer. Show top N patterns. Default: 10.

- style:

  Character. `"auto"` (default) draws the back-to-back pyramid for
  exactly 2 groups and the heatmap for any other number; `"pyramid"` and
  `"heatmap"` force a specific style.

- sort:

  Character. `"statistic"` (default) ranks patterns by test statistic or
  residual magnitude. `"frequency"` ranks by total occurrence count
  across all groups.

- alpha:

  Numeric. Significance threshold for p-value display in the pyramid:
  patterns with `p_value < alpha` are starred and drawn in bold dark
  text, the rest stay plain grey. Default: 0.05.

- show_residuals:

  Logical. If `TRUE`, print the standardized residual value inside each
  pyramid bar. Default: `FALSE`. Ignored for the heatmap (which always
  shows residuals).

## Value

An object of class `"net_sequence_comparison"` containing:

- patterns:

  Tidy data.frame, one row per retained k-gram pattern. Always present:
  `pattern`, `length`, and one `freq_<group>`, `prop_<group>` and
  `resid_<group>` column per group. If `test = "permutation"`:
  `effect_size`, `p_value`. If `test = "chisq"`: `statistic`, `p_value`.
  Rows are ordered by ascending adjusted `p_value` when a test was run,
  and by descending maximum absolute residual otherwise.

- groups:

  Character vector of group names, sorted.

- n_patterns:

  Integer. Number of rows in `patterns`, i.e. the patterns meeting
  `min_freq` in every group.

- params:

  List of sub, min_freq, test, iter, adjust.

In `print.net_sequence_comparison()`: The input object, invisibly.

In `summary.net_sequence_comparison()`: The `patterns` data.frame: tidy,
one row per k-gram pattern, with a frequency, proportion and
standardized-residual column per group, and the test columns when `test`
was not `"none"`.

In `plot.net_sequence_comparison()`: The drawn `ggplot` object,
invisibly (the plot is also printed). `NULL`, invisibly, when the object
holds no patterns.

## Details

Standardized residuals are always computed from a 2xG contingency table
of (this pattern vs. everything else) using the textbook formula
`(o - e) / sqrt(e * (1 - r/N) * (1 - c/N))`. They describe how much each
group's count for a given pattern deviates from expectation under
independence, scaled to be approximately N(0,1) under the null.

The optional `test` argument chooses an inference method:

- `"permutation"`:

  Shuffles group labels across sequences and recomputes a per-pattern
  statistic (row-wise Euclidean residual norm). Answers: "is this
  pattern's distribution associated with group membership at the *actor*
  level?" Respects the sequence as the unit of analysis; can be
  underpowered when the number of sequences is small.

- `"chisq"`:

  Runs `chisq.test` on the 2xG table per pattern. Answers: "do the group
  *streams* generate this pattern at different rates?" Treats each
  k-gram occurrence as an event; fast and powerful even with few
  sequences, but the iid assumption it makes is optimistic when
  sequences are strongly autocorrelated.

- `"none"`:

  Skip inference. Only residuals, frequencies, and proportions are
  returned.

P-values are adjusted once across all patterns (not per-pattern) using
any method supported by
[`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). The default is
`"fdr"` (Benjamini-Hochberg).

## References

Haberman, S. J. (1973). The analysis of residuals in cross-classified
tables. *Biometrics*, 29(1), 205–220. (standardized residuals)

Benjamini, Y. & Hochberg, Y. (1995). Controlling the false discovery
rate. *Journal of the Royal Statistical Society B*, 57(1), 289–300. (the
default `adjust = "fdr"`)

## Examples

``` r
set.seed(1)
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 60, TRUE),
  V2 = sample(LETTERS[1:4], 60, TRUE),
  V3 = sample(LETTERS[1:4], 60, TRUE),
  V4 = sample(LETTERS[1:4], 60, TRUE)
)
grp <- rep(c("X", "Y"), 30)
net <- build_network(seqs, method = "relative")
res <- sequence_compare(net, group = grp, sub = 2:3, test = "chisq")
```
