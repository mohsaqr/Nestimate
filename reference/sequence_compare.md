# Compare Subsequence Patterns Between Groups

`sequence_compare()` compares how often short subsequences occur in two
or more groups of sequences. For each pattern length in `sub` it counts
every k-gram (a run of k consecutive states, such as `"A->B->C"`) in
each group, keeps the patterns that are frequent enough in every group,
and returns a `net_sequence_comparison` whose pattern table has one row
per pattern with its count, proportion and standardised residual in each
group and, optionally, a permutation or chi-square test of whether its
frequency differs between the groups. It takes a grouped network, a
single network or wide sequences with a grouping.

A residual says how far a group's count of a pattern lies from the count
expected if pattern and group were independent, in approximately
standard-normal units (Haberman, 1973); a positive value means the
pattern is more frequent in that group than independence predicts. The
p-values are adjusted across every pattern tested. The permutation test
treats the sequence as the unit, while the chi-square test treats every
occurrence as independent, which overstates the evidence when
occurrences within a sequence depend on each other (see *Details*).
Patterns of different lengths overlap (each 3-gram contains two
2-grams), so their tests are not independent of each other.
[`sequence_patterns()`](https://pak.dynasite.org/Nestimate/reference/sequence_patterns.md)
with `group` asks the related question of whether a pattern is present
in more sequences of one group.

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

  The sequences: a `netobject_group` (from
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  with `group`, or from a clustering), whose groups are compared; a
  `netobject`, `tna` model or `cograph_network` built from sequences; a
  wide data frame with one row per sequence and one column per position,
  whose position columns are named as a numbered run, as in
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  (any column beside them other than the `group` column raises an error
  of class `nestimate_bad_input`); or a character matrix with one row
  per sequence and one column per position, whatever its column names. A
  numeric or logical matrix of only 0 and 1, and one-hot or association
  networks, raise `nestimate_bad_input`. Every input except a
  `netobject_group` needs `group`. For the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: a
  `net_sequence_comparison`.

- group:

  The group of each sequence: the name of a column of a data frame or
  matrix `x` (that column is then not read as a position), or a vector
  with one label per sequence. Ignored for a `netobject_group`; required
  otherwise. At least two groups are needed. Default `NULL`.

- sub:

  The pattern lengths to count, positive whole numbers. Lengths longer
  than the sequences are dropped, and an error is raised when none
  remains. Default `3:5`.

- min_freq:

  The smallest count a pattern must reach in every group to be kept, a
  positive whole number; when no pattern of any length reaches it, an
  error is raised. Default `5`.

- test:

  The inference: `"permutation"` (default), `"chisq"` or `"none"`; see
  *Details*.

- iter:

  The number of permutations for `test = "permutation"`, a positive
  whole number; ignored otherwise. Default `1000`.

- adjust:

  The p-value adjustment, any method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)
  (`"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"`, `"BH"`, `"BY"`,
  `"fdr"` or `"none"`). Default `"fdr"` (Benjamini-Hochberg).

- ...:

  Not used.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `net_sequence_comparison`.

- top_n:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): the
  number of patterns drawn, the first ones after ranking by `sort`.
  Default `10`.

- style:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): `"auto"`
  (default) draws the back-to-back pyramid for exactly 2 groups and the
  heatmap for any other number; `"pyramid"` and `"heatmap"` force one.
  `"pyramid"` with other than 2 groups raises an error.

- sort:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): how to
  pick and order the patterns. `"statistic"` (default) ranks them by the
  absolute effect size after a permutation test, by the chi-square
  statistic after a chi-square test, and by the largest absolute
  residual without a test. `"frequency"` ranks them by their total count
  across the groups.

- alpha:

  For the pyramid: the significance level for the p-value labels.
  Patterns with `p_value < alpha` are starred and drawn in bold dark
  text, the rest in plain grey. Default `0.05`. The heatmap and a
  comparison without a test ignore it.

- show_residuals:

  For the pyramid: whether to print each bar's standardised residual
  inside it. Default `FALSE`. The heatmap always prints the residuals.

## Value

A `net_sequence_comparison` object.
[`print()`](https://rdrr.io/r/base/print.html) shows the numbers of
patterns and groups, the settings and the first 10 rows of the pattern
table, and returns its input invisibly;
[`summary()`](https://rdrr.io/r/base/summary.html) returns the pattern
table, a data frame with one row per pattern;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
residuals of the top patterns (see *Plot styles*). The object holds:

- patterns:

  The pattern table, one row per k-gram kept. Always: `pattern`,
  `length`, and per group `freq_<group>` (the count), `prop_<group>`
  (the count as a share of all k-grams of that length in the group, rare
  ones included) and `resid_<group>` (the standardised residual). With
  `test = "permutation"` also `effect_size` and `p_value`; with
  `test = "chisq"`, `statistic` and `p_value`. The p-values are
  adjusted. Rows are ordered by ascending `p_value` (ties by the larger
  absolute effect size, or the larger statistic) when a test was run,
  and by descending largest absolute residual otherwise.

- groups:

  The group names, sorted.

- n_patterns:

  The number of rows of the pattern table: the patterns that reach
  `min_freq` in every group.

- params:

  The settings used: `sub` (after dropping lengths longer than the
  sequences), `min_freq`, `test`, `iter` and `adjust`.

## Details

The standardised residuals come from the 2 x G table of the pattern's
counts against the counts of all other k-grams of the same length in
each group, with the formula
`(o - e) / sqrt(e * (1 - r/N) * (1 - c/N))`. They are always computed,
whatever `test` is. A k-gram that spans a missing position is not
counted.

`test` chooses the inference:

- `"permutation"`:

  Shuffles the group labels across sequences `iter` times and
  recomputes, for each pattern, the Euclidean norm of its row of
  differences between observed and expected counts in the
  patterns-by-groups table of the patterns kept at that length (expected
  under independence). It asks whether the pattern's distribution is
  associated with group membership at the level of the sequence, which
  respects the sequence as the unit of analysis but has little power
  when there are few sequences. The p-value is `(b + 1) / (iter + 1)`,
  where `b` counts the permutations with a statistic at least as large
  as the observed one, so it cannot fall below `1 / (iter + 1)`. The
  effect size is the observed statistic minus the mean of the permuted
  ones, divided by their standard deviation (taken as 1 when the
  permuted statistics do not vary; with `iter = 1` the effect size is
  `NA`, since at least two permutations are needed). The shuffles use
  R's random number generator: call
  [`base::set.seed()`](https://rdrr.io/r/base/Random.html) first for
  reproducible p-values.

- `"chisq"`:

  Runs [`stats::chisq.test()`](https://rdrr.io/r/stats/chisq.test.html)
  on each pattern's 2 x G table (with Yates' continuity correction when
  there are two groups). It asks whether the groups' streams of k-grams
  produce the pattern at different rates, treating every occurrence as
  an independent event: it is fast and has power with few sequences, but
  the independence it assumes is optimistic when sequences are strongly
  autocorrelated.

- `"none"`:

  No test; the counts, proportions and residuals only.

The p-values are adjusted once across all patterns of all lengths, with
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

## Plot styles

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
standardised residuals of the top `top_n` patterns, prints the plot and
returns the `ggplot` invisibly (`NULL`, with a message, when there are
no patterns). `style = "auto"` (the default) picks the style from the
number of groups:

- `"pyramid"`:

  Back-to-back bars of each pattern's proportion in the two groups, each
  bar shaded by that group's standardised residual (the colour saturates
  at plus or minus 3), with the (adjusted) p-value in the gap between
  the bars when a test was run. Needs exactly 2 groups.

- `"heatmap"`:

  One tile per pattern and group, coloured by the standardised residual
  (the colour saturates at plus or minus 3) and labelled with its value.
  Works for any number of groups.

The residuals are filled whatever `test` was chosen, so both styles work
after `test = "none"`.

## References

Haberman, S. J. (1973). The analysis of residuals in cross-classified
tables. *Biometrics*, 29(1), 205–220. (standardised residuals)

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
