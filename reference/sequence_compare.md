# Compare subsequence patterns between groups

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
  adjust = "fdr",
  actor = NULL
)

# S3 method for class 'net_sequence_comparison'
print(x, n = 10L, digits = 2L, ...)

# S3 method for class 'net_sequence_comparison'
summary(object, format = c("wide", "long"), ...)

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

  The sequences:

  - a `netobject_group`, from
    [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
    with `group` or from a clustering, whose groups are compared;

  - a `netobject`, `tna` model or `cograph_network` built from
    sequences;

  - a wide data frame with numbered position columns, as in
    [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
    or a character matrix, one row per sequence;

  - a `group_tna`, read as the group of networks it converts to (see
    [`as_netobject()`](https://pak.dynasite.org/Nestimate/reference/as_netobject.md)).

  Every input except a `netobject_group` needs `group`. One-hot data,
  association networks and data frames with extra columns raise
  `nestimate_bad_input`. With `actor`, a tna model must keep its
  metadata, as one built with
  [`tna::prepare_data()`](https://sonsoles.me/tna/reference/prepare_data.html)
  does. For [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html): a
  `sequence_compare()` result.

- group:

  The group of each sequence: the name of a column of a data frame or
  matrix `x` (that column is then not read as a position), or a vector
  with one label per sequence. Ignored for a `netobject_group`; required
  otherwise. At least two groups are needed. Sequences with a missing
  group are left out of the counts with a warning of class
  `nestimate_missing_group`. Default `NULL`.

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
  whole number; ignored otherwise. With `iter = 1` the effect size is
  undefined, and a warning says so. Default `1000`.

- adjust:

  The p-value adjustment, any method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)
  (`"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"`, `"BH"`, `"BY"`,
  `"fdr"` or `"none"`). Default `"fdr"` (Benjamini-Hochberg).

- actor:

  The column identifying the actor each sequence belongs to, or several
  columns for nested levels, outermost first (`c("Group", "Student")`);
  see *Nested data and actor*. With several levels, labels are shuffled
  at the level where they are assigned, within the units above it, and
  `n_actors` counts the innermost level. For a `netobject_group`,
  `netobject` or `tna` model the columns are looked up in the network's
  `$metadata` (for example the `Group` column kept by
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  when students are nested in teams); for a data frame or matrix `x`,
  they are columns of `x`, which are then not read as positions. Only
  used with `test = "permutation"`; with `"chisq"` or `"none"` it is
  ignored with a warning of class `nestimate_ignored_argument`.
  Duplicated, missing or empty names, and levels that are not nested
  outermost first (a level that adds no units below the one before it),
  raise an error of class `nestimate_bad_actor`; a column that is absent
  or has missing values, `nestimate_actor_missing`; a column that does
  not give one id per sequence, `nestimate_actor_misaligned`. Default
  `NULL`: every sequence is shuffled on its own.

- n:

  For [`print()`](https://rdrr.io/r/base/print.html): the number of
  patterns shown, the first ones of the pattern table. Default `10`.

- digits:

  For [`print()`](https://rdrr.io/r/base/print.html): the decimals of
  the counts' shares, residuals and effect sizes. Default `2`; p-values
  show 3 significant digits.

- ...:

  Not used.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `net_sequence_comparison`.

- format:

  For [`summary()`](https://rdrr.io/r/base/summary.html): `"wide"`
  (default), the pattern table, one row per pattern with a `freq_`,
  `prop_` and `resid_` column per group; `"long"`, one row per pattern
  and group, with columns `pattern`, `length`, `group`, `count`,
  `proportion`, `residual` and the pattern's `effect_size` (or
  `statistic`) and `p_value`.

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
patterns and groups, the settings and the first `n` patterns, one row
per pattern and group, rounded, and returns its input invisibly;
[`summary()`](https://rdrr.io/r/base/summary.html) returns the pattern
table, a data frame with one row per pattern (`format = "long"`: one row
per pattern and group);
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
  sequences), `min_freq`, `test`, `iter` and `adjust`; with `actor` and
  `test = "permutation"` also `actor` (the column name), `n_actors` (the
  number of distinct actors among the sequences with a group) and
  `design` (`"between"` when every actor lies in one group, `"within"`
  when every actor lies in several, `"mixed"` otherwise). With several
  `actor` levels, `n_actors` counts the innermost level, `design` reads
  `"<level>: <design> (<n> actors)"` for every level, and `actor_levels`
  is a data frame with one row per level: `level`, `n_actors`, `n_pure`
  and `n_crossed` (the units the shuffle met at that level, inside
  crossed units above, in one group or in several) and `design`
  (`"not opened"` when no unit above was crossed).

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
  patterns-by-groups table of every pattern observed at that length
  (expected under independence), whether or not `min_freq` keeps it, so
  the filter does not change the null distribution of a reported
  pattern. It asks whether the pattern's distribution is associated with
  group membership at the level of the sequence, which respects the
  sequence as the unit of analysis but has little power when there are
  few sequences. The p-value is `(b + 1) / (iter + 1)`, where `b` counts
  the permutations with a statistic at least as large as the observed
  one, so it cannot fall below `1 / (iter + 1)`. The effect size is the
  observed statistic minus the mean of the permuted ones, divided by
  their standard deviation (taken as 1 when the permuted statistics do
  not vary; with `iter = 1` the effect size is `NA`, since at least two
  permutations are needed). The shuffles use R's random number
  generator: call
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

## Nested data and `actor`

The permutation test treats every sequence as exchangeable. When several
sequences come from the same actor (sessions of one student, students of
one team), they resemble each other, and shuffling single sequences
understates how much the pattern counts vary by chance: when the groups
are made of different actors, the test then rejects too often. `actor`
names the column identifying the actor of each sequence, and the shuffle
respects it, as in
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
(Good, 2005; Anderson & ter Braak, 2003): an actor whose sequences all
lie in one group keeps them together and moves as a whole, its group
label shuffled among the other such actors; an actor with sequences in
several groups (a within-actor design, such as pre and post sessions of
each student) has its labels shuffled among its own sequences only, so
it keeps its number of sequences in each group. Mixed designs combine
the two. Use it whenever the groups are made of different actors; for
groups within actors it costs little. The observed counts, proportions
and residuals do not change; only the p-values and effect sizes do.
Sequences with a missing group keep it. With few actors there are few
distinct ways to shuffle them, and a warning of class
`nestimate_few_actors` is raised when no p-value could fall below 0.05.

Several columns, outermost first (`actor = c("Group", "Student")`), are
nested levels, the unit at each level being the combination of its
column with the ones before it. The labels are then shuffled at the
level where they are assigned, within the units above it, as in
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md):
units whose sequences all lie in one group exchange their labels with
the other such units of the same parent (all top-level units at the
first level), each moving whole; a unit with sequences in several groups
is opened at the next level, and inside such a unit of the last level
single sequences exchange labels. One column runs exactly the one-level
scheme above.

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

Anderson, M. J., & ter Braak, C. J. F. (2003). Permutation tests for
multi-factorial analysis of variance. *Journal of Statistical
Computation and Simulation*, 73(2), 85-113. (the `actor` shuffle)

Good, P. (2005). *Permutation, Parametric, and Bootstrap Tests of
Hypotheses* (3rd ed.). Springer. (the `actor` shuffle)

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
# \donttest{
# Students are nested in teams, and Achiever is a team-level label:
# shuffle whole teams, not single students
net <- build_network(group_regulation_long, method = "relative",
                     actor = "Actor", action = "Action", time = "Time",
                     group = "Achiever")
set.seed(1)
sequence_compare(net, sub = 2, iter = 20, actor = "Group")
#> Sequence Comparison  [67 patterns | 2 groups: High, Low]
#>   Lengths: 2  |  min_freq: 5  |  permutation: 20 iter (fdr), actor: Group (200 actors, between)
#> 
#>  pattern               group count prop resid  effect p     
#>  cohesion->consensus   High   503  0.04   5.78 14.65  0.0778
#>  cohesion->consensus   Low    341  0.03  -5.78 14.65  0.0778
#>  discuss->consensus    High   851  0.07  12.60 13.93  0.0778
#>  discuss->consensus    Low    418  0.03 -12.60 13.93  0.0778
#>  plan->plan            High   948  0.07  -8.73 13.21  0.0778
#>  plan->plan            Low   1356  0.11   8.73 13.21  0.0778
#>  consensus->discuss    High   789  0.06  11.65 12.24  0.0778
#>  consensus->discuss    Low    401  0.03 -11.65 12.24  0.0778
#>  discuss->adapt        High    48  0.00 -11.08 12.05  0.0778
#>  discuss->adapt        Low    234  0.02  11.08 12.05  0.0778
#>  discuss->synthesis    High   213  0.02  -5.53  9.72  0.0778
#>  discuss->synthesis    Low    344  0.03   5.53  9.72  0.0778
#>  cohesion->cohesion    High    41  0.00   5.34  9.03  0.0778
#>  cohesion->cohesion    Low      5  0.00  -5.34  9.03  0.0778
#>  emotion->consensus    High   521  0.04   4.60  8.32  0.0778
#>  emotion->consensus    Low    388  0.03  -4.60  8.32  0.0778
#>  coregulate->consensus High    97  0.01  -4.33  8.21  0.0778
#>  coregulate->consensus Low    168  0.01   4.33  8.21  0.0778
#>  adapt->consensus      High    73  0.01  -6.20  8.13  0.0778
#>  adapt->consensus      Low    170  0.01   6.20  8.13  0.0778
#>   ... and 57 more patterns (114 rows); summary() returns them all
# }
```
