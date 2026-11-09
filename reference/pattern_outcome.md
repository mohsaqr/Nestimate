# Model a Binary Outcome From Sequential Patterns

`pattern_outcome()` selects the `n` sequential patterns most associated
with a two-valued outcome and fits one logistic regression of the
outcome on their presence (or their counts), returning one row per
selected pattern with its log odds ratio, odds ratio and Wald interval.
It reads any sequence data
[`extract_sequences()`](https://pak.dynasite.org/Nestimate/reference/extract_sequences.md)
reads, takes the outcome from a stored covariate, a vector or each
sequence's last state, and can adjust for covariates and add a random
intercept per actor.

Each estimate is the association of one pattern with the outcome,
adjusted for the other selected patterns and the covariates; it is not a
causal effect. The patterns are chosen by their association with the
outcome in the same data the model is fitted to, so the p-values do not
account for that selection and overstate the evidence: treat them as a
ranking, or confirm the selected patterns on new data. A pattern and a
longer pattern containing it often occur together, and each coefficient
is then adjusted for the other;
[`sequence_patterns()`](https://pak.dynasite.org/Nestimate/reference/sequence_patterns.md)
shows how often each pattern occurs, overall and by group.

## Usage

``` r
pattern_outcome(
  x,
  outcome = "last_obs",
  type = c("ngram", "gapped", "repeated"),
  len = 1:2,
  gap = 1,
  min_freq = 5,
  min_support = 0.01,
  start = NULL,
  end = NULL,
  contain = NULL,
  n = 10,
  priority = "chisq",
  desc = TRUE,
  occurrences = FALSE,
  adjust = NULL,
  adjust_length = FALSE,
  group = NULL,
  reference = NULL,
  ci_level = 0.95,
  cols = NULL,
  actor = NULL,
  action = NULL,
  time = NULL,
  ...
)
```

## Arguments

- x:

  The sequence data, in any form
  [`extract_sequences()`](https://pak.dynasite.org/Nestimate/reference/extract_sequences.md)
  accepts: prepared data, a long event log, a wide data frame, a
  `netobject`, `netobject_group` or `tna` model.

- outcome:

  The outcome, with exactly two distinct values: the name of a column
  stored with the sequences (a covariate of a wide data frame or long
  event log, or of a network's metadata), or a vector with one value per
  sequence; use a vector for a model that stores no covariates, such as
  a `tna` model. The default, `"last_obs"`, uses each sequence's last
  observed state as the outcome and removes the states that end some
  sequence from all sequences, as codyna does; it suits data whose final
  state codes the outcome, such as an appended `"pass"` or `"fail"`. An
  outcome without exactly two values, a vector of the wrong length or a
  column that is not stored with the sequences raises an error of class
  `nestimate_bad_input`.

- type, len, gap, min_freq, min_support, start, end, contain:

  The pattern definition and filters; see
  [`sequence_patterns()`](https://pak.dynasite.org/Nestimate/reference/sequence_patterns.md).
  The defaults are codyna's: `type = "ngram"`, `len = 1:2`, `gap = 1`,
  `min_freq = 5` and `min_support = 0.01`. When no pattern passes the
  filters under both outcome values, an error of class
  `nestimate_no_patterns` is raised.

- n:

  The number of patterns entered into the model, a positive whole
  number. Default `10`.

- priority:

  The columns ranking the patterns for selection, applied in order:
  `"chisq"` (default), `"p_value"`, `"support"`, `"n_sequences"`,
  `"occurrences"`, `"proportion"`, `"lift"` and `"length"`.

- desc:

  Whether each `priority` column ranks in descending order, recycled to
  the length of `priority`. Default `TRUE`, so `"p_value"` needs
  `desc = FALSE` to select the smallest p-values.

- occurrences:

  Whether each pattern enters the model as its number of occurrences in
  the sequence instead of its presence (0/1). Default `FALSE`.

- adjust:

  Covariate columns added to the model. Default `NULL`. A name that is
  not a column stored with the sequences raises an error of class
  `nestimate_bad_input`.

- adjust_length:

  Whether to add the number of observed states per sequence as a
  covariate. Default `FALSE`, as in codyna.

- group:

  The column identifying the actor (or other cluster) each sequence
  belongs to; it adds a random intercept per group. It needs lme4,
  otherwise an error of class `nestimate_no_lme4` is raised. Default
  `NULL`. A name that is not a column stored with the sequences raises
  an error of class `nestimate_bad_input`.

- reference:

  The reference outcome value; the odds ratios are for the other value.
  Default `NULL`, the first value in sorted order; a value that does not
  occur raises an error of class `nestimate_bad_input`.

- ci_level:

  The confidence level of the Wald intervals, a number between 0 and 1.
  Default `0.95`.

- cols:

  The state columns of a wide data frame; see
  [`extract_sequences()`](https://pak.dynasite.org/Nestimate/reference/extract_sequences.md).
  Default `NULL`.

- actor, action, time:

  The actor, action and timestamp columns of a long event log, as in
  [`extract_sequences()`](https://pak.dynasite.org/Nestimate/reference/extract_sequences.md):
  setting `action` declares `x` an event log. Default `NULL`.

- ...:

  Further arguments passed to
  [`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md)
  for a long event log, such as `time_threshold` or `order`.

## Value

A data frame with one row per selected pattern, ordered by pattern name,
with the columns `pattern`, `length`, `occurrences`, `n_sequences`,
`support` and `lift` (as in
[`sequence_patterns()`](https://pak.dynasite.org/Nestimate/reference/sequence_patterns.md));
the share of the sequences of each outcome value that contain the
pattern (`support_<value>`); `chisq` (the selection statistic);
`estimate` (the log odds ratio, adjusted for the other terms),
`std_error`, `ci_lower`, `ci_upper` and `p_value` (Wald z); and
`odds_ratio`, `or_lower` and `or_upper`.

Sequences with a missing outcome, covariate or group are dropped with a
warning of class `nestimate_missing_outcome`. A selected pattern
collinear with the other terms gets `NA` estimates (and `NA` standard
error, intervals, p-value and odds ratios), its descriptive columns
kept, and a warning of class `nestimate_aliased_pattern`. Warnings from
the model fit (such as separation) are passed on unchanged.

## Details

`pattern_outcome()` reproduces `codyna::analyze_outcome()`: the same
pattern discovery and filters, the same ranking and the same model, so
the coefficients agree. The steps:

1.  Patterns are found as in
    [`sequence_patterns()`](https://pak.dynasite.org/Nestimate/reference/sequence_patterns.md)
    and filtered by `min_freq`, `min_support`, `start`, `end` and
    `contain`. Patterns absent from all sequences of one outcome value
    are dropped.

2.  Each pattern gets the chi-square statistic of the number of
    sequences containing it under each outcome value, against an equal
    split across the values (codyna's statistic; it does not weight by
    the size of each outcome group).

3.  The top `n` patterns by `priority` enter one model,
    `outcome ~ adjust + pattern_1 + ... + pattern_n`, as presence
    indicators (or occurrence counts with `occurrences = TRUE`). Ties
    keep codyna's order: more occurrences first, then first appearance.
    With `group`, a random intercept per group is added
    ([`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html));
    otherwise the model is a plain
    [`glm()`](https://rdrr.io/r/stats/glm.html).

## References

Saqr, M., Tikka, S., & López-Pernas, S. (2026). codyna: Complex Dynamic
Systems. R package version 1.0.0.

## See also

[`sequence_patterns`](https://pak.dynasite.org/Nestimate/reference/sequence_patterns.md)

## Examples

``` r
pattern_outcome(group_regulation_long, outcome = "Achiever",
                actor = "Actor", action = "Action", time = "Time")
#>                pattern length occurrences n_sequences support       lift
#> 1                adapt      1         554         475  0.2375 1.00000000
#> 2      adapt->cohesion      2         139         137  0.0685 0.50203839
#> 3     adapt->consensus      2         243         222  0.1110 0.51218457
#> 4   cohesion->cohesion      2          46          42  0.0210 0.06362668
#> 5   consensus->discuss      2        1190         828  0.4140 0.57539458
#> 6  coregulate->discuss      2         539         452  0.2260 0.47258065
#> 7       discuss->adapt      2         282         258  0.1290 0.68884958
#> 8   discuss->consensus      2        1269         880  0.4400 0.61153048
#> 9     emotion->discuss      2         289         267  0.1335 0.22848693
#> 10    synthesis->adapt      2         153         146  0.0730 1.01946408
#>    support_High support_Low     chisq    estimate std_error   ci_lower
#> 1         0.145       0.330  72.05263 -0.35794847 0.2848304 -0.9162057
#> 2         0.037       0.100  28.97080  0.14226023 0.2879495 -0.4221105
#> 3         0.068       0.154  33.31532 -0.06048539 0.2541657 -0.5586410
#> 4         0.037       0.005  24.38095 -1.95703140 0.5129015 -2.9622999
#> 5         0.513       0.315  47.34783 -0.84235929 0.1129553 -1.0637476
#> 6         0.179       0.273  19.54867  0.56362598 0.1273555  0.3140139
#> 7         0.047       0.211 104.24806  2.24802414 0.2818028  1.6957007
#> 8         0.545       0.335  50.11364 -0.66382982 0.1135967 -0.8864753
#> 9         0.171       0.096  21.06742 -0.62744477 0.1563031 -0.9337931
#> 10        0.039       0.107  31.67123  1.57967906 0.2981577  0.9953006
#>      ci_upper      p_value odds_ratio   or_lower   or_upper
#> 1   0.2003088 2.088596e-01  0.6991091 0.40003400  1.2217800
#> 2   0.7066309 6.212739e-01  1.1528766 0.65566161  2.0271501
#> 3   0.4376702 8.118996e-01  0.9413075 0.57198585  1.5490940
#> 4  -0.9517629 1.358477e-04  0.1412772 0.05169988  0.3860598
#> 5  -0.6209710 8.820688e-14  0.4306932 0.34515988  0.5374223
#> 6   0.8132381 9.616864e-06  1.7570319 1.36890874  2.2551987
#> 7   2.8003476 1.495760e-15  9.4690079 5.45046393 16.4503631
#> 8  -0.4411844 5.104103e-09  0.5148757 0.41210576  0.6432741
#> 9  -0.3210964 5.962671e-05  0.5339544 0.39305994  0.7253533
#> 10  2.1640575 1.169935e-07  4.8533979 2.70553754  8.7063924
```
