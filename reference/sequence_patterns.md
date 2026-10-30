# Discover Sequential Patterns

Finds the ordered patterns that occur in a set of sequences and reports,
for each one, how often it occurs, how many sequences contain it, and
its lift. Three pattern types are supported:

- `"ngram"`:

  Contiguous runs of `len` states, e.g. `A->B->C`.

- `"gapped"`:

  Two states separated by `gap` positions of any state, e.g. `A->*->B`
  (gap 1) or `A->**->B` (gap 2).

- `"repeated"`:

  Runs of `len` identical states, e.g. `A->A->A`.

A window that contains a missing value is skipped, so a pattern never
spans a gap or the end of a sequence.

## Usage

``` r
sequence_patterns(
  x,
  type = c("ngram", "gapped", "repeated"),
  len = 2:5,
  gap = 1:3,
  min_freq = 2,
  min_support = 0.01,
  start = NULL,
  end = NULL,
  contain = NULL,
  pattern = NULL,
  group = NULL,
  reference = NULL,
  adjust = "BH",
  top = NULL,
  cols = NULL,
  ...
)
```

## Arguments

- x:

  Sequence data, in any form accepted by
  [`extract_sequences()`](https://pak.dynasite.org/Nestimate/reference/extract_sequences.md):
  prepared data, a long event log (with `action`), a wide data frame or
  character matrix, a `netobject`, `netobject_group` or `tna` model.

- type:

  Pattern type: `"ngram"` (default), `"gapped"` or `"repeated"`.

- len:

  Integer vector of pattern lengths for `"ngram"` and `"repeated"`.
  Default `2:5`.

- gap:

  Integer vector of gap sizes for `"gapped"`. Default `1:3`.

- min_freq:

  Minimum number of occurrences. Default `2`.

- min_support:

  Minimum share of sequences containing the pattern. Default `0.01`.

- start, end:

  Optional character vectors of states. Keep only patterns that begin
  (`start`) or end (`end`) with one of them.

- contain:

  Optional character vector of states. Keep only patterns that include
  at least one of them.

- pattern:

  Optional template to search for instead of discovering all patterns,
  e.g. `"discuss->*->consensus"`: each `*` is one position holding any
  state (`**` is two). Every concrete instance is returned as its own
  row (`"discuss->plan->consensus"`, ...). `type`, `len` and `gap` are
  then ignored.

- group:

  Optional grouping of the sequences: the name of a column stored with
  them (a covariate of a wide data frame or long event log, or of a
  network's metadata), or a vector with one value per sequence. For a
  `netobject_group` its own grouping is used when `group` is `NULL`.
  Adds a test of pattern presence by group; see Details.
  `group = "last_obs"` uses each sequence's last observed state as its
  group and removes the states that end some sequence from all sequences
  (for data whose final state codes the outcome, e.g. an appended
  `"pass"`/`"fail"`; codyna's `last_obs`).

- reference:

  For two groups, the reference group; the odds ratio is for the other
  group. Default: the first group in sorted order.

- adjust:

  Multiplicity correction for the group test, passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). Default `"BH"`.

- top:

  Optional number of rows to return (the first `top` in the output
  order). Default `NULL`, all rows.

- cols:

  State columns of a wide data frame; see
  [`extract_sequences()`](https://pak.dynasite.org/Nestimate/reference/extract_sequences.md).

- ...:

  Passed to
  [`extract_sequences()`](https://pak.dynasite.org/Nestimate/reference/extract_sequences.md),
  e.g. `actor`, `action` and `time` for a long event log.

## Value

A `data.frame` with one row per pattern, ordered by `support` (then
`occurrences`, then `pattern`):

- pattern:

  States joined by `"->"`; `*` marks a wildcard position.

- length:

  Number of positions the pattern spans.

- occurrences:

  Total number of times it occurs.

- proportion:

  `occurrences` as a share of all patterns of the same length.

- n_sequences:

  Number of sequences that contain it.

- support:

  `n_sequences` divided by the number of sequences.

- lift:

  Support relative to independent presence of its states.

With `group`, also `support_<group>` (the support within each group);
for two groups `odds_ratio`, `or_lower`, `or_upper`; for more, `chisq`
and `cramers_v`; then `p_value` and `p_adj`. Rows are ordered by
`p_value` (then `support`).

Conditions: `nestimate_bad_input` for a state name containing `"->"` or
consisting only of asterisks (reserved for wildcards), a missing group
column, a group vector of the wrong length, or a group with fewer than
two values; warning `nestimate_missing_group` when sequences with a
missing group are dropped.

## Details

`support` is the share of sequences that contain the pattern at least
once. `lift` compares it with the support expected if the states in the
pattern appeared in sequences independently:
`support / prod(state_support)`, where the support of a state is the
share of sequences in which it occurs and wildcard positions are left
out. A lift above 1 means the states co-occur in this order more often
than independent presence predicts. It does not account for the
transition structure: a pattern that the transition probabilities make
likely can have a lift above 1.

`proportion` is computed over all patterns of the same length before
`min_freq`, `min_support`, `start`, `end` and `contain` are applied, so
it does not change with the filters.

## Group test

With `group`, each pattern's presence (in a sequence or not) is
cross-tabulated against the group. For two groups, Fisher's exact test
gives the p-value, the conditional maximum-likelihood odds ratio of
containing the pattern in the non-reference group, and its 95\\
confidence interval. A pattern that occurs in one group only has an
infinite (or zero) odds ratio, with a finite lower (or upper) limit; a
pattern present in every sequence has none (`NA`). For three or more
groups, the Pearson chi-square of the presence-by-group table gives the
p-value (G - 1 degrees of freedom), with Cramér's V as the effect size.
Unlike a comparison of counts against an equal split, both tests account
for unequal group sizes. p-values are adjusted across the returned
patterns' family (all patterns passing the filters) with `adjust`. The
rows are then ordered by `p_value`.

## References

Agrawal, R., & Srikant, R. (1995). Mining sequential patterns.
*Proceedings of the Eleventh International Conference on Data
Engineering*, 3–14.

Fisher, R. A. (1935). The logic of inductive inference. *Journal of the
Royal Statistical Society*, 98(1), 39–82.

## See also

[`pattern_outcome`](https://pak.dynasite.org/Nestimate/reference/pattern_outcome.md),
[`sequence_compare`](https://pak.dynasite.org/Nestimate/reference/sequence_compare.md)

## Examples

``` r
sequence_patterns(trajectories, len = 2:3, top = 10)
#>                      pattern length occurrences proportion n_sequences
#> 1           Average->Average      2         458 0.26412918          99
#> 2            Active->Average      2         176 0.10149942          82
#> 3             Active->Active      2         459 0.26470588          79
#> 4            Average->Active      2         153 0.08823529          77
#> 5  Average->Average->Average      3         284 0.17716781          68
#> 6        Disengaged->Average      2         129 0.07439446          67
#> 7        Average->Disengaged      2         140 0.08073818          66
#> 8   Active->Average->Average      3          83 0.05177792          64
#> 9    Active->Active->Average      3          84 0.05240175          59
#> 10    Active->Active->Active      3         340 0.21210231          57
#>      support      lift
#> 1  0.7173913 0.9487500
#> 2  0.5942029 0.9155340
#> 3  0.5724638 1.0276181
#> 4  0.5579710 0.8597087
#> 5  0.4927536 0.7494167
#> 6  0.4855072 0.8375000
#> 7  0.4782609 0.8250000
#> 8  0.4637681 0.8217476
#> 9  0.4275362 0.8825808
#> 10 0.4130435 0.9933936

sequence_patterns(trajectories, type = "gapped", gap = 1)
#>                     pattern length occurrences proportion n_sequences   support
#> 1       Average->*->Average      3         420 0.26200873          93 0.6739130
#> 2        Active->*->Average      3         174 0.10854648          83 0.6014493
#> 3        Average->*->Active      3         159 0.09918902          75 0.5434783
#> 4         Active->*->Active      3         409 0.25514660          68 0.4927536
#> 5    Disengaged->*->Average      3         109 0.06799750          59 0.4275362
#> 6    Average->*->Disengaged      3         118 0.07361198          57 0.4130435
#> 7 Disengaged->*->Disengaged      3         144 0.08983157          49 0.3550725
#> 8     Disengaged->*->Active      3          37 0.02308172          30 0.2173913
#> 9     Active->*->Disengaged      3          33 0.02058640          28 0.2028986
#>        lift
#> 1 0.8912500
#> 2 0.9266990
#> 3 0.8373786
#> 4 0.8845320
#> 5 0.7375000
#> 6 0.7125000
#> 7 0.7989130
#> 8 0.4368932
#> 9 0.4077670

sequence_patterns(trajectories, type = "repeated", len = 3:5,
                  start = "Disengaged")
#>                                                      pattern length occurrences
#> 1                         Disengaged->Disengaged->Disengaged      3          94
#> 2             Disengaged->Disengaged->Disengaged->Disengaged      4          57
#> 3 Disengaged->Disengaged->Disengaged->Disengaged->Disengaged      5          35
#>   proportion n_sequences    support      lift
#> 1 0.13091922          28 0.20289855 0.6847826
#> 2 0.11422846          21 0.15217391 0.7703804
#> 3 0.09887006          10 0.07246377 0.5502717

sequence_patterns(group_regulation_long, actor = "Actor",
                  action = "Action", time = "Time", len = 2,
                  group = "Achiever", top = 10)
#>                pattern length occurrences  proportion n_sequences support
#> 1       discuss->adapt      2         282 0.011044531         258  0.1290
#> 2   discuss->consensus      2        1269 0.049700388         880  0.4400
#> 3   consensus->discuss      2        1190 0.046606353         828  0.4140
#> 4     adapt->consensus      2         243 0.009517096         222  0.1110
#> 5     synthesis->adapt      2         153 0.005992245         146  0.0730
#> 6      adapt->cohesion      2         139 0.005443935         137  0.0685
#> 7   cohesion->cohesion      2          46 0.001801590          42  0.0210
#> 8  coregulate->discuss      2         539 0.021109936         452  0.2260
#> 9  cohesion->consensus      2         844 0.033055262         688  0.3440
#> 10    emotion->discuss      2         289 0.011318686         267  0.1335
#>          lift support_High support_Low odds_ratio   or_lower  or_upper
#> 1  0.68884958        0.047       0.211  5.4182988 3.87371327 7.7074607
#> 2  0.61153048        0.545       0.335  0.4207572 0.34950934 0.5060687
#> 3  0.57539458        0.513       0.315  0.4367325 0.36230806 0.5259653
#> 4  0.51218457        0.068       0.154  2.4938250 1.83313387 3.4206047
#> 5  1.01946408        0.039       0.107  2.9509995 2.00246467 4.4251357
#> 6  0.50203839        0.037       0.100  2.8904423 1.94036332 4.3856532
#> 7  0.06362668        0.037       0.005  0.1308943 0.03998606 0.3357926
#> 8  0.47258065        0.179       0.273  1.7218367 1.38415246 2.1454614
#> 9  0.65619896        0.397       0.291  0.6235728 0.51515543 0.7542389
#> 10 0.22848693        0.171       0.096  0.5149905 0.38979180 0.6774734
#>         p_value        p_adj
#> 1  2.456517e-29 1.645866e-27
#> 2  3.328919e-21 1.115188e-19
#> 3  2.782011e-19 6.213158e-18
#> 4  9.461193e-10 1.584750e-08
#> 5  4.688100e-09 6.282055e-08
#> 6  2.408664e-08 2.689675e-07
#> 7  3.434169e-07 3.286990e-06
#> 8  6.173873e-07 5.170618e-06
#> 9  7.443520e-07 5.541287e-06
#> 10 9.992119e-07 6.694720e-06

sequence_patterns(trajectories, pattern = "Active->*->Disengaged")
#>                          pattern length occurrences proportion n_sequences
#> 1    Active->Average->Disengaged      3          17  0.5151515          17
#> 2     Active->Active->Disengaged      3           9  0.2727273           9
#> 3 Active->Disengaged->Disengaged      3           7  0.2121212           7
#>      support      lift
#> 1 0.12318841 0.2847087
#> 2 0.06521739 0.1756056
#> 3 0.05072464 0.1529126
```
