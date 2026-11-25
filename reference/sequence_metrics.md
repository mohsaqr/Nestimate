# Compute descriptive indices for each sequence

`sequence_metrics()` computes descriptive indices for every sequence and
returns one row per sequence: length and validity, composition and
diversity, spells, transitions, the fate of the initial state, the
dominant and emergent states, cyclic recurrence and a composite
complexity index, plus the integrative potential of a set of favourable
states when `favorable` is given. It reads any input
[`extract_sequences()`](https://pak.dynasite.org/Nestimate/reference/extract_sequences.md)
reads: a wide table of positions, an event log, the output of
[`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md),
a sequence network or group of networks, or a `tna` model. The indices
are those of `codyna::sequence_indices()`, with three corrections listed
under *Details*.

The indices describe each sequence as it was observed, and several of
them (spells, transitions, recurrence) rest on few positions in a short
sequence; read them beside `valid_n`. Because there is one row per
sequence, the result can be summarised, plotted or used as predictors;
[`sequence_patterns()`](https://pak.dynasite.org/Nestimate/reference/sequence_patterns.md)
and
[`pattern_outcome()`](https://pak.dynasite.org/Nestimate/reference/pattern_outcome.md)
look at which subsequences recur and how they relate to an outcome.

## Usage

``` r
sequence_metrics(
  x,
  favorable = NULL,
  omega = 1,
  meta = FALSE,
  cols = NULL,
  actor = NULL,
  action = NULL,
  time = NULL,
  group = NULL,
  indices = NULL,
  iter = 1000L,
  adjust = "BH",
  format = c("long", "wide"),
  ...
)

# S3 method for class 'net_sequence_metrics'
print(x, n = 10L, digits = 2L, ...)

# S3 method for class 'net_sequence_metrics'
summary(object, ...)

# S3 method for class 'net_sequence_metrics'
x[i, j, ...]
```

## Arguments

- x:

  The sequence data, in any form
  [`extract_sequences()`](https://pak.dynasite.org/Nestimate/reference/extract_sequences.md)
  accepts.

- favorable:

  The favourable states; giving them adds `integrative_potential`.
  States that do not occur in the data are ignored, and an error of
  class `nestimate_bad_input` is raised when none of them occur. Default
  `NULL`.

- omega:

  The exponent of the position weights in `integrative_potential`; a
  positive value gives later positions more weight, `0` weights every
  position equally and a negative value gives earlier positions more
  weight; any single finite number is accepted. Ignored without
  `favorable`. Default `1`.

- meta:

  Whether to put each sequence's covariates (actor, session, group,
  metadata) before the indices. Default `FALSE`.

- cols:

  The state columns of a wide data frame, as in
  [`extract_sequences()`](https://pak.dynasite.org/Nestimate/reference/extract_sequences.md).
  Default `NULL`.

- actor, action, time:

  The actor, action and timestamp columns of a long event log; setting
  `action` marks `x` as one. Default `NULL`. With `group` and sequences
  that are not an event log (a network, a wide table), `actor` instead
  names the covariate that nests the sequences (the session of a
  half-session, the student of a course), and the permutation test
  shuffles the group labels within it (see *Comparing groups*). For a
  network built with `actor`, the default nests the sequences in that
  actor (its `.session_label`), so the two halves of one session are
  compared within the session.

- group:

  A covariate of the sequences holding each sequence's group: for a
  grouped network its group column (for example `"phase"` after
  [`mark_phases()`](https://pak.dynasite.org/Nestimate/reference/mark_phases.md)),
  for a wide table one of its non-state columns. With `group`, the
  result compares the indices between the groups instead of listing them
  per sequence. A grouped network is compared by its own groups by
  default; `meta = TRUE` lists its sequences instead. A name that is not
  a covariate of the sequences, or fewer than two groups, raises an
  error of class `nestimate_bad_input`. Default `NULL`.

- indices:

  The indices to return, by name (for example
  `c("valid_n", "transition_rate")`). Default `NULL` returns them all. A
  name that is not an index raises an error of class
  `nestimate_bad_input`.

- iter:

  The number of permutations of the group comparison, a whole number of
  at least 1. Default `1000`.

- adjust:

  The method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) applied
  to the p-values of the group comparison. Default `"BH"`.

- format:

  The layout of a group comparison: `"long"` (default), one row per
  index and group, or `"wide"`, one row per index with a column per
  group. Ignored without `group`.

- ...:

  Further arguments passed to
  [`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md)
  for a long event log, such as `time_threshold` or `order`.

- n:

  For [`print()`](https://rdrr.io/r/base/print.html) of a group
  comparison: the number of indices shown. Default `10`.

- digits:

  For [`print()`](https://rdrr.io/r/base/print.html): the decimals of
  the means, SDs and effect sizes. Default `2`; p-values show 3
  significant digits.

- object:

  For [`summary()`](https://rdrr.io/r/base/summary.html): a group
  comparison from `sequence_metrics(group = )`.

- i, j:

  For `[`: rows and columns, as for a data frame. A subset of a group
  comparison is a plain data frame.

## Value

A data frame with one row per sequence, in the order of
[`extract_sequences()`](https://pak.dynasite.org/Nestimate/reference/extract_sequences.md),
and the columns described in *Details* (`integrative_potential` only
with `favorable`). With `meta = TRUE`, the input's covariate columns, if
it has any, come first.

With `group`, a data frame of class
`c("net_sequence_metrics", "data.frame")` with one row per numeric index
and group: `index`, `group`, `n` (the group's sequences with a value of
the index), `mean` and `sd`; then the index's comparison, which repeats
on each of its rows: with two groups `difference` (the second group's
mean minus the first's, in sorted order), `hedges_g` and its 95%
interval `g_lower`, `g_upper`; for any number of groups `eta_squared`,
the one-way `statistic` (F), the permutation `p_value` and the adjusted
`p_adj`. With `format = "wide"`, one row per index with `n_<group>`,
`mean_<group>` and `sd_<group>` columns instead. The settings are
attributes (`groups`, `group`, `actor`, `iter`, `adjust`, `format`,
`n_indices`). [`print()`](https://rdrr.io/r/base/print.html) shows the
first `n` indices, rounded;
[`summary()`](https://rdrr.io/r/base/summary.html) returns the whole
table as a plain data frame, and a subset is a plain data frame.

## Details

Positions after a sequence's last observed state are not part of it; `p`
below is the position of the last observed state.

- valid_n, valid_proportion:

  The number of observed states, and their share of the first `p`
  positions.

- unique_states:

  The number of distinct states.

- mean_spell_duration, max_spell_duration:

  The mean and the longest run of one state.

- longitudinal_entropy:

  The Shannon entropy of the state proportions, divided by the log of
  the number of states in the data (0 to 1; 0 when the data hold one
  state).

- simpson_diversity:

  One minus the sum of squared proportions.

- self_loop_tendency, transition_rate:

  The share of transitions (pairs of adjacent observed states) that
  repeat the state, and the share that change it; the two sum to 1.

- transition_complexity:

  The number of distinct state-changing transitions, as a share of the
  \\a(a - 1)\\ possible ones, where \\a\\ is the number of states in the
  data; `NA` when the data hold one state.

- initial_state_persistence:

  The length of the opening spell as a share of `p`.

- initial_state_proportion:

  The share of the observed states that are the first state.

- initial_state_influence_decay:

  The share of the first third spent in the first state minus the share
  of the last third.

- cyclic_feedback_strength:

  The largest share, over lags of 2 or more, of positions whose state
  recurs that many steps later.

- first_state, last_state, dominant_state, dominant_proportion,
  dominant_max_spell:

  The opening and closing states; the most frequent state (ties go to
  the state that sorts first), its share and its longest spell.

- emergent_state, emergent_state_persistence, emergent_state_proportion:

  A state that settles in later, its spell and its share; `NA` when
  there is none. The candidates are a spell of at least 3 that outlasts
  the dominant spell, the dominant state when it is not the first and
  has a spell of at least 3, and a return to the first state in a spell
  of at least 3 and more than half the opening one; the longest of them
  is taken.

- integrative_potential:

  With `favorable`: the position-weighted share of favourable states
  over the observed positions, with weights `position^omega`
  (Brzinsky-Fay 2018).

- complexity_index:

  `0.4 * longitudinal_entropy + 0.4 * transition_rate + 0.2 * min(CV of spell durations, 1)`;
  0 for a sequence of one spell (one state, no change). This is codyna's
  composite, not the complexity index of Gabadinho et al. (2011) that
  `TraMineR::seqici()` computes.

Missing values: a missing position is not an observation. The first
state is the first observed one, a missing position inside a sequence
ends a spell, and a pair of states across it is not a transition.
Indices that are undefined for a sequence (for example a rate with no
transition) are `NA`; a sequence with no observed state has `valid_n`
and `unique_states` 0 and every other index `NA`. On sequences without
internal gaps the indices shared with TraMineR (`seqindic()`,
`seqistatd()`, `seqdur()`, `seqintegr()`) agree exactly; at an internal
gap TraMineR joins the states on either side into one sequence, so
spells and transitions differ there.

Corrections relative to `codyna::sequence_indices()` 1.0.0, which
otherwise agrees value for value: `dominant_max_spell` is computed per
sequence (codyna reports the last sequence's value for every row);
`longitudinal_entropy` is normalised by the log of the number of states
(codyna counts `NA` as one when the data have missing values); and
`complexity_index` uses the normalised entropy once (codyna divides it
by the log of the number of states a second time) and is 0, not `NA`,
for a single spell. Leading and internal missing values follow the rule
above, where codyna reads position 1 as the first state and counts
missing positions as unfavourable and as possible transitions.

## Comparing groups

For each numeric index, the groups are compared with the one-way F
statistic, which for two groups is the square of the t statistic. Its
p-value is `(b + 1) / (iter + 1)`, where `b` counts the permutations of
the group labels whose F is at least the observed one (Good, 2005); call
[`base::set.seed()`](https://rdrr.io/r/base/Random.html) first for
reproducible p-values. A sequence with a missing value of an index (a
one-state sequence has no transition rate) takes no part in that index.
With `actor`, the labels are shuffled within each actor whose sequences
are in several groups and between the actors whose sequences are all in
one group, as in
[`sequence_compare()`](https://pak.dynasite.org/Nestimate/reference/sequence_compare.md).
With two groups, `hedges_g` is the standardised mean difference
corrected for small-sample bias, with an interval from its normal
approximation (Hedges, 1981). `eta_squared` is the share of the index's
variance between the groups.

## References

Saqr, M., Tikka, S., & López-Pernas, S. (2026). codyna: Complex Dynamic
Systems. R package version 1.0.0.

Brzinsky-Fay, C. (2018). Unused resources: Sequence and trajectory
indicators. International Symposium on Sequence Analysis and Related
Methods, Monte Verità.

Good, P. (2005). *Permutation, parametric and bootstrap tests of
hypotheses*. Springer.
[doi:10.1007/b138696](https://doi.org/10.1007/b138696)

Hedges, L. V. (1981). Distribution theory for Glass's estimator of
effect size and related estimators. *Journal of Educational Statistics*,
6(2), 107-128.
[doi:10.3102/10769986006002107](https://doi.org/10.3102/10769986006002107)

Gabadinho, A., Ritschard, G., Müller, N. S., & Studer, M. (2011).
Analyzing and visualizing state sequences in R with TraMineR. Journal of
Statistical Software, 40(4), 1-37.
[doi:10.18637/jss.v040.i04](https://doi.org/10.18637/jss.v040.i04)

## See also

[`sequence_patterns`](https://pak.dynasite.org/Nestimate/reference/sequence_patterns.md),
[`extract_sequences`](https://pak.dynasite.org/Nestimate/reference/extract_sequences.md)

## Examples

``` r
head(sequence_metrics(trajectories))
#>   valid_n valid_proportion unique_states mean_spell_duration max_spell_duration
#> 1      12                1             3            1.714286                  4
#> 2      15                1             2            5.000000                  9
#> 3      15                1             2            3.750000                  8
#> 4      15                1             2            1.875000                  5
#> 5      15                1             3            1.875000                  3
#> 6      15                1             2            1.666667                  3
#>   longitudinal_entropy simpson_diversity self_loop_tendency transition_rate
#> 1            0.8734790         0.5694444          0.4545455       0.5454545
#> 2            0.2229449         0.1244444          0.8571429       0.1428571
#> 3            0.3574277         0.2311111          0.7857143       0.2142857
#> 4            0.5793802         0.4444444          0.5000000       0.5000000
#> 5            0.7312088         0.4800000          0.5000000       0.5000000
#> 6            0.5793802         0.4444444          0.4285714       0.5714286
#>   transition_complexity initial_state_persistence initial_state_proportion
#> 1             0.6666667                0.08333333                0.5833333
#> 2             0.3333333                0.60000000                0.9333333
#> 3             0.3333333                0.06666667                0.1333333
#> 4             0.3333333                0.33333333                0.6666667
#> 5             0.6666667                0.20000000                0.6666667
#> 6             0.3333333                0.13333333                0.6666667
#>   initial_state_influence_decay cyclic_feedback_strength first_state last_state
#> 1                   -0.35000000                1.0000000      Active     Active
#> 2                    0.16666667                1.0000000     Average    Average
#> 3                    0.03333333                1.0000000     Average     Active
#> 4                    0.50000000                0.6666667      Active    Average
#> 5                    0.30000000                0.8000000      Active    Average
#> 6                   -0.06666667                1.0000000     Average    Average
#>   dominant_state dominant_proportion dominant_max_spell emergent_state
#> 1         Active           0.5833333                  4         Active
#> 2        Average           0.9333333                  9        Average
#> 3         Active           0.8666667                  8         Active
#> 4         Active           0.6666667                  5         Active
#> 5         Active           0.6666667                  3         Active
#> 6        Average           0.6666667                  3        Average
#>   emergent_state_persistence emergent_state_proportion complexity_index
#> 1                          4                 0.5833333        0.7138228
#> 2                          5                 0.9333333        0.3063208
#> 3                          8                 0.8666667        0.4102016
#> 4                          3                 0.6666667        0.5872441
#> 5                          3                 0.6666667        0.5981935
#> 6                          3                 0.6666667        0.5642465

head(sequence_metrics(trajectories, favorable = "Active"))
#>   valid_n valid_proportion unique_states mean_spell_duration max_spell_duration
#> 1      12                1             3            1.714286                  4
#> 2      15                1             2            5.000000                  9
#> 3      15                1             2            3.750000                  8
#> 4      15                1             2            1.875000                  5
#> 5      15                1             3            1.875000                  3
#> 6      15                1             2            1.666667                  3
#>   longitudinal_entropy simpson_diversity self_loop_tendency transition_rate
#> 1            0.8734790         0.5694444          0.4545455       0.5454545
#> 2            0.2229449         0.1244444          0.8571429       0.1428571
#> 3            0.3574277         0.2311111          0.7857143       0.2142857
#> 4            0.5793802         0.4444444          0.5000000       0.5000000
#> 5            0.7312088         0.4800000          0.5000000       0.5000000
#> 6            0.5793802         0.4444444          0.4285714       0.5714286
#>   transition_complexity initial_state_persistence initial_state_proportion
#> 1             0.6666667                0.08333333                0.5833333
#> 2             0.3333333                0.60000000                0.9333333
#> 3             0.3333333                0.06666667                0.1333333
#> 4             0.3333333                0.33333333                0.6666667
#> 5             0.6666667                0.20000000                0.6666667
#> 6             0.3333333                0.13333333                0.6666667
#>   initial_state_influence_decay cyclic_feedback_strength first_state last_state
#> 1                   -0.35000000                1.0000000      Active     Active
#> 2                    0.16666667                1.0000000     Average    Average
#> 3                    0.03333333                1.0000000     Average     Active
#> 4                    0.50000000                0.6666667      Active    Average
#> 5                    0.30000000                0.8000000      Active    Average
#> 6                   -0.06666667                1.0000000     Average    Average
#>   dominant_state dominant_proportion dominant_max_spell emergent_state
#> 1         Active           0.5833333                  4         Active
#> 2        Average           0.9333333                  9        Average
#> 3         Active           0.8666667                  8         Active
#> 4         Active           0.6666667                  5         Active
#> 5         Active           0.6666667                  3         Active
#> 6        Average           0.6666667                  3        Average
#>   emergent_state_persistence emergent_state_proportion integrative_potential
#> 1                          4                 0.5833333            0.62820513
#> 2                          5                 0.9333333            0.08333333
#> 3                          8                 0.8666667            0.90833333
#> 4                          3                 0.6666667            0.56666667
#> 5                          3                 0.6666667            0.59166667
#> 6                          3                 0.6666667            0.00000000
#>   complexity_index
#> 1        0.7138228
#> 2        0.3063208
#> 3        0.4102016
#> 4        0.5872441
#> 5        0.5981935
#> 6        0.5642465

# compare the indices of the first and second halves of each session
halves <- mark_phases(human_long, actor = "session_id",
                      time = "order_in_session",
                      labels = c("early", "late"))
net <- build_network(halves, actor = "session_id", action = "cluster",
                     group = "phase", method = "relative")
#> Metadata aggregated per session: ties resolved by first occurrence in 'session_date' (2 sessions), 'code' (125 sessions)
#> Metadata aggregated per session: ties resolved by first occurrence in 'session_date' (2 sessions), 'code' (114 sessions)
set.seed(1)
sequence_metrics(net, iter = 99)
#> Sequence metrics by phase  [19 indices | 2 groups: early, late]
#>   permutation of the one-way F: 99 iter (BH), within session
#> 
#>  index                 group n   mean  sd    diff  g     g_lo  g_hi  eta_sq
#>  valid_n               early 425 12.44 12.00  0.40  0.03 -0.10  0.17 0.00  
#>  valid_n               late  429 12.84 11.98  0.40  0.03 -0.10  0.17 0.00  
#>  valid_proportion      early 425  1.00  0.00  0.00    NA    NA    NA   NA  
#>  valid_proportion      late  429  1.00  0.00  0.00    NA    NA    NA   NA  
#>  unique_states         early 425  2.42  0.74  0.20  0.31  0.17  0.44 0.02  
#>  unique_states         late  429  2.62  0.58  0.20  0.31  0.17  0.44 0.02  
#>  mean_spell_duration   early 425  2.21  0.93 -0.36 -0.41 -0.54 -0.27 0.04  
#>  mean_spell_duration   late  429  1.85  0.84 -0.36 -0.41 -0.54 -0.27 0.04  
#>  max_spell_duration    early 425  3.97  1.75 -0.28 -0.14 -0.28 -0.01 0.01  
#>  max_spell_duration    late  429  3.69  2.16 -0.28 -0.14 -0.28 -0.01 0.01  
#>  longitudinal_entropy  early 425  0.61  0.30  0.13  0.47  0.33  0.60 0.05  
#>  longitudinal_entropy  late  429  0.73  0.24  0.13  0.47  0.33  0.60 0.05  
#>  simpson_diversity     early 425  0.41  0.20  0.09  0.50  0.36  0.63 0.06  
#>  simpson_diversity     late  429  0.50  0.16  0.09  0.50  0.36  0.63 0.06  
#>  self_loop_tendency    early 418  0.59  0.23 -0.14 -0.60 -0.74 -0.46 0.08  
#>  self_loop_tendency    late  423  0.45  0.22 -0.14 -0.60 -0.74 -0.46 0.08  
#>  transition_rate       early 418  0.41  0.23  0.14  0.60  0.46  0.74 0.08  
#>  transition_rate       late  423  0.55  0.22  0.14  0.60  0.46  0.74 0.08  
#>  transition_complexity early 425  0.47  0.31  0.06  0.20  0.07  0.34 0.01  
#>  transition_complexity late  429  0.53  0.27  0.06  0.20  0.07  0.34 0.01  
#>  F     p    p_adj 
#>   0.24 0.01 0.0129
#>   0.24 0.01 0.0129
#>     NA NA   NA    
#>     NA NA   NA    
#>  20.03 0.01 0.0129
#>  20.03 0.01 0.0129
#>  35.24 0.01 0.0129
#>  35.24 0.01 0.0129
#>   4.38 0.01 0.0129
#>   4.38 0.01 0.0129
#>  46.76 0.01 0.0129
#>  46.76 0.01 0.0129
#>  52.75 0.01 0.0129
#>  52.75 0.01 0.0129
#>  75.70 0.01 0.0129
#>  75.70 0.01 0.0129
#>  75.70 0.01 0.0129
#>  75.70 0.01 0.0129
#>   8.65 0.01 0.0129
#>   8.65 0.01 0.0129
#>   ... and 9 more indices (18 rows); summary() returns them all
```
