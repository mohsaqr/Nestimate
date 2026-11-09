# Per-Sequence Indices

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
  ...
)
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

  The actor, action and timestamp columns of a long event log, as in
  [`extract_sequences()`](https://pak.dynasite.org/Nestimate/reference/extract_sequences.md):
  setting `action` declares `x` an event log. Default `NULL`.

- ...:

  Further arguments passed to
  [`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md)
  for a long event log, such as `time_threshold` or `order`.

## Value

A data frame with one row per sequence, in the order of
[`extract_sequences()`](https://pak.dynasite.org/Nestimate/reference/extract_sequences.md),
and the columns described in *Details* (`integrative_potential` only
with `favorable`). With `meta = TRUE`, the input's covariate columns, if
it has any, come first.

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

## References

Saqr, M., Tikka, S., & López-Pernas, S. (2026). codyna: Complex Dynamic
Systems. R package version 1.0.0.

Brzinsky-Fay, C. (2018). Unused resources: Sequence and trajectory
indicators. International Symposium on Sequence Analysis and Related
Methods, Monte Verità.

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
```
