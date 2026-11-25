# Estimate a gap transition network analysis model

`build_network(method = "gap")`, or for short `build_gap()`, estimates a
gap transition network analysis model, which links each state to the
states up to `max_gap + 1` positions later, each pair weighted by one
over its distance, as `tna::build_model(type = "gap")` does.

The weights are raw weighted counts; add `scaling = "normalize"` for row
probabilities.

## Usage

``` r
build_gap(
  data,
  start = FALSE,
  end = FALSE,
  actor = NULL,
  action = NULL,
  time = NULL,
  session = NULL,
  order = NULL,
  time_threshold = 900,
  timezone = "UTC",
  state_cols = NULL,
  metadata_cols = NULL,
  max_gap = 1L,
  group = NULL,
  scaling = NULL,
  threshold = 0,
  labels = NULL,
  params = list(),
  ...
)
```

## Arguments

- data:

  The sequences: a data frame with one row per sequence and one column
  per position, or an event log (one row per event, with `action` set).

- start:

  A start state put before every sequence, so the network shows where
  sequences begin (an edge from it to each sequence's first state).
  `FALSE` (default) adds none, `TRUE` adds one labelled `"Start"`, and a
  single string adds one with that label. One-hot data, whose rows are
  time points rather than sequences, raises an error.

- end:

  An end state put in the cell after each sequence's last observed
  state, so the network shows where sequences stop. It has no outgoing
  edges and no self-loop, unlike
  [`mark_terminal_state()`](https://pak.dynasite.org/Nestimate/reference/mark_terminal_state.md),
  which fills all trailing missing cells with an absorbing state.
  `FALSE` (default), `TRUE` (`"End"`) or a single string, as `start`.

- actor:

  The column naming who each row belongs to (a person, a team). In an
  event log, each actor's events form a sequence, split further by
  `session` or `time_threshold`; in one-hot data, each actor's rows.
  Default `NULL`.

- action:

  The column holding each event's state; setting it declares the data an
  event log. Column names give no role: `action`, `actor`, `time` and
  `session` are read only when set. A column called `Action` with
  nothing declared raises an error of class `nestimate_bad_input`, and
  an undeclared time or session column beside `action` is named in a
  message of class `nestimate_column_role`. Default `NULL`.

- time:

  The timestamp column of an event log. Events are ordered by it, and a
  gap longer than `time_threshold` starts a new sequence. Default
  `NULL`.

- session:

  The column identifying sessions in an event log; each session is its
  own sequence, also without `actor`. Default `NULL`.

- order:

  The column that orders events sharing a timestamp; without it, those
  events keep their row order. Default `NULL`.

- time_threshold:

  The longest gap, in seconds, between two events of the same sequence
  in an event log; a longer gap starts a new sequence. `FALSE` turns the
  split off, so each actor (or actor-session) is one sequence. Default
  `900`.

- timezone:

  The Olson time zone for timestamps without an offset (passed to
  [`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md));
  timestamps with an offset, such as `...Z` or `+02:00`, are converted
  from it. Default `"UTC"`.

- state_cols:

  The columns of a wide table that are sequence positions; every other
  column goes to the network's metadata. Declared columns are positions
  even when they hold only 0 and 1. With nothing declared, the frame
  must be one column or hold position columns named as a numbered run (a
  shared prefix and a number, `T1`, `T2`, ..., or `1`, `2`, ...);
  numeric columns beside the run and character or factor columns outside
  it are then covariates. Any other frame, and columns of only 0 and 1,
  raise an error of class `nestimate_bad_input`: nothing says whether
  they are positions or one-hot codes. Default `NULL`.

- metadata_cols:

  Columns to keep as the network's metadata (covariates), never as
  positions; they cannot overlap `state_cols`. A column left out of a
  numbered run of positions is named in a message of class
  `nestimate_column_role`, and listing it here or in `state_cols`
  settles it. Default `NULL`.

- max_gap:

  How many states may lie between the two states of a pair, a whole
  number of at least 0. Default `1`.

- group:

  A column to split the data by: one network per value, returned as a
  `netobject_group`. Default `NULL`.

- scaling:

  Rescaling applied to the weights after estimation, in the order given:
  `"minmax"` rescales the weights to run from 0 to 1, over every cell
  for the transition methods (tna's rule, so with no zero cell the
  smallest weight becomes 0) and over the non-zero weights, into (0,
  1\], for the other methods, `"max"` divides by the largest absolute
  weight, `"rank"` replaces the non-zero weights by their ranks, and
  `"normalize"` makes each row's absolute weights sum to 1; an
  undirected network whose matrix it makes asymmetric becomes directed,
  so both of a pair's values are kept as edges. Zero weights stay zero
  under every option. Combine them as in `c("rank", "minmax")`. Default
  `NULL` (none).

- threshold:

  Weights whose absolute value is below this are set to zero. Default
  `0` (none).

- labels:

  New display labels for the nodes, applied after estimation: a
  two-column data frame (`name`, `label`), a named character vector
  `c(name = "label")` or a named list. It rewrites the node labels and
  the weight matrix's dimnames; nodes not listed keep their names.
  Default `NULL`.

- params:

  Settings for the estimator, as a named list, as for
  [`build_tna()`](https://pak.dynasite.org/Nestimate/reference/build_tna.md);
  also `duration` (one value per position, accumulated into times) or
  `time_matrix` (one row per sequence, one column per position), which
  measure the distance in time instead of positions. Default
  [`list()`](https://rdrr.io/r/base/list.html).

- ...:

  Further arguments for
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md):
  `format` (`"auto"`, `"wide"`, `"long"` or `"onehot"`, to force a
  reading of `data`), and for event logs
  [`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md)'s
  time options `unix_time_unit`, `is_unix_time` and `custom_format`.
  `method` cannot be given: each alias has its own (error of class
  `nestimate_bad_input`).

## Value

A `netobject` (a `netobject_group` with `group`), as
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
returns it.

## References

Tikka, S., López-Pernas, S., & Saqr, M. (2025). tna: An R package for
transition network analysis. *Applied Psychological Measurement*.
[doi:10.1177/01466216251348840](https://doi.org/10.1177/01466216251348840)

## See also

[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
for every method;
[`build_ngram()`](https://pak.dynasite.org/Nestimate/reference/build_ngram.md)
and
[`build_reverse()`](https://pak.dynasite.org/Nestimate/reference/build_reverse.md)
for the other tna model types.

## Examples

``` r
build_network(trajectories, method = "gap", max_gap = 2)
#> Network (method: gap) [directed]
#>   Weights: [51.167, 783.333]  |  mean: 336.315
#> 
#>   Strongest edges (9 of 9):
#>    from       to         weight 
#>    Average    Average    783.333
#>    Active     Active     783.167
#>    Active     Average    321.333
#>    Average    Active     284.500
#>    Disengaged Disengaged 261.667
#>    Average    Disengaged 243.333
#>    Disengaged Average    227.833
#>    Disengaged Active      70.500
#>    Active     Disengaged  51.167 
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix;
#>   extract_initial_probs() the initial probabilities.
build_gap(trajectories, max_gap = 2)
#> Network (method: gap) [directed]
#>   Weights: [51.167, 783.333]  |  mean: 336.315
#> 
#>   Strongest edges (9 of 9):
#>    from       to         weight 
#>    Average    Average    783.333
#>    Active     Active     783.167
#>    Active     Average    321.333
#>    Average    Active     284.500
#>    Disengaged Disengaged 261.667
#>    Average    Disengaged 243.333
#>    Disengaged Average    227.833
#>    Disengaged Active      70.500
#>    Active     Disengaged  51.167 
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix;
#>   extract_initial_probs() the initial probabilities.
```
