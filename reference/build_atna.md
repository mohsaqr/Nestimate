# Estimate an attention transition network analysis (ATNA) model

`build_network(method = "attention")`, or for short `build_atna()`,
estimates an attention transition network analysis (ATNA) model, in
which every later state in a sequence counts, not only the next one. A
pair of positions d steps apart adds `exp(-d / lambda)` to the edge from
the earlier state to the later one, so near successors weigh most (Saqr
et al., 2025).

The weights are summed decay weights, not counts or probabilities, so
they compare within one network; add `scaling = "normalize"` for rows
that sum to 1. Attention networks read sequences and event logs, not
one-hot data.

## Usage

``` r
build_atna(
  data,
  start = FALSE,
  end = FALSE,
  actor = NULL,
  action = NULL,
  time = NULL,
  session = NULL,
  ...,
  order = NULL,
  time_threshold = 900,
  timezone = "UTC",
  state_cols = NULL,
  metadata_cols = NULL,
  lambda = 1,
  direction = c("forward", "backward", "both"),
  decay = NULL,
  group = NULL,
  scaling = NULL,
  threshold = 0,
  labels = NULL,
  params = list()
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

- ...:

  Further arguments for
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md):
  `format` (`"auto"`, `"wide"`, `"long"` or `"onehot"`, to force a
  reading of `data`), and for event logs
  [`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md)'s
  time options `unix_time_unit`, `is_unix_time` and `custom_format`.
  `method` cannot be given: each alias has its own (error of class
  `nestimate_bad_input`).

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

- lambda:

  The decay length: larger values let distant states count more. Default
  `1`.

- direction:

  `"forward"` (default) links each state to the later ones, `"backward"`
  to the earlier ones, and `"both"` to all.

- decay:

  A function of `(ti, tj, lambda)` returning one weight per pair of
  positions, in place of the exponential decay. Default `NULL`.

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

## Value

A `netobject` (a `netobject_group` with `group`), as
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
returns it.

## References

Saqr, M., López-Pernas, S., Törmänen, T., Kaliisa, R., Misiejuk, K., &
Tikka, S. (2025). Transition network analysis: A novel framework for
modeling, visualizing, and identifying the temporal patterns of learners
and learning processes. In *Proceedings of the 15th International
Learning Analytics and Knowledge Conference (LAK25)* (pp. 351–361). ACM.
[doi:10.1145/3706468.3706513](https://doi.org/10.1145/3706468.3706513)

## See also

[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
for every method;
[`build_tna()`](https://pak.dynasite.org/Nestimate/reference/build_tna.md)
for the transitions to the next state only.

## Examples

``` r
build_network(trajectories, method = "attention")
#> Attention Network (decay-weighted transitions) [directed]
#>   Weights: [15.638, 251.913]  |  mean: 107.233
#> 
#>   Strongest edges (9 of 9):
#>    from       to         weight 
#>    Average    Average    251.913
#>    Active     Active     251.098
#>    Active     Average    101.100
#>    Average    Active      89.399
#>    Disengaged Disengaged  85.083
#>    Average    Disengaged  76.872
#>    Disengaged Average     71.697
#>    Disengaged Active      22.296
#>    Active     Disengaged  15.638 
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix;
#>   extract_initial_probs() the initial probabilities.
build_atna(trajectories)
#> Attention Network (decay-weighted transitions) [directed]
#>   Weights: [15.638, 251.913]  |  mean: 107.233
#> 
#>   Strongest edges (9 of 9):
#>    from       to         weight 
#>    Average    Average    251.913
#>    Active     Active     251.098
#>    Active     Average    101.100
#>    Average    Active      89.399
#>    Disengaged Disengaged  85.083
#>    Average    Disengaged  76.872
#>    Disengaged Average     71.697
#>    Disengaged Active      22.296
#>    Active     Disengaged  15.638 
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix;
#>   extract_initial_probs() the initial probabilities.

build_atna(trajectories, lambda = 2, direction = "both")
#> Attention Network (decay-weighted transitions) [directed]
#>   Weights: [101.161, 1217.841]  |  mean: 526.043
#> 
#>   Strongest edges (9 of 9):
#>    from       to         weight  
#>    Average    Average    1217.841
#>    Active     Active     1211.021
#>    Average    Active      482.735
#>    Active     Average     482.735
#>    Disengaged Disengaged  401.810
#>    Disengaged Average     367.961
#>    Average    Disengaged  367.961
#>    Disengaged Active      101.161
#>    Active     Disengaged  101.161 
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix;
#>   extract_initial_probs() the initial probabilities.
```
