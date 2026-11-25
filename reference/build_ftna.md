# Estimate a frequency transition network analysis (FTNA) model

`build_network(method = "frequency")`, or for short `build_ftna()`,
estimates a frequency transition network analysis (FTNA) model from
sequences of states. The edge from state A to state B is the number of
times B comes right after A (Saqr et al., 2025).

Counts grow with the number and length of sequences, so they compare
within one network; for transition probabilities use
[`build_tna()`](https://pak.dynasite.org/Nestimate/reference/build_tna.md).
A matrix given as `data` must hold non-negative whole counts. Check how
far the counts move under resampling with
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md).

## Usage

``` r
build_ftna(
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
  codes = NULL,
  window_size = 3L,
  mode = c("non-overlapping", "overlapping"),
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
  per position, an event log (one row per event, with `action` set),
  one-hot code columns (declared with `codes`), or a square matrix of
  transition counts or weights, which is row-normalised (every row needs
  a positive sum).

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

- codes:

  The one-hot code columns (one per state, 1 when the state is active);
  setting them declares one-hot data. Default `NULL`.

- window_size:

  For one-hot data: how many consecutive rows form one window. Default
  `3`. Windows on sequence data are set through `params` instead.

- mode:

  For one-hot data: `"non-overlapping"` (default) or `"overlapping"`
  windows;
  [`build_wtna()`](https://pak.dynasite.org/Nestimate/reference/build_wtna.md)
  describes how overlapping windows are counted.

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

  Settings for the estimator, as a named list. Transitions are counted
  with tna's `weighted` and `concat` options, and `window_size` (with
  `window_span`, or `windowed = TRUE` for a size of 2) counts them
  between consecutive non-overlapping blocks of that many positions, as
  [`tna::tna()`](https://sonsoles.me/tna/reference/build_model.html)
  does; `mode = "overlapping"` then raises a warning of class
  `nestimate_ignored_mode`. The network keeps these settings, and the
  resampling verbs such as
  [`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  re-estimate with them. Default
  [`list()`](https://rdrr.io/r/base/list.html).

## Value

A `netobject` (a `netobject_group` with `group`), as
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
returns it: print it for an overview, read the edges with
[`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md)
and the counts with
[`extract_transition_matrix()`](https://pak.dynasite.org/Nestimate/reference/extract_transition_matrix.md),
and draw it with
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html).

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
[`build_tna()`](https://pak.dynasite.org/Nestimate/reference/build_tna.md),
[`build_atna()`](https://pak.dynasite.org/Nestimate/reference/build_atna.md)
and
[`build_cna()`](https://pak.dynasite.org/Nestimate/reference/build_cna.md)
for the other transition networks.

## Examples

``` r
build_network(trajectories, method = "frequency")
#> Transition Network (frequency counts) [directed]
#>   Weights: [23.000, 459.000]  |  mean: 192.667
#> 
#>   Strongest edges (9 of 9):
#>    from       to         weight
#>    Active     Active     459   
#>    Average    Average    458   
#>    Active     Average    176   
#>    Disengaged Disengaged 157   
#>    Average    Active     153   
#>    Average    Disengaged 140   
#>    Disengaged Average    129   
#>    Disengaged Active      39   
#>    Active     Disengaged  23    
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix;
#>   extract_initial_probs() the initial probabilities.
build_ftna(trajectories)
#> Transition Network (frequency counts) [directed]
#>   Weights: [23.000, 459.000]  |  mean: 192.667
#> 
#>   Strongest edges (9 of 9):
#>    from       to         weight
#>    Active     Active     459   
#>    Average    Average    458   
#>    Active     Average    176   
#>    Disengaged Disengaged 157   
#>    Average    Active     153   
#>    Average    Disengaged 140   
#>    Disengaged Average    129   
#>    Disengaged Active      39   
#>    Active     Disengaged  23    
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix;
#>   extract_initial_probs() the initial probabilities.

build_ftna(group_regulation_long, actor = "Actor", action = "Action",
           time = "Time", group = "Achiever")
#> Group Networks (2 groups, group_col: Achiever)
#> 
#>   Group  Nodes  Edges  Weights
#>   High   9      76     [2.000, 1236.000]
#>   Low    9      75     [1.000, 1356.000]
```
