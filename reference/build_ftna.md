# Build a Frequency Transition Network (FTNA)

`build_ftna()` estimates a frequency transition network (FTNA): how many
times each state is followed by each state, itself included, as raw
counts. It is
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
with `method = "frequency"`, and takes the same input as
[`build_tna()`](https://pak.dynasite.org/Nestimate/reference/build_tna.md);
a matrix must hold non-negative whole counts.

Counts grow with the number and length of sequences, so they compare
within one network; for transition probabilities use
[`build_tna()`](https://pak.dynasite.org/Nestimate/reference/build_tna.md).
Check how far the counts move under resampling with
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
  ...
)
```

## Arguments

- data:

  The sequences: a data frame with one row per sequence and one column
  per position, an event log (one row per event, with `action` set),
  one-hot code columns, or a square matrix of transition counts
  (non-negative whole numbers). Column roles are declared as in
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md).

- start:

  A start state put before every sequence, so the network shows where
  sequences begin. `FALSE` (default) adds none, `TRUE` adds one labelled
  `"Start"`, a single string adds one with that label.

- end:

  An end state put after each sequence's last observed state, so the
  network shows where sequences stop; it has no outgoing edges. `FALSE`
  (default), `TRUE` (`"End"`) or a single string, as `start`.

- actor:

  The column naming who each row belongs to (a person, a team). In an
  event log, each actor's events form a sequence, split further by
  `session` or `time_threshold`. Default `NULL`.

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
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
  such as `group`, `params`, `scaling`, `threshold` or `state_cols`.
  `method` cannot be given: it is always `"frequency"` (error of class
  `nestimate_bad_input`).

## Value

A `netobject` (a `netobject_group` with `group`), exactly as
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
returns it: print it for an overview, read the edges with
[`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md)
and the weights with
[`extract_transition_matrix()`](https://pak.dynasite.org/Nestimate/reference/extract_transition_matrix.md),
and draw it with
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html).

## See also

[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md)

## Examples

``` r
net <- build_ftna(trajectories)
net
#> Transition Network (frequency counts) [directed]
#>   Weights: [23.000, 459.000]  |  mean: 192.667
#> 
#>   Weight matrix:
#>              Active Average Disengaged
#>   Active        459     176         23
#>   Average       153     458        140
#>   Disengaged     39     129        157 
#> 
#>   Initial probabilities:
#>   Active        0.382  ████████████████████████████████████████
#>   Average       0.368  ██████████████████████████████████████
#>   Disengaged    0.250  ██████████████████████████
```
