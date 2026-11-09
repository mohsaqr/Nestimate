# Build a Transition Network (TNA)

`build_tna()` estimates a transition network (TNA): the probability that
each state is followed by each state, itself included. It is
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
with `method = "relative"`: transition counts normalised to sum to 1 in
each row.

An edge from A to B is the share of A's outgoing transitions that go to
B, so a state's edges compare within its row, not across rows. The
weights are point estimates: check how far they move under resampling
with
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md),
and compare groups with
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md).

## Usage

``` r
build_tna(
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
  one-hot code columns (declared with `codes` or `format = "onehot"`),
  or a square matrix of transition counts or weights, which is
  row-normalised (every row needs a positive sum).

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
  `method` cannot be given: it is always `"relative"` (error of class
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
net <- build_tna(trajectories)
net
#> Transition Network (relative probabilities) [directed]
#>   Weights: [0.035, 0.698]  |  mean: 0.333
#> 
#>   Weight matrix:
#>              Active Average Disengaged
#>   Active      0.698   0.267      0.035
#>   Average     0.204   0.610      0.186
#>   Disengaged  0.120   0.397      0.483 
#> 
#>   Initial probabilities:
#>   Active        0.382  ████████████████████████████████████████
#>   Average       0.368  ██████████████████████████████████████
#>   Disengaged    0.250  ██████████████████████████
```
