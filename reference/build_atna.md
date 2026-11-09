# Build an Attention-Weighted Transition Network (ATNA)

`build_atna()` estimates an attention-weighted transition network
(ATNA): every later state in a sequence, not only the next one, adds to
the edge from an earlier state, with a weight that decays with the
number of positions between them, `exp(-distance / lambda)`. It is
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
with `method = "attention"`. The decay is set in `params`: `lambda`
(default `1`; larger values decay more slowly) and `direction`
(`"forward"`, the default, `"backward"` or `"both"`). `duration` (one
value per position, accumulated into times) or `time_matrix` (one row
per sequence, one column per position) measure the distance in time
instead of positions, and `decay`, a function of `(ti, tj, lambda)`
returning one weight per pair, replaces the decay function.

The weights are summed decay weights, not counts or probabilities, so
they compare within one network; add `scaling = "normalize"` for rows
that sum to 1. Check how far they move under resampling with
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md).

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
  ...
)
```

## Arguments

- data:

  The sequences: a data frame with one row per sequence and one column
  per position, or an event log (one row per event, with `action` set).
  Column roles are declared as in
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md).
  One-hot data and matrices are refused.

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
  `method` cannot be given: it is always `"attention"` (error of class
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
net <- build_atna(trajectories)
net
#> Attention Network (decay-weighted transitions) [directed]
#>   Weights: [15.638, 251.913]  |  mean: 107.233
#> 
#>   Weight matrix:
#>               Active Average Disengaged
#>   Active     251.098 101.100     15.638
#>   Average     89.399 251.913     76.872
#>   Disengaged  22.296  71.697     85.083 
#> 
#>   Initial probabilities:
#>   Active        0.382  ████████████████████████████████████████
#>   Average       0.368  ██████████████████████████████████████
#>   Disengaged    0.250  ██████████████████████████
```
