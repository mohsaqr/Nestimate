# Build a Co-occurrence Network (CNA)

`build_cna()` estimates an undirected co-occurrence network (CNA): for
every pair of positions in a sequence where both are observed, the two
states count as co-occurring once, in either order. It is
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
with `method = "co_occurrence"`, and takes sequences in a wide table, an
event log, one-hot codes (declared with `codes` or `format = "onehot"`,
and counted in windows as
[`wtna()`](https://pak.dynasite.org/Nestimate/reference/wtna.md) does)
or a square matrix of counts.

An edge counts how often two states occur together within sequences
(over pairs of positions, so a sequence can add to one edge several
times), not which came first; for order, use
[`build_tna()`](https://pak.dynasite.org/Nestimate/reference/build_tna.md).
Check how far the counts move under resampling with
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md).

## Usage

``` r
build_cna(
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
  one-hot code columns, or a square matrix of co-occurrence counts
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
  `method` cannot be given: it is always `"co_occurrence"` (error of
  class `nestimate_bad_input`).

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

[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
[`cooccurrence`](https://pak.dynasite.org/Nestimate/reference/cooccurrence.md)
for delimited-field, bipartite, and other non-sequence co-occurrence
formats.

## Examples

``` r
net <- build_cna(trajectories)
net
#> Co-occurrence Network [undirected]
#>   Weights: [663.000, 3179.000]  |  mean: 2134.000
#> 
#>   Weight matrix:
#>              Active Average Disengaged
#>   Active       3156    2783        663
#>   Average      2783    3179       2075
#>   Disengaged    663    2075        948 
```
