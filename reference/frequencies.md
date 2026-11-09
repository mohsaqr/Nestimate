# Build a Transition Frequency Matrix

`frequencies()` counts how often each state is followed by each other
state across all sequences and returns the counts as a square integer
matrix. It reads an event log (one row per event, `action` set) or a
wide table (one row per sequence, one column per position).

## Usage

``` r
frequencies(
  data,
  action = NULL,
  id = NULL,
  time = NULL,
  cols = NULL,
  format = c("auto", "long", "wide")
)

# S3 method for class 'nest_transition_counts'
summary(object, ...)
```

## Arguments

- data:

  A data frame of sequences: an event log or a wide table.

- action:

  The column holding each event's state; setting it declares an event
  log. Default `NULL`. Column names give no role: a column called
  `Action` with nothing declared raises an error of class
  `nestimate_bad_input`, and so does an `action` not in `data`.

- id:

  The column(s) identifying sequences. Default `NULL`. In an event log,
  each unique combination of their values is one sequence, and without
  `id` all rows form a single sequence. In a wide table, they are left
  out of the state columns. A column not in `data` raises an error of
  class `nestimate_bad_input`.

- time:

  The column that orders the events of each sequence in an event log.
  Default `NULL`: events keep their row order. A wide table does not use
  it and leaves it out of the positions. A column not in `data`, or one
  also listed in `cols`, raises an error of class `nestimate_bad_input`.

- cols:

  The state columns of a wide table; setting them declares a wide table.
  Default `NULL`. With neither `action` nor `cols`, the frame must be
  one column or hold position columns named as a numbered run (`T1`,
  `T2`, ..., or `1`, `2`, ...), with nothing beside them but the `id`
  (and `time`) columns; anything else raises an error of class
  `nestimate_bad_input`, as does a listed column not in `data`. A single
  column is read as one position per row, so it counts no transitions;
  declare an event log with `action`.

- format:

  How to read `data`. `"auto"` (default) reads an event log when
  `action` is set and a wide table otherwise. `"long"` forces an event
  log and needs `action`. `"wide"` forces a wide table and reads every
  column outside `id` and `time` (or the columns in `cols`) as a
  position, without the check on undeclared columns; it ignores
  `action`, so that column is counted as a position unless `cols` leaves
  it out. `"long"` ignores `cols`.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `nest_transition_counts` matrix from `frequencies()`.

- ...:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method:
  ignored.

## Value

A square integer matrix of transition counts, of class
`c("nest_transition_counts", "matrix", "array")`, where `mat[i, j]` is
the number of times state i was followed by state j; the row and column
names are the sorted unique states. It behaves as an ordinary matrix and
can be passed directly to
[`tna::tna()`](https://sonsoles.me/tna/reference/build_model.html). Its
[`summary()`](https://rdrr.io/r/base/summary.html) returns the same
counts as a data frame with one row per non-zero transition and columns
`from`, `to` and `count`.

## Details

In an **event log**, the sequences are defined by the `id` column(s),
and the events of each sequence are ordered by `time` (events sharing a
time keep their row order) or, without `time`, by row order. Each pair
of consecutive events is one transition.

In a **wide table**, each row is a sequence and the state columns are
consecutive positions; each pair of adjacent positions is one
transition.

In both, a missing state (`NA` or an empty string) breaks the sequence:
no transition is counted into or out of it, and the states on either
side are not joined.
[`convert_sequence_format()`](https://pak.dynasite.org/Nestimate/reference/convert_sequence_format.md)
closes such gaps instead.

## See also

[`convert_sequence_format`](https://pak.dynasite.org/Nestimate/reference/convert_sequence_format.md)
for converting to other representations (frequency counts, one-hot, edge
lists).

## Examples

``` r
# Wide format
seqs <- data.frame(V1 = c("A","B","A"), V2 = c("B","A","C"), V3 = c("A","C","B"))
freq <- frequencies(seqs, format = "wide")

# Long format
long <- data.frame(
  Actor = rep(1:2, each = 3), Time = rep(1:3, 2),
  Action = c("A","B","C","B","A","C")
)
freq <- frequencies(long, action = "Action", id = "Actor")
```
