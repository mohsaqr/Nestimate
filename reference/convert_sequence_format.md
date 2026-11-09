# Convert Sequence Data to Different Formats

`convert_sequence_format()` turns wide or long sequence data into one of
four data frames: per-sequence state counts, per-sequence state presence
(one-hot), the list of consecutive transitions, or each event paired
with the one before it. It reads an event log (`action` set) or a wide
table (one row per sequence, one column per position). The counts can
serve as per-sequence behaviour profiles for the association estimators
of
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
and for
[`build_mcml_pc()`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md).
The id columns (those that `id_col` names, and `.rid` for a wide table
without `id_col`) are not variables and must be set aside there
(`id_col` in
[`build_mcml_pc()`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md),
`actor` in
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md));
the `rid` column is dropped automatically.

## Usage

``` r
convert_sequence_format(
  data,
  seq_cols = NULL,
  id_col = NULL,
  action = NULL,
  time = NULL,
  format = c("frequency", "onehot", "edgelist", "follows")
)
```

## Arguments

- data:

  A data frame of sequences: an event log or a wide table.

- seq_cols:

  The state columns of a wide table; setting them declares a wide table.
  Default `NULL`: every column except `id_col` and `time`. With neither
  `action` nor `seq_cols`, the frame must be one column or hold position
  columns named as a numbered run (`T1`, `T2`, ..., or `1`, `2`, ...),
  with nothing beside them but the `id_col` (and `time`) columns;
  anything else raises an error of class `nestimate_bad_input`, as does
  a column listed both here and in `id_col` or `time`.

- id_col:

  The identifier column(s), copied to every output row. In an event log,
  each unique combination of their values is one sequence, and it is
  required: an event log without `id_col` raises an error of class
  `nestimate_bad_input` (no column is taken by position). In a wide
  table, each row is one sequence and no id is inferred from the data:
  without `id_col`, a row-index column `.rid` is added. Default `NULL`.

- action:

  The column holding each event's state; setting it declares an event
  log. Default `NULL`.

- time:

  The column that orders the events within each sequence of an event
  log; without it, events keep their row order. A wide table leaves it
  out of the state columns. Default `NULL`.

- format:

  The output format: `"frequency"` (default), the count of each state
  per sequence; `"onehot"`, whether each state occurs in the sequence
  (0/1); `"edgelist"`, each pair of consecutive states (`from`, `to`);
  `"follows"`, each event (`act`) paired with the event before it
  (`follows`).

## Value

A data frame with the id column(s) first. For `"frequency"` and
`"onehot"`: one row per sequence, a `rid` column identifying the
sequence, and one integer column per state, holding counts or 0/1. For
`"edgelist"`: one row per transition, with columns `from` and `to`. For
`"follows"`: one row per event after the first of its sequence, with
columns `act` and `follows`. An event log's rows come out sorted by
`id_col`.

## Details

A missing state (`NA` or an empty string) is dropped and the states on
either side are joined, so a sequence `A, NA, B` gives the transition
`A -> B`.
[`frequencies()`](https://pak.dynasite.org/Nestimate/reference/frequencies.md)
counts no transition across a missing state instead.

## See also

[`frequencies`](https://pak.dynasite.org/Nestimate/reference/frequencies.md)
for building transition frequency matrices.

## Examples

``` r
# Wide format input
seqs <- data.frame(V1 = c("A","B","A"), V2 = c("B","A","C"), V3 = c("A","C","B"))
convert_sequence_format(seqs, format = "frequency")
#>   .rid rid A B C
#> 1    1   1 2 1 0
#> 2    2   2 1 1 1
#> 3    3   3 1 1 1
convert_sequence_format(seqs, format = "edgelist")
#>   .rid from to
#> 1    1    A  B
#> 2    1    B  A
#> 3    2    B  A
#> 4    2    A  C
#> 5    3    A  C
#> 6    3    C  B
```
