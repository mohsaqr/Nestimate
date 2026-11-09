# Convert Wide Sequences to Long Format

`wide_to_long()` turns sequences in wide format (one row per sequence,
one column per time point) into long format, with one row per sequence
and time point. The time point columns are those named `time_prefix`
followed by a number (`V1`, `V2`, ...); every other column is kept and
repeated on each row of its sequence.

## Usage

``` r
wide_to_long(
  data,
  id_col = NULL,
  time_prefix = "V",
  action_col = "Action",
  time_col = "Time",
  drop_na = TRUE
)
```

## Arguments

- data:

  A data frame in wide format, one row per sequence.

- id_col:

  The column identifying each sequence. Default `NULL`, which adds a
  column `id` holding the row number of the wide input (replacing any
  existing column of that name). A name that is not a column raises an
  error.

- time_prefix:

  The prefix of the time point columns: with `"V"` (default), the
  columns `V1`, `V2`, ... are read. No column matching it raises an
  error.

- action_col:

  The name of the state column in the output. Default `"Action"`.

- time_col:

  The name of the time column in the output. Default `"Time"`.

- drop_na:

  Whether to drop the rows whose state is missing, such as the padding
  after a shorter sequence ends. Default `TRUE`.

## Value

A data frame in long format, one row per sequence and time point, sorted
by identifier and then time. Its columns are the identifier, the other
kept columns, the time and the state:

- id:

  The sequence identifier, named by `id_col`; with `id_col = NULL` it is
  called `id` and holds the row number of the wide input (integer).

- Time:

  The time point within the sequence (integer), taken from the number in
  the wide column name; named by `time_col`.

- Action:

  The state at that time point; named by `action_col`.

## Details

The long format is the event-log layout read by
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
with `action` (and
[`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md)):
give `action = "Action"`, `actor = "id"` and `time = "Time"` (or the
names chosen here).

## See also

[`long_to_wide`](https://pak.dynasite.org/Nestimate/reference/long_to_wide.md)
for the reverse conversion,
[`prepare_for_tna`](https://pak.dynasite.org/Nestimate/reference/prepare_for_tna.md)
for preparing data for TNA analysis.

## Examples

``` r
wide_data <- data.frame(
  V1 = c("A", "B", "C"), V2 = c("B", "C", "A"), V3 = c("C", "A", "B")
)
long_data <- wide_to_long(wide_data)
head(long_data)
#>   id Time Action
#> 1  1    1      A
#> 2  1    2      B
#> 3  1    3      C
#> 4  2    1      B
#> 5  2    2      C
#> 6  2    3      A
```
