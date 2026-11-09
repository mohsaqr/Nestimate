# Convert Long Format to Wide Sequences

`long_to_wide()` turns a long event table (one row per action) into wide
sequences, with one row per sequence and one column per position (`V1`,
`V2`, ...). Each sequence's actions are ordered by `time_col`, and the
k-th action fills column k.

## Usage

``` r
long_to_wide(
  data,
  id_col = "Actor",
  time_col = "Time",
  action_col = "Action",
  time_prefix = "V",
  fill_na = TRUE
)
```

## Arguments

- data:

  A data frame in long format, one row per action.

- id_col:

  The column identifying each sequence. Default `"Actor"`; a missing
  column raises an error.

- time_col:

  The column that orders the actions within a sequence. Default
  `"Time"`. Its values only order the actions; they do not set the
  column numbers, so gaps between time points are not kept. Actions with
  the same time keep their row order. When the column is absent, the row
  order is used.

- action_col:

  The column holding the actions (states). Default `"Action"`; a missing
  column raises an error.

- time_prefix:

  The prefix of the position columns in the output. Default `"V"`.

- fill_na:

  How to handle sequences of different lengths. With `TRUE` (default),
  shorter sequences are padded with `NA`, so every row has the columns
  of the longest sequence. With `FALSE`, only the positions present in
  every sequence are kept, so longer sequences are cut to the length of
  the shortest.

## Value

A data frame in wide format, one row per sequence in the sorted order of
`id_col`: the `id_col` column followed by the position columns (`V1`,
`V2`, ..., named with `time_prefix`) holding the action at each
position. Other columns of `data` are dropped. A missing action stays a
missing cell.

## Details

The wide format is the one
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
[`tna::tna()`](https://sonsoles.me/tna/reference/build_model.html) and
related functions read without column declarations, since the position
columns form a numbered run;
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
sets the `id_col` column aside as metadata, with a message.

## See also

[`wide_to_long`](https://pak.dynasite.org/Nestimate/reference/wide_to_long.md)
for the reverse conversion,
[`prepare_for_tna`](https://pak.dynasite.org/Nestimate/reference/prepare_for_tna.md)
for preparing data for TNA analysis.

## Examples

``` r
long_data <- data.frame(
  Actor = rep(1:3, each = 4),
  Time = rep(1:4, 3),
  Action = sample(c("A", "B", "C"), 12, replace = TRUE)
)
wide_data <- long_to_wide(long_data, id_col = "Actor")
head(wide_data)
#>   Actor V1 V2 V3 V4
#> 1     1  C  C  A  C
#> 2     2  B  A  B  A
#> 3     3  B  C  A  B
```
