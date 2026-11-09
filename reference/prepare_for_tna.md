# Prepare Data for TNA Analysis

`prepare_for_tna()` returns the sequence columns of a data set as a wide
data frame for
[`tna::tna()`](https://sonsoles.me/tna/reference/build_model.html) and
related functions: one row per sequence and the position columns `V1`,
`V2`, ... in numeric order, with factors turned into character. Long
input is converted with
[`long_to_wide()`](https://pak.dynasite.org/Nestimate/reference/long_to_wide.md)
first, padding shorter sequences with `NA`. Every other column (an id,
metadata) is dropped.

## Usage

``` r
prepare_for_tna(
  data,
  type = c("sequences", "long", "auto"),
  state_names = NULL,
  id_col = "Actor",
  time_col = "Time",
  action_col = "Action",
  validate = TRUE
)
```

## Arguments

- data:

  A data frame of sequences: wide, with position columns named `V1`,
  `V2`, ..., or long, with one row per action.

- type:

  How to read `data`. `"sequences"` (default) reads wide position
  columns; `"long"` converts a long table with
  [`long_to_wide()`](https://pak.dynasite.org/Nestimate/reference/long_to_wide.md),
  which needs `id_col` and `action_col`; `"auto"` reads long data when
  `action_col` is a column and no `V<number>` column exists, wide data
  in the opposite case, and, when both exist, long data only if `id_col`
  is a column and there are more than twice as many rows as ids. Data
  with neither raises an error. Whatever the type, the result needs
  `V<number>` columns; others, such as `T1`, `T2`, are not renamed and
  raise an error.

- state_names:

  The expected states, used with `validate`. Default `NULL`, which
  accepts the states found in the data.

- id_col:

  The column identifying sequences in long data. Default `"Actor"`.

- time_col:

  The column ordering the actions in long data. Default `"Time"`.

- action_col:

  The column holding the actions in long data. Default `"Action"`.

- validate:

  Whether to check the states against `state_names`; states not listed
  are named in a warning. Ignored when `state_names` is `NULL`. Default
  `TRUE`.

## Value

A data frame with one row per sequence and only the position columns
`V1`, `V2`, ..., in numeric order, holding character states (or the
original type for non-factor columns), with no row names.

## See also

[`wide_to_long`](https://pak.dynasite.org/Nestimate/reference/wide_to_long.md),
[`long_to_wide`](https://pak.dynasite.org/Nestimate/reference/long_to_wide.md)
for format conversions.

## Examples

``` r
# From wide format sequences
sequences <- data.frame(
  V1 = c("A","B","C","A"), V2 = c("B","C","A","B"),
  V3 = c("C","A","B","C"), V4 = c("A","B","A","B")
)
tna_data <- prepare_for_tna(sequences, type = "sequences")
```
