# Convert Long Format to Wide Sequences

Convert sequence data from long format (one row per action) to wide
format (one row per sequence, columns as time points).

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

  Data frame in long format.

- id_col:

  Character. Name of the column identifying sequences. Default: "Actor".

- time_col:

  Character. Name of the column identifying time points. Default:
  "Time".

- action_col:

  Character. Name of the column containing actions/states. Default:
  "Action".

- time_prefix:

  Character. Prefix for time point columns in output. Default: "V".

- fill_na:

  Logical. Whether to fill missing time points with NA. Default: TRUE.

## Value

A data frame in wide format, one row per sequence: the `id_col` column
followed by the time point columns `V1`, `V2`, ... (named with
`time_prefix`) holding the action at each time point. With
`fill_na = TRUE` short sequences are padded with `NA` so every row
shares the same columns; with `fill_na = FALSE` only the time points
present in every sequence are kept.

## Details

Converts long format data (one row per action) to the wide format
expected by
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md),
[`tna::tna()`](http://sonsoles.me/tna/reference/build_model.md) and
related functions.

If `time_col` contains non-integer values (e.g., timestamps), the
function will use the ordering within each sequence to create time
indices.

## See also

[`wide_to_long`](https://saqr.me/Nestimate/reference/wide_to_long.md)
for the reverse conversion,
[`prepare_for_tna`](https://saqr.me/Nestimate/reference/prepare_for_tna.md)
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
#> 1     1  C  B  A  C
#> 2     2  C  C  B  A
#> 3     3  B  A  B  C
```
