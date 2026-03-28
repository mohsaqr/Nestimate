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

A data frame in wide format where each row is a sequence and columns V1,
V2, ... contain the actions at each time point.

## Details

This function converts long format data (like that from
`simulate_long_data()`) to the wide format expected by
[`tna::tna()`](http://sonsoles.me/tna/reference/build_model.md) and
related functions.

If `time_col` contains non-integer values (e.g., timestamps), the
function will use the ordering within each sequence to create time
indices.

## See also

[`wide_to_long`](https://mohsaqr.github.io/Nestimate/reference/wide_to_long.md)
for the reverse conversion,
[`prepare_for_tna`](https://mohsaqr.github.io/Nestimate/reference/prepare_for_tna.md)
for preparing data for TNA analysis.

## Examples

``` r
# \donttest{
long_data <- data.frame(
  Actor = rep(1:3, each = 4),
  Time = rep(1:4, 3),
  Action = sample(c("A", "B", "C"), 12, replace = TRUE)
)
wide_data <- long_to_wide(long_data, id_col = "Actor")
head(wide_data)
#>   Actor V1 V2 V3 V4
#> 1     1  B  A  C  B
#> 2     2  A  B  A  C
#> 3     3  C  A  C  A
# }
```
