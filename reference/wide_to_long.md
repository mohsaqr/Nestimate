# Convert Wide Sequences to Long Format

Convert sequence data from wide format (one row per sequence, columns as
time points) to long format (one row per action).

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

  Data frame in wide format with sequences in rows.

- id_col:

  Character. Name of the ID column, or NULL to auto-generate IDs.
  Default: NULL.

- time_prefix:

  Character. Prefix for time point columns (e.g., "V" for V1, V2, ...).
  Default: "V".

- action_col:

  Character. Name of the action column in output. Default: "Action".

- time_col:

  Character. Name of the time column in output. Default: "Time".

- drop_na:

  Logical. Whether to drop NA values. Default: TRUE.

## Value

A data frame in long format with columns:

- id:

  Sequence identifier (integer).

- Time:

  Time point within the sequence (integer).

- Action:

  The action/state at that time point (character).

Any additional columns from the original data are preserved.

## Details

This function converts data from the format produced by
`simulate_sequences()` to the long format used by many TNA functions and
analyses.

## See also

[`long_to_wide`](https://mohsaqr.github.io/Nestimate/reference/long_to_wide.md)
for the reverse conversion,
[`prepare_for_tna`](https://mohsaqr.github.io/Nestimate/reference/prepare_for_tna.md)
for preparing data for TNA analysis.

## Examples

``` r
# \donttest{
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
# }
```
