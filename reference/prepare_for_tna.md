# Prepare Data for TNA Analysis

Prepare simulated or real data for use with
[`tna::tna()`](http://sonsoles.me/tna/reference/build_model.md) and
related functions. Handles various input formats and ensures the output
is compatible with TNA models.

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

  Data frame containing sequence data.

- type:

  Character. Type of input data:

  "sequences"

  :   Wide format with one row per sequence (default).

  "long"

  :   Long format with one row per action.

  "auto"

  :   Automatically detect format based on column names.

- state_names:

  Character vector. Expected state names, or NULL to extract from data.
  Default: NULL.

- id_col:

  Character. Name of ID column for long format data. Default: "Actor".

- time_col:

  Character. Name of time column for long format data. Default: "Time".

- action_col:

  Character. Name of action column for long format data. Default:
  "Action".

- validate:

  Logical. Whether to validate that all actions are in state_names.
  Default: TRUE.

## Value

A data frame ready for use with TNA functions. For "sequences" type,
returns a data frame where each row is a sequence and columns are time
points (V1, V2, ...). For "long" type, converts to wide format first.

## Details

This function performs several preparations:

1.  Converts long format to wide format if needed.

2.  Validates that all actions/states are recognized.

3.  Removes any non-sequence columns (e.g., id, metadata).

4.  Converts factors to characters.

5.  Ensures consistent column naming (V1, V2, ...).

## See also

[`wide_to_long`](https://mohsaqr.github.io/Nestimate/reference/wide_to_long.md),
[`long_to_wide`](https://mohsaqr.github.io/Nestimate/reference/long_to_wide.md)
for format conversions.

## Examples

``` r
# \donttest{
# From wide format sequences
sequences <- data.frame(
  V1 = c("A","B","C","A"), V2 = c("B","C","A","B"),
  V3 = c("C","A","B","C"), V4 = c("A","B","A","B")
)
tna_data <- prepare_for_tna(sequences, type = "sequences")
# }
```
