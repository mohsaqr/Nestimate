# Convert Action Column to One-Hot Encoding

Convert a categorical Action column to one-hot (binary indicator)
columns.

## Usage

``` r
action_to_onehot(
  data,
  action_col = "Action",
  states = NULL,
  drop_action = TRUE,
  sort_states = FALSE,
  prefix = ""
)
```

## Arguments

- data:

  Data frame containing an action column.

- action_col:

  Character. Name of the action column. Default: "Action".

- states:

  Character vector or NULL. States to include as columns. If NULL, uses
  all unique values. Default: NULL.

- drop_action:

  Logical. Remove the original action column. Default: TRUE.

- sort_states:

  Logical. Sort state columns alphabetically. Default: FALSE.

- prefix:

  Character. Prefix for state column names. Default: "".

## Value

Data frame with one-hot encoded columns (0/1 integers).

## Examples

``` r
long_data <- data.frame(
  Actor = rep(1:3, each = 4),
  Time = rep(1:4, 3),
  Action = sample(c("A", "B", "C"), 12, replace = TRUE)
)
onehot_data <- action_to_onehot(long_data)
head(onehot_data)
#>   Actor Time A C B
#> 1     1    1 1 0 0
#> 2     1    2 0 1 0
#> 3     1    3 0 1 0
#> 4     1    4 1 0 0
#> 5     2    1 0 0 1
#> 6     2    2 0 1 0
```
