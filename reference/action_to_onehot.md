# Convert Action Column to One-Hot Encoding

`action_to_onehot()` adds one 0/1 indicator column per state to a data
frame with a categorical action column: the column of a row's state is 1
and the others are 0. The result is the one-hot layout that
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
reads with `codes` or `format = "onehot"`, and
[`wtna()`](https://pak.dynasite.org/Nestimate/reference/wtna.md) reads
directly.

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

  A data frame with an action column.

- action_col:

  The action column. Default `"Action"`; a name that is not a column
  raises an error.

- states:

  The states to make columns for, in that order. Default `NULL`, every
  non-missing value of the action column in order of first appearance. A
  listed state that never occurs gets a column of zeros, and a row whose
  state is not listed is 0 in every column.

- drop_action:

  Whether to remove the action column. Default `TRUE`.

- sort_states:

  Whether to sort the state columns alphabetically (in level order for a
  factor action column). Default `FALSE`.

- prefix:

  A prefix for the state column names. Without one, a state named like
  an existing column replaces that column. Default `""`.

## Value

The input data frame with one integer column per state appended, named
`paste0(prefix, state)`. All other columns are kept; the action column
is removed unless `drop_action = FALSE`. A row with a missing action is
`NA` in every state column.

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
