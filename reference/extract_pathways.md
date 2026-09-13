# Cut an Event Log into Pathways

Splits a long event log into pathways and returns them as a tidy
`data.frame`, one row per pathway, ready to tabulate against an outcome
or to feed to a sequence model.

## Usage

``` r
extract_pathways(
  data,
  action,
  group,
  order = NULL,
  type = c("unit", "segments", "anchored"),
  anchor = NULL,
  terminal = NULL,
  resolve = NULL,
  sep = " -> "
)
```

## Arguments

- data:

  A long-format `data.frame`, one row per event, already in
  chronological order within each group (or with `order` supplied).

- action:

  Name of the column holding the state / event label.

- group:

  Character vector of column names defining the unit a pathway may not
  cross – typically the actor, or actor and session.

- order:

  Optional column giving the within-group event order. When `NULL`
  (default) the existing row order is used.

- type:

  The cut to make:

  `"unit"`

  :   (default) one pathway per group. With `terminal`, the group is
      truncated at its last terminal state so the pathway ends on an
      outcome rather than on whatever followed.

  `"segments"`

  :   consecutive, non-overlapping pathways: the group is cut after
      every terminal state, so each pathway ends on one. Every event
      belongs to exactly one pathway.

  `"anchored"`

  :   one pathway per occurrence of `anchor`, closing at the next
      terminal state. Pathways may overlap and events between a terminal
      and the next anchor belong to none.

- anchor:

  State that opens a pathway when `type = "anchored"`. Ignored by the
  other cuts. Default `NULL`.

- terminal:

  States that close a pathway, the closing state included. Required for
  `"segments"` and `"anchored"`; optional for `"unit"`, where it
  truncates.

- resolve:

  Optional named list giving a resolution label to append as the
  pathway's final state. Each element is a character vector of states; a
  pathway takes the label of the **first** element whose states all
  occur in it, so order the list from most specific to least. Pathways
  matching nothing are labelled `"Unresolved"`. The label becomes the
  last element of `path` and the value of `closes`, and the raw final
  state is kept in `ends`.

- sep:

  Separator used when pasting the path. Default `" -> "`.

## Value

A `data.frame` with one row per pathway: the `group` columns, `pathway`
(an integer id within group), `path` (the states pasted with `sep`),
`length`, `opens` (first state) and `closes` (last state). With
`resolve`, `closes` holds the resolution label instead and a further
`ends` column carries the raw final state. Pathways are returned in the
order they occur.

## See also

[`sequence_compare`](https://saqr.me/Nestimate/reference/sequence_compare.md),
[`outcome_model`](https://saqr.me/Nestimate/reference/outcome_model.md)

## Examples

``` r
log <- data.frame(
  id  = c(1, 1, 1, 1, 1, 1, 2, 2, 2),
  act = c("Try", "Wrong", "Hint", "Retry", "Right", "Praise",
          "Try", "Wrong", "Retry"),
  stringsAsFactors = FALSE
)
# Each failure and its consequence:
extract_pathways(log, action = "act", group = "id", type = "anchored",
                 anchor = "Wrong", terminal = c("Right", "Wrong"))
#>   id pathway                            path length opens closes
#> 1  1       1 Wrong -> Hint -> Retry -> Right      4 Wrong  Right

# Each id as one episode:
extract_pathways(log, action = "act", group = "id")
#>   id pathway                                             path length opens
#> 1  1       1 Try -> Wrong -> Hint -> Retry -> Right -> Praise      6   Try
#> 2  2       1                            Try -> Wrong -> Retry      3   Try
#>   closes
#> 1 Praise
#> 2  Retry
```
