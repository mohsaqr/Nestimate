# Cut an Event Log into Pathways

`extract_pathways()` cuts a long event log into pathways, runs of
consecutive events within one unit (an actor, or an actor's session),
and returns them as a data frame with one row per pathway, ready to
tabulate against an outcome or to pass to a sequence model. Three cuts
are available through `type`: each unit as one pathway, consecutive
segments that each end on a terminal state, or one pathway from every
occurrence of an `anchor` state to the next terminal state.

With `terminal`, the cuts keep only what reaches a terminal state: a
unit without one, events after a unit's last terminal state, and anchors
with no later terminal state are left out without a message, so count
the units and events kept before comparing pathway frequencies across
groups.

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

  A long data frame with one row per event, in time order within each
  unit unless `order` is given. A unit's rows need not be contiguous:
  rows of other units in between are not part of its pathways.

- action:

  The name of the column holding each event's state.

- group:

  The names of the columns defining the unit a pathway may not cross,
  typically the actor, or the actor and the session.

- order:

  The name of a column giving the order of events within a unit; the
  rows are sorted by `group` and then by it. Default `NULL`, which keeps
  the row order.

- type:

  The cut to make:

  `"unit"`

  :   (default) one pathway per unit. With `terminal`, each unit is cut
      after its last terminal state, so the pathway ends on an outcome
      rather than on whatever followed; units without a terminal state
      are left out.

  `"segments"`

  :   consecutive, non-overlapping pathways: each unit is cut after
      every terminal state, so each pathway ends on one. Events after a
      unit's last terminal state belong to no pathway.

  `"anchored"`

  :   one pathway per occurrence of `anchor`, closing at the first
      terminal state after it. Pathways may overlap, and events between
      a terminal state and the next anchor belong to none.

- anchor:

  The single state that opens a pathway when `type = "anchored"`, where
  it is required; the other cuts ignore it. It may also be one of the
  `terminal` states: a pathway then closes at the next terminal state
  after the anchor, not at the anchor itself. Default `NULL`.

- terminal:

  The states that close a pathway, the closing state included. Required
  for `"segments"` and `"anchored"`; optional for `"unit"`, where it
  cuts each unit after its last terminal state. Default `NULL`.

- resolve:

  A named list that appends a resolution label as each pathway's final
  state. Each element is a character vector of states; a pathway takes
  the name of the first element whose states all occur in it, so order
  the list from most to least specific, and pathways matching none are
  labelled `"Unresolved"`. The label becomes the last state of `path`,
  counts in `length`, and is the value of `closes`; the observed last
  state is kept in `ends`. Default `NULL` (no label).

- sep:

  The separator placed between states in `path`. Default `" -> "`.

## Value

A data frame with one row per pathway, in the order the pathways occur
within each unit. Units appear in the order they first occur in `data`,
or sorted by `group` when `order` is given and when `type = "unit"` is
combined with `terminal`. Columns: the `group` columns, `pathway` (the
pathway's number within its unit, from 1), `path` (the states joined
with `sep`), `length` (the number of states), `opens` (the first state)
and `closes` (the last state). With `resolve`, `closes` holds the
resolution label and an `ends` column the observed last state. When no
pathway is found (for example, the anchor or the terminal states never
occur), an error of class `nestimate_no_pathway` is raised.

## See also

[`sequence_compare`](https://pak.dynasite.org/Nestimate/reference/sequence_compare.md),
[`outcome_model`](https://pak.dynasite.org/Nestimate/reference/outcome_model.md)

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
