# Per-Class State Distribution as a Tidy Data Frame

`state_distribution()` counts how often each state occurs in the
sequences behind a network and returns one row per group and state, with
the count and the share within the group. The groups follow the class: a
single `netobject` is one group, `"all"`, unless it carries node groups,
which then form the groups (states in none are `"(ungrouped)"`); an
`htna` network has one group per actor; a `netobject_group` one per
network, named as in the group (or `"Group 1"`, `"Group 2"`, ...); an
`mcml` one per cluster, with the node labels it was given in place of
the raw codes. Missing and empty cells are not counted.

The table is the one behind
[`plot_state_frequencies()`](https://pak.dynasite.org/Nestimate/reference/plot_state_frequencies.md),
and the per-group counterpart of
[`state_frequencies()`](https://pak.dynasite.org/Nestimate/reference/state_frequencies.md),
which counts the states of raw sequence data.

## Usage

``` r
state_distribution(x, ...)

# S3 method for class 'netobject'
state_distribution(x, ...)

# S3 method for class 'htna'
state_distribution(x, ...)

# S3 method for class 'mcml'
state_distribution(x, include_macro = FALSE, ...)

# S3 method for class 'netobject_group'
state_distribution(x, ...)

# S3 method for class 'state_freq'
state_distribution(x, ...)

# Default S3 method
state_distribution(x, ...)
```

## Arguments

- x:

  A `netobject`, `netobject_group`, `mcml` or `htna` network built from
  sequences, or the `state_freq` result of
  [`plot_state_frequencies()`](https://pak.dynasite.org/Nestimate/reference/plot_state_frequencies.md),
  whose table is returned as drawn. A network without sequence data
  raises an error; in a `netobject_group` or `mcml`, groups without it
  are left out, and an error is raised only when none has any. Other
  classes raise an error.

- ...:

  Not used. The `mcml` method raises an error for any argument it does
  not take; the other methods ignore it.

- include_macro:

  For an `mcml`: whether to add a `"macro"` group in front, holding the
  state counts over all clusters. Default `FALSE`.

## Value

A data frame with one row per group and state and the columns `group`
(character), `state` (character), `count` (integer) and `proportion`
(the state's share of its group's events). Within a group the rows run
from the most to the least frequent state; the `"macro"` group is in
alphabetical order.

## Examples

``` r
data(group_regulation_long, package = "Nestimate")
net <- build_network(group_regulation_long, method = "frequency",
                     format = "long", actor = "Actor", action = "Action",
                     order = "Time", group = "Course")
state_distribution(net)
#>    group      state count proportion
#> 1      A  consensus  3298 0.26618241
#> 2      A       plan  2805 0.22639225
#> 3      A    discuss  1960 0.15819209
#> 4      A    emotion  1517 0.12243745
#> 5      A   cohesion   923 0.07449556
#> 6      A coregulate   855 0.06900726
#> 7      A    monitor   602 0.04858757
#> 8      A  synthesis   290 0.02340597
#> 9      A      adapt   140 0.01129944
#> 10     B       plan  2445 0.25399958
#> 11     B  consensus  2226 0.23124870
#> 12     B    discuss  1453 0.15094536
#> 13     B    emotion   995 0.10336588
#> 14     B coregulate   817 0.08487430
#> 15     B   cohesion   590 0.06129233
#> 16     B    monitor   565 0.05869520
#> 17     B  synthesis   274 0.02846458
#> 18     B      adapt   261 0.02711407
#> 19     C       plan  1373 0.24886714
#> 20     C  consensus  1273 0.23074134
#> 21     C    discuss   854 0.15479427
#> 22     C    emotion   563 0.10204821
#> 23     C coregulate   461 0.08355991
#> 24     C    monitor   349 0.06325902
#> 25     C   cohesion   326 0.05909009
#> 26     C  synthesis   165 0.02990756
#> 27     C      adapt   153 0.02773246
```
