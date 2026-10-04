# The state colours an object will draw with

Reads back the palette an object resolves to: the colours set with
[`set_state_colors`](https://pak.dynasite.org/Nestimate/reference/set_state_colors.md)
plus the defaults filled in for everything else, so the table is what
the figures actually use.

## Usage

``` r
state_colors(x, ...)

# Default S3 method
state_colors(x, ...)

# S3 method for class 'netobject'
state_colors(x, ...)

# S3 method for class 'htna'
state_colors(x, ...)

# S3 method for class 'mcml'
state_colors(x, ...)

# S3 method for class 'netobject_group'
state_colors(x, ...)
```

## Arguments

- x:

  A `netobject`, `netobject_group`, `mcml` or `htna`.

- ...:

  Ignored, for method consistency.

## Value

A `data.frame`, one row per colour key the object carries, with columns
`state` (the key), `color` (the hex colour it draws with) and `source`
(`"set"` when the palette named it, `"default"` when it fell back to
Okabe-Ito). For an `mcml` the cluster names appear after the states.

## See also

[`set_state_colors`](https://pak.dynasite.org/Nestimate/reference/set_state_colors.md).

## Examples

``` r
net <- build_network(group_regulation_long, method = "relative",
                     actor = "Actor", action = "Action", time = "Time")
state_colors(set_state_colors(net, c(plan = "#0072B2")))
#>        state   color  source
#> 1      adapt #E69F00 default
#> 2   cohesion #56B4E9 default
#> 3  consensus #009E73 default
#> 4 coregulate #F0E442 default
#> 5    discuss #D55E00 default
#> 6    emotion #CC79A7 default
#> 7    monitor #999999 default
#> 8       plan #0072B2     set
#> 9  synthesis #000000 default
```
