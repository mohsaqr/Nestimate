# Set the state colours carried by a network object

`set_state_colors()` attaches a palette to a network object, so every
figure drawn from it uses the same colours:
[`sequence_plot()`](https://pak.dynasite.org/Nestimate/reference/sequence_plot.md),
[`distribution_plot()`](https://pak.dynasite.org/Nestimate/reference/distribution_plot.md),
[`plot_state_frequencies()`](https://pak.dynasite.org/Nestimate/reference/plot_state_frequencies.md)
and
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html).
A colour argument given to one of those plot verbs still wins for that
figure, and the object keeps its palette. `state_colors(x) <- value` is
the replacement form, and
[`state_colors()`](https://pak.dynasite.org/Nestimate/reference/state_colors.md)
reads back the colours in use.

## Usage

``` r
set_state_colors(x, colors)

# Default S3 method
set_state_colors(x, colors)

# S3 method for class 'netobject'
set_state_colors(x, colors)

# S3 method for class 'htna'
set_state_colors(x, colors)

# S3 method for class 'mcml'
set_state_colors(x, colors)

# S3 method for class 'netobject_group'
set_state_colors(x, colors)

state_colors(x) <- value
```

## Arguments

- x:

  A `netobject`, `netobject_group`, `mcml` or `htna`; any other object
  raises an error.

- colors:

  The palette: a named character vector of colours, such as
  `c(plan = "#0072B2", monitor = "#D55E00")`. Names are states, and for
  an `mcml` may also be cluster names. A name the object does not carry
  is reported in a message and draws nothing, so one project-wide
  palette can be attached to every object. States left unnamed are dealt
  the Okabe-Ito colours the palette does not already use, in the
  object's state order (recycled when there are more such states than
  colours). The palette replaces any palette set earlier, and `NULL`
  removes it. An unnamed vector raises an error.

- value:

  The palette, as for `colors`, in the replacement form
  `state_colors(x) <- value`.

## Value

`x` with the palette stored on it; the class is unchanged.
[`state_colors()`](https://pak.dynasite.org/Nestimate/reference/state_colors.md)
shows the colours it resolves to. For a `netobject` or `htna`, the
palette is also written, in node order, into the `meta$splot` producer
contract (`defaults$node_fill`), so
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
draws the same colours without being passed any. A `netobject_group`
stores the palette on the group and on each member network, so a member
drawn on its own carries it too.

## See also

[`state_colors`](https://pak.dynasite.org/Nestimate/reference/state_colors.md)
to read the resolved palette back.

## Examples

``` r
net <- build_network(group_regulation_long, method = "relative",
                     actor = "Actor", action = "Action", time = "Time")
net <- set_state_colors(net, c(plan = "#0072B2", monitor = "#D55E00"))
state_colors(net)
#>        state   color  source
#> 1      adapt #E69F00 default
#> 2   cohesion #56B4E9 default
#> 3  consensus #009E73 default
#> 4 coregulate #F0E442 default
#> 5    discuss #CC79A7 default
#> 6    emotion #999999 default
#> 7    monitor #D55E00     set
#> 8       plan #0072B2     set
#> 9  synthesis #000000 default
```
