# Set the state colours carried by a network object

Attaches a palette to the object so every figure drawn from it uses the
same colours:
[`sequence_plot`](https://saqr.me/Nestimate/reference/sequence_plot.md),
[`distribution_plot`](https://saqr.me/Nestimate/reference/distribution_plot.md),
[`plot_state_frequencies`](https://saqr.me/Nestimate/reference/plot_state_frequencies.md)
and
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html).

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

  A `netobject`, `netobject_group`, `mcml` or `htna`.

- colors:

  A named character vector of colours, e.g.
  `c(plan = "#0072B2", monitor = "#D55E00")`. Names are states, and for
  an `mcml` may also be cluster names. Names the object does not carry
  are dropped with a message, so one project-wide palette can be
  attached to every object. States you do not name keep the default
  Okabe-Ito colour. `NULL` removes a palette set earlier.

- value:

  The palette, as for `colors`, in the replacement form
  `state_colors(x) <- value`.

## Value

`x`, with the palette stored in `x$state_colors` and, for an object
carrying `$nodes`, mirrored into `x$meta$splot$defaults$node_fill` in
node order so
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
honours it. The class is unchanged.

## See also

[`state_colors`](https://saqr.me/Nestimate/reference/state_colors.md) to
read the resolved palette back.

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
