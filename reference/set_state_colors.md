# Set the state colours of a network object

`set_state_colors()` attaches a palette to a network object, so every
figure drawn from it uses the same colours. `state_colors(x) <- value`
is the replacement form, and
[`state_colors()`](https://pak.dynasite.org/Nestimate/reference/state_colors.md)
reads back the colours in use.

One rule holds for colours across the package:

- With no colours set, the states are drawn in the Okabe-Ito colours,
  dealt in the alphabetical order of the states, so a state has the same
  colour in every figure whatever the order of the network's nodes.

- A palette set here replaces those colours for the states it names.

- Objects computed from the network carry its colours: groups and their
  members,
  [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md),
  [`build_mmm()`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
  and the networks
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  and
  [`cluster_network()`](https://pak.dynasite.org/Nestimate/reference/cluster_network.md)
  build from them,
  [`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md),
  [`certainty()`](https://pak.dynasite.org/Nestimate/reference/certainty.md),
  [`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md),
  [`bayes_compare()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md),
  [`net_edge_betweenness()`](https://pak.dynasite.org/Nestimate/reference/net_edge_betweenness.md),
  [`entropy_network()`](https://pak.dynasite.org/Nestimate/reference/entropy_network.md),
  [`net_prune()`](https://pak.dynasite.org/Nestimate/reference/net_prune.md),
  [`stability_landscape()`](https://pak.dynasite.org/Nestimate/reference/stability_landscape.md),
  [`landscape_trend()`](https://pak.dynasite.org/Nestimate/reference/landscape_trend.md)
  and
  [`landscape_surface()`](https://pak.dynasite.org/Nestimate/reference/landscape_surface.md).
  A cluster or group that lacks a state keeps the colours of the others.

- Every figure that colours states follows them:
  [`sequence_plot()`](https://pak.dynasite.org/Nestimate/reference/sequence_plot.md),
  [`distribution_plot()`](https://pak.dynasite.org/Nestimate/reference/distribution_plot.md),
  [`plot_state_frequencies()`](https://pak.dynasite.org/Nestimate/reference/plot_state_frequencies.md),
  the [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods
  of
  [`stability_landscape()`](https://pak.dynasite.org/Nestimate/reference/stability_landscape.md)
  and
  [`landscape_trend()`](https://pak.dynasite.org/Nestimate/reference/landscape_trend.md),
  and
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
  (through the network's `meta$splot` contract).

- `state_colors =` on one of those figures changes that figure only.

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
  `c(plan = "#0072B2", monitor = "#D55E00")`. Names are states (for a
  network built with `labels`, the labels it draws), and for an `mcml`
  may also be cluster names. A name the object does not carry is
  reported in a message and draws nothing, so one project-wide palette
  can be attached to every object. States left unnamed are dealt the
  Okabe-Ito colours the palette does not already use, in the object's
  state order (recycled when there are more such states than colours).
  The palette replaces any palette set earlier, and `NULL` removes it.
  An unnamed vector raises an error.

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
