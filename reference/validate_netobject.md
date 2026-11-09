# Validate a netobject / cograph_network against the shared schema

`validate_netobject()` checks that an object meets the structure
Nestimate `netobject`s and psychnets objects share, so either can be
passed across the package boundary. It is the one place that says what a
network object is, so a change on either side (a renamed field, a
mistyped edge column) fails here instead of further down, for instance
when the network is drawn.

The checks are the shared part only: the object inherits
`cograph_network`; `weights` is a square numeric matrix with node names
as dimnames; `nodes` is a data frame with columns `id`, `label` and
`name`, one row per node; `edges` is a data frame with columns `from`,
`to` and a numeric `weight`, whose endpoints are integer node indices
(Nestimate) or node labels (psychnets); `directed` is a single logical
value; and `method` is a single string. The `x`/`y` layout columns of
`nodes` and Nestimate's own fields (`data`, `level`, ...) are not
required.

## Usage

``` r
validate_netobject(x)
```

## Arguments

- x:

  The object to check, expected to be a `cograph_network`.

## Value

`TRUE`, invisibly, when `x` passes every check. Otherwise an error that
lists every check it failed.

## See also

[`as_netobject`](https://pak.dynasite.org/Nestimate/reference/as_netobject.md)

## Examples

``` r
net <- build_cor(data.frame(a = rnorm(50), b = rnorm(50), c = rnorm(50)))
validate_netobject(net)
```
