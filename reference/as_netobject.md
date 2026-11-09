# Coerce a network object to a Nestimate netobject

`as_netobject()` converts a psychnets result (class
`c("psychnet", "cograph_network")`) or any bare `cograph_network` into
the `c("netobject", "cograph_network")` object the rest of Nestimate
reads, so it can be passed to
[`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
and the other verbs that read a network. The resampling verbs
([`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md),
[`network_reliability()`](https://pak.dynasite.org/Nestimate/reference/network_reliability.md))
need the source data, which a psychnets result does not carry, so build
the network with
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
for those. A `netobject` is returned unchanged.

## Usage

``` r
as_netobject(x)

# S3 method for class 'netobject'
as_netobject(x)

# S3 method for class 'psychnet'
as_netobject(x)

# S3 method for class 'cograph_network'
as_netobject(x)

# Default S3 method
as_netobject(x)
```

## Arguments

- x:

  The network to convert: a `psychnet` object, a `cograph_network` or a
  `netobject`. Any other object raises an error.

## Value

A `netobject` (class `c("netobject", "cograph_network")`), read by the
same verbs as one from
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md);
[`validate_netobject()`](https://pak.dynasite.org/Nestimate/reference/validate_netobject.md)
checks its structure.

## Details

For a `psychnet`, the edge table is rebuilt with integer node indices (a
`psychnet` stores character-labelled edges), the estimator name is kept
in `method`, and every field specific to the `psychnet`, including the
graphical-lasso `kkt` optimality certificate, is kept in the `psychnet`
entry of `meta`, so nothing is lost in the conversion.

For a bare `cograph_network`, the method is read from the `tna` entry of
its `meta`; without one, a symmetric weight matrix is taken as
`"co_occurrence"` and any other as `"relative"`. A missing `directed`
field is read as `TRUE`. Numeric state codes in the stored data of a
converted tna model are decoded to their labels.

## Note

The idiographic package exports a function with the same name. With both
packages attached, the one attached last masks the other; write
`Nestimate::as_netobject()` to call this one.

## See also

[`validate_netobject`](https://pak.dynasite.org/Nestimate/reference/validate_netobject.md)

## Examples

``` r
net <- build_cor(data.frame(a = rnorm(50), b = rnorm(50), c = rnorm(50)))
identical(as_netobject(net), net) # netobjects pass through unchanged
#> [1] TRUE
```
