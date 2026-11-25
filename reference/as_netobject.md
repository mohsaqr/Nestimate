# Convert a network from another package to a netobject

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

  The network to convert: a `psychnet` object, a `cograph_network`, a
  `netobject`, or a `tna` model or `group_tna` (package tna). Any other
  object raises an error.

## Value

A `netobject` (class `c("netobject", "cograph_network")`), read by the
same verbs as one from
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md);
[`validate_netobject()`](https://pak.dynasite.org/Nestimate/reference/validate_netobject.md)
checks its structure. For a `group_tna`: a `netobject_group`, a named
list with one such `netobject` per group.

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

For a `tna` model (package tna), the network is rebuilt with
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
from the sequences the model carries, with the method and the parameters
stored in its attributes (`type`, `params`), and used only when the
rebuilt weights equal the model's own. It then behaves as any network
built from those sequences, in the resampling verbs too. A model that
cannot be reproduced that way (scaled weights, no sequence data, or a
type or parameter
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
does not reproduce) becomes a network of its weights without sequence
data: verbs that only read weights accept it, and verbs that resample or
read sequences raise an error of class `nestimate_tna_unrebuilt` that
names the reason. A pruned model keeps its data, and Nestimate reads the
data: it converts to the full network, with a warning of class
`nestimate_tna_pruned`. A `group_tna` gives a `netobject_group`, one
network per group, named as the groups are. The covariates the model
keeps become the network's metadata: `actor` is read from the metadata
the model keeps (`attr(model$data, "meta")`, kept when the model is
built from
[`tna::prepare_data()`](https://sonsoles.me/tna/reference/prepare_data.html));
a model without it raises `nestimate_actor_missing`. Every verb that
takes a network applies this conversion itself, so calling
`as_netobject()` first is not needed.

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
