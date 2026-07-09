# Coerce a network object to a Nestimate netobject

Promotes a psychnet result (class `c("psychnet", "cograph_network")`) or
any bare `cograph_network` to the dual-class
`c("netobject", "cograph_network")` used throughout Nestimate, so it
dispatches to the package's verbs (`centrality()`,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html), bootstrap,
reliability, ...). A `netobject` is returned unchanged.

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

  A `psychnet` object, a `cograph_network`, or a `netobject`.

## Value

A `c("netobject", "cograph_network")` object.

## Details

The psychnet method re-derives the integer-indexed edge table that
Nestimate expects (psychnet stores character-labelled edges), preserves
the estimator name in `$method`, and parks every psychnet-specific
field - including the graphical-lasso `$kkt` optimality certificate -
under `$meta$psychnet` so nothing is lost in translation.

## See also

[`validate_netobject`](https://saqr.me/Nestimate/reference/validate_netobject.md)

## Examples

``` r
net <- build_cor(data.frame(a = rnorm(50), b = rnorm(50), c = rnorm(50)))
identical(as_netobject(net), net) # netobjects pass through unchanged
#> [1] TRUE
```
