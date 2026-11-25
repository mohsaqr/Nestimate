# Estimate a network (deprecated)

`estimate_network()` is deprecated: use
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
instead. It passes its arguments to
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
unchanged and returns what that returns, after a warning of class
`deprecatedWarning`.

## Usage

``` r
estimate_network(
  data,
  method = "relative",
  params = list(),
  scaling = NULL,
  threshold = 0,
  level = NULL,
  ...
)
```

## Arguments

- data:

  The data, as for
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md).

- method:

  The estimator, as in
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md).
  Default `"relative"`, kept for backward compatibility;
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  has no default.

- params:

  Settings for the estimator, as a named list, as on the method's page
  (for example
  [`build_glasso()`](https://pak.dynasite.org/Nestimate/reference/build_glasso.md)).
  Default [`list()`](https://rdrr.io/r/base/list.html).

- scaling, threshold:

  As in
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md).

- level:

  Multilevel decomposition for the association methods, as in
  [`build_cor()`](https://pak.dynasite.org/Nestimate/reference/build_cor.md).
  Default `NULL`.

- ...:

  Further arguments for
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
  such as `actor`, `action`, `time` or `group`.

## Value

A `netobject`, exactly as
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
returns it (a `netobject_group` or `netobject_ml` when `...` or `level`
asks for one).

## See also

[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md)

## Examples

``` r
data <- data.frame(T1 = c("x","y","z","x"), T2 = c("y","x","z","y"))
net <- estimate_network(data, method = "relative")
#> Warning: 'estimate_network' is deprecated.
#> Use 'build_network' instead.
#> See help("Deprecated")
```
