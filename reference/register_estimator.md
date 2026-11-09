# Register a Network Estimator

`register_estimator()` adds an estimator function to the registry that
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
reads, so that `build_network(data, method = name)` calls it the way it
calls a built-in method; the package registers its fourteen built-in
methods with this function when it loads. A name already in the registry
is replaced, built-ins included, without a message.
[`list_estimators()`](https://pak.dynasite.org/Nestimate/reference/list_estimators.md)
shows what is registered and
[`remove_estimator()`](https://pak.dynasite.org/Nestimate/reference/remove_estimator.md)
takes an entry out.

A name that
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
maps to a different built-in method (an alias such as `"tna"` or
`"ebicglasso"`) is resolved to that method before the registry is read,
so an estimator registered under it is never called; pick a name that is
not an alias.

## Usage

``` r
register_estimator(name, fn, description, directed)
```

## Arguments

- name:

  The name the estimator is called by, as `method` in
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md):
  a single non-empty string. A name already registered is replaced.

- fn:

  The estimator function. Its first argument is `data`, and it must
  accept `...`, because
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  passes every setting by name. It returns a list with `matrix` (the
  square numeric weight matrix), `nodes` (the node names, in the
  matrix's order) and `directed` (`TRUE` or `FALSE`); see *Details* for
  the optional elements.

- description:

  A short description of the estimator, a single string, shown by
  [`list_estimators()`](https://pak.dynasite.org/Nestimate/reference/list_estimators.md).

- directed:

  Whether the estimator produces directed networks, `TRUE` or `FALSE`.
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  reads it to decide whether `level` (the multilevel decomposition) is
  allowed: a directed estimator refuses it. The edges of the network
  follow the `directed` element `fn` returns, but
  [`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  reads this registered flag to build its edge tables, so the two should
  agree.

## Value

`NULL`, invisibly. The estimator is added to the registry as a side
effect.

## Details

[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
calls the function as `fn(data = data, ...)`, where `...` holds the
`params` list, the further arguments given to
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
and a `format` setting, which is why `fn` must accept `...`. When
`actor` names a column of wide or association data (long data is
prepared by
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
itself) and `fn` has an `id_col` argument (or, failing that and for data
that is not one-hot, an `id` argument), the actor column is passed under
that name unless `params` already sets it. Scaling, thresholding and the
edge list are then applied to the returned `matrix`, as for a built-in
method.

Besides `matrix`, `nodes` and `directed`, the returned list may hold
`cleaned_data`, the data frame the network was estimated from, which
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
splits into the network's `data` and `metadata` fields; any other named
element is copied onto the `netobject`, replacing a field of the same
name. A return value without `matrix`, `nodes` or `directed` makes
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
stop with an error.

## See also

[`get_estimator`](https://pak.dynasite.org/Nestimate/reference/get_estimator.md),
[`list_estimators`](https://pak.dynasite.org/Nestimate/reference/list_estimators.md),
[`remove_estimator`](https://pak.dynasite.org/Nestimate/reference/remove_estimator.md),
[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md)

## Examples

``` r
my_fn <- function(data, ...) {
  m <- cor(data)
  diag(m) <- 0
  list(matrix = m, nodes = colnames(m), directed = FALSE)
}
register_estimator("my_cor", my_fn, "Custom correlation", directed = FALSE)
df <- data.frame(A = rnorm(20), B = rnorm(20), C = rnorm(20))
net <- build_network(df, method = "my_cor")
remove_estimator("my_cor")
```
