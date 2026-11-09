# Retrieve a Registered Estimator

`get_estimator()` returns the registry entry of one estimator: its
function, description and directedness, as
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
uses them. A name that is not registered raises an error that lists the
registered names. Method aliases such as `"tna"` are not resolved here;
give the registered name, as
[`list_estimators()`](https://pak.dynasite.org/Nestimate/reference/list_estimators.md)
shows it.

## Usage

``` r
get_estimator(name)
```

## Arguments

- name:

  The registered name of the estimator, a single string.

## Value

A list with `fn` (the estimator function), `description` (a string) and
`directed` (`TRUE` or `FALSE`).

## Note

The idiographic package exports a function with the same name. With both
packages attached, the one attached last masks the other; write
`Nestimate::get_estimator()` to call this one.

## See also

[`register_estimator`](https://pak.dynasite.org/Nestimate/reference/register_estimator.md),
[`list_estimators`](https://pak.dynasite.org/Nestimate/reference/list_estimators.md)

## Examples

``` r
est <- get_estimator("relative")
```
