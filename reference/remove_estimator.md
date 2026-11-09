# Remove a Registered Estimator

`remove_estimator()` takes an estimator out of the registry, so
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
no longer finds it by name. Built-ins can be removed too; they come back
when the package is loaded again. A name that is not registered raises
an error.

## Usage

``` r
remove_estimator(name)
```

## Arguments

- name:

  The registered name of the estimator to remove, a single string.

## Value

`NULL`, invisibly. The entry is removed from the registry as a side
effect.

## See also

[`register_estimator`](https://pak.dynasite.org/Nestimate/reference/register_estimator.md),
[`list_estimators`](https://pak.dynasite.org/Nestimate/reference/list_estimators.md)

## Examples

``` r
register_estimator("test_est", function(data, ...) diag(3),
  description = "test", directed = FALSE)
remove_estimator("test_est")
```
