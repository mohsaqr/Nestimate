# Remove a Registered Estimator

Remove a network estimator from the registry.

## Usage

``` r
remove_estimator(name)
```

## Arguments

- name:

  Character. Name of the estimator to remove.

## Value

Invisible `NULL`.

## See also

[`register_estimator`](https://mohsaqr.github.io/Nestimate/reference/register_estimator.md),
[`list_estimators`](https://mohsaqr.github.io/Nestimate/reference/list_estimators.md)

## Examples

``` r
# \donttest{
register_estimator("test_est", function(data, ...) diag(3),
  description = "test", directed = FALSE)
remove_estimator("test_est")
# }
```
