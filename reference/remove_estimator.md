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

[`register_estimator`](https://saqr.me/Nestimate/reference/register_estimator.md),
[`list_estimators`](https://saqr.me/Nestimate/reference/list_estimators.md)

## Examples

``` r
register_estimator("test_est", function(data, ...) diag(3),
  description = "test", directed = FALSE)
remove_estimator("test_est")
```
