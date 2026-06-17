# Retrieve a Registered Estimator

Retrieve a registered network estimator by name.

## Usage

``` r
get_estimator(name)
```

## Arguments

- name:

  Character. Name of the estimator to retrieve.

## Value

A list with elements `fn`, `description`, `directed`.

## See also

[`register_estimator`](https://saqr.me/Nestimate/reference/register_estimator.md),
[`list_estimators`](https://saqr.me/Nestimate/reference/list_estimators.md)

## Examples

``` r
est <- get_estimator("relative")
```
