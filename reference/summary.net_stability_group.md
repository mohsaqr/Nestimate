# Summary Method for net_stability_group

Per-network stability as a tidy data frame. Stacks
[`summary()`](https://rdrr.io/r/base/summary.html) results for each
network with a `group` column.

## Usage

``` r
# S3 method for class 'net_stability_group'
summary(object, ...)
```

## Arguments

- object:

  A `net_stability_group` object.

- ...:

  Additional arguments (ignored).

## Value

A data frame with columns `group`, `measure`, `drop_prop`, `mean_cor`,
`sd_cor`, `prop_above`.
