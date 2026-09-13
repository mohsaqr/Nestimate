# Summary method for mosaic_analysis objects

Summary method for mosaic_analysis objects

## Usage

``` r
# S3 method for class 'mosaic_analysis'
summary(object, ...)
```

## Arguments

- object:

  A `mosaic_analysis` object.

- ...:

  Ignored.

## Value

The tidy per-cell `data.frame`: one row per (var1, var2) cell, with the
two variable columns (named after `var1` / `var2`) plus `observed`,
`expected`, `residual` and `pct`. The one-row test summary is attached
as the `"stats"` attribute.
