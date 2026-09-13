# Plot method for mosaic_analysis objects

Re-renders the flat mosaic from the stored contingency table and
residuals, so styling can be changed without re-running the test. Any
flat-mosaic styling argument (`tile_label`, `pct_base`,
`col_label_side`, `legend_size`, ...) may be overridden via `...`.

## Usage

``` r
# S3 method for class 'mosaic_analysis'
plot(x, ...)
```

## Arguments

- x:

  A `mosaic_analysis` object.

- ...:

  Styling overrides forwarded to the flat renderer.

## Value

The re-rendered flat mosaic `ggplot` object, invisibly; the plot is
drawn on the active device as a side effect.

## Examples

``` r
# \donttest{
data(group_regulation_long, package = "Nestimate")
res <- mosaic_analysis(group_regulation_long, "Course", "Action",
                       min_count = 20)
plot(res, tile_label = "percent", legend_position = "bottom")

# }
```
