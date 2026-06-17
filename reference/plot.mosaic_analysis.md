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

The `ggplot` object, invisibly (drawn as a side effect).

## Examples

``` r
df <- data.frame(
  a = sample(c("X", "Y", "Z"), 200, replace = TRUE),
  b = sample(c("P", "Q"), 200, replace = TRUE)
)
res <- mosaic_analysis(df, "a", "b", min_count = 5)
# \donttest{
plot(res, tile_label = "percent", legend_position = "bottom")

# }
```
