# Compare MMM fits across different k

Compare MMM fits across different k

## Usage

``` r
compare_mmm(data, k = 2:5, return_fits = FALSE, ...)

# S3 method for class 'mmm_compare'
print(x, ...)

# S3 method for class 'mmm_compare'
summary(object, ...)

# S3 method for class 'mmm_compare'
plot(x, ...)
```

## Arguments

- data:

  Data frame, netobject, or tna model.

- k:

  Integer vector of component counts. Values must be whole finite
  numbers \>= 2. Default: 2:5.

- return_fits:

  Logical. When `TRUE` the fitted models are retained on the result via
  `attr(result, "fits")` (a list of `net_mmm` objects, named by `k`), so
  the user can pick the chosen model without re-running the EM. Default
  `FALSE` keeps the historical lightweight return shape – only the
  comparison table is allocated.

- ...:

  Arguments passed to
  [`build_mmm`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md).
  In `plot.mmm_compare()`, `print.mmm_compare()` and
  `summary.mmm_compare()`: Unsupported. Supplying unused arguments
  raises an error.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: an
  object of class `mmm_compare`.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: an
  object of class `mmm_compare`.

## Value

A `mmm_compare` data frame, one row per requested `k`, with columns `k`,
`log_likelihood`, `AIC`, `BIC`, `ICL`, `AvePP`, `Entropy` and
`converged`. When `return_fits = TRUE`, the fitted `net_mmm` models are
attached as `attr(result, "fits")`.

In `print.mmm_compare()`: The comparison table, invisibly, with the
printed `best` marker column (`"<-- BIC"` / `"<-- ICL"`) added.

In `summary.mmm_compare()`: A tidy data frame with one row per `k`, plus
a `best` character column flagging the minimum-BIC and minimum-ICL
solutions.

In `plot.mmm_compare()`: A `ggplot` object, invisibly.

## Examples

``` r
seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
                   V2 = sample(c("A","B","C"), 30, TRUE))
comp <- compare_mmm(seqs, k = 2:3, n_starts = 1, max_iter = 10, seed = 1)
comp
#> MMM Model Comparison
#> 
#>  k log_likelihood AIC      BIC      ICL      AvePP     Entropy   converged
#>  2 -62.32972      158.6594 182.4798 184.6518 0.9647836 0.2121121  TRUE    
#>  3 -62.32987      176.6597 213.0909 217.2174 0.9343265 0.2487206 FALSE    
#>  best   
#>  <-- BIC
#>         
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:3], 30, TRUE), V2 = sample(LETTERS[1:3], 30, TRUE),
  V3 = sample(LETTERS[1:3], 30, TRUE), V4 = sample(LETTERS[1:3], 30, TRUE)
)
comp <- compare_mmm(seqs, k = 2:3, seed = 42)
print(comp)
#> MMM Model Comparison
#> 
#>  k log_likelihood AIC      BIC      ICL      AvePP     Entropy   converged
#>  2 -122.4385      278.8771 302.6974 310.2070 0.8974786 0.3160632 TRUE     
#>  3 -118.1764      288.3527 324.7838 332.3711 0.9006835 0.1882551 TRUE     
#>  best   
#>  <-- BIC
#>         

# Retain the fits so the chosen model needs no re-run; summary() marks
# the minimum-BIC and minimum-ICL rows in its `best` column.
comp_with_fits <- compare_mmm(seqs, k = 2:3, seed = 42, return_fits = TRUE)
summary(comp_with_fits)
#>   k log_likelihood      AIC      BIC      ICL     AvePP   Entropy converged
#> 1 2      -122.4385 278.8771 302.6974 310.2070 0.8974786 0.3160632      TRUE
#> 2 3      -118.1764 288.3527 324.7838 332.3711 0.9006835 0.1882551      TRUE
#>   best
#> 1  BIC
#> 2     
# }
```
