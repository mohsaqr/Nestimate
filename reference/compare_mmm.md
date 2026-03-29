# Compare MMM fits across different k

Compare MMM fits across different k

## Usage

``` r
compare_mmm(data, k = 2:5, ...)
```

## Arguments

- data:

  Data frame, netobject, or tna model.

- k:

  Integer vector of component counts. Default: 2:5.

- ...:

  Arguments passed to
  [`build_mmm`](https://mohsaqr.github.io/Nestimate/reference/build_mmm.md).

## Value

A `mmm_compare` data frame with BIC, AIC, ICL, AvePP, entropy per k.

## Examples

``` r
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:3], 30, TRUE), V2 = sample(LETTERS[1:3], 30, TRUE),
  V3 = sample(LETTERS[1:3], 30, TRUE), V4 = sample(LETTERS[1:3], 30, TRUE)
)
comp <- compare_mmm(seqs, k = 2:3, seed = 42)
print(comp)
#> MMM Model Comparison
#> 
#>  k log_likelihood AIC      BIC      ICL      AvePP     Entropy    converged
#>  2 -124.1824      282.3648 306.1851 309.0453 0.9598347 0.13863520 TRUE     
#>  3 -120.6159      293.2317 329.6628 331.8871 0.9697914 0.06640355 TRUE     
#>  best   
#>  <-- BIC
#>         
# }
```
