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
#>  2 -116.3993      266.7986 290.6190 293.2251 0.9654796 0.09341584 TRUE     
#>  3 -110.9706      273.9411 310.3723 312.0948 0.9751276 0.05463341 TRUE     
#>  best   
#>  <-- BIC
#>         
# }
```
