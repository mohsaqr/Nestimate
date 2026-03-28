# Print Method for mmm_compare

Print Method for mmm_compare

## Usage

``` r
# S3 method for class 'mmm_compare'
print(x, ...)
```

## Arguments

- x:

  An `mmm_compare` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(c("A","B","C"), 30, TRUE),
  V2 = sample(c("A","B","C"), 30, TRUE),
  V3 = sample(c("A","B","C"), 30, TRUE)
)
cmp <- compare_mmm(seqs, k = 2:3, n_starts = 5, seed = 1)
print(cmp)
#> MMM Model Comparison
#> 
#>  k log_likelihood AIC      BIC      ICL      AvePP     Entropy   converged
#>  2 -89.29212      212.5842 236.4046 241.9893 0.9284601 0.1726947 TRUE     
#>  3 -87.61855      227.2371 263.6682 270.9108 0.9094695 0.1893851 TRUE     
#>  best   
#>  <-- BIC
#>         
# }
```
