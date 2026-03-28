# Print Method for net_mmm

Print Method for net_mmm

## Usage

``` r
# S3 method for class 'net_mmm'
print(x, ...)
```

## Arguments

- x:

  A `net_mmm` object.

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
mmm <- build_mmm(seqs, k = 2, n_starts = 5, seed = 1)
print(mmm)
#> Mixed Markov Model
#>   k = 2 | 30 sequences | 3 states
#>   LL = -89.3 | BIC = 236.4 | ICL = 242.0
#> 
#>   Cluster  Size  Mix%%   AvePP
#>   ------------------------------
#>         1    24  72.9%  0.911
#>         2     6  27.1%  0.999
#> 
#>   Overall AvePP = 0.928 | Entropy = 0.173 | Class.Err = 0.0%
# }
```
