# Summary Method for net_mmm

Summary Method for net_mmm

## Usage

``` r
# S3 method for class 'net_mmm'
summary(object, ...)
```

## Arguments

- object:

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
summary(mmm)
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
#> 
#> --- Cluster 1 (72.9%, n=24) ---
#>       A     B     C
#> A 0.001 0.624 0.375
#> B 0.430 0.219 0.351
#> C 0.599 0.001 0.400
#> 
#> --- Cluster 2 (27.1%, n=6) ---
#>       A     B     C
#> A 0.988 0.008 0.004
#> B 0.256 0.445 0.299
#> C 0.003 0.994 0.003
#> 
# }
```
