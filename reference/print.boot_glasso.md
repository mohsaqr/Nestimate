# Print Method for boot_glasso

Print Method for boot_glasso

## Usage

``` r
# S3 method for class 'boot_glasso'
print(x, ...)
```

## Arguments

- x:

  A `boot_glasso` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
# \donttest{
set.seed(42)
mat <- matrix(rnorm(60), ncol = 4)
colnames(mat) <- LETTERS[1:4]
boot <- boot_glasso(as.data.frame(mat), iter = 20, cs_iter = 10,
  centrality = "strength", seed = 42)
print(boot)
#> GLASSO Bootstrap (20 iterations, 10 case-drop per proportion)
#>   Data: 15 x 4  |  Alpha: 0.05  |  Gamma: 0.50
#>   Edges: 0/6 significant (CI excludes zero)
#>   Mean inclusion probability: 0.34
#> 
#>   Centrality Stability (CS-coefficient):
#>     strength:              0.00 [Unstable]
#> 
#>   Edge differences: 1/15 pairs significantly different
#>   Timing: 0.1s (bootstrap: 0.1s, case-drop: 0.0s)
# }
```
