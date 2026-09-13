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

  Unsupported. Supplying unused arguments raises an error.

## Value

A per-component summary `data.frame`. The class and visibility depend on
whether the model was fitted with covariates:

- No covariates:

  A plain `data.frame` with one row per component and columns
  `component`, `prior`, `n_assigned`, `mean_posterior`, `avepp`,
  returned *visibly* (so it auto-prints after the printed summary
  block).

- With covariates:

  A `tidy_covariates`/`data.frame` (the tidied covariate table, with the
  per-component stats attached), returned *invisibly*.

In both cases the printed summary (model fit, per-cluster transition
matrices, optional covariate profiles) is emitted as a side effect.

## Examples

``` r
seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
                   V2 = sample(c("A","B","C"), 30, TRUE))
mmm <- build_mmm(seqs, k = 2, n_starts = 1, max_iter = 10, seed = 1)
summary(mmm)
#> Mixed Markov Model
#>   Sequences: 30  |  Clusters: 2  |  States: 3
#>   ICs: LL = -62.330  |  BIC = 182.480  |  AIC = 158.659  |  ICL = 184.652
#>   Quality: AvePP = 0.965  |  Entropy = 0.212  |  Class.Err = 0.0%
#> 
#>   Cluster  N           Mix%   AvePP
#>   1        24 (80.0%)  78.1%  0.966
#>   2        6 (20.0%)   21.9%  0.961
#> 
#> --- Cluster 1 (78.1%, n=24) ---
#>       A     B     C
#> A 0.248 0.376 0.376
#> B 0.447 0.333 0.220
#> C 0.842 0.034 0.124
#> 
#> --- Cluster 2 (21.9%, n=6) ---
#>       A     B     C
#> A 0.330 0.335 0.335
#> B 0.337 0.334 0.329
#> C 0.035 0.939 0.027
#> 
#>   component     prior n_assigned mean_posterior     avepp
#> 1         1 0.7809305         24      0.9656098 0.9656098
#> 2         2 0.2190695          6      0.9614788 0.9614788
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
#>   Sequences: 30  |  Clusters: 2  |  States: 3
#>   ICs: LL = -89.292  |  BIC = 236.405  |  AIC = 212.584  |  ICL = 241.989
#>   Quality: AvePP = 0.928  |  Entropy = 0.173  |  Class.Err = 0.0%
#> 
#>   Cluster  N           Mix%   AvePP
#>   1        24 (80.0%)  72.9%  0.911
#>   2        6 (20.0%)   27.1%  0.999
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
#>   component    prior n_assigned mean_posterior     avepp
#> 1         1 0.728925         24      0.9108722 0.9108722
#> 2         2 0.271075          6      0.9988113 0.9988113
# }
```
