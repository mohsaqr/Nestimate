# Fit a Mixed Markov Model

Discovers latent subgroups with different transition dynamics using
Expectation-Maximization. Each mixture component has its own transition
matrix. Sequences are probabilistically assigned to components.

## Usage

``` r
build_mmm(
  data,
  k = 2L,
  n_starts = 50L,
  max_iter = 200L,
  tol = 1e-06,
  smooth = 0.01,
  seed = NULL,
  covariates = NULL
)
```

## Arguments

- data:

  A data.frame (wide format), `netobject`, or `tna` model. For tna
  objects, extracts the stored data.

- k:

  Integer. Number of mixture components. Default: 2.

- n_starts:

  Integer. Number of random restarts. Default: 50.

- max_iter:

  Integer. Maximum EM iterations per start. Default: 200.

- tol:

  Numeric. Convergence tolerance. Default: 1e-6.

- smooth:

  Numeric. Laplace smoothing constant. Default: 0.01.

- seed:

  Integer or NULL. Random seed.

- covariates:

  Optional. Covariates integrated into the EM algorithm to model
  covariate-dependent mixing proportions. Accepts formula, character
  vector, string, or data.frame (same forms as
  [`build_clusters`](https://mohsaqr.github.io/Nestimate/reference/build_clusters.md)).
  Unlike the post-hoc analysis in
  [`build_clusters()`](https://mohsaqr.github.io/Nestimate/reference/build_clusters.md),
  these covariates directly influence cluster membership during
  estimation. Requires the nnet package.

## Value

An object of class `net_mmm` with components:

- models:

  List of `netobject`s, one per component.

- k:

  Number of components.

- mixing:

  Numeric vector of mixing proportions.

- posterior:

  N x k matrix of posterior probabilities.

- assignments:

  Integer vector of hard assignments (1..k).

- quality:

  List: `avepp` (per-class), `avepp_overall`, `entropy`,
  `relative_entropy`, `classification_error`.

- log_likelihood, BIC, AIC, ICL:

  Model fit statistics.

- states:

  Character vector of state names.

## See also

[`compare_mmm`](https://mohsaqr.github.io/Nestimate/reference/compare_mmm.md),
[`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md)

## Examples

``` r
seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
                   V2 = sample(c("A","B","C"), 30, TRUE))
mmm <- build_mmm(seqs, k = 2, n_starts = 1, max_iter = 10, seed = 1)
mmm
#> Mixed Markov Model
#>   Sequences: 30  |  Clusters: 2  |  States: 3
#>   ICs: LL = -57.596  |  BIC = 173.013  |  AIC = 149.192  |  ICL = 174.931
#>   Quality: AvePP = 0.969  |  Entropy = 0.192  |  Class.Err = 0.0%
#>   Status: did not converge in 10 iterations
#> 
#>   Cluster  N           Mix%   AvePP
#>   1        24 (80.0%)  78.5%  0.971
#>   2        6 (20.0%)   21.5%  0.961
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:3], 30, TRUE), V2 = sample(LETTERS[1:3], 30, TRUE),
  V3 = sample(LETTERS[1:3], 30, TRUE), V4 = sample(LETTERS[1:3], 30, TRUE)
)
mmm <- build_mmm(seqs, k = 2, seed = 42)
print(mmm)
#> Mixed Markov Model
#>   Sequences: 30  |  Clusters: 2  |  States: 3
#>   ICs: LL = -124.400  |  BIC = 306.621  |  AIC = 282.800  |  ICL = 313.645
#>   Quality: AvePP = 0.904  |  Entropy = 0.304  |  Class.Err = 0.0%
#> 
#>   Cluster  N           Mix%   AvePP
#>   1        20 (66.7%)  62.4%  0.896
#>   2        10 (33.3%)  37.6%  0.919
summary(mmm)
#> Mixed Markov Model
#>   Sequences: 30  |  Clusters: 2  |  States: 3
#>   ICs: LL = -124.400  |  BIC = 306.621  |  AIC = 282.800  |  ICL = 313.645
#>   Quality: AvePP = 0.904  |  Entropy = 0.304  |  Class.Err = 0.0%
#> 
#>   Cluster  N           Mix%   AvePP
#>   1        20 (66.7%)  62.4%  0.896
#>   2        10 (33.3%)  37.6%  0.919
#> 
#> --- Cluster 1 (62.4%, n=20) ---
#>       A     B     C
#> A 0.001 0.568 0.431
#> B 0.485 0.203 0.313
#> C 0.189 0.629 0.181
#> 
#> --- Cluster 2 (37.6%, n=10) ---
#>       A     B     C
#> A 0.593 0.261 0.145
#> B 0.149 0.397 0.454
#> C 0.394 0.002 0.604
#> 
#>   component     prior n_assigned mean_posterior     avepp
#> 1         1 0.6241614         20      0.8959195 0.8959195
#> 2         2 0.3758386         10      0.9193350 0.9193350
# }
```
