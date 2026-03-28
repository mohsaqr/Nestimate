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
  [`cluster_data`](https://mohsaqr.github.io/Nestimate/reference/cluster_data.md)).
  Unlike the post-hoc analysis in
  [`cluster_data()`](https://mohsaqr.github.io/Nestimate/reference/cluster_data.md),
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
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:3], 30, TRUE), V2 = sample(LETTERS[1:3], 30, TRUE),
  V3 = sample(LETTERS[1:3], 30, TRUE), V4 = sample(LETTERS[1:3], 30, TRUE)
)
mmm <- build_mmm(seqs, k = 2, seed = 42)
print(mmm)
#> Mixed Markov Model
#>   k = 2 | 30 sequences | 3 states
#>   LL = -122.4 | BIC = 302.7 | ICL = 302.8
#> 
#>   Cluster  Size  Mix%%   AvePP
#>   ------------------------------
#>         1    24  80.0%  0.998
#>         2     6  20.0%  0.995
#> 
#>   Overall AvePP = 0.998 | Entropy = 0.019 | Class.Err = 0.0%
summary(mmm)
#> Mixed Markov Model
#>   k = 2 | 30 sequences | 3 states
#>   LL = -122.4 | BIC = 302.7 | ICL = 302.8
#> 
#>   Cluster  Size  Mix%%   AvePP
#>   ------------------------------
#>         1    24  80.0%  0.998
#>         2     6  20.0%  0.995
#> 
#>   Overall AvePP = 0.998 | Entropy = 0.019 | Class.Err = 0.0%
#> 
#> --- Cluster 1 (80.0%, n=24) ---
#>       A     B     C
#> A 0.250 0.250 0.500
#> B 0.302 0.333 0.364
#> C 0.334 0.466 0.200
#> 
#> --- Cluster 2 (20.0%, n=6) ---
#>       A     B     C
#> A 0.002 0.798 0.200
#> B 0.746 0.251 0.003
#> C 0.333 0.444 0.223
#> 
# }
```
