# Cluster sequences using Mixed Markov Models

Convenience alias for
[`build_mmm`](https://mohsaqr.github.io/Nestimate/reference/build_mmm.md).
Fits a mixture of Markov chains to sequence data and returns
per-component transition networks with EM-fitted initial state
probabilities.

## Usage

``` r
cluster_mmm(
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

A `net_mmm` object. See
[`build_mmm`](https://mohsaqr.github.io/Nestimate/reference/build_mmm.md)
for details.

## Details

Use
[`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md)
on the result to extract per-cluster networks with any estimation
method, or use
[`cluster_network`](https://mohsaqr.github.io/Nestimate/reference/cluster_network.md)
for a one-shot clustering + network call.

## See also

[`build_mmm`](https://mohsaqr.github.io/Nestimate/reference/build_mmm.md),
[`cluster_network`](https://mohsaqr.github.io/Nestimate/reference/cluster_network.md)

## Examples

``` r
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:3], 40, TRUE),
  V2 = sample(LETTERS[1:3], 40, TRUE),
  V3 = sample(LETTERS[1:3], 40, TRUE)
)
mmm <- cluster_mmm(seqs, k = 2)
print(mmm)
#> Mixed Markov Model
#>   k = 2 | 40 sequences | 3 states
#>   LL = -122.9 | BIC = 308.6 | ICL = 314.2
#> 
#>   Cluster  Size  Mix%%   AvePP
#>   ------------------------------
#>         1    23  56.6%  0.944
#>         2    17  43.4%  0.945
#> 
#>   Overall AvePP = 0.944 | Entropy = 0.153 | Class.Err = 0.0%
# }
```
