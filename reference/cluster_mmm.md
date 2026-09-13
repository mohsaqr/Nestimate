# Cluster sequences using Mixed Markov Models

Fits a mixture of Markov chains to sequence data and returns the fitted
`net_mmm` clustering object. The fit retains assignments, posterior
probabilities, mixing proportions, information criteria, and the fitted
component models.

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
  covariates = NULL,
  covariate_effect = c("em", "posthoc"),
  estimator = c("auto", "firth", "multinom", "chisq"),
  cluster_by = "mmm",
  ...
)
```

## Arguments

- data:

  A data.frame (wide format), `netobject`, `cograph_network`, or `tna`
  model. For tna and cograph_network objects the stored
  (integer-encoded) data is extracted and decoded to state labels.

- k:

  Integer. Whole finite number of mixture components, \>= 2. Default: 2.

- n_starts:

  Integer. Positive whole finite number of random restarts. Default: 50.

- max_iter:

  Integer. Positive whole finite maximum EM iterations per start.
  Default: 200.

- tol:

  Numeric. Finite positive convergence tolerance. Default: 1e-6.

- smooth:

  Numeric. Finite non-negative Laplace smoothing constant. Default:
  0.01.

- seed:

  Integer or NULL. Random seed.

- covariates:

  Optional. Covariates integrated into the EM algorithm to model
  covariate-dependent mixing proportions. Accepts a string, character
  vector, formula, or data.frame (same forms as
  [`build_clusters`](https://saqr.me/Nestimate/reference/build_clusters.md)).
  For `netobject` or `cograph_network` input, names are resolved against
  `$metadata` first, so a typical call is
  `build_mmm(net, k = 3, covariates = "session_label")`. Unlike the
  post-hoc analysis in
  [`build_clusters()`](https://saqr.me/Nestimate/reference/build_clusters.md),
  these covariates directly influence cluster membership during EM
  estimation (see `covariate_effect`).

- covariate_effect:

  How `covariates` enter the model. `"em"` (default) folds them into the
  EM as covariate-dependent mixing proportions, so they shape the
  cluster fit itself (and rows with missing covariates are dropped
  before fitting). `"posthoc"` fits a plain mixture on every sequence
  and uses the covariates only for the after-fit multinomial logit, so
  covariate values — and their missingness — never change which clusters
  are found. Ignored when `covariates` is `NULL`.

- estimator:

  Multinomial fitter for the post-hoc covariate analysis (does not
  affect EM): `"auto"` (default) inspects the cluster x covariate
  cross-tab and falls back to `"firth"` only when any cell has fewer
  than 5 observations (separation risk), otherwise the much faster
  `"multinom"`; `"firth"` forces Firth's penalised likelihood via
  [`brglm2::brmultinom`](https://rdrr.io/pkg/brglm2/man/brmultinom.html)
  (finite under separation); `"multinom"` forces
  [`nnet::multinom`](https://rdrr.io/pkg/nnet/man/multinom.html) (warns
  about separation risk); `"chisq"` runs descriptive tests (no logit).
  See
  [`build_clusters`](https://saqr.me/Nestimate/reference/build_clusters.md)
  for full details.

- cluster_by:

  Character. Accepted only as `"mmm"` (the default). Present so
  `cluster_mmm()` and
  [`cluster_network()`](https://saqr.me/Nestimate/reference/cluster_network.md)
  share the same call shape; any other value raises an error pointing at
  [`cluster_network`](https://saqr.me/Nestimate/reference/cluster_network.md).

- ...:

  Unsupported. Supplying unused arguments raises an error.

## Value

A fitted `net_mmm` clustering object. This is the same object contract
returned by
[`build_mmm`](https://saqr.me/Nestimate/reference/build_mmm.md). For
HTNA input, its preserved actor partition is restored when the fit is
materialized with
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md)
or
[`Nestimate::as_htna()`](https://saqr.me/Nestimate/reference/as_htna.md).

## Details

To materialize one network per fitted cluster, pass the result to
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md)
or use `cluster_network(..., cluster_by = "mmm")` for fitting and
network construction in one call.

## See also

[`build_mmm`](https://saqr.me/Nestimate/reference/build_mmm.md),
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md),
and
[`cluster_network`](https://saqr.me/Nestimate/reference/cluster_network.md)
for fitting and immediately materializing per-cluster networks

## Examples

``` r
seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
                   V2 = sample(c("A","B","C"), 30, TRUE))
fit <- cluster_mmm(seqs, k = 2, n_starts = 1, max_iter = 10, seed = 1)
fit
#> Mixed Markov Model
#>   Sequences: 30  |  Clusters: 2  |  States: 3
#>   ICs: LL = -59.844  |  BIC = 177.508  |  AIC = 153.688  |  ICL = 179.621
#>   Quality: AvePP = 0.966  |  Entropy = 0.202  |  Class.Err = 0.0%
#>   Status: did not converge in 10 iterations
#> 
#>   Cluster  N           Mix%   AvePP
#>   1        26 (86.7%)  84.8%  0.969
#>   2        4 (13.3%)   15.2%  0.943
cluster_diagnostics(fit)
#> Cluster Diagnostics (mmm) [k = 2]
#>   Sequences: 30  |  Clusters: 2  |  States: 3
#>   Quality: AvePP = 0.966  |  Entropy = 0.202  |  Class.Err = 0.0%
#>   ICs: LL = -59.844  |  BIC = 177.508  |  AIC = 153.688  |  ICL = 179.621
#> 
#>   Cluster  N           Mix%   AvePP  Class.Err%
#>   1        26 (86.7%)  84.8%  0.969   0.0%
#>   2        4 (13.3%)   15.2%  0.943   0.0%
# \donttest{
# Visualise with sequence_plot
seqs <- data.frame(
  V1 = sample(LETTERS[1:3], 40, TRUE),
  V2 = sample(LETTERS[1:3], 40, TRUE),
  V3 = sample(LETTERS[1:3], 40, TRUE)
)
fit <- cluster_mmm(seqs, k = 2)
sequence_plot(fit, type = "index")

# }
```
