# Cluster Sequences by Dissimilarity

Clusters wide-format sequences using pairwise string dissimilarity and
either PAM (Partitioning Around Medoids) or hierarchical clustering.
Supports 9 distance metrics including temporal weighting for Hamming
distance. When the stringdist package is available, uses C-level
distance computation for 100-1000x speedup on edit distances.

## Usage

``` r
cluster_data(
  data,
  k,
  dissimilarity = "hamming",
  method = "pam",
  na_syms = c("*", "%"),
  weighted = FALSE,
  lambda = 1,
  seed = NULL,
  q = 2L,
  p = 0.1,
  covariates = NULL,
  ...
)

cluster_sequences(
  data,
  k,
  dissimilarity = "hamming",
  method = "pam",
  na_syms = c("*", "%"),
  weighted = FALSE,
  lambda = 1,
  seed = NULL,
  q = 2L,
  p = 0.1,
  covariates = NULL,
  ...
)
```

## Arguments

- data:

  Input data. Accepts multiple formats:

  data.frame / matrix

  :   Wide-format sequences (rows = sequences, columns = time points,
      values = state names).

  netobject

  :   A network object from
      [`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md).
      Extracts the stored sequence data. Only valid for sequence-based
      methods (relative, frequency, co_occurrence, attention).

  tna

  :   A tna model from the tna package. Decodes the integer-encoded
      sequence data using stored labels.

  cograph_network

  :   A cograph network object. Extracts the stored sequence data.

- k:

  Integer. Number of clusters (must be between 2 and `nrow(data) - 1`).

- dissimilarity:

  Character. Distance metric. One of `"hamming"`, `"osa"` (optimal
  string alignment), `"lv"` (Levenshtein), `"dl"` (Damerau-Levenshtein),
  `"lcs"` (longest common subsequence), `"qgram"`, `"cosine"`,
  `"jaccard"`, `"jw"` (Jaro-Winkler). Default: `"hamming"`.

- method:

  Character. Clustering method. `"pam"` for Partitioning Around Medoids,
  or a hierarchical method: `"ward.D2"`, `"ward.D"`, `"complete"`,
  `"average"`, `"single"`, `"mcquitty"`, `"median"`, `"centroid"`.
  Default: `"pam"`.

- na_syms:

  Character vector. Symbols treated as missing values. Default:
  `c("*", "%")`.

- weighted:

  Logical. Apply exponential decay weighting to Hamming distance
  positions? Only valid when `dissimilarity = "hamming"`. Default:
  `FALSE`.

- lambda:

  Numeric. Decay rate for weighted Hamming. Higher values weight earlier
  positions more strongly. Default: 1.

- seed:

  Integer or NULL. Random seed for reproducibility. Default: `NULL`.

- q:

  Integer. Size of q-grams for `"qgram"`, `"cosine"`, and `"jaccard"`
  distances. Default: `2L`.

- p:

  Numeric. Winkler prefix penalty for Jaro-Winkler distance (clamped to
  0–0.25). Default: `0.1`.

- covariates:

  Optional. Post-hoc covariate analysis of cluster membership via
  multinomial logistic regression. Accepts:

  formula

  :   `~ Age + Gender`

  character vector

  :   `c("Age", "Gender")`

  string

  :   `"Age + Gender"`

  data.frame

  :   All columns used as covariates

  NULL

  :   No covariate analysis (default)

  Covariates are looked up in `netobject$metadata` or non-sequence
  columns of the input data. For `tna` and `cograph_network` inputs,
  pass covariates as a data.frame. Results stored in `$covariates`.
  Requires the nnet package.

- ...:

  Additional arguments (currently unused).

## Value

An object of class `"net_clustering"` containing:

- data:

  The original input data.

- k:

  Number of clusters.

- assignments:

  Named integer vector of cluster assignments.

- silhouette:

  Overall average silhouette width.

- sizes:

  Named integer vector of cluster sizes.

- method:

  Clustering method used.

- dissimilarity:

  Distance metric used.

- distance:

  The computed dissimilarity matrix (`dist` object).

- medoids:

  Integer vector of medoid row indices (PAM only; NULL for hierarchical
  methods).

- seed:

  Seed used (or NULL).

- weighted:

  Logical, whether weighted Hamming was used.

- lambda:

  Lambda value used (0 if not weighted).

## Examples

``` r
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:3], 20, TRUE), V2 = sample(LETTERS[1:3], 20, TRUE),
  V3 = sample(LETTERS[1:3], 20, TRUE), V4 = sample(LETTERS[1:3], 20, TRUE)
)
cl <- cluster_data(seqs, k = 2)
print(cl)
#> Sequence Clustering
#>   Method:        pam 
#>   Dissimilarity: hamming  
#>   Clusters:      2 
#>   Silhouette:    0.2704 
#>   Cluster sizes: 8, 12 
#>   Medoids:       18, 16 
summary(cl)
#> Sequence Clustering Summary
#>   Method:        pam 
#>   Dissimilarity: hamming 
#>   Silhouette:    0.2704 
#> 
#> Per-cluster statistics:
#>  cluster size mean_within_dist
#>        1    8         2.428571
#>        2   12         2.136364
# }
```
