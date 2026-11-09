# Cluster Sequences by Dissimilarity

`build_clusters()` splits a set of sequences into `k` groups of similar
sequences. It computes the dissimilarity between every pair of sequences
(the Hamming distance, an edit distance, a q-gram distance or the
Jaro-Winkler distance) and partitions that matrix with PAM (partitioning
around medoids) or a hierarchical method. It takes wide sequence data, a
`netobject` or bare `cograph_network` built from sequences, or a `tna`
model. A one-hot network built with `"relative"`, `"frequency"` or
`"co_occurrence"` is clustered by unit, each actor's or session's
windows forming one sequence; a `wtna` network is refused (see `data`).
Covariates can be related to the clusters after the fit (see
`covariates`); they never change the assignments.

Every dissimilarity and method returns `k` groups, whether or not the
sequences fall into groups. The average silhouette width (Rousseeuw,
1987), which [`print()`](https://rdrr.io/r/base/print.html) reports,
compares each sequence's mean dissimilarity to its own cluster with that
to the nearest other cluster: it runs from -1 to 1, values near 1
meaning well-separated clusters and values near 0 overlapping ones. It
describes one solution and is not a test. Compare values of `k`,
dissimilarities and methods with
[`cluster_choice()`](https://pak.dynasite.org/Nestimate/reference/cluster_choice.md),
check the classification with
[`cluster_diagnostics()`](https://pak.dynasite.org/Nestimate/reference/cluster_diagnostics.md),
and build one network per cluster with
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
(or
[`cluster_network()`](https://pak.dynasite.org/Nestimate/reference/cluster_network.md)
in one call) to see how the groups differ.

## Usage

``` r
build_clusters(
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
  estimator = c("auto", "firth", "multinom", "chisq"),
  ...
)

# S3 method for class 'net_clustering'
print(x, digits = 3L, ...)

# S3 method for class 'net_clustering'
summary(object, ...)

# S3 method for class 'net_clustering'
plot(
  x,
  type = c("silhouette", "mds", "heatmap", "predictors"),
  combined = TRUE,
  ...
)

# S3 method for class 'tidy_covariates'
print(x, ...)
```

## Arguments

- data:

  The sequences: a data frame or matrix with one row per sequence and
  one column per position, holding state names; a `netobject` or bare
  `cograph_network` built from sequences with method `"relative"`,
  `"frequency"`, `"co_occurrence"` or `"attention"`, whose stored
  sequences are used; or a `tna` model, whose integer-coded sequences
  are decoded to their labels. Columns of a data frame take the same
  roles as in
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md):
  an id or a covariate is not a position, and the columns set aside are
  named in a message and kept as metadata. A one-hot network needs an
  `actor` or `session` to define its units. A one-hot network without
  one, and a network of any other method, raise an error of class
  `nestimate_bad_input`.

- k:

  The number of clusters, a whole number from 2 to one less than the
  number of sequences. Required.

- dissimilarity:

  How two sequences are compared. `"hamming"` (default) counts the
  positions where they differ. `"lv"` (Levenshtein) counts the
  insertions, deletions and substitutions that turn one into the other;
  `"osa"` (optimal string alignment) and `"dl"` (Damerau-Levenshtein)
  also count a swap of two adjacent states as one edit; `"lcs"` counts
  the states left unpaired by their longest common subsequence.
  `"qgram"`, `"cosine"` and `"jaccard"` compare the runs of `q`
  consecutive states (q-grams) the sequences contain: the summed
  absolute difference of their counts, one minus the cosine similarity
  of the counts, and one minus the share of distinct q-grams they have
  in common. `"jw"` is the Jaro-Winkler distance, which weighs agreement
  in the first states more (see `p`).

- method:

  The clustering algorithm. `"pam"` (default) partitions around medoids:
  it chooses `k` sequences as centres and assigns every sequence to the
  nearest one (Kaufman & Rousseeuw, 1990). `"ward.D2"`, `"ward.D"`,
  `"complete"`, `"average"`, `"single"`, `"mcquitty"`, `"median"` and
  `"centroid"` build a tree with
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html) using that
  linkage and cut it into `k` groups.

- na_syms:

  The symbols read as missing values, besides `NA`. Default
  `c("*", "%")`. A missing cell is not dropped from a comparison: it is
  treated as one more state, so two missing cells at the same position
  match and a missing cell against an observed state is a mismatch (a
  distance of 1 for Hamming). This is the usual rule for aligned
  sequences, because dropping missing cells pair by pair would change
  each pair's effective length and break the metric. For every
  dissimilarity except Hamming, a sequence ends at its last observed
  state, so trailing missing cells are not compared. To treat missing
  cells differently, drop or recode them before the call.

- weighted:

  Whether to weight the Hamming distance by position, so that a
  difference early in the sequences counts more than a late one (see
  `lambda`). Default `FALSE`. With any dissimilarity other than
  `"hamming"`, `TRUE` raises an error.

- lambda:

  The decay rate of the position weights when `weighted = TRUE`:
  position \\t\\ gets weight \\e^{-\lambda (t - 1)}\\, so a larger value
  puts more of the distance on the earlier positions. A finite
  non-negative number; default `1`. Ignored, and stored as `0`, when
  `weighted = FALSE`.

- seed:

  An integer used as the random seed for clustering; the session's own
  random stream is left as it was. Default `NULL`. PAM and the
  hierarchical methods are deterministic, so it does not change the
  clusters; it is stored with the result.

- q:

  The length of the q-grams compared by `"qgram"`, `"cosine"` and
  `"jaccard"`, a positive whole number. Default `2`. Ignored by the
  other dissimilarities.

- p:

  The Winkler prefix weight of `"jw"`: how much agreement in the first
  states (up to four) lowers the Jaro distance, from 0 (the plain Jaro
  distance) to 0.25. Default `0.1`. Ignored by the other
  dissimilarities.

- covariates:

  Variables describing each sequence (such as an achievement level or a
  course), related to the clusters after they are found; they never
  change the assignments. Give a column name, a character vector of
  names, a one-sided formula such as `~ Age + Gender` (or the string
  `"Age + Gender"`), or a data frame with one row per sequence, all of
  whose columns are used. For a `netobject` or `cograph_network`, names
  are looked up in its metadata first and then among the columns of its
  data, so a typical call is
  `build_clusters(net, k = 3, covariates = "Achiever")`; for a data
  frame or matrix, among its columns. A `tna` model takes only the
  data-frame form. Named covariate columns are taken out of the sequence
  data before the distances are computed, and so are the columns of a
  supplied data frame that match a sequence column in name and values; a
  sequence column that shares only a name is kept. For a network, a name
  found in its metadata leaves a same-named data column in place as a
  sequence position. Rows with a missing covariate are left out of the
  analysis with a warning, and a constant covariate raises an error. How
  the analysis is run is set by `estimator`. Default `NULL` (none).

- estimator:

  How the relation between the covariates and the clusters is estimated;
  it never affects the clustering. `"auto"` (default) fits a multinomial
  logit of cluster membership on the covariates: with Firth's
  bias-reduced penalised likelihood
  ([`brglm2::brmultinom()`](https://rdrr.io/pkg/brglm2/man/brmultinom.html),
  finite under quasi-complete separation) when a cell of the
  cluster-by-covariate table of a categorical covariate holds fewer than
  5 sequences, and with
  [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html)
  otherwise. `"firth"` forces the Firth fit, which is about 200 times
  slower on well-conditioned data. `"multinom"` forces classical maximum
  likelihood and warns that rare cells can give extreme odds ratios with
  degenerate confidence intervals. `"chisq"` runs descriptive tests
  instead of a logit, as WeightedCluster does: for each categorical
  covariate a chi-square test with Cramer's V and standardised adjusted
  residuals, for each numeric one a Kruskal-Wallis test with
  eta-squared. Ignored without `covariates`.

- ...:

  Not used: an extra argument raises an error, here and in the
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods for a
  `net_clustering`. [`print()`](https://rdrr.io/r/base/print.html) of a
  `tidy_covariates` table ignores it.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: a
  `net_clustering`; for [`print()`](https://rdrr.io/r/base/print.html),
  also the `tidy_covariates` table that
  [`summary()`](https://rdrr.io/r/base/summary.html) returns.

- digits:

  For [`print()`](https://rdrr.io/r/base/print.html) of a
  `net_clustering`: decimal places of the silhouette width and the
  within-cluster distances, a non-negative whole number. Default `3`.

- object:

  For [`summary()`](https://rdrr.io/r/base/summary.html): a
  `net_clustering`.

- type:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): the
  view. `"silhouette"` (default) draws each sequence's silhouette width
  as a bar, grouped by cluster, with the average as a dashed line;
  `"mds"` places the sequences on the first two dimensions of a
  classical multidimensional scaling of the distances, coloured by
  cluster; `"heatmap"` draws the distance matrix with the sequences
  ordered by cluster; `"predictors"` draws a forest plot of the
  covariates' odds ratios (log scale) for each cluster against the
  reference cluster. `"predictors"` needs a fit with `covariates` and a
  logit estimator; with `estimator = "chisq"`, which has no odds ratios,
  it raises an error.

- combined:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html) with
  `type = "predictors"`: `TRUE` (default) draws one plot with a panel
  per non-reference cluster; `FALSE` returns a list of plots, one per
  non-reference cluster. Ignored by the other views.

## Value

A `net_clustering` object.
[`print()`](https://rdrr.io/r/base/print.html) shows the method,
dissimilarity and average silhouette width, then one row per cluster
with its size (count and percentage), mean within-cluster distance and,
for PAM, its medoid; after a logit covariate analysis it adds a line
naming the predictors.
[`summary()`](https://rdrr.io/r/base/summary.html) prints the method,
dissimilarity and average silhouette width and a per-cluster table of
sizes and mean within-cluster distances, followed by the covariate
analysis when there is one, and returns a per-cluster table, or the
covariate table after a covariate analysis (see the methods' return
value below). [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
draws the silhouette widths, a map of the distances, the distance matrix
or the covariate effects.
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
turns the clustering into one network per cluster (a `netobject_group`),
and
[`cluster_diagnostics()`](https://pak.dynasite.org/Nestimate/reference/cluster_diagnostics.md)
and
[`session_ids()`](https://pak.dynasite.org/Nestimate/reference/session_ids.md)
give the classification quality and each sequence's cluster. The
object's fields:

- data:

  The sequence positions that were clustered, one row per sequence:
  covariates and set-aside columns are removed, and a one-hot network
  gives one row per unit.

- k:

  Number of clusters.

- assignments:

  Integer vector of cluster memberships (1 to `k`), one per sequence in
  row order.

- silhouette:

  Average silhouette width over all sequences.

- sizes:

  Named integer vector of cluster sizes.

- method:

  Clustering method used.

- dissimilarity:

  Dissimilarity used.

- distance:

  The pairwise dissimilarities, a `dist` object.

- medoids:

  For PAM, the row numbers of the `k` medoid sequences; NULL for the
  hierarchical methods.

- seed:

  Seed given (or NULL).

- weighted:

  Whether weighted Hamming was used.

- lambda:

  Decay rate used (0 when not weighted).

- covariates:

  The covariate analysis (a list; see `estimator`), or NULL when
  `covariates = NULL`.

- network_method, build_args:

  For `netobject` input, the source network's method and stored build
  arguments, so per-cluster networks can be rebuilt the same way. NULL
  otherwise, and `build_args` is NULL for a one-hot network.

- metadata:

  For `netobject` input, its per-sequence metadata, one row per
  clustered sequence, so
  [`session_ids`](https://pak.dynasite.org/Nestimate/reference/session_ids.md)
  can name each sequence. For a data frame, the columns set aside as
  covariates (see `data`). NULL otherwise, and for a one-hot network.

- htna_partition:

  For HTNA input, the preserved node-to-actor partition used to restore
  HTNA children when networks are built.

For the [`print()`](https://rdrr.io/r/base/print.html) methods: the
input, invisibly. [`print()`](https://rdrr.io/r/base/print.html) of a
`tidy_covariates` table prints a header naming the estimator, then the
table. For [`summary()`](https://rdrr.io/r/base/summary.html): a data
frame with one row per cluster and columns `cluster`, `size` and
`mean_within_dist`, returned visibly after the printed summary. For a
clustering fitted with `covariates`, the printed summary adds the
covariate profiles and the logit or test tables, and a `tidy_covariates`
data frame is returned invisibly instead: one row per cluster and
predictor (`cluster`, `predictor`, `OR`, `CI_low`, `CI_high`, `p`,
`sig`) for a logit, or one row per covariate (`variable`, `type`,
`statistic`, `df`, `p`, `effect`, `effect_kind`, `sig`) for
`estimator = "chisq"`, with the per-cluster table, fit statistics and
profiles attached as attributes. For
[`plot()`](https://rdrr.io/r/graphics/plot.default.html): a `ggplot`
object, drawn and returned invisibly; with `type = "predictors"` and
`combined = FALSE`, a list of `ggplot` objects, one per non-reference
cluster and named by it, returned invisibly without drawing.

## Details

With the stringdist package installed, the distances are computed by
[`stringdist::stringdistmatrix()`](https://rdrr.io/pkg/stringdist/man/stringdist.html)
in C, many times faster than the base R implementation for the edit
distances. The base R implementation is used when stringdist is absent,
for weighted Hamming, and when the data hold 62 or more states.

## References

Kaufman, L., & Rousseeuw, P. J. (1990). *Finding Groups in Data: An
Introduction to Cluster Analysis*. Wiley.

Rousseeuw, P. J. (1987). Silhouettes: A graphical aid to the
interpretation and validation of cluster analysis. *Journal of
Computational and Applied Mathematics*, 20, 53–65.
[doi:10.1016/0377-0427(87)90125-7](https://doi.org/10.1016/0377-0427%2887%2990125-7)

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","C","A","B"), V2 = c("B","C","A","B","A"),
                   V3 = c("C","A","B","C","B"))
cl <- build_clusters(seqs, k = 2)
cl
#> Sequence Clustering [pam]
#>   Sequences: 5  |  Clusters: 2
#>   Dissimilarity: hamming
#>   Quality: silhouette = 0.600
#> 
#>   Cluster  N          Mean within-dist  Medoid
#>   1        2 (40.0%)  0.000             4
#>   2        3 (60.0%)  2.000             5
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:3], 20, TRUE), V2 = sample(LETTERS[1:3], 20, TRUE),
  V3 = sample(LETTERS[1:3], 20, TRUE), V4 = sample(LETTERS[1:3], 20, TRUE)
)
cl <- build_clusters(seqs, k = 2)
print(cl)
#> Sequence Clustering [pam]
#>   Sequences: 20  |  Clusters: 2
#>   Dissimilarity: hamming
#>   Quality: silhouette = 0.315
#> 
#>   Cluster  N           Mean within-dist  Medoid
#>   1        9 (45.0%)   2.333             5
#>   2        11 (55.0%)  2.000             15
summary(cl)
#> Sequence Clustering Summary
#>   Method:        pam 
#>   Dissimilarity: hamming 
#>   Silhouette:    0.315 
#> 
#> Per-cluster statistics:
#>  cluster size mean_within_dist
#>        1    9         2.333333
#>        2   11         2.000000
#>   cluster size mean_within_dist
#> 1       1    9         2.333333
#> 2       2   11         2.000000
# }
```
