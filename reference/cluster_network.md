# Cluster data and build per-cluster networks in one step

`cluster_network()` clusters a set of sequences and estimates one
network per cluster in a single call. It runs
[`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
(dissimilarity clustering with PAM or a hierarchical method) or, with
`cluster_by = "mmm"`,
[`build_mmm()`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
(a mixture of Markov chains), then
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
on each cluster's sequences. The network method is the `method` given in
`...`; without one, a `netobject` input passes on its own method and
stored build arguments, so the per-cluster networks are of the same type
as the input, and other input gives `"relative"` transition networks.

The clusters were formed to make the sequences within each one alike, so
their networks differ partly by construction. Check the classification
with
[`cluster_diagnostics()`](https://pak.dynasite.org/Nestimate/reference/cluster_diagnostics.md)
and the choice of `k` with
[`cluster_choice()`](https://pak.dynasite.org/Nestimate/reference/cluster_choice.md)
or
[`compare_mmm()`](https://pak.dynasite.org/Nestimate/reference/compare_mmm.md),
and see how far each cluster's edges move under resampling with
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
before reading differences between the networks.

## Usage

``` r
cluster_network(data, k, cluster_by = "pam", dissimilarity = "hamming", ...)
```

## Arguments

- data:

  The sequences: a data frame, matrix, `netobject`, bare
  `cograph_network` or `tna` model, as accepted by
  [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  (or by
  [`build_mmm()`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
  when `cluster_by = "mmm"`). A `wtna` network is accepted only with
  `cluster_by = "mmm"`; dissimilarity clustering refuses it with an
  error of class `nestimate_bad_input`.

- k:

  The number of clusters, a whole number from 2 to one less than the
  number of sequences. Required.

- cluster_by:

  The clustering algorithm. `"pam"` (default) and the hierarchical
  methods `"ward.D2"`, `"ward.D"`, `"complete"`, `"average"`,
  `"single"`, `"mcquitty"`, `"median"` and `"centroid"` are passed to
  [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  as its `method`; `"mmm"` fits a mixed Markov model with
  [`build_mmm()`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
  instead.

- dissimilarity:

  The distance between sequences for dissimilarity clustering, one of
  those
  [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  accepts. Default `"hamming"`. Giving it with `cluster_by = "mmm"`
  raises an error.

- ...:

  Arguments for the two stages, routed by name. For dissimilarity
  clustering,
  [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)'s
  `na_syms`, `weighted`, `lambda`, `seed`, `q`, `p`, `covariates` and
  `estimator` go to the clustering; with `cluster_by = "mmm"`,
  [`build_mmm()`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)'s
  `n_starts`, `max_iter`, `tol`, `smooth`, `seed`, `covariates`,
  `covariate_effect` and `estimator` go to the mixture fit. Everything
  else is meant for the network stage and must be an argument of
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  or of a registered estimator; any other name raises an error before
  either stage runs. A name both stages accept goes to the clustering:
  under dissimilarity clustering a `lambda` given here is the Hamming
  decay of the clustering and never reaches an attention network, and
  the same holds for `weighted`. When `data` is a `netobject`, its
  stored build arguments are reused for the per-cluster networks, and a
  value given in the call replaces the stored one. The *Details* section
  lists which network arguments take effect.

## Value

A `netobject_group`: one `netobject` per cluster, named `"Cluster 1"`,
`"Cluster 2"`, and so on, with the clustering (a `net_clustering`, or
for `"mmm"` a `net_mmm_clustering`) kept as its `"clustering"`
attribute. [`print()`](https://rdrr.io/r/base/print.html) names the
clustering method and shows each cluster's network size, weight range
and share of the sequences;
[`cluster_diagnostics()`](https://pak.dynasite.org/Nestimate/reference/cluster_diagnostics.md)
reads the clustering, and
[`summary()`](https://rdrr.io/r/base/summary.html) (the networks'
metrics, a column per cluster),
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
and
[`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
take the group as it is.

## Details

`method`, the estimators' own arguments (such as `n_gram` or
`direction`) and
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)'s
`threshold`, `scaling`, `start`, `end` and `params` reach the
per-cluster networks. A `netobject` input of the same method also passes
on its estimator settings (`$params`, such as an attention network's
`lambda`) unless `params` is given. With `cluster_by = "mmm"` and the
default relative networks, which are the fitted chains, `scaling` and
`threshold` are applied to the fitted weights, and an estimator
argument, `start`, `end` or a new `params` raises an error of class
`nestimate_bad_input`. For a one-hot `netobject` the stored build
arguments are not used. Under `"mmm"`, rebuilding a network built from
long-format data with any method other than `"relative"` currently
raises an error, because its column names are passed on to the wide
per-cluster data.

## See also

[`build_clusters`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md),
[`cluster_mmm`](https://pak.dynasite.org/Nestimate/reference/cluster_mmm.md),
[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md)

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","C","A","B"), V2 = c("B","C","A","B","A"),
                   V3 = c("C","A","B","C","B"))
grp <- cluster_network(seqs, k = 2)
grp
#> Group Networks (2 clusters via pam / hamming)
#> 
#>   Group      Nodes  Edges  Weights         N
#>   Cluster 1  3      2      [1.000, 1.000]  2 (40.0%)
#>   Cluster 2  3      4      [0.500, 1.000]  3 (60.0%)
# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 50, TRUE), V2 = sample(LETTERS[1:4], 50, TRUE),
  V3 = sample(LETTERS[1:4], 50, TRUE), V4 = sample(LETTERS[1:4], 50, TRUE)
)
# Default: PAM clustering, relative (transition) networks
grp <- cluster_network(seqs, k = 3)


# MMM-based clustering
grp <- cluster_network(seqs, k = 2, cluster_by = "mmm")
# }
```
