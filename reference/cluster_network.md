# Cluster data and build per-cluster networks in one step

Combines sequence clustering and network estimation into a single call.
Clusters the data using the specified algorithm, then calls
[`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md)
on each cluster subset.

## Usage

``` r
cluster_network(data, k, cluster_by = "pam", dissimilarity = "hamming", ...)
```

## Arguments

- data:

  Sequence data. Accepts a data frame, matrix, or `netobject`. See
  [`cluster_data`](https://mohsaqr.github.io/Nestimate/reference/cluster_data.md)
  for supported formats.

- k:

  Integer. Number of clusters.

- cluster_by:

  Character. Clustering algorithm passed to
  [`cluster_data`](https://mohsaqr.github.io/Nestimate/reference/cluster_data.md)'s
  `method` parameter (`"pam"`, `"ward.D2"`, `"ward.D"`, `"complete"`,
  `"average"`, `"single"`, `"mcquitty"`, `"median"`, `"centroid"`), or
  `"mmm"` for Mixed Markov Model clustering. Default: `"pam"`.

- dissimilarity:

  Character. Distance metric for sequence clustering (ignored when
  `cluster_by = "mmm"`). Default: `"hamming"`.

- ...:

  Passed directly to
  [`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md).
  Use `method` to specify the network type; `threshold`, `scaling`, and
  all other `build_network` arguments are supported.

## Value

A `netobject_group`.

## Details

If `data` is a `netobject` and `method` is not provided in `...`, the
original network method is inherited automatically so the per-cluster
networks match the type of the input network.

## See also

[`cluster_data`](https://mohsaqr.github.io/Nestimate/reference/cluster_data.md),
[`cluster_mmm`](https://mohsaqr.github.io/Nestimate/reference/cluster_mmm.md),
[`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md)

## Examples

``` r
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 50, TRUE), V2 = sample(LETTERS[1:4], 50, TRUE),
  V3 = sample(LETTERS[1:4], 50, TRUE), V4 = sample(LETTERS[1:4], 50, TRUE)
)
# Default: PAM clustering, relative (transition) networks
grp <- cluster_network(seqs, k = 3)

# Specify network method (cor requires numeric panel data)
if (FALSE) { # \dontrun{
panel <- as.data.frame(matrix(rnorm(1500), nrow = 300, ncol = 5))
grp <- cluster_network(panel, k = 2, method = "cor")
} # }

# MMM-based clustering
grp <- cluster_network(seqs, k = 2, cluster_by = "mmm")
# }
```
