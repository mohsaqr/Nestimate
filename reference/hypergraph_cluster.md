# Spectral clustering of hypergraph vertices

Partitions the nodes of a hypergraph into `k` clusters with the
Laplacian-eigenmap + k-means algorithm of Hayashi et al. (2020,
"RDC-Spec"): the eigenvectors of the `k` smallest eigenvalues of the
normalized hypergraph Laplacian are row-normalized to unit length and
clustered with k-means. With `type = "random_walk"` and a weighted
incidence (e.g. from
[`bipartite_groups()`](https://pak.dynasite.org/Nestimate/reference/bipartite_groups.md)
with `weight =`), the edge-dependent vertex weights genuinely change the
partition - with edge-independent weights the walk collapses to a graph
random walk (Chitra & Raphael 2019).

## Usage

``` r
hypergraph_cluster(
  hg,
  k,
  type = c("zhou", "random_walk"),
  edge_weights = NULL,
  nstart = 25L,
  seed = NULL
)

# S3 method for class 'net_hypergraph_cluster'
print(x, ...)

# S3 method for class 'net_hypergraph_cluster'
summary(object, ...)

# S3 method for class 'net_hypergraph_cluster'
as.data.frame(x, ...)

# S3 method for class 'net_hypergraph_cluster'
plot(x, what = c("both", "spectrum", "embedding"), n_values = NULL, ...)
```

## Arguments

- hg:

  A connected `net_hypergraph`.

- k:

  Integer number of clusters, `2 <= k <= n_nodes - 1`.

- type, edge_weights:

  Passed to
  [`hypergraph_laplacian()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_laplacian.md).

- nstart:

  Integer. k-means random restarts (default 25).

- seed:

  Optional integer seed for the k-means initialization.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html),
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: an
  object of class `net_hypergraph_cluster`.

- ...:

  For the S3 methods: further arguments passed to or from other methods.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: an
  object of class `net_hypergraph_cluster`.

- what:

  Character. `"both"` (default), `"spectrum"`, or `"embedding"`.

- n_values:

  Integer. How many smallest eigenvalues to show in the spectrum panel
  (default: `min(3 * k, n_nodes)`).

## Value

An object of class `net_hypergraph_cluster`: a list with `$clusters`
(data.frame, one row per node: `node`, `cluster` - labels `"Cluster 1"`,
`"Cluster 2"`, ... ordered by first appearance), `$embedding` (node x k
row-normalized spectral embedding used by k-means, dims `dim1..dimk`),
`$k`, `$type`, `$eigenvalues` (full Laplacian spectrum, increasing),
`$eigengap` (gap after the k-th eigenvalue), `$sizes` (data.frame
`cluster`/`size`), `$pi` (named stationary distribution), `$n_nodes`,
`$n_hyperedges` and `$params` (the `edge_weights` used, `nstart`,
`seed`, `tot_withinss`). Has `print`, `summary`, `plot` and
`as.data.frame` methods;
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
one row per node with `node`, `cluster`, the stationary probability
`pi`, and the embedding coordinates.

In `print.net_hypergraph_cluster()`: The input object, invisibly.

In `summary.net_hypergraph_cluster()`: A data.frame, one row per
cluster: `cluster`, `size`, `share`.

In `as.data.frame.net_hypergraph_cluster()`: The tidy assignment table:
one row per node, columns `node`, `cluster`, `pi` (stationary
probability of the node under the Laplacian's random walk) and the
spectral-embedding coordinates `dim1..dimk`.

In `plot.net_hypergraph_cluster()`: For `"spectrum"`/`"embedding"`, the
ggplot object. For `"both"`, the arranged gtable when gridExtra is
installed (drawn on the current device), otherwise the two panels are
drawn via grid viewports and the list of the two ggplots is returned
invisibly.

## Details

k-means is stochastic: `nstart` restarts are used and a `seed` fixes the
result. Report stability across seeds for consequential results.

## References

Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
random walks, Laplacians, and clustering. *CIKM 2020*, 495-504.
[doi:10.1145/3340531.3412034](https://doi.org/10.1145/3340531.3412034)

Chitra, U., & Raphael, B. J. (2019). Random walks on hypergraphs with
edge-dependent vertex weights. *ICML 2019*.

## Examples

``` r
events <- data.frame(
  person = c("a", "b", "c", "a", "b", "c", "d", "e", "f",
             "d", "e", "f", "c", "d"),
  meeting = c("m1", "m1", "m1", "m2", "m2", "m2", "m3", "m3", "m3",
              "m4", "m4", "m4", "m5", "m5")
)
hg <- bipartite_groups(events, player = "person", group = "meeting")
cl <- hypergraph_cluster(hg, k = 2, seed = 1)
cl
#> Hypergraph spectral clustering (zhou Laplacian)
#>   Nodes: 6 | Hyperedges: 5 | k: 2
#>   Cluster sizes: Cluster 1 = 3, Cluster 2 = 3
#>   Eigengap after k: 0.6667
as.data.frame(cl)
#>   node   cluster        pi       dim1       dim2
#> 1    a Cluster 1 0.1428571 -0.6575959  0.7533708
#> 2    b Cluster 1 0.1428571 -0.6575959  0.7533708
#> 3    c Cluster 1 0.2142857 -0.7947194  0.6069770
#> 4    d Cluster 2 0.2142857 -0.7947194 -0.6069770
#> 5    e Cluster 2 0.1428571 -0.6575959 -0.7533708
#> 6    f Cluster 2 0.1428571 -0.6575959 -0.7533708
```
