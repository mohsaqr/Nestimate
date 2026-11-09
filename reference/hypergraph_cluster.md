# Spectral clustering of hypergraph vertices

`hypergraph_cluster()` partitions the nodes of a hypergraph into `k`
clusters with the spectral algorithm of Hayashi et al. (2020,
"RDC-Spec"): it takes the eigenvectors of the `k` smallest eigenvalues
of the normalised Laplacian from
[`hypergraph_laplacian()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_laplacian.md),
scales each node's row of them to unit length, and clusters the rows
with k-means. It works with either Laplacian. With
`type = "random_walk"` and a weighted incidence (for example from
[`bipartite_groups()`](https://pak.dynasite.org/Nestimate/reference/bipartite_groups.md)
with `weight`), the edge-dependent vertex weights can change the
partition; with edge-independent weights the walk reduces to a random
walk on a graph (Chitra & Raphael 2019).

The eigengap, the distance between the `k`-th and the next eigenvalue,
shows how well the spectrum separates `k` groups: a large gap supports
`k`, a small one means the spectrum does not single it out.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
spectrum, so neighbouring values of `k` can be compared. k-means is
stochastic, so a partition that matters is worth refitting with other
seeds to see whether it holds.

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

  The hypergraph: a connected `net_hypergraph`, as for
  [`hypergraph_laplacian()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_laplacian.md).

- k:

  The number of clusters, a whole number from 2 to one less than the
  number of nodes. Required; a value that is not a single whole number
  in that range raises an error.

- type, edge_weights:

  The Laplacian and its hyperedge weights, passed to
  [`hypergraph_laplacian()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_laplacian.md).
  Defaults `"zhou"` and `NULL`.

- nstart:

  The number of random starts k-means tries, keeping the best. Default
  `25`.

- seed:

  A whole number that makes the k-means starts reproducible. The
  session's own random stream is left as it was. Default `NULL`.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html),
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: a
  `net_hypergraph_cluster` object.

- ...:

  For the methods: ignored.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `net_hypergraph_cluster` object.

- what:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: the panel to draw. `"spectrum"` is a scree plot of the
  smallest eigenvalues with the `k` used marked, so the eigengap after
  `k` can be read; `"embedding"` places the nodes on the first two
  embedding dimensions, labelled, coloured and shaped by cluster and
  sized by stationary probability, which is the geometry k-means
  clustered; `"both"` (default) draws the two side by side.

- n_values:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: how many of the smallest eigenvalues the spectrum panel shows.
  Default `NULL`, which shows `min(3 * k, n_nodes)`.

## Value

A `net_hypergraph_cluster` object.
[`print()`](https://rdrr.io/r/base/print.html) reports the Laplacian,
the cluster sizes and the eigengap;
[`summary()`](https://rdrr.io/r/base/summary.html) returns one row per
cluster; [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws
the spectrum and the embedding. The object's fields:

- clusters:

  A data frame with one row per node: `node` and `cluster`. Clusters are
  labelled `"Cluster 1"`, `"Cluster 2"`, ... in the order they first
  appear among the nodes.

- embedding:

  The node-by-`k` spectral embedding k-means clustered, each row of unit
  length, with columns `dim1` to `dimk`.

- k, type:

  The number of clusters and the Laplacian type.

- eigenvalues:

  The full spectrum of the Laplacian, increasing.

- eigengap:

  The `(k + 1)`-th minus the `k`-th eigenvalue.

- sizes:

  A data frame of `cluster` and `size`.

- pi:

  The stationary distribution of the walk, named by node.

- n_nodes, n_hyperedges:

  The size of the hypergraph.

- params:

  The `edge_weights` used, `nstart`, `seed` and the k-means total
  within-cluster sum of squares (`tot_withinss`).

The [`print()`](https://rdrr.io/r/base/print.html) method returns its
input invisibly. The [`summary()`](https://rdrr.io/r/base/summary.html)
method returns a data frame with one row per cluster: `cluster`, `size`
and `share` (the proportion of nodes). The
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) method
returns a data frame with one row per node: `node`, `cluster`, `pi` (the
node's stationary probability under the Laplacian's walk) and the
embedding coordinates `dim1` to `dimk`. The
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method returns
a `ggplot` for `"spectrum"` or `"embedding"`. For `"both"` it draws the
two panels and returns the arranged `gtable` invisibly when gridExtra is
installed, and otherwise draws them in grid viewports and returns the
two `ggplot`s as a list, invisibly.

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
