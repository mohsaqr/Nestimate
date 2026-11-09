# Higher-order hypergraph from a network's clique structure

`build_hypergraph()` turns a network into a hypergraph by promoting its
cliques to hyperedges: every set of 3 to `max_size` nodes that are all
tied to each other becomes a hyperedge, kept with probability `p`, and
the pairwise edges can be kept as hyperedges of size 2. It takes a
`netobject`, a bare `cograph_network`, a `simplicial_complex` or an
adjacency or weight matrix, and returns a `net_hypergraph`, the input of
the hypergraph verbs such as
[`hypergraph_measures()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_measures.md),
[`hypergraph_centrality()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_centrality.md)
and
[`clique_expansion()`](https://pak.dynasite.org/Nestimate/reference/clique_expansion.md).

A hyperedge built this way records that every pair of its nodes is tied
in the thresholded network; a pairwise network cannot show whether those
nodes also act together as a group. When group memberships are observed
directly,
[`bipartite_groups()`](https://pak.dynasite.org/Nestimate/reference/bipartite_groups.md)
builds the hypergraph from them instead. Each clique also contributes
its smaller cliques: with `max_size = 4`, a 4-clique gives one hyperedge
of size 4 and four of size 3.

## Usage

``` r
build_hypergraph(
  net,
  p = 1,
  method = c("clique", "vr", "rips"),
  include_pairwise = TRUE,
  max_size = 3L,
  threshold = 0,
  seed = NULL
)

# S3 method for class 'net_hypergraph'
print(x, ...)

# S3 method for class 'net_hypergraph'
summary(object, ...)
```

## Arguments

- net:

  The network: a `netobject`, bare `cograph_network`,
  `simplicial_complex` (its edges are read as ties of weight 1 and the
  cliques are found again up to `max_size`), or a numeric or logical
  adjacency or weight matrix (nodes are named `V1`, `V2`, ... when it
  has no row names). A tie counts in either direction, so a directed
  network is symmetrised, and self-loops are ignored.

- p:

  The probability that each clique of 3 or more nodes becomes a
  hyperedge, a number in `[0, 1]`; each clique is kept or dropped
  independently. Default `1` (every clique is kept). Pairwise edges are
  not sampled.

- method:

  How hyperedges are found. `"clique"` (default) promotes the cliques of
  the thresholded network. `"vr"` and `"rips"` (a metric Vietoris-Rips
  construction) are not implemented and raise an error instead of
  falling back to `"clique"`.

- include_pairwise:

  Whether the pairwise edges of the thresholded network become
  hyperedges of size 2. `FALSE` gives a hypergraph of hyperedges of size
  3 or more only. Default `TRUE`.

- max_size:

  The largest hyperedge size, at least 2. Default `3` (triangles); `4`
  adds the 4-cliques, and so on. With `2`, only the pairwise edges
  remain (none, with `include_pairwise = FALSE`).

- threshold:

  The cut-off that turns weights into ties: an edge is kept when its
  absolute weight is non-zero and at least `threshold`, so negative
  weights count by their size. A negative value raises an error. Default
  `0` (every non-zero weight is a tie).

- seed:

  An integer that makes the sampling reproducible when `0 < p < 1`. The
  session's own random stream is left as it was. Default `NULL`.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) method: a
  `net_hypergraph`.

- ...:

  Not used.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `net_hypergraph`.

## Value

A `net_hypergraph` object.
[`print()`](https://rdrr.io/r/base/print.html) shows the numbers of
nodes and hyperedges, the size distribution and the settings, and
returns its input invisibly.
[`summary()`](https://rdrr.io/r/base/summary.html) prints the counts
with the mean and largest hyperedge size, and returns a data frame with
one row per node: `node` and `degree`, the number of hyperedges that
contain it. The object holds:

- `hyperedges`:

  A list of integer vectors, one per hyperedge, each the sorted indices
  of the nodes it spans; pairwise edges, when kept, come first.

- `incidence`:

  An integer matrix of nodes by hyperedges, 1 where a node belongs to a
  hyperedge and 0 elsewhere. Its row names are the node names and its
  column names `h1`, `h2`, ...

- `nodes`:

  The node names.

- `n_nodes`, `n_hyperedges`:

  The numbers of nodes and hyperedges.

- `size_distribution`:

  A named integer vector with the number of hyperedges of each size,
  named `size_2`, `size_3`, ...

- `params`:

  The settings used: `method`, `p`, `include_pairwise`, `max_size`,
  `threshold` and `seed`.

## Details

The construction follows Burgio, Matamalas, Gomez and Arenas (2020) on
cooperation in networks and hypergraphs. The cliques are found with
[`build_simplicial()`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md),
and each clique of 3 or more nodes then goes through an independent
Bernoulli(`p`) trial that decides whether it becomes a hyperedge. The
pairwise edges are not sampled: with `include_pairwise = TRUE` all of
them are kept, so the hypergraph holds the pairwise structure and the
sampled higher-order interactions together.

At the limits:

- `p = 0` with `include_pairwise = TRUE` gives the thresholded pairwise
  network as a hypergraph of size-2 edges.

- `p = 1` with `include_pairwise = FALSE` gives only the hyperedges of
  size 3 or more, every clique of the network's clique complex up to
  `max_size`.

## References

Burgio, G., Matamalas, J. T., Gomez, S., & Arenas, A. (2020). Evolution
of cooperation in the presence of higher-order interactions: from
networks to hypergraphs. *Entropy* 22(7), 744.
[doi:10.3390/e22070744](https://doi.org/10.3390/e22070744)

## See also

[`build_simplicial()`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md)
(the clique enumeration),
[`bipartite_groups()`](https://pak.dynasite.org/Nestimate/reference/bipartite_groups.md),
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md).

## Examples

``` r
set.seed(1)
n <- 8
adj <- matrix(stats::rbinom(n * n, 1, 0.5), n, n)
diag(adj) <- 0
adj <- (adj + t(adj)) > 0
rownames(adj) <- colnames(adj) <- LETTERS[seq_len(n)]
hg <- build_hypergraph(adj, p = 1, max_size = 3L)
print(hg)
#> Hypergraph: 8 nodes, 51 hyperedges
#> Size distribution:
#>   size_2   : 23
#>   size_3   : 28
#> Method: clique, p = 1.00, include_pairwise = TRUE, max_size = 3
summary(hg)
#> Hypergraph summary
#>   Nodes:         8
#>   Hyperedges:    51
#>   Mean size:     2.55
#>   Max size:      3
#>   node degree
#> 1    A     17
#> 2    B     17
#> 3    C     23
#> 4    D     17
#> 5    E     17
#> 6    F     13
#> 7    G     13
#> 8    H     13
```
