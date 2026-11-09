# Structural measures for a hypergraph

`hypergraph_measures()` computes the structural statistics of a
hypergraph at three levels: for each node (how many hyperedges it
belongs to, how large they are, and how many it shares with each other
node), for each hyperedge and pair of hyperedges (their sizes and how
much they overlap), and for the hypergraph as a whole (density, mean
hyperedge size, the intersection profile). It reads the `net_hypergraph`
objects that
[`build_hypergraph()`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md)
and
[`bipartite_groups()`](https://pak.dynasite.org/Nestimate/reference/bipartite_groups.md)
return.

The measures describe the hypergraph as observed: they carry no
uncertainty and test nothing. The intersection profile shows whether
hyperedges mostly share no or few nodes or share substantial cores (Do
et al., 2020; Lee, Choe & Shin, 2021). To rank nodes by their position
in the hypergraph, use
[`hypergraph_centrality()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_centrality.md);
to read it as an ordinary pairwise network, project it with
[`clique_expansion()`](https://pak.dynasite.org/Nestimate/reference/clique_expansion.md).

## Usage

``` r
hypergraph_measures(hg)

# S3 method for class 'hypergraph_measures'
print(x, ...)
```

## Arguments

- hg:

  A `net_hypergraph`, from
  [`build_hypergraph()`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md)
  or
  [`bipartite_groups()`](https://pak.dynasite.org/Nestimate/reference/bipartite_groups.md).
  Any other object raises an error.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) method: a
  `hypergraph_measures` object.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html) method: not
  used.

## Value

A `hypergraph_measures` object, a named list;
[`print()`](https://rdrr.io/r/base/print.html) shows the density, mean
hyperedge size and pairwise participation, the minimum, median and
maximum of hyperdegree, node strength and hyperedge size, and the
intersection profile, and returns its input invisibly. Its components,
node-level first:

- **Node-level (one value per node)**:

- `hyperdegree`:

  The number of hyperedges containing each node.

- `node_strength`:

  Total participation: for node \\i\\, \\\sum\_{e \ni i} \|e\|\\. A node
  in many large hyperedges has high strength.

- `max_edge_size`:

  The size of the largest hyperedge containing each node (0 for a node
  in none).

- `co_degree`:

  An `n_nodes` x `n_nodes` matrix:
  `co_degree[i, j] = |{e : i, j in e}|`, the number of hyperedges
  containing both nodes \\i\\ and \\j\\. The diagonal is zero.

- **Hyperedge-level (one value per hyperedge, or m x m)**:

- `edge_sizes`:

  The hyperedge sizes \\\|e\|\\.

- `edge_pairwise_overlap`:

  An `m` x `m` matrix of `|e_i intersect e_j|`. The diagonal is zero.

- `overlap_coefficient`:

  An `m` x `m` matrix of `|e_i intersect e_j| / min(|e_i|, |e_j|)`: how
  much of the smaller hyperedge lies in the larger. The diagonal is
  zero.

- `jaccard`:

  An `m` x `m` matrix of the symmetric overlap index
  `|e_i intersect e_j| / |e_i union e_j|`. The diagonal is zero.

- **Global (scalars)**:

- `density`:

  For a hypergraph whose hyperedges all have size `k`:
  `m / choose(n, k)`. For mixed sizes: `sum(|e|) / (n * m)`, the mean
  fraction of nodes per hyperedge.

- `avg_edge_size`:

  The mean of `edge_sizes`.

- `size_distribution`:

  The count of hyperedges of each size, passed through from `hg`.

- `intersection_profile`:

  The number of hyperedge pairs at each intersection size that occurs,
  named `overlap_<size>` (sizes with no pairs are left out); empty with
  fewer than two hyperedges.

- `pairwise_participation`:

  The fraction of node pairs that appear together in at least one
  hyperedge (0 with fewer than two nodes).

- `n_nodes`, `n_hyperedges`:

  The numbers of nodes and hyperedges.

## Details

All measures are matrix operations on the binary incidence matrix \\B =
(b\_{ij})\\, where \\b\_{ij} = 1\\ if node \\i\\ is in hyperedge \\j\\;
a weighted incidence is binarised first:

- `hyperdegree = rowSums(B)`, `edge_sizes = colSums(B)`

- `co_degree = tcrossprod(B)` (with zero diagonal)

- `edge_pairwise_overlap = crossprod(B)` (with zero diagonal)

- `overlap_coefficient[i, j] = overlap[i, j] / min(edge_sizes[i], edge_sizes[j])`

- `jaccard[i, j] = overlap[i, j] / (edge_sizes[i] + edge_sizes[j] - overlap[i, j])`

A ratio whose denominator is 0 is set to 0.

A hypergraph without hyperedges (`n_hyperedges == 0`) gives zeros for
the node measures, empty hyperedge matrices, a density and pairwise
participation of 0, and an `avg_edge_size` of `NA`.

## References

Lee, G., Choe, M., & Shin, K. (2021). How do hyperedges overlap in
real-world hypergraphs? Patterns, measures, and generators. *Proceedings
of the Web Conference 2021*, 3396-3407.

Do, M. T., Yoon, S., Hooi, B., & Shin, K. (2020). Structural patterns
and generative models of real-world hypergraphs. arXiv:2006.07060.

## See also

[`build_hypergraph()`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md),
[`bipartite_groups()`](https://pak.dynasite.org/Nestimate/reference/bipartite_groups.md),
[`hypergraph_centrality()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_centrality.md),
[`clique_expansion()`](https://pak.dynasite.org/Nestimate/reference/clique_expansion.md).

## Examples

``` r
df <- data.frame(
  player  = c("A", "B", "C", "A", "B", "D", "C", "D"),
  session = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S3")
)
hg <- bipartite_groups(df, "player", "session")
m  <- hypergraph_measures(hg)
print(m)   # node-level, hyperedge-level and global measures
#> Hypergraph measures: 4 nodes, 3 hyperedges
#>   Density:                0.5833
#>   Mean edge size:         2.33
#>   Pairwise participation: 0.6667
#> 
#>   Hyperdegree:    min=1  med=2.0  max=2
#>   Node strength:  min=2.0  med=5.0  max=5.0
#>   Edge sizes:     min=2  med=2.0  max=3
#> 
#> Intersection profile (pair counts by overlap size):
#>   overlap_0   : 1
#>   overlap_1   : 1
#>   overlap_2   : 1
```
