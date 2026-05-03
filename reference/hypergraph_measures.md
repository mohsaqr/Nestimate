# Structural measures for a hypergraph

Computes a comprehensive structural-statistics suite for a
[net_hypergraph](https://mohsaqr.github.io/Nestimate/reference/build_hypergraph.md):
node-level, hyperedge-level, and global measures. All measures are
derived in a few BLAS calls on the incidence matrix.

## Usage

``` r
hypergraph_measures(hg)

# S3 method for class 'hypergraph_measures'
print(x, ...)
```

## Arguments

- hg:

  A `net_hypergraph` (from
  [`build_hypergraph()`](https://mohsaqr.github.io/Nestimate/reference/build_hypergraph.md)
  or
  [`bipartite_groups()`](https://mohsaqr.github.io/Nestimate/reference/bipartite_groups.md)).

- x:

  A `hypergraph_measures` object.

- ...:

  Additional arguments (ignored).

## Value

An object of class `hypergraph_measures` (a named list) with components:

- **Node-level (length `n_nodes`)**:

- `hyperdegree`:

  Number of hyperedges containing each node.

- `node_strength`:

  Total participation: for node \\i\\, \\\sum\_{e \ni i} \|e\|\\. A node
  in many large hyperedges has high strength.

- `max_edge_size`:

  Size of the largest hyperedge containing each node.

- `co_degree`:

  `n_nodes` x `n_nodes` matrix: `co_degree[i, j] = |{e : i, j in e}|`
  (number of hyperedges co-containing nodes \\i\\ and \\j\\). Diagonal
  is zero.

- **Hyperedge-level (length `n_hyperedges` or m x m)**:

- `edge_sizes`:

  Hyperedge sizes \\\|e\|\\.

- `edge_pairwise_overlap`:

  `m` x `m` matrix: `|e_i intersect e_j|`. Diagonal is zero.

- `overlap_coefficient`:

  `m` x `m`: `|e_i and e_j| / min(|e_i|, |e_j|)`. Measures how much the
  smaller hyperedge is contained in the larger.

- `jaccard`:

  `m` x `m`: symmetric overlap index `|e_i and e_j| / |e_i union e_j|`.

- **Global (scalars)**:

- `density`:

  For `k`-uniform hypergraphs: `m / choose(n, k)`. For mixed sizes:
  `sum(|e|) / (n * m)` (mean fraction of nodes per hyperedge).

- `avg_edge_size`:

  Mean of `edge_sizes`.

- `size_distribution`:

  Tabulation of hyperedge sizes (passed through from `hg`).

- `intersection_profile`:

  Distribution of pairwise hyperedge intersection sizes — useful for
  spotting whether hyperedges overlap mostly trivially or share
  substantial cores (Do et al. 2020).

- `pairwise_participation`:

  Fraction of node pairs co-appearing in at least one hyperedge.

- `n_nodes`, `n_hyperedges`:

  Convenience scalars.

The input `x` invisibly.

## Details

All measures are computed via standard matrix operations on the binary
incidence \\B = (b\_{ij})\\ where \\b\_{ij} = 1\\ iff node \\i\\ is in
hyperedge \\j\\:

- `hyperdegree = rowSums(B)`, `edge_sizes = colSums(B)`

- `co_degree = tcrossprod(B)` (with zero diagonal)

- `edge_pairwise_overlap = crossprod(B)` (with zero diagonal)

- `overlap_coefficient[i, j] = overlap[i, j] / min(edge_sizes[i], edge_sizes[j])`

- `jaccard[i, j] = overlap[i, j] / (edge_sizes[i] + edge_sizes[j] - overlap[i, j])`

Empty hypergraph (`n_hyperedges == 0`) returns trivial zeros and empty
matrices.

## References

Lee, G., Choe, M., & Shin, K. (2024). A survey on hypergraph
representation, learning and mining. *Data Mining & Knowledge Discovery*
37, 1-39.

Do, M. T., Yoon, S., Hooi, B., & Shin, K. (2020). Structural patterns
and generative models of real-world hypergraphs. arXiv:2006.07060.

## See also

[`build_hypergraph()`](https://mohsaqr.github.io/Nestimate/reference/build_hypergraph.md),
[`bipartite_groups()`](https://mohsaqr.github.io/Nestimate/reference/bipartite_groups.md),
[`clique_expansion()`](https://mohsaqr.github.io/Nestimate/reference/clique_expansion.md).

## Examples

``` r
df <- data.frame(
  player  = c("A", "B", "C", "A", "B", "D", "C", "D"),
  session = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S3")
)
hg <- bipartite_groups(df, "player", "session")
m  <- hypergraph_measures(hg)
print(m)
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
m$hyperdegree         # how many sessions each player joined
#> A B C D 
#> 2 2 2 1 
m$co_degree           # pairwise co-membership counts
#>   A B C D
#> A 0 2 1 0
#> B 2 0 1 0
#> C 1 1 0 1
#> D 0 0 1 0
m$jaccard             # symmetric overlap between sessions
#>           S1        S2   S3
#> S1 0.0000000 0.6666667 0.25
#> S2 0.6666667 0.0000000 0.00
#> S3 0.2500000 0.0000000 0.00
```
