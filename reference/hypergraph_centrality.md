# Hypergraph eigenvector centralities

Computes one or more eigenvector-style centralities on a
[net_hypergraph](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md):
*clique-motif* (CEC), *Z-eigenvector* (ZEC), and *H-eigenvector* (HEC).
Each variant captures influence differently - CEC flattens group
structure via clique expansion, while ZEC and HEC propagate through the
higher-order groups directly.

## Usage

``` r
hypergraph_centrality(
  hg,
  type = "clique",
  max_iter = 1000L,
  tol = 1e-08,
  normalize = TRUE
)
```

## Arguments

- hg:

  A `net_hypergraph` (from
  [`build_hypergraph()`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md)
  or
  [`bipartite_groups()`](https://pak.dynasite.org/Nestimate/reference/bipartite_groups.md)).

- type:

  Character vector, any subset of `c("clique", "Z", "H")`. Default
  `"clique"`. `"Z"` and `"H"` need a uniform hypergraph (every hyperedge
  with two or more members has the same size).

- max_iter:

  Maximum number of power-iteration steps of the `"Z"` and `"H"`
  iterations. Default `1000`. `"clique"` is solved directly.

- tol:

  Convergence tolerance on the L1 change between successive iterates of
  `"Z"` and `"H"`. Default `1e-8`. An iteration that does not reach it
  within `max_iter` steps warns with class `nestimate_no_converge`.

- normalize:

  Logical. If `TRUE` (default), each returned centrality vector is
  L2-normalized to unit norm (compatible with
  [`igraph::eigen_centrality()`](https://r.igraph.org/reference/eigen_centrality.html)'s
  scale for type `"clique"`). If `FALSE`, the vector is rescaled so its
  largest absolute entry is 1.

## Value

A named list, one component per requested `type`, in the order given by
`type`. Each component is a numeric vector of length `hg$n_nodes` named
by `hg$nodes`. A hypergraph with no hyperedges yields all-zero vectors.
Requesting `"Z"` or `"H"` on a non-uniform hypergraph raises an error of
class `nestimate_nonuniform_hypergraph`.

## Details

**Clique-motif eigenvector centrality (CEC)**: forms the clique-expanded
pairwise graph \\W\\ where \\W\_{ij} = \|\\e : i, j \in e\\\|\\ and
returns the leading eigenvector of \\W\\. Equivalent to running
[`igraph::eigen_centrality()`](https://r.igraph.org/reference/eigen_centrality.html)
on
[`clique_expansion()`](https://pak.dynasite.org/Nestimate/reference/clique_expansion.md)
output. The vector is solved directly from a symmetric
eigendecomposition of each connected component (Perron-Frobenius makes
the leading eigenvalue of a connected component simple), not by power
iteration, which oscillates when the clique expansion is bipartite (for
example a star of pairs). When several components share the largest
eigenvalue, the result is the projection of the uniform vector onto
their joint leading eigenspace; components with a smaller leading
eigenvalue score zero.

**Z-eigenvector centrality (ZEC)**: solves the eigen-equation on the
hyperedge tensor, \$\$\lambda\\ x_i \\=\\ \sum\_{e \ni i}\\ \prod\_{j
\in e,\\ j \neq i} x_j,\$\$ by shifted power iteration. Tensor products
and normalization are evaluated in log space to retain finite,
normalized centralities for large hyperedges.

**H-eigenvector centrality (HEC)**: solves \$\$\lambda\\ x_i^{k-1} \\=\\
\sum\_{e \ni i}\\ \prod\_{j \in e,\\ j \neq i} x_j.\$\$

**Uniform hypergraphs only.** Benson (2019) defines ZEC and HEC for
\\k\\-uniform hypergraphs. With mixed sizes the right-hand side mixes
polynomial degrees, so the equation has no scale-free solution: the
vector satisfies it only at the normalisation chosen inside the
iteration. Hyperedges with one member contribute nothing and are ignored
by the uniformity check. For mixed sizes use `"clique"`, or keep one
size with
[`build_hypergraph()`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md)
(`include_pairwise = FALSE`).

## Note

The `"clique"` (CEC) variant is validated against
[`igraph::eigen_centrality`](https://r.igraph.org/reference/eigen_centrality.html)
(cosine ~ 1). The `"Z"` and `"H"` variants are **(experimental)** -
validated only against a clean-room list-based tensor power iteration
(same operator, different loop structure); no R package exposes tensor
eigenvectors as a primitive for independent comparison.

## References

Benson, A. R. (2019). Three hypergraph eigenvector centralities. *SIAM
Journal on Mathematics of Data Science* 1(2), 293-312. arXiv:1807.09644.

## See also

[`build_hypergraph()`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md),
[`clique_expansion()`](https://pak.dynasite.org/Nestimate/reference/clique_expansion.md),
[`hypergraph_measures()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_measures.md).

## Examples

``` r
df <- data.frame(
  player  = c("A", "B", "C", "A", "B", "D", "C", "D", "E"),
  session = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S3", "S3")
)
hg <- bipartite_groups(df, "player", "session")
hypergraph_centrality(hg)
#> $clique
#>         A         B         C         D         E 
#> 0.5345225 0.5345225 0.5345225 0.2672612 0.2672612 
#> 

# Z and H need a uniform hypergraph: here every session has three players
df3 <- data.frame(
  player  = c("A", "B", "C", "A", "B", "D", "C", "D", "E"),
  session = c("S1", "S1", "S1", "S2", "S2", "S2", "S3", "S3", "S3")
)
hypergraph_centrality(bipartite_groups(df3, "player", "session"),
                      type = c("clique", "Z", "H"))
#> $clique
#>         A         B         C         D         E 
#> 0.5127870 0.5127870 0.4559848 0.4559848 0.2413603 
#> 
#> $Z
#>         A         B         C         D         E 
#> 0.5345225 0.5345225 0.4364358 0.4364358 0.2182179 
#> 
#> $H
#>         A         B         C         D         E 
#> 0.4867656 0.4867656 0.4556438 0.4556438 0.3330103 
#> 
```
