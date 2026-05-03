# Hypergraph eigenvector centralities

Computes one or more eigenvector-style centralities on a
[net_hypergraph](https://mohsaqr.github.io/Nestimate/reference/build_hypergraph.md):
*clique-motif* (CEC), *Z-eigenvector* (ZEC), and *H-eigenvector* (HEC).
Each variant captures influence differently — CEC flattens group
structure via clique expansion, while ZEC and HEC propagate through the
higher-order groups directly.

## Usage

``` r
hypergraph_centrality(
  hg,
  type = c("clique", "Z", "H"),
  max_iter = 1000L,
  tol = 1e-08,
  normalize = TRUE
)
```

## Arguments

- hg:

  A `net_hypergraph` (from
  [`build_hypergraph()`](https://mohsaqr.github.io/Nestimate/reference/build_hypergraph.md)
  or
  [`bipartite_groups()`](https://mohsaqr.github.io/Nestimate/reference/bipartite_groups.md)).

- type:

  Character vector, any subset of `c("clique", "Z", "H")`. Default
  computes all three.

- max_iter:

  Maximum number of power-iteration steps. Default `1000`.

- tol:

  Convergence tolerance on the L1 change between successive iterates.
  Default `1e-8`.

- normalize:

  Logical. If `TRUE` (default), each returned centrality vector is
  L2-normalized to unit norm (compatible with
  [`igraph::eigen_centrality()`](https://r.igraph.org/reference/eigen_centrality.html)'s
  scale for type `"clique"`).

## Value

A named list; one component per requested `type`. Each component is a
named numeric vector of length `hg$n_nodes`.

## Details

**Clique-motif eigenvector centrality (CEC)**: forms the clique-expanded
pairwise graph \\W\\ where \\W\_{ij} = \|\\e : i, j \in e\\\|\\ and
returns the leading eigenvector of \\W\\. Equivalent to running
[`igraph::eigen_centrality()`](https://r.igraph.org/reference/eigen_centrality.html)
on
[`clique_expansion()`](https://mohsaqr.github.io/Nestimate/reference/clique_expansion.md)
output.

**Z-eigenvector centrality (ZEC)**: solves the linear eigen-equation on
the hyperedge tensor, \$\$\lambda\\ x_i \\=\\ \sum\_{e \ni i}\\
\prod\_{j \in e,\\ j \neq i} x_j,\$\$ via power iteration. Works for
hypergraphs with mixed edge sizes.

**H-eigenvector centrality (HEC)**: solves the power-k-1 eigen-equation,
\$\$\lambda\\ x_i^{k-1} \\=\\ \sum\_{e \ni i}\\ \prod\_{j \in e,\\ j
\neq i} x_j.\$\$ For uniform hypergraphs (all hyperedges of size \\k\\),
this is equivalent to normalizing the ZEC update by the geometric-mean
exponent \\1/(k-1)\\. For mixed sizes, the effective exponent is taken
from the largest hyperedge; expect slightly different rankings from ZEC
in the mixed case.

## Note

The `"clique"` (CEC) variant is validated against
[`igraph::eigen_centrality`](https://r.igraph.org/reference/eigen_centrality.html)
(cosine ~ 1). The `"Z"` and `"H"` variants are **(experimental)** —
validated only against a clean-room list-based tensor power iteration
(same operator, different loop structure); no R package exposes tensor
eigenvectors as a primitive for independent comparison.

## References

Benson, A. R. (2019). Three hypergraph eigenvector centralities. *SIAM
Journal on Mathematics of Data Science* 1(2), 293-312. arXiv:1807.09644.

## See also

[`build_hypergraph()`](https://mohsaqr.github.io/Nestimate/reference/build_hypergraph.md),
[`clique_expansion()`](https://mohsaqr.github.io/Nestimate/reference/clique_expansion.md),
[`hypergraph_measures()`](https://mohsaqr.github.io/Nestimate/reference/hypergraph_measures.md).

## Examples

``` r
df <- data.frame(
  player  = c("A", "B", "C", "A", "B", "D", "C", "D", "E"),
  session = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S3", "S3")
)
hg <- bipartite_groups(df, "player", "session")
cent <- hypergraph_centrality(hg)
# Compare rankings across the three variants
do.call(cbind, cent)
#>      clique            Z         H
#> A 0.5345225 6.666667e-01 0.6147176
#> B 0.5345225 6.666667e-01 0.6147176
#> C 0.5345225 3.333333e-01 0.4216259
#> D 0.2672612 3.983626e-09 0.1823130
#> E 0.2672612 3.983626e-09 0.1823130
```
