# Hypergraph eigenvector centralities

`hypergraph_centrality()` scores how central each node of a hypergraph
is, with up to three eigenvector centralities from Benson (2019):
clique-motif (CEC, `"clique"`), Z-eigenvector (ZEC, `"Z"`) and
H-eigenvector (HEC, `"H"`). CEC flattens the groups into pairs by clique
expansion and takes the eigenvector centrality of that pairwise graph.
ZEC and HEC work on the hyperedges directly: a node's score grows with
the product of the scores of the other members of each hyperedge it
belongs to. It reads a `net_hypergraph` from
[`build_hypergraph()`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md)
or
[`bipartite_groups()`](https://pak.dynasite.org/Nestimate/reference/bipartite_groups.md).

The variants can rank the nodes differently, since CEC counts shared
memberships pair by pair while ZEC and HEC weigh each group as a whole.
ZEC and HEC are defined for uniform hypergraphs only, where every
hyperedge of two or more members has the same size. On a disconnected
hypergraph, compare scores within one component: CEC gives zero to every
component whose leading eigenvalue is below the largest, and the
Z-eigenvector is not unique.
[`hypergraph_measures()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_measures.md)
describes the structure the scores come from.

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

  The hypergraph: a `net_hypergraph` from
  [`build_hypergraph()`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md)
  or
  [`bipartite_groups()`](https://pak.dynasite.org/Nestimate/reference/bipartite_groups.md).

- type:

  The centralities to compute, any of `"clique"`, `"Z"` and `"H"`,
  returned in the order given. Default `"clique"`. `"Z"` and `"H"` need
  a uniform hypergraph (every hyperedge with two or more members has the
  same size); on any other they raise an error of class
  `nestimate_nonuniform_hypergraph`.

- max_iter:

  The largest number of power-iteration steps for `"Z"` and `"H"`.
  Default `1000`. `"clique"` is solved directly and ignores it.

- tol:

  The convergence tolerance for `"Z"` and `"H"`: the iteration stops
  once the L1 change between successive iterates is below it. Default
  `1e-8`. An iteration that does not get there within `max_iter` steps,
  or that stops without solving the eigen-equation, warns with class
  `nestimate_no_converge`.

- normalize:

  How each vector is scaled. `TRUE` (default) gives it unit Euclidean
  (L2) norm; `FALSE` divides it by its largest absolute entry, so the
  top node scores 1, the scale of
  [`igraph::eigen_centrality()`](https://r.igraph.org/reference/eigen_centrality.html).

## Value

A named list with one numeric vector per requested `type`, in the order
given, each with one value per node of `hg`, named by the nodes. A
hypergraph with no hyperedges, or with only one-member hyperedges, gives
all-zero vectors, and one without nodes empty vectors.

## Details

**Clique-motif eigenvector centrality (CEC)** is the leading eigenvector
of the clique-expanded pairwise graph \\W\\ with \\W\_{ij} = \|\\e : i,
j \in e\\\|\\, the number of hyperedges two nodes share. It reads
membership only, so incidence weights are ignored: it equals
[`igraph::eigen_centrality()`](https://r.igraph.org/reference/eigen_centrality.html)
on `clique_expansion(hg, weighted = FALSE)`, and on
[`clique_expansion()`](https://pak.dynasite.org/Nestimate/reference/clique_expansion.md)'s
default output only for an unweighted hypergraph. The vector is solved
directly from a symmetric eigendecomposition of each connected component
(Perron-Frobenius makes the leading eigenvalue of a connected component
simple), not by power iteration, which oscillates when the clique
expansion is bipartite (for example a star of pairs). When several
components share the largest eigenvalue, the result is the projection of
the uniform vector onto their joint leading eigenspace; components with
a smaller leading eigenvalue score zero.

**Z-eigenvector centrality (ZEC)** solves the eigen-equation on the
hyperedge tensor, \$\$\lambda\\ x_i \\=\\ \sum\_{e \ni i}\\ \prod\_{j
\in e,\\ j \neq i} x_j,\$\$ by a shifted power iteration (a variant of
SS-HOPM, Kolda & Mayo, 2011, that normalises the tensor step before
adding the shift). Tensor products and normalisation are evaluated in
log space, so the centralities stay finite and normalised for large
hyperedges. Besides the step tolerance, the iteration checks that the
eigen-equation residual is small (relative residual at most
`max(sqrt(tol), 1e-6)`); otherwise it warns.

**H-eigenvector centrality (HEC)** solves, by the same iteration,
\$\$\lambda\\ x_i^{k-1} \\=\\ \sum\_{e \ni i}\\ \prod\_{j \in e,\\ j
\neq i} x_j.\$\$

**Uniform hypergraphs only.** Benson (2019) defines ZEC and HEC for
\\k\\-uniform hypergraphs. With mixed sizes the right-hand side mixes
polynomial degrees, so the equation has no scale-free solution: the
vector satisfies it only at the normalisation chosen inside the
iteration. Hyperedges with one member contribute nothing and are ignored
by the uniformity check. For mixed sizes use `"clique"`, or keep one
size with
[`build_hypergraph()`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md)
(`include_pairwise = FALSE`).

**Disconnected hypergraphs.** The Z-eigenvector is not unique when the
hypergraph has several components: different iterations can solve the
equation with different vectors, so compare `"Z"` values within one
connected hypergraph only.

## Note

The `"clique"` (CEC) variant is validated against
[`igraph::eigen_centrality()`](https://r.igraph.org/reference/eigen_centrality.html)
(cosine ~ 1). The `"Z"` and `"H"` variants are **experimental**: they
are validated only against a clean-room list-based tensor power
iteration (the same operator with a different loop structure), as no R
package exposes tensor eigenvectors for an independent comparison.

## References

Benson, A. R. (2019). Three hypergraph eigenvector centralities. *SIAM
Journal on Mathematics of Data Science* 1(2), 293-312. arXiv:1807.09644.

Kolda, T. G., & Mayo, J. R. (2011). Shifted power method for computing
tensor eigenpairs. *SIAM Journal on Matrix Analysis and Applications*
32(4), 1095-1124.
[doi:10.1137/100801482](https://doi.org/10.1137/100801482)

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
