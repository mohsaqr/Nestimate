# Normalized hypergraph Laplacian

Computes the normalized Laplacian of a hypergraph, either the classic
Zhou-Huang-Scholkopf form on the binary incidence pattern
(`type = "zhou"`) or the random-walk form with edge-dependent vertex
weights (`type = "random_walk"`), in which the weighted incidence cells
(e.g. the summed weights produced by
[`bipartite_groups()`](https://saqr.me/Nestimate/reference/bipartite_groups.md))
determine where a random walker lands inside a hyperedge, and the
resulting non-reversible walk is symmetrized through its stationary
distribution (Chung 2005). Both Laplacians are symmetric positive
semi-definite with eigenvalues in `[0, 2]`; for a binary incidence with
unit hyperedge weights they coincide.

## Usage

``` r
hypergraph_laplacian(hg, type = c("zhou", "random_walk"), edge_weights = NULL)
```

## Arguments

- hg:

  A `net_hypergraph` from
  [`build_hypergraph()`](https://saqr.me/Nestimate/reference/build_hypergraph.md)
  or
  [`bipartite_groups()`](https://saqr.me/Nestimate/reference/bipartite_groups.md).
  Must be connected and have at least one hyperedge.

- type:

  Character. `"zhou"` (default) for the Zhou et al. (2006) normalized
  Laplacian on the binary incidence pattern, or `"random_walk"` for the
  Hayashi et al. (2020) EDVW random-walk Laplacian on the weighted
  incidence.

- edge_weights:

  Numeric vector of positive hyperedge weights (length
  `hg$n_hyperedges`), or `NULL` for the type-specific default: unit
  weights for `"zhou"`; for `"random_walk"` the Hayashi et al.
  heuristic - the population standard deviation of each hyperedge's
  (non-zero) vertex weights plus one - which reduces to unit weights on
  a binary incidence.

## Value

A symmetric `n_nodes` x `n_nodes` numeric matrix (node names as
dimnames) with attributes `type` (the Laplacian type), `pi` (named
stationary distribution of the underlying random walk) and
`edge_weights` (the hyperedge weights actually used).

## References

Zhou, D., Huang, J., & Scholkopf, B. (2006). Learning with hypergraphs:
Clustering, classification, and embedding. *NeurIPS 19*.

Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
random walks, Laplacians, and clustering. *CIKM 2020*, 495-504.
[doi:10.1145/3340531.3412034](https://doi.org/10.1145/3340531.3412034)

Chung, F. (2005). Laplacians and the Cheeger inequality for directed
graphs. *Annals of Combinatorics*, 9(1), 1-19.

## Examples

``` r
events <- data.frame(
  person = c("a", "b", "c", "a", "b", "d", "c", "d", "e", "e", "a"),
  meeting = c("m1", "m1", "m1", "m2", "m2", "m2", "m3", "m3", "m3",
              "m4", "m4"),
  hours = c(2, 1, 1, 3, 2, 1, 2, 2, 4, 1, 1)
)
hg <- bipartite_groups(events, player = "person", group = "meeting",
                       weight = "hours")
L <- hypergraph_laplacian(hg, type = "random_walk")
range(eigen(L, symmetric = TRUE, only.values = TRUE)$values)
#> [1] -2.899699e-16  1.001279e+00
```
