# Clique expansion of a hypergraph

Projects a
[net_hypergraph](https://saqr.me/Nestimate/reference/build_hypergraph.md)
to a standard pairwise
[netobject](https://saqr.me/Nestimate/reference/build_network.md) (the
*clique expansion* - also called the "downgrade" of a hypergraph to a
dyadic graph). Each hyperedge of size k contributes 1 (or its weight) to
every pair of its members. The resulting edge weight `W[i, j]` equals
the number of hyperedges containing both `i` and `j` (binary incidence)
or the sum of incidence products (weighted incidence).

## Usage

``` r
clique_expansion(hg, weighted = TRUE)
```

## Arguments

- hg:

  A `net_hypergraph` object as returned by
  [`build_hypergraph()`](https://saqr.me/Nestimate/reference/build_hypergraph.md)
  or
  [`bipartite_groups()`](https://saqr.me/Nestimate/reference/bipartite_groups.md).

- weighted:

  Logical. If `TRUE` (default), use the hypergraph's incidence values
  directly (so weighted hypergraphs from
  [`bipartite_groups()`](https://saqr.me/Nestimate/reference/bipartite_groups.md)
  produce weighted projections). If `FALSE`, binarise the incidence
  first so `W[i, j]` is just the count of shared hyperedges.

## Value

A `netobject` (also `cograph_network`) with
`method = "clique_expansion"`, undirected, with weighted symmetric
adjacency `W = incidence %*% t(incidence)` and zero diagonal. The
standard `netobject` fields are present (`$weights`, `$nodes`,
`$edges` - one row per non-zero upper-triangle cell with integer
`from`/`to` node indices and `weight` - `$n_nodes`, `$n_edges`,
`$meta`); `$params` records `source`, `weighted`, `n_hyperedges` and
`hypergraph_size_distribution`.

## Details

The clique expansion is the standard "loss-y but lossless-on-pairwise"
projection: it preserves *which pairs co-occurred* and *how often* but
discards the higher-order grouping. Comparing `clique_expansion(hg)` to
a directly-estimated pairwise network (e.g. via
[`cooccurrence()`](https://saqr.me/Nestimate/reference/cooccurrence.md)
on the same data) quantifies how much information was carried by the
hyperedge structure.

Computed in one BLAS call via `tcrossprod(incidence)`; runs in
`O(n_nodes^2 * n_hyperedges)` time, fast for typical sizes.

Closes the I/O cycle: event data -\>
[`bipartite_groups()`](https://saqr.me/Nestimate/reference/bipartite_groups.md)
-\> `clique_expansion()` -\> any function that accepts a `netobject`
(centrality, bootstrap, clustering, plotting via cograph).

## Note

(experimental) Validated against `tcrossprod(incidence)` with zero
diagonal. No external R package exposes clique expansion as a primitive;
the implementation is a direct one-line restatement of the definition.

## References

Tian, H., & Zafarani, R. (2024). Higher-order networks representation
and learning: A survey. *ACM SIGKDD Explorations Newsletter* 26(1),
1-18.

## See also

[`build_hypergraph()`](https://saqr.me/Nestimate/reference/build_hypergraph.md),
[`bipartite_groups()`](https://saqr.me/Nestimate/reference/bipartite_groups.md),
[`build_network()`](https://saqr.me/Nestimate/reference/build_network.md).

## Examples

``` r
df <- data.frame(
  player  = c("A", "B", "C", "A", "B", "D", "C", "D", "E"),
  session = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S3", "S3")
)
hg  <- bipartite_groups(df, player = "player", group = "session")
net <- clique_expansion(hg)
extract_edges(net, threshold = 1)
#>    from to weight
#> 1     B  A      2
#> 2     A  B      2
#> 3     C  A      1
#> 4     C  B      1
#> 5     A  C      1
#> 6     B  C      1
#> 7     D  C      1
#> 8     E  C      1
#> 9     C  D      1
#> 10    E  D      1
#> 11    C  E      1
#> 12    D  E      1
```
