# Hypergraph from bipartite group / event data

Constructs a
[net_hypergraph](https://mohsaqr.github.io/Nestimate/reference/build_hypergraph.md)
from long-format event data in which each row records a `player`
participating in a `group` (a session, team, project, transaction, or
any group context). Each unique group becomes one hyperedge spanning the
players that appeared in it. Optional `weight` column produces a
weighted incidence matrix.

## Usage

``` r
bipartite_groups(data, player, group, weight = NULL)
```

## Arguments

- data:

  Data frame in long format. Must contain `player` and `group` columns;
  optionally a `weight` column.

- player:

  Character. Name of the column whose values become the hypergraph's
  nodes (players, participants, actors).

- group:

  Character. Name of the column whose values become the hypergraph's
  hyperedges (groups, sessions, teams).

- weight:

  Character or `NULL`. If supplied, the column is summed per
  `(player, group)` pair to produce a weighted incidence matrix. Default
  `NULL` produces a 0/1 binary incidence matrix.

## Value

A `net_hypergraph` object with the same structure produced by
[`build_hypergraph()`](https://mohsaqr.github.io/Nestimate/reference/build_hypergraph.md)
(`hyperedges`, `incidence`, `nodes`, `n_nodes`, `n_hyperedges`,
`size_distribution`, `params`). The `params` list records
`source = "bipartite_groups"` and the original column names.

## Details

The bipartite representation preserves the full group structure without
projecting to a pairwise network. A group of three players A, B, C
produces a single 3-hyperedge containing all three, not three pairwise
edges AB, AC, BC. This avoids information loss when group interactions
are the primary unit of analysis (Perc et al. 2013).

Unlike
[`build_hypergraph()`](https://mohsaqr.github.io/Nestimate/reference/build_hypergraph.md)
(which derives hyperedges from a network's clique structure),
`bipartite_groups()` takes group memberships directly. The two functions
are complementary:

- `bipartite_groups()` — when group membership is observed (sessions,
  transactions, co-authorships).

- [`build_hypergraph()`](https://mohsaqr.github.io/Nestimate/reference/build_hypergraph.md)
  — when only pairwise interactions are observed and triadic structure
  must be inferred from triangles.

Rows with `NA` in either the `player` or `group` column (or, when
supplied, the `weight` column) are dropped silently.

## Note

(experimental) Validated against a hand-computed
[`table()`](https://rdrr.io/r/base/table.html) incidence reference only;
no independent R package exposes the long-format-to-binary-incidence
primitive, because the operation is definitionally
[`table()`](https://rdrr.io/r/base/table.html). The code path is a
direct one-to-one restatement of its definition.

## References

Perc, M., Gomez-Gardenes, J., Szolnoki, A., Floria, L. M., & Moreno, Y.
(2013). Evolutionary dynamics of group interactions on structured
populations: a review. *Journal of the Royal Society Interface* 10(80),
20120997.
[doi:10.1098/rsif.2012.0997](https://doi.org/10.1098/rsif.2012.0997)

## See also

[`build_hypergraph()`](https://mohsaqr.github.io/Nestimate/reference/build_hypergraph.md)
for the clique-based constructor.

## Examples

``` r
df <- data.frame(
  player = c("Alice", "Bob", "Carol", "Alice", "Bob",
             "Dave", "Carol", "Dave", "Eve"),
  session = c("S1", "S1", "S1", "S2", "S2",
              "S3", "S3", "S3", "S3")
)
hg <- bipartite_groups(df, player = "player", group = "session")
print(hg)
#> Hypergraph: 5 nodes, 3 hyperedges
#> Size distribution:
#>   size_2   : 1
#>   size_3   : 2
summary(hg)
#> Hypergraph summary
#>   Nodes:         5
#>   Hyperedges:    3
#>   Mean size:     2.67
#>   Max size:      3
#>    node degree
#> 1 Alice      2
#> 2   Bob      2
#> 3 Carol      2
#> 4  Dave      1
#> 5   Eve      1
```
