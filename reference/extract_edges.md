# Extract Edge List with Weights

`extract_edges()` returns the edges of a network as a data frame, one
row per edge with its weight, read from the weight matrix (rows are the
`from` states, columns the `to` states). It reads every input
[`extract_transition_matrix()`](https://pak.dynasite.org/Nestimate/reference/extract_transition_matrix.md)
reads; an `mcml` gives one edge list per layer.

An edge is a cell whose absolute weight is above `threshold`, so with
the default `threshold = 0` every non-zero cell is an edge, negative
weights (in correlation, partial-correlation and mlVAR networks)
included, and zero cells are not. A directed network lists every ordered
pair; an undirected one lists each edge once, with `from` the earlier
state in the matrix order. A network's own `directed` field decides
which; a `tna` model is directed unless it is a co-occurrence model, and
a bare matrix is undirected when it is symmetric.

## Usage

``` r
extract_edges(model, threshold = 0, include_self = FALSE, sort_by = "weight")
```

## Arguments

- model:

  The network: a `netobject`, a `tna`-family model, an `mcml` or a
  weight matrix, as for
  [`extract_transition_matrix()`](https://pak.dynasite.org/Nestimate/reference/extract_transition_matrix.md).

- threshold:

  Cells whose absolute weight is at most `threshold` are dropped. A
  single finite number of at least 0, otherwise an error is raised
  (class `nestimate_bad_input` for a negative value); default `0`.

- include_self:

  Whether to keep self-loops (a state followed by itself), `TRUE` or
  `FALSE`. Default `FALSE`.

- sort_by:

  The row order: `"weight"` (default) by decreasing weight, `"from"` by
  source then target state, `"to"` by target then source state, or
  `NULL` for the order of the weight matrix (column by column). Any
  other value raises an error.

## Value

For a single network: a data frame with one row per edge and the
columns:

- from:

  The source state.

- to:

  The target state.

- weight:

  The edge weight (for a `"relative"` network, the transition
  probability).

For an `mcml`: a named list of such data frames, `macro` first and then
one per cluster.

## Details

An edge list is the form most other network tools (igraph, export to
Gephi or a spreadsheet) take.

## Note

The idiographic package exports a function with the same name. With both
packages attached, the one attached last masks the other; write
`Nestimate::extract_edges()` to call this one.

## See also

[`extract_transition_matrix`](https://pak.dynasite.org/Nestimate/reference/extract_transition_matrix.md)
for the full matrix,
[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
for network estimation.

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","A"), V2 = c("B","A","C"), V3 = c("A","C","B"))
net <- build_network(seqs, method = "relative")
edges <- extract_edges(net, threshold = 0.05)
head(edges)
#>   from to    weight
#> 1    B  A 1.0000000
#> 2    C  B 1.0000000
#> 3    A  C 0.6666667
#> 4    A  B 0.3333333
```
