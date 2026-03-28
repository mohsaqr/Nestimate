# Extract Edge List with Weights

Extract an edge list from a TNA model, representing the network as a
data frame of from-to-weight tuples.

## Usage

``` r
extract_edges(model, threshold = 0, include_self = FALSE, sort_by = "weight")
```

## Arguments

- model:

  A TNA model object or a matrix of weights.

- threshold:

  Numeric. Minimum weight to include an edge. Default: 0.

- include_self:

  Logical. Whether to include self-loops. Default: FALSE.

- sort_by:

  Character. Column to sort by: "weight" (descending), "from", "to", or
  NULL for no sorting. Default: "weight".

## Value

A data frame with columns:

- from:

  Source state name.

- to:

  Target state name.

- weight:

  Edge weight (transition probability).

## Details

This function converts the transition matrix into an edge list format,
which is useful for visualization, analysis with igraph, or export to
other network tools.

## See also

[`extract_transition_matrix`](https://mohsaqr.github.io/Nestimate/reference/extract_transition_matrix.md)
for the full matrix,
[`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md)
for network estimation.

## Examples

``` r
# \donttest{
seqs <- data.frame(V1 = c("A","B","A"), V2 = c("B","A","C"), V3 = c("A","C","B"))
net <- build_network(seqs, method = "relative")
edges <- extract_edges(net, threshold = 0.05)
head(edges)
#>   from to    weight
#> 1    B  A 1.0000000
#> 2    C  B 1.0000000
#> 3    A  C 0.6666667
#> 4    A  B 0.3333333
# }
```
