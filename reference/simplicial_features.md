# Tidy Topological Features for One or Many Networks

Builds a simplicial complex per network and returns its topological
summaries as a tidy `data.frame` – one row per network, one column per
feature – ready to use as regression predictors or to join onto
unit-level outcomes.

## Usage

``` r
simplicial_features(
  x,
  threshold = 0,
  max_dim = 4L,
  normalize = FALSE,
  type = "clique"
)
```

## Arguments

- x:

  A `netobject`, `netobject_group`, `mcml`, a square weight matrix, or a
  named list of any of these. A group or list yields one row per member,
  an `mcml` one row per cluster, and a single network one row.

- threshold:

  Minimum absolute edge weight for an edge to exist (passed to
  [`build_simplicial`](https://saqr.me/Nestimate/reference/build_simplicial.md)).
  Topology is a step function of this value, so a single threshold is a
  choice, not a result – pass a vector to sweep it and get one row per
  network per threshold.

- max_dim:

  Maximum simplex dimension retained. Default `4`.

- normalize:

  Divide simplex counts by the number of nodes, so networks of different
  size are comparable. Default `FALSE`.

- type:

  Complex type passed to
  [`build_simplicial`](https://saqr.me/Nestimate/reference/build_simplicial.md).
  Default `"clique"`.

## Value

A `data.frame` with one row per network (per threshold), and columns
`network`, `threshold`, `n_nodes`, `n_edges`, `b0`, `b1` (Betti
numbers), `euler`, `max_q`, `d1` ... `d<max_dim>` (simplex counts by
dimension), and `higher_order` (the total of `d2` upward).

## Details

Higher-order structure is reported as `d2`, `d3`, ... : the number of
simplices of that dimension. A 2-simplex is a triangle of three mutually
connected states, a 3-simplex a tetrahedron of four. These count
*co-participation in a dense region*, not statistical interaction.

## See also

[`build_simplicial`](https://saqr.me/Nestimate/reference/build_simplicial.md),
[`betti_numbers`](https://saqr.me/Nestimate/reference/betti_numbers.md),
[`q_analysis`](https://saqr.me/Nestimate/reference/q_analysis.md),
[`outcome_model`](https://saqr.me/Nestimate/reference/outcome_model.md)

## Examples

``` r
m1 <- matrix(c(0, .6, .5, .6, 0, .4, .5, .4, 0), 3, 3,
             dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
m2 <- matrix(c(0, .2, 0, .2, 0, .1, 0, .1, 0), 3, 3,
             dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
simplicial_features(list(dense = m1, sparse = m2), threshold = 0.3)
#>   network threshold n_nodes n_edges b0 b1 euler max_q d1 d2 d3 d4 higher_order
#> 1   dense       0.3       3       3  1  0     1     2  3  1  0  0            1
#> 2  sparse       0.3       3       0  3  0     3     0  0  0  0  0            0

# Sweep the threshold rather than committing to one.
simplicial_features(list(dense = m1), threshold = c(0.1, 0.3, 0.5))
#>   network threshold n_nodes n_edges b0 b1 euler max_q d1 d2 d3 d4 higher_order
#> 1   dense       0.1       3       3  1  0     1     2  3  1  0  0            1
#> 2   dense       0.3       3       3  1  0     1     2  3  1  0  0            1
#> 3   dense       0.5       3       2  1  0     1     1  2  0  0  0            0
```
