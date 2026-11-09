# Tidy Topological Features for One or Many Networks

`simplicial_features()` builds a simplicial complex for each network
with
[`build_simplicial()`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md)
and returns its topological summaries as a tidy data frame, one row per
network and one column per feature, ready to use as regression
predictors or to join onto unit-level outcomes. It takes one network, a
group of networks, the within-cluster networks of an `mcml`, or a list
of networks, at one threshold or a sweep of them.

Higher-order structure is reported as `d2`, `d3`, ...: the number of
simplices of that dimension. A 2-simplex is a triangle of three mutually
connected states, a 3-simplex a tetrahedron of four. These count
co-participation in a dense region, not statistical interaction. The
counts are a step function of `threshold`, so one threshold is a choice,
not a result: sweep several and check that a conclusion holds across
them.

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

  The networks: a `netobject`, `cograph_network`, `tna` model or square
  weight matrix (one row); a `netobject_group`, or a list of single
  networks or matrices (one row per member, named by the list, or
  `network_1`, `network_2`, ... when unnamed); or an `mcml` object (one
  row per cluster, from its within-cluster networks; one built without
  them raises an error).

- threshold:

  The smallest absolute edge weight that makes an edge, passed to
  [`build_simplicial()`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md).
  Default `0` (every non-zero edge). A vector sweeps it, giving one row
  per network per threshold. Values must be finite and not negative.
  `type = "vr"` and `type = "pathway"` ignore it, although the
  `threshold` column still records it.

- max_dim:

  The highest simplex dimension built, a whole number of at least 1 (a
  fractional value is truncated). Default `4`. When a network holds a
  larger clique, the complex is truncated and
  [`betti_numbers()`](https://pak.dynasite.org/Nestimate/reference/betti_numbers.md)
  warns with class `nestimate_truncated_complex`; the warning concerns
  the Betti number of dimension `max_dim`, so it affects `b1` only when
  `max_dim = 1`.

- normalize:

  Whether to divide the simplex counts (`d1` onwards and `higher_order`)
  by the number of nodes, so networks of different size are comparable;
  `n_edges` stays a count. Default `FALSE`.

- type:

  The complex construction passed to
  [`build_simplicial()`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md):
  `"clique"` (default), `"pathway"`, `"vr"` or its alias `"rips"`.
  `"pathway"` needs networks of class `netobject` or `tna`, from which a
  higher-order network is built; a matrix, a bare `cograph_network` or
  an `mcml` (whose cluster layers are matrices) raises an error.

## Value

A data frame with one row per network and threshold, ordered by network
name and then threshold, and columns `network`, `threshold`, `n_nodes`,
`n_edges`, `b0` and `b1` (the Betti numbers: connected components and
independent cycles, 0 when the complex has no such dimension), `euler`
(the Euler characteristic of the complex as built), `max_q` (the top
q-level of
[`q_analysis()`](https://pak.dynasite.org/Nestimate/reference/q_analysis.md),
the dimension of the largest maximal simplex), `d1`, `d2`, ... up to
dimension `max_dim` (simplex counts by dimension) and `higher_order`
(the total from `d2` upward).

## See also

[`build_simplicial`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md),
[`betti_numbers`](https://pak.dynasite.org/Nestimate/reference/betti_numbers.md),
[`q_analysis`](https://pak.dynasite.org/Nestimate/reference/q_analysis.md),
[`outcome_model`](https://pak.dynasite.org/Nestimate/reference/outcome_model.md)

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
