# Verify Simplicial Complex Against igraph

`verify_simplicial()` builds the clique complex of `mat` with
[`build_simplicial()`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md)
and checks it against two independent references: its simplices against
the cliques
[`igraph::cliques()`](https://r.igraph.org/reference/cliques.html) finds
in the same thresholded graph, and its Euler characteristic against the
alternating sum of its Betti numbers (the Euler-Poincare theorem). It
prints the result of each check and the Betti numbers, and is meant for
testing. It needs igraph, and raises an error without it.

## Usage

``` r
verify_simplicial(mat, threshold = 0)
```

## Arguments

- mat:

  A square weight matrix, or a network object whose weight matrix
  [`build_simplicial()`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md)
  reads (a `netobject`, `tna` model, `cograph_network` or `net_hon`).

- threshold:

  The smallest absolute edge weight an edge needs to enter, as in
  [`build_simplicial()`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md).
  Default `0`.

## Value

A list, invisibly: `cliques_match` (whether the simplices equal the
cliques from
[`igraph::cliques()`](https://r.igraph.org/reference/cliques.html)
exactly), `n_simplices_ours`, `n_simplices_igraph`, `betti`, `euler` and
`f_vector`. The two checks are printed as a side effect.

## Details

The complex is built with the default `max_dim = 10`, while
[`igraph::cliques()`](https://r.igraph.org/reference/cliques.html) has
no size limit, so a graph with a clique of more than 11 nodes reports a
mismatch, with a warning of class `nestimate_truncated_complex` from
[`betti_numbers()`](https://pak.dynasite.org/Nestimate/reference/betti_numbers.md).

## Examples

``` r
mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
colnames(mat) <- rownames(mat) <- c("A","B","C")
verify_simplicial(mat, threshold = 0.3)
#>   Cliques:  MATCH (7 simplices)
#>   Betti:    b0=1 b1=0 b2=0
#>   Euler:    1 (Euler-Poincare: VERIFIED)
```
