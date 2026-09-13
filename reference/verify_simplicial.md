# Verify Simplicial Complex Against igraph

Cross-validates clique finding and Betti numbers against igraph and
known topological invariants. Useful for testing.

## Usage

``` r
verify_simplicial(mat, threshold = 0)
```

## Arguments

- mat:

  A square adjacency matrix.

- threshold:

  Edge weight threshold.

## Value

Invisibly, a list with `cliques_match` (logical: do the simplices match
[`igraph::cliques()`](https://r.igraph.org/reference/cliques.html)
exactly), `n_simplices_ours`, `n_simplices_igraph`, `betti`, `euler`,
and `f_vector`. The comparison is also printed to the console.

## Examples

``` r
mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
colnames(mat) <- rownames(mat) <- c("A","B","C")
verify_simplicial(mat, threshold = 0.3)
#>   Cliques:  MATCH (7 simplices)
#>   Betti:    b0=1 b1=0 b2=0
#>   Euler:    1 (Euler-Poincare: VERIFIED)
```
