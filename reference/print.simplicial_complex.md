# Print a simplicial complex

Print a simplicial complex

## Usage

``` r
# S3 method for class 'simplicial_complex'
print(x, ...)
```

## Arguments

- x:

  A `simplicial_complex` object.

- ...:

  Additional arguments (unused).

## Value

The input object, invisibly.

## Examples

``` r
# \donttest{
seqs <- data.frame(
  V1 = c("A","B","C","A","B"),
  V2 = c("B","C","A","B","C"),
  V3 = c("C","A","B","C","A")
)
net <- build_network(seqs, method = "relative")
sc  <- build_simplicial(net, type = "clique")
print(sc)
#> Clique Complex 
#>   3 nodes, 7 simplices, dimension 2
#>   Density: 100.0%  |  Mean dim: 0.71  |  Euler: 1
#>   f-vector: (f0=3 f1=3 f2=1)
#>   Betti: b0=1
#>   Nodes: A, B, C 
# }
```
