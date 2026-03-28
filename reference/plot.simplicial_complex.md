# Plot a Simplicial Complex

Produces a multi-panel summary: f-vector, simplicial degree ranking, and
degree-by-dimension heatmap.

## Usage

``` r
# S3 method for class 'simplicial_complex'
plot(x, ...)
```

## Arguments

- x:

  A `simplicial_complex` object.

- ...:

  Ignored.

## Value

A grid grob (invisibly).

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
plot(sc)

# }
```
