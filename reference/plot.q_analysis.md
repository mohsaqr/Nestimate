# Plot Q-Analysis

Two panels: Q-vector (components at each connectivity level) and
structure vector (max simplex dimension per node).

## Usage

``` r
# S3 method for class 'q_analysis'
plot(x, ...)
```

## Arguments

- x:

  A `q_analysis` object.

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
qa  <- q_analysis(sc)
plot(qa)

# }
```
