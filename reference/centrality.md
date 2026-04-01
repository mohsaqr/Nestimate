# Compute Centrality Measures for a Network

Computes centrality measures from a `netobject` or `netobject_group`.
For directed networks the default measures are InStrength, OutStrength,
and Betweenness. For undirected networks the defaults are Strength
(column sums) and Betweenness.

## Usage

``` r
centrality(x, ...)

# S3 method for class 'netobject'
centrality(x, measures = NULL, loops = FALSE, centrality_fn = NULL, ...)

# S3 method for class 'netobject_group'
centrality(x, measures = NULL, loops = FALSE, centrality_fn = NULL, ...)

# S3 method for class 'cograph_network'
centrality(x, measures = NULL, loops = FALSE, centrality_fn = NULL, ...)

# S3 method for class 'mcml'
centrality(x, measures = NULL, loops = FALSE, centrality_fn = NULL, ...)
```

## Arguments

- x:

  A `netobject` or `netobject_group`.

- ...:

  Additional arguments (ignored).

- measures:

  Character vector. Centrality measures to compute. Built-in:
  `"InStrength"`, `"OutStrength"`, `"Betweenness"`, `"InCloseness"`,
  `"OutCloseness"`, `"Closeness"`. Default depends on directedness.

- loops:

  Logical. Include self-loops (diagonal) in computation? Default:
  `FALSE`.

- centrality_fn:

  Optional function. Custom centrality function that takes a weight
  matrix and returns a named list of centrality vectors.

## Value

For a `netobject`: a data frame with node names as rows and centrality
measures as columns. For a `netobject_group`: a named list of such data
frames (one per group).

## Examples

``` r
# \donttest{
seqs <- data.frame(
  V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
  V3 = c("C","A","C","B"))
net <- build_network(seqs, method = "relative")
centrality(net)
#>   InStrength OutStrength Betweenness
#> A          1           1           2
#> B          1           1           2
#> C          1           1           2
# }
```
