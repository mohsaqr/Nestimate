# Print Q-analysis results

Print Q-analysis results

## Usage

``` r
# S3 method for class 'q_analysis'
print(x, ...)
```

## Arguments

- x:

  A `q_analysis` object.

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
qa  <- q_analysis(sc)
print(qa)
#> Q-Analysis (max q = 2)
#>   Components: q2:1 q1:1 q0:1
#>   Fully connected at all q levels
#>   Structure: A:2 B:2 C:2
# }
```
