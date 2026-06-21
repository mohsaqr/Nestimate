# Print method for net_bayes

Print method for net_bayes

## Usage

``` r
# S3 method for class 'net_bayes'
print(x, ...)
```

## Arguments

- x:

  A `net_bayes` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
s1 <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
s2 <- data.frame(V1 = c("A","C","B"), V2 = c("C","B","A"))
b <- bayes_compare(build_network(s1, method = "relative"),
                   build_network(s2, method = "relative"),
                   draws = 500, seed = 1)
print(b)
#> Bayesian Dirichlet-Multinomial Comparison: Transition Network (relative probabilities) 
#>   Prior: Dirichlet(0.50)  |  Draws: 500  |  CI: 95%
#>   Thresholds: |mean diff| > 0.010, nearest CI bound > 0.001
#>   Nodes: 3  |  Edges compared: 6  |  Credibly different: 0
```
