# Print method for net_bayes_group

Print method for net_bayes_group

## Usage

``` r
# S3 method for class 'net_bayes_group'
print(x, ...)
```

## Arguments

- x:

  A `net_bayes_group` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
s <- data.frame(V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
                grp = c("X","X","Y","Y"))
nets <- build_network(s, method = "relative", group = "grp")
print(bayes_compare(nets, draws = 200, seed = 1))
#> Grouped Bayesian Dirichlet-Multinomial Comparison
#> Comparisons: X vs Y 
#> Use summary() for combined edge-level results.
```
