# Summary Method for net_association_rules

Summary Method for net_association_rules

## Usage

``` r
# S3 method for class 'net_association_rules'
summary(object, ...)
```

## Arguments

- object:

  A `net_association_rules` object.

- ...:

  Additional arguments (ignored).

## Value

A data frame summarizing the rules, invisibly.

## Examples

``` r
trans <- list(c("A","B","C"), c("A","B"), c("B","C","D"), c("A","C","D"))
rules <- association_rules(trans, min_support = 0.3, min_confidence = 0.5,
                           min_lift = 0)
summary(rules)
#>   antecedent consequent support confidence      lift conviction count
#> 1          D          C     0.5  1.0000000 1.3333333        Inf     2
#> 2          C          D     0.5  0.6666667 1.3333333       1.50     2
#> 3          A          B     0.5  0.6666667 0.8888889       0.75     2
#> 4          B          A     0.5  0.6666667 0.8888889       0.75     2
#> 5          A          C     0.5  0.6666667 0.8888889       0.75     2
#> 6          C          A     0.5  0.6666667 0.8888889       0.75     2
#> 7          B          C     0.5  0.6666667 0.8888889       0.75     2
#> 8          C          B     0.5  0.6666667 0.8888889       0.75     2
#>   n_transactions
#> 1              4
#> 2              4
#> 3              4
#> 4              4
#> 5              4
#> 6              4
#> 7              4
#> 8              4
```
