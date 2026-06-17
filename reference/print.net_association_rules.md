# Print Method for net_association_rules

Print Method for net_association_rules

## Usage

``` r
# S3 method for class 'net_association_rules'
print(x, ...)
```

## Arguments

- x:

  A `net_association_rules` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
trans <- list(c("A","B","C"), c("A","B"), c("B","C","D"), c("A","C","D"))
rules <- association_rules(trans, min_support = 0.3, min_confidence = 0.5,
                           min_lift = 0)
print(rules)
#> Association Rules  [8 rules | 4 items | 4 transactions]
#>   Support >= 0.30  |  Confidence >= 0.50  |  Lift >= 0.00
#> 
#>   Top rules (by lift):
#>     1. D -> C  (sup=0.500 conf=1.000 lift=1.33)
#>     2. C -> D  (sup=0.500 conf=0.667 lift=1.33)
#>     3. A -> B  (sup=0.500 conf=0.667 lift=0.89)
#>     4. B -> A  (sup=0.500 conf=0.667 lift=0.89)
#>     5. A -> C  (sup=0.500 conf=0.667 lift=0.89)
#>     6. C -> A  (sup=0.500 conf=0.667 lift=0.89)
#>     7. B -> C  (sup=0.500 conf=0.667 lift=0.89)
#>     8. C -> B  (sup=0.500 conf=0.667 lift=0.89)
```
