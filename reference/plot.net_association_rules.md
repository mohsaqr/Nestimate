# Plot Method for net_association_rules

Scatter plot of association rules: support vs confidence, with point
size proportional to lift.

## Usage

``` r
# S3 method for class 'net_association_rules'
plot(x, ...)
```

## Arguments

- x:

  A `net_association_rules` object.

- ...:

  Additional arguments passed to `ggplot2` functions.

## Value

A `ggplot` object, invisibly.

## Examples

``` r
trans <- list(c("A","B","C"), c("A","B"), c("B","C","D"),
              c("A","C","D"), c("A","B","D"), c("B","C"))
rules <- association_rules(trans, min_support = 0.3, min_confidence = 0.3,
                           min_lift = 0)
plot(rules)

```
