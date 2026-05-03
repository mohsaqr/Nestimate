# List All Registered Estimators

Return a data frame summarising all registered network estimators.

## Usage

``` r
list_estimators()
```

## Value

A data frame with columns `name`, `description`, `directed`.

## See also

[`register_estimator`](https://mohsaqr.github.io/Nestimate/reference/register_estimator.md),
[`get_estimator`](https://mohsaqr.github.io/Nestimate/reference/get_estimator.md)

## Examples

``` r
list_estimators()
#>                 name                                                description
#> 1          attention                       Decay-weighted attention transitions
#> 2      co_occurrence                             Co-occurrence within sequences
#> 3                cor                               Pairwise correlation network
#> 4          frequency                            Raw transition frequency counts
#> 5             glasso                EBICglasso regularized partial correlations
#> 6              ising             Ising model (L1-penalized logistic regression)
#> 7                mgm Mixed Graphical Model (nodewise lasso, EBIC, LW threshold)
#> 8               pcor                         Unregularized partial correlations
#> 9           relative                    Row-normalized transition probabilities
#> 10              wtna                     Window-based TNA transitions (one-hot)
#> 11 wtna_cooccurrence                   Window-based TNA co-occurrence (one-hot)
#>    directed
#> 1      TRUE
#> 2     FALSE
#> 3     FALSE
#> 4      TRUE
#> 5     FALSE
#> 6     FALSE
#> 7     FALSE
#> 8     FALSE
#> 9      TRUE
#> 10     TRUE
#> 11    FALSE
```
