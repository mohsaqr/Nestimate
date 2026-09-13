# List All Registered Estimators

Return a data frame summarising all registered network estimators.

## Usage

``` r
list_estimators()
```

## Value

A data frame with columns `name`, `description`, `directed`.

## See also

[`register_estimator`](https://saqr.me/Nestimate/reference/register_estimator.md),
[`get_estimator`](https://saqr.me/Nestimate/reference/get_estimator.md)

## Examples

``` r
list_estimators()
#>                 name                                                description
#> 1          attention                       Decay-weighted attention transitions
#> 2      co_occurrence                             Co-occurrence within sequences
#> 3                cor                               Pairwise correlation network
#> 4          frequency                            Raw transition frequency counts
#> 5                gap             Gap-allowed transitions weighted by 1/distance
#> 6             glasso                EBICglasso regularized partial correlations
#> 7              ising             Ising model (L1-penalized logistic regression)
#> 8                mgm Mixed Graphical Model (nodewise lasso, EBIC, LW threshold)
#> 9              ngram      n-gram transitions (adjacent pairs per n-gram window)
#> 10              pcor                         Unregularized partial correlations
#> 11          relative                    Row-normalized transition probabilities
#> 12           reverse        Reverse (reply) transitions: transpose of frequency
#> 13              wtna                     Window-based TNA transitions (one-hot)
#> 14 wtna_cooccurrence                   Window-based TNA co-occurrence (one-hot)
#>    directed
#> 1      TRUE
#> 2     FALSE
#> 3     FALSE
#> 4      TRUE
#> 5      TRUE
#> 6     FALSE
#> 7     FALSE
#> 8     FALSE
#> 9      TRUE
#> 10    FALSE
#> 11     TRUE
#> 12     TRUE
#> 13     TRUE
#> 14    FALSE
```
