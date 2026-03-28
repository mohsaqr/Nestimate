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
#>                 name                                    description directed
#> 1          attention           Decay-weighted attention transitions     TRUE
#> 2      co_occurrence                 Co-occurrence within sequences    FALSE
#> 3                cor                   Pairwise correlation network    FALSE
#> 4          frequency                Raw transition frequency counts     TRUE
#> 5             glasso    EBICglasso regularized partial correlations    FALSE
#> 6              ising Ising model (L1-penalized logistic regression)    FALSE
#> 7               pcor             Unregularized partial correlations    FALSE
#> 8           relative        Row-normalized transition probabilities     TRUE
#> 9               wtna         Window-based TNA transitions (one-hot)     TRUE
#> 10 wtna_cooccurrence       Window-based TNA co-occurrence (one-hot)    FALSE
```
