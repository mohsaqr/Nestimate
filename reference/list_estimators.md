# List the registered estimators

`list_estimators()` returns the estimators
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
can call by name: the fifteen built-ins and any added with
[`register_estimator()`](https://pak.dynasite.org/Nestimate/reference/register_estimator.md).
Aliases such as `"tna"` or `"ebicglasso"` are not listed; the
*Estimators* section of
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
gives them.

## Usage

``` r
list_estimators()
```

## Value

A data frame with one row per registered estimator, sorted by name, and
columns `name`, `description` and `directed` (`TRUE` for directed
networks). With an empty registry it has no rows.

## Note

The idiographic package exports a function with the same name. With both
packages attached, the one attached last masks the other; write
`Nestimate::list_estimators()` to call this one.

## See also

[`register_estimator`](https://pak.dynasite.org/Nestimate/reference/register_estimator.md),
[`get_estimator`](https://pak.dynasite.org/Nestimate/reference/get_estimator.md)

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
#> 8            latency         Reply network weighted by response speed (latency)
#> 9                mgm Mixed Graphical Model (nodewise lasso, EBIC, LW threshold)
#> 10             ngram      n-gram transitions (adjacent pairs per n-gram window)
#> 11              pcor                         Unregularized partial correlations
#> 12          relative                    Row-normalized transition probabilities
#> 13           reverse        Reverse (reply) transitions: transpose of frequency
#> 14              wtna                     Window-based TNA transitions (one-hot)
#> 15 wtna_cooccurrence                   Window-based TNA co-occurrence (one-hot)
#>    directed
#> 1      TRUE
#> 2     FALSE
#> 3     FALSE
#> 4      TRUE
#> 5      TRUE
#> 6     FALSE
#> 7     FALSE
#> 8      TRUE
#> 9     FALSE
#> 10     TRUE
#> 11    FALSE
#> 12     TRUE
#> 13     TRUE
#> 14     TRUE
#> 15    FALSE
```
