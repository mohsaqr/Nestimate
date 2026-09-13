# Centrality Stability Coefficient (CS-coefficient)

Estimates the stability of centrality indices under case-dropping. For
each drop proportion, sequences are randomly removed and the network is
re-estimated. The correlation between the original and subset centrality
values is computed. The CS-coefficient is the maximum proportion of
cases that can be dropped while maintaining a correlation above
`threshold` in at least `certainty` of bootstrap samples.

For transition methods, uses pre-computed per-sequence count matrices
for fast resampling. Strength centralities (InStrength, OutStrength) are
computed directly from the matrix without igraph.

## Usage

``` r
centrality_stability(
  x,
  measures = c("InStrength", "OutStrength", "Betweenness"),
  iter = 1000L,
  drop_prop = seq(0.1, 0.9, by = 0.1),
  threshold = 0.7,
  certainty = 0.95,
  method = "pearson",
  centrality_fn = NULL,
  loops = FALSE,
  normalize = FALSE,
  invert = TRUE,
  normalize_diffusion = TRUE,
  seed = NULL
)
```

## Arguments

- x:

  A `netobject` from
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md),
  a `cograph_network`, or a `netobject_group` / `mcml` (each constituent
  network is assessed and a `net_stability_group` is returned).

- measures:

  Character vector. Centrality measures to assess. Defaults to
  `c("InStrength", "OutStrength", "Betweenness")`. Pass `"all"` for
  every built-in measure: `"OutStrength"`, `"InStrength"`,
  `"ClosenessIn"`, `"ClosenessOut"`, `"Closeness"`, `"Betweenness"`,
  `"BetweennessRSP"`, `"Diffusion"`, and `"Clustering"`. The legacy
  aliases `"InCloseness"` and `"OutCloseness"` are also accepted. Custom
  measures beyond these are valid only when a `centrality_fn` is
  supplied to resolve them.

- iter:

  Integer. Number of bootstrap iterations per drop proportion (default:
  1000).

- drop_prop:

  Numeric vector. Proportions of cases to drop (default:
  `seq(0.1, 0.9, by = 0.1)`).

- threshold:

  Numeric. Minimum correlation to consider stable (default: 0.7).

- certainty:

  Numeric. Required proportion of iterations above threshold (default:
  0.95).

- method:

  Character. Correlation method: `"pearson"`, `"spearman"`, or
  `"kendall"` (default: `"pearson"`).

- centrality_fn:

  Optional function. A custom centrality function that takes a weight
  matrix and returns a named list of centrality vectors. When `NULL`
  (default), all built-in measures are computed internally:
  `"InStrength"`/`"OutStrength"` via `colSums`/`rowSums`,
  `"Betweenness"`/ `"ClosenessIn"`/`"ClosenessOut"`/`"Closeness"` via an
  internal Floyd-Warshall shortest-path routine, and `"BetweennessRSP"`,
  `"Diffusion"` and `"Clustering"` from the weight matrix directly. When
  provided, the function is called as `centrality_fn(mat)` and is used
  only for requested measures that are not one of the built-ins; it
  should return a named list (e.g., `list(my_metric = ...)`).

- loops:

  Logical. If `FALSE` (default), self-loops (diagonal) are excluded from
  centrality computation. This does not modify the stored matrix.

- normalize:

  Logical. Range-normalize all requested measures using the same
  transformation as `tna::centralities(normalize = TRUE)`. Default:
  `FALSE`.

- invert:

  Logical. Invert weights for shortest-path measures? Default: `TRUE`,
  matching `tna`.

- normalize_diffusion:

  Logical. Range-normalize `Diffusion` even when `normalize = FALSE`.
  Default: `TRUE`.

- seed:

  Integer or NULL. RNG seed for reproducibility.

## Value

An object of class `"net_stability"`: a list with

- cs:

  Named numeric vector of CS-coefficients, one per retained measure.

- correlations:

  Named list of `iter` x `length(drop_prop)` matrices of correlation
  values, one per retained measure.

- measures:

  Character vector of the measures actually assessed (see the
  zero-variance rule below).

- drop_prop:

  Drop proportions used.

- threshold:

  Stability threshold.

- certainty:

  Required certainty level.

- iter:

  Number of iterations.

- method:

  Correlation method.

A `netobject_group` or `mcml` input instead returns a
`"net_stability_group"`: a named list of one `net_stability` per
constituent network.

Zero-variance measures are handled by two different rules, both
long-standing behaviour. When *some* requested measures have zero
variance on the original network (for example `"OutStrength"` on a
row-normalised transition network), those measures are **dropped**:
`$cs`, `$correlations` and `$measures` cover only the retained ones.
When *every* requested measure has zero variance a warning is issued and
**all** requested names are returned with `cs = 0` and all-`NA`
correlation matrices.

## References

Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
psychological networks and their accuracy: A tutorial paper. *Behavior
Research Methods* 50(1), 195-212.
[doi:10.3758/s13428-017-0862-1](https://doi.org/10.3758/s13428-017-0862-1)

## See also

[`build_network`](https://saqr.me/Nestimate/reference/build_network.md),
[`network_reliability`](https://saqr.me/Nestimate/reference/network_reliability.md)

## Examples

``` r
seqs <- data.frame(
  T1 = c("plan", "code", "debug", "plan", "test", "code"),
  T2 = c("code", "debug", "code", "plan", "code", "test"),
  T3 = c("debug", "code", "plan", "code", "debug", "plan"),
  T4 = c("test", "plan", "test", "debug", "plan", "code")
)
net <- build_network(seqs, method = "relative")
cs <- centrality_stability(net, iter = 10, drop_prop = 0.3, seed = 1)
# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
  V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
)
net <- build_network(seqs, method = "relative")
cs <- centrality_stability(net, iter = 100, seed = 42,
  measures = c("InStrength", "OutStrength"))
print(cs)
#> Centrality Stability (100 iterations, threshold = 0.7)
#>   Drop proportions: 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9
#> 
#>   CS-coefficients:
#>     InStrength       0.00
#>     OutStrength      0.30
# }
```
