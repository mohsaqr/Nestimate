# Persistent Homology

Computes persistent homology by building simplicial complexes at
decreasing weight thresholds and tracking the birth/death of topological
features.

## Usage

``` r
persistent_homology(x, n_steps = 20L, max_dim = 3L)
```

## Arguments

- x:

  A square matrix, `tna`, or `netobject`.

- n_steps:

  Number of filtration steps (default 20).

- max_dim:

  Maximum simplex dimension to track (default 3).

## Value

A `persistent_homology` object with:

- betti_curve:

  Data frame: `threshold`, `dimension`, `betti` at each filtration step.

- persistence:

  Data frame of birth-death pairs: `dimension`, `birth`, `death`,
  `persistence`.

- thresholds:

  Numeric vector of filtration thresholds.

## Examples

``` r
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
  V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
)
net <- build_network(seqs, method = "relative")
ph <- persistent_homology(net, n_steps = 15)
print(ph)
#> Persistent Homology
#>   15 filtration steps [0.3636 → 0.0036]
#>   Features: b0: 4 (1 persistent) 
#>   Longest-lived:
#>     b0: 0.3636 → 0.0000 (life: 0.3636)
#>     b0: 0.3636 → 0.3122 (life: 0.0514)
#>     b0: 0.3636 → 0.3122 (life: 0.0514)
# }
```
