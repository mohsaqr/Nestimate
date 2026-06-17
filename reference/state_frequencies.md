# Compute State Frequencies from Trajectory Data

Counts how often each state appears across all trajectories. Returns a
data frame sorted by frequency (descending).

## Usage

``` r
state_frequencies(data)
```

## Arguments

- data:

  A list of character vectors (trajectories) or a data.frame.

## Value

A data frame with columns: `state`, `count`, `proportion`.

## Examples

``` r
trajs <- list(c("A","B","C"), c("A","B","A"))
state_frequencies(trajs)
#>   state count proportion
#> 1     A     3     0.5000
#> 2     B     2     0.3333
#> 3     C     1     0.1667
```
