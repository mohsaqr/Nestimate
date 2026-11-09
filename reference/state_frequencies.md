# Compute State Frequencies from Trajectory Data

`state_frequencies()` counts how often each state occurs across all
trajectories, and returns one row per state, most frequent first.
Missing values are not counted.

## Usage

``` r
state_frequencies(data)
```

## Arguments

- data:

  The trajectories: a list of state vectors, a data frame or character
  matrix with one row per trajectory and one column per time point, a
  `tna` model or a `netobject` built from sequences. A data frame needs
  position columns named as a numbered run (`T1`, `T2`, ...); any other
  frame, a one-hot or association network and a 0/1 matrix raise an
  error of class `nestimate_bad_input`.

## Value

A data frame with one row per distinct state, sorted by `count`
(descending), with columns `state`, `count` and `proportion` (the
state's share of all observed states, rounded to 4 decimal places).

## Examples

``` r
trajs <- list(c("A","B","C"), c("A","B","A"))
state_frequencies(trajs)
#>   state count proportion
#> 1     A     3     0.5000
#> 2     B     2     0.3333
#> 3     C     1     0.1667
```
