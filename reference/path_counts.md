# Count Path Frequencies in Trajectory Data

`path_counts()` counts how often each path of `k` consecutive states (a
k-gram) occurs across all trajectories, and returns one row per path,
most frequent first. Every position of a trajectory starts a path, so
the paths of one trajectory overlap. The counts are descriptive: they
show which runs of states dominate the data before a model such as
[`build_mogen()`](https://pak.dynasite.org/Nestimate/reference/build_mogen.md)
or
[`build_hon()`](https://pak.dynasite.org/Nestimate/reference/build_hon.md)
is fitted to them.

## Usage

``` r
path_counts(data, k = 2L, top = NULL)
```

## Arguments

- data:

  The trajectories: a list of state vectors, a data frame or character
  matrix with one row per trajectory and one column per time point, a
  `tna` model or a `netobject` built from sequences. A data frame needs
  position columns named as a numbered run (`T1`, `T2`, ...); any other
  frame, a one-hot or association network and a 0/1 matrix raise an
  error of class `nestimate_bad_input`. A missing value ends a run of
  observed states, and no path spans it.

- k:

  The number of states in each path, a whole number of at least 2; a
  non-integer value is an error, not truncated. Default `2`, which
  counts single transitions; `3` counts two-step paths, and so on.

- top:

  The number of most frequent paths to keep. Default `NULL`, every path.

## Value

A data frame with one row per distinct path, sorted by `count`
(descending), with columns `path` (the k states in arrow notation, such
as `"A -> B"`), `count` and `proportion` (the path's share of all paths
counted, rounded to 4 decimal places). When no trajectory holds `k`
consecutive observed states, the data frame has no rows.

## Examples

``` r
trajs <- list(c("A","B","C","D"), c("A","B","D","C"))
path_counts(trajs, k = 2)
#>     path count proportion
#> 1 A -> B     2     0.3333
#> 2 B -> C     1     0.1667
#> 3 B -> D     1     0.1667
#> 4 C -> D     1     0.1667
#> 5 D -> C     1     0.1667

# \donttest{
path_counts(trajs, k = 3, top = 10)
#>          path count proportion
#> 1 A -> B -> C     1       0.25
#> 2 A -> B -> D     1       0.25
#> 3 B -> C -> D     1       0.25
#> 4 B -> D -> C     1       0.25
# }
```
