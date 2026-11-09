# Rename the models of a `netobject_group`

`rename_models()` replaces the names of the networks in a
`netobject_group` and returns the group with nothing else changed. The
names label every later report (print methods, plot panels, bootstrap
and permutation output), so this is the step that turns the generic
labels
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
gives clusters (`"Cluster 1"`, `"Cluster 2"`) into meaningful ones
(`"High engagement"`, `"Low engagement"`) once it is known what each
cluster represents. It is the Nestimate counterpart of
[`tna::rename_groups()`](https://sonsoles.me/tna/reference/rename_groups.html).

## Usage

``` r
rename_models(x, new_names)

# S3 method for class 'netobject_group'
rename_models(x, new_names)

# Default S3 method
rename_models(x, new_names)
```

## Arguments

- x:

  A `netobject_group`, or an object inheriting from it. Any other object
  raises an error. A `net_mlvar` is accepted too, but its
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) find the three
  layers by their names (`temporal`, `contemporaneous`, `between`):
  after renaming, [`print()`](https://rdrr.io/r/base/print.html) shows
  none of the model's dimensions and
  [`summary()`](https://rdrr.io/r/base/summary.html) fails.

- new_names:

  The new names, a character vector in the order of the networks in `x`:
  one per network, none missing or empty, and no two the same. Anything
  else raises an error.

## Value

An object of the same class and length as `x`, its networks and
attributes (such as `group_col`) unchanged and its names replaced by
`new_names`.

## Examples

``` r
grp <- build_network(group_regulation_long, method = "tna",
                     actor = "Actor", action = "Action", time = "Time",
                     group = "Achiever")
names(grp)
#> [1] "High" "Low" 
names(rename_models(grp, c("High achievers", "Low achievers")))
#> [1] "High achievers" "Low achievers" 
```
