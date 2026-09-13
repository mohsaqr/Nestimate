# Rename the models of a `netobject_group`

Replaces the names of the constituent networks in a `netobject_group`
(or any object inheriting from it). Useful when
[`build_network()`](https://saqr.me/Nestimate/reference/build_network.md)
produced generic labels (e.g. `"Cluster 1"`, `"Cluster 2"`) and you want
to substitute meaningful ones (e.g. `"High engagement"`,
`"Low engagement"`).

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

  A `netobject_group` (or any object inheriting from it, such as
  `net_mlvar`).

- new_names:

  A character vector of new names. Must have the same length as `x`,
  contain no `NA` or empty strings, and be unique.

## Value

A `netobject_group` of the same class and length as `x`, with
[`names()`](https://rdrr.io/r/base/names.html) replaced by `new_names`.
The constituent networks are returned unchanged.

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
