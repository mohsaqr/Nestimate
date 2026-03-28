# Print Method for nestimate_data

Print Method for nestimate_data

## Usage

``` r
# S3 method for class 'nestimate_data'
print(x, ...)
```

## Arguments

- x:

  A `nestimate_data` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
# \donttest{
events <- data.frame(
  actor  = c("u1","u1","u1","u2","u2","u2"),
  action = c("A","B","C","B","A","C"),
  time   = c(1,2,3,1,2,3)
)
nd <- prepare_data(events, action = "action",
                   actor = "actor", time = "time")
print(nd)
#> Prepared Data for Network Estimation
#>   Sessions: 2  |  Actions: 6  |  Max length: 3
#>   Actors: 2
#>   Time data: available
# }
```
