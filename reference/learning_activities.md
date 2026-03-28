# Online Learning Activity Indicators

Simulated binary time-series data for 200 students across 30 time
points. At each time point, one or more learning activities may be
active (1) or inactive (0). Activities: Reading, Video, Forum, Quiz,
Coding, Review. Includes temporal persistence (activities tend to
continue across adjacent time points).

## Usage

``` r
learning_activities
```

## Format

A data frame with 6,000 rows and 7 columns:

- student:

  Integer. Student identifier (1–200).

- Reading:

  Integer (0/1). Reading activity indicator.

- Video:

  Integer (0/1). Video watching indicator.

- Forum:

  Integer (0/1). Discussion forum indicator.

- Quiz:

  Integer (0/1). Quiz/assessment indicator.

- Coding:

  Integer (0/1). Coding practice indicator.

- Review:

  Integer (0/1). Review/revision indicator.

## Examples

``` r
# \donttest{
net <- build_network(learning_activities, method = "cna",
                     actor = "student", codes = c("Reading", "Video",
                     "Forum", "Quiz", "Coding", "Review"), window_size = 3)
net
#> Co-occurrence Network [undirected]
#>   Weights: [2681.000, 3290.000]  |  mean: 3047.333
#> 
#>   Weight matrix:
#>           Reading Video Forum Quiz Coding Review
#>   Reading    6100  3169  3211 2891   3138   3113
#>   Video      3169  6262  3183 2903   3278   3290
#>   Forum      3211  3183  5692 2942   2958   3045
#>   Quiz       2891  2903  2942 5379   2681   2725
#>   Coding     3138  3278  2958 2681   5886   3183
#>   Review     3113  3290  3045 2725   3183   6186 
# }
```
