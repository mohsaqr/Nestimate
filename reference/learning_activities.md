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
head(learning_activities)
#>   student Reading Video Forum Quiz Coding Review
#> 1       1       1     0     0    1      0      0
#> 2       1       1     0     0    1      1      0
#> 3       1       1     0     0    0      1      1
#> 4       1       1     0     0    1      1      1
#> 5       1       1     1     0    0      1      0
#> 6       1       1     1     1    0      1      1
```
