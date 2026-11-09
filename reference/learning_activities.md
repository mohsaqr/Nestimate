# Online Learning Activity Indicators

Simulated binary (one-hot) indicator time series of six online learning
activities (Reading, Video, Forum, Quiz, Coding, Review) for 200
students over 30 time points. Each row is one time point of one student
and marks every activity active then (1) or not (0); at least one
activity is active in every row, and several often are. Activities tend
to continue from one time point to the next. The rows of each student
are consecutive and in time order, so the data are ready for
[`wtna()`](https://pak.dynasite.org/Nestimate/reference/wtna.md) with
`actor = "student"`.

## Usage

``` r
learning_activities
```

## Format

A data frame with 6,000 rows and 7 columns:

- student:

  Integer student identifier, 1 to 200, with 30 rows each.

- Reading:

  Integer reading indicator (0/1).

- Video:

  Integer video-watching indicator (0/1).

- Forum:

  Integer discussion-forum indicator (0/1).

- Quiz:

  Integer quiz or assessment indicator (0/1).

- Coding:

  Integer coding-practice indicator (0/1).

- Review:

  Integer review or revision indicator (0/1).

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
