# Co-occurrence Networks with Nestimate

## Introduction

Co-occurrence networks analyze **binary indicator data**—situations
where multiple states can be active (1) or inactive (0) simultaneously.
Co-occurrence networks differ from transition networks in that the
former methods capture **contemporaneous relationships** (i.e., which
states tend to occur together?) whereas the latter capture sequences
(i.e., which states occur after one another?)

Nestimate provides two methods for binary data:

- Simple **Co-occurrence networks** (`method = "co_occurrence"` or
  `"cna"`) count how often pairs of states are both active at the same
  time point
- **Ising networks** (`method = "ising"`) use L1-regularized logistic
  regression to estimate conditional dependencies, producing sparse
  networks

This vignette demonstrates both methods using the `learning_activities`
dataset—binary indicators of 6 learning activities across 200 students
and 30 time points.

## Data

The `learning_activities` dataset contains 6,000 observations (200
students x 30 time points). At each time point, each of 6 learning
activities is either active (1) or inactive (0).

``` r
library(Nestimate)
data(learning_activities)
head(learning_activities, 10)
#>    student Reading Video Forum Quiz Coding Review
#> 1        1       1     0     0    1      0      0
#> 2        1       1     0     0    1      1      0
#> 3        1       1     0     0    0      1      1
#> 4        1       1     0     0    1      1      1
#> 5        1       1     1     0    0      1      0
#> 6        1       1     1     1    0      1      1
#> 7        1       1     1     1    0      1      1
#> 8        1       1     1     1    0      1      1
#> 9        1       1     1     1    0      1      1
#> 10       1       0     1     1    0      1      1
```

The 6 activities are:

``` r
activities <- c("Reading", "Video", "Forum", "Quiz", "Coding", "Review")
colSums(learning_activities[, activities])
#> Reading   Video   Forum    Quiz  Coding  Review 
#>    2542    2622    2424    2273    2488    2594
```

Multiple activities can be active simultaneously:

``` r
# How many activities are active at each time point?
n_active <- rowSums(learning_activities[, activities])
table(n_active)
#> n_active
#>    1    2    3    4    5    6 
#> 1367 1766 1712  884  254   17
```

## Co-occurrence Network

Co-occurrence networks count how often pairs of states are **both
active** at the same time point. An edge between A and B indicates they
frequently co-occur.

### Basic Co-occurrence

``` r
net_cna <- build_network(learning_activities,
                         method = "co_occurrence",
                         codes = activities,
                         actor = "student")
net_cna
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
```

**Interpretation**: Edge weights represent raw co-occurrence counts
summed across all students and time points. Higher weights mean those
activities frequently happen together.

### Normalized Co-occurrence

Raw counts can be misleading—frequent activities will have high
co-occurrence simply because they’re common. Normalizing by expected
co-occurrence (under independence) reveals associations beyond base
rates.

``` r
# View the co-occurrence matrix
round(net_cna$weights, 0)
#>         Reading Video Forum Quiz Coding Review
#> Reading    6100  3169  3211 2891   3138   3113
#> Video      3169  6262  3183 2903   3278   3290
#> Forum      3211  3183  5692 2942   2958   3045
#> Quiz       2891  2903  2942 5379   2681   2725
#> Coding     3138  3278  2958 2681   5886   3183
#> Review     3113  3290  3045 2725   3183   6186
```

The diagonal shows how often each activity occurs (self-co-occurrence =
frequency).

### Windowed Co-occurrence

You can aggregate across time windows before computing co-occurrence.
This captures activities that occur in the same temporal neighborhood,
not just the exact same time point:

``` r
net_cna_windowed <- build_network(learning_activities,
                                   method = "co_occurrence",
                                   codes = activities,
                                   actor = "student",
                                   window_size = 10)
net_cna_windowed
#> Co-occurrence Network [undirected]
#>   Weights: [9096.000, 11181.000]  |  mean: 10230.667
#> 
#>   Weight matrix:
#>           Reading Video Forum  Quiz Coding Review
#>   Reading   16894 10639 10603  9735  10493  10699
#>   Video     10639 17012 10511  9710  10981  11181
#>   Forum     10603 10511 15090  9777   9992  10353
#>   Quiz       9735  9710  9777 14495   9123   9096
#>   Coding    10493 10981  9992  9123  15956  10567
#>   Review    10699 11181 10353  9096  10567  16696
```

Larger windows capture broader temporal associations at the cost of
temporal precision.

## Ising Network

Ising networks use **L1-regularized logistic regression** to estimate
conditional dependencies between binary variables. Each variable is
regressed on all others, and the resulting coefficients form the network
edges.

Key advantages over simple co-occurrence:

- **Sparsity**: L1 regularization shrinks weak associations to exactly
  zero
- **Conditional**: Edges represent direct relationships controlling for
  other variables
- **Interpretable**: Edge weights reflect log-odds of co-activation

### Basic Ising Network

``` r
# Aggregate to student-level summaries for Ising (requires cross-sectional data)
student_summary <- aggregate(learning_activities[, activities],
                              by = list(student = learning_activities$student),
                              FUN = function(x) as.integer(mean(x) > 0.5))
student_summary <- student_summary[, -1]  # Remove student column

net_ising <- build_network(student_summary,
                           method = "ising",
                           params = list(gamma = 0.25))
net_ising
#> Ising Model Network [undirected]
#>   Sample size: 200
#> 
#>   Weight matrix:
#>           Reading Video Forum Quiz Coding Review
#>   Reading       0     0     0    0      0      0
#>   Video         0     0     0    0      0      0
#>   Forum         0     0     0    0      0      0
#>   Quiz          0     0     0    0      0      0
#>   Coding        0     0     0    0      0      0
#>   Review        0     0     0    0      0      0 
#> 
#>   Gamma: 0.25  |  Rule: AND
#>   Thresholds: [-0.895, -0.532]
```

### Ising Parameters

The `gamma` parameter controls sparsity via EBIC model selection:

- `gamma = 0`: Less sparse (BIC-like selection)
- `gamma = 0.25`: Moderate sparsity (default)
- `gamma = 0.5`: More sparse

``` r
net_ising_sparse <- build_network(student_summary,
                                   method = "ising",
                                   params = list(gamma = 0.5))

# Compare edge counts
cat("Gamma = 0.25:", sum(net_ising$weights != 0), "edges\n")
#> Gamma = 0.25: 0 edges
cat("Gamma = 0.50:", sum(net_ising_sparse$weights != 0), "edges\n")
#> Gamma = 0.50: 0 edges
```

### Symmetrization Rules

Ising estimation produces asymmetric coefficients (A predicting B may
differ from B predicting A). The `rule` parameter controls
symmetrization:

- `"AND"` (default): Keep edge only if both directions are non-zero
- `"OR"`: Keep edge if either direction is non-zero

``` r
net_and <- build_network(student_summary, method = "ising",
                         params = list(gamma = 0.25, rule = "AND"))
net_or <- build_network(student_summary, method = "ising",
                        params = list(gamma = 0.25, rule = "OR"))

cat("AND rule edges:", sum(net_and$weights != 0), "\n")
#> AND rule edges: 0
cat("OR rule edges:", sum(net_or$weights != 0), "\n")
#> OR rule edges: 0
```

`"AND"` is more conservative; `"OR"` retains more edges.
