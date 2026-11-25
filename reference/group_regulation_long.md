# Group regulation actions in collaborative learning

`group_regulation_long` is a coded dataset of students' regulation
actions during collaborative learning: 27,533 timestamped actions by
2,000 students working in 200 groups across three courses. It is in long
format, one row per action.

## Usage

``` r
group_regulation_long
```

## Format

A data frame with 27,533 rows and 6 columns:

- Actor:

  Integer student identifier, 1 to 2,000.

- Achiever:

  Character achievement level, `"High"` or `"Low"`.

- Group:

  Numeric collaboration group identifier, 1 to 200; each student belongs
  to one group.

- Course:

  Character course identifier, `"A"`, `"B"` or `"C"`.

- Time:

  POSIXct timestamp of the action.

- Action:

  Character regulation action, one of adapt, cohesion, consensus,
  coregulate, discuss, emotion, monitor, plan and synthesis.

## Source

Synthetically generated from the `group_regulation` dataset in the tna
package.

## See also

[`learning_activities`](https://pak.dynasite.org/Nestimate/reference/learning_activities.md),
[`srl_strategies`](https://pak.dynasite.org/Nestimate/reference/srl_strategies.md)

## Examples

``` r
net <- build_network(group_regulation_long,
                     method = "relative",
                     actor = "Actor", action = "Action", time = "Time")
net
#> Transition Network (relative probabilities) [directed]
#>   Weights: [0.001, 0.498]  |  mean: 0.115
#> 
#>   Strongest edges (10 of 78):
#>    from      to        weight
#>    cohesion  consensus 0.498 
#>    adapt     consensus 0.477 
#>    synthesis consensus 0.466 
#>    consensus plan      0.396 
#>    monitor   discuss   0.375 
#>    plan      plan      0.374 
#>    emotion   cohesion  0.325 
#>    discuss   consensus 0.321 
#>    emotion   consensus 0.320 
#>    plan      consensus 0.290  
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix;
#>   extract_initial_probs() the initial probabilities.
```
