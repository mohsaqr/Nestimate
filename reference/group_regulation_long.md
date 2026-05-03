# Group Regulation in Collaborative Learning (Long Format)

Students' regulation strategies during collaborative learning, in long
format. Contains 27,533 timestamped action records from multiple
students working in groups across two courses.

## Usage

``` r
group_regulation_long
```

## Format

A data frame with 27,533 rows and 6 columns:

- Actor:

  Integer. Student identifier.

- Achiever:

  Character. Achievement level: `"High"` or `"Low"`.

- Group:

  Numeric. Collaboration group identifier.

- Course:

  Character. Course identifier (`"A"`, `"B"`, or `"C"`).

- Time:

  POSIXct. Timestamp of the action.

- Action:

  Character. Regulation action (e.g., cohesion, consensus, discuss,
  synthesis).

## Source

Synthetically generated from the `group_regulation` dataset in the tna
package.

## See also

[`learning_activities`](https://mohsaqr.github.io/Nestimate/reference/learning_activities.md),
[`srl_strategies`](https://mohsaqr.github.io/Nestimate/reference/srl_strategies.md)

## Examples

``` r
net <- build_network(group_regulation_long,
                     method = "relative",
                     actor = "Actor", action = "Action", time = "Time")
net
#> Transition Network (relative probabilities) [directed]
#>   Weights: [0.001, 0.498]  |  mean: 0.115
#> 
#>   Weight matrix:
#>              adapt cohesion consensus coregulate discuss emotion monitor  plan
#>   adapt      0.000    0.273     0.477      0.022   0.059   0.120   0.033 0.016
#>   cohesion   0.003    0.027     0.498      0.119   0.060   0.116   0.033 0.141
#>   consensus  0.005    0.015     0.082      0.188   0.188   0.073   0.047 0.396
#>   coregulate 0.016    0.036     0.135      0.023   0.274   0.172   0.086 0.239
#>   discuss    0.071    0.048     0.321      0.084   0.195   0.106   0.022 0.012
#>   emotion    0.002    0.325     0.320      0.034   0.102   0.077   0.036 0.100
#>   monitor    0.011    0.056     0.159      0.058   0.375   0.091   0.018 0.216
#>   plan       0.001    0.025     0.290      0.017   0.068   0.147   0.076 0.374
#>   synthesis  0.235    0.034     0.466      0.044   0.063   0.071   0.012 0.075
#>              synthesis
#>   adapt          0.000
#>   cohesion       0.004
#>   consensus      0.008
#>   coregulate     0.019
#>   discuss        0.141
#>   emotion        0.003
#>   monitor        0.016
#>   plan           0.002
#>   synthesis      0.000 
#> 
#>   Initial probabilities:
#>   consensus     0.214  ████████████████████████████████████████
#>   plan          0.204  ██████████████████████████████████████
#>   discuss       0.175  █████████████████████████████████
#>   emotion       0.151  ████████████████████████████
#>   monitor       0.144  ███████████████████████████
#>   cohesion      0.060  ███████████
#>   synthesis     0.019  ████
#>   coregulate    0.019  ████
#>   adapt         0.011  ██
```
