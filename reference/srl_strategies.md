# Self-Regulated Learning Strategy Frequencies

Simulated frequency counts of 9 self-regulated learning (SRL) strategies
for 250 university students. Strategies are grouped into three clusters:
metacognitive (Planning, Monitoring, Evaluating), cognitive
(Elaboration, Organization, Rehearsal), and resource management
(Help_Seeking, Time_Mgmt, Effort_Reg). Within-cluster correlations are
moderate (0.3–0.6), cross-cluster correlations are weaker.

## Usage

``` r
srl_strategies
```

## Format

A data frame with 250 rows and 9 columns. Each column is an integer
count of how often the student used that strategy.

## Examples

``` r
# \donttest{
net <- build_network(srl_strategies, method = "glasso",
                     params = list(gamma = 0.5))
net
#> Partial Correlation Network (EBICglasso) [undirected]
#>   Sample size: 250
#>   Weights: [0.089, 0.413]  |  +13 / -0 edges
#> 
#>   Weight matrix:
#>                Planning Monitoring Evaluating Elaboration Organization Rehearsal
#>   Planning        0.000      0.295      0.161       0.000        0.000     0.000
#>   Monitoring      0.295      0.000      0.361       0.105        0.000     0.000
#>   Evaluating      0.161      0.361      0.000       0.000        0.221     0.000
#>   Elaboration     0.000      0.105      0.000       0.000        0.329     0.228
#>   Organization    0.000      0.000      0.221       0.329        0.000     0.218
#>   Rehearsal       0.000      0.000      0.000       0.228        0.218     0.000
#>   Help_Seeking    0.000      0.000      0.000       0.000        0.000     0.000
#>   Time_Mgmt       0.205      0.000      0.000       0.000        0.000     0.000
#>   Effort_Reg      0.000      0.000      0.161       0.000        0.000     0.000
#>                Help_Seeking Time_Mgmt Effort_Reg
#>   Planning            0.000     0.205      0.000
#>   Monitoring          0.000     0.000      0.000
#>   Evaluating          0.000     0.000      0.161
#>   Elaboration         0.000     0.000      0.000
#>   Organization        0.000     0.000      0.000
#>   Rehearsal           0.000     0.000      0.000
#>   Help_Seeking        0.000     0.141      0.089
#>   Time_Mgmt           0.141     0.000      0.413
#>   Effort_Reg          0.089     0.413      0.000 
#> 
#>   Gamma: 0.50  |  Lambda: 0.1319
# }
```
