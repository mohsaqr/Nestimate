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
net <- build_network(srl_strategies, method = "glasso",
                     params = list(gamma = 0.5))
net
#> Partial Correlation Network (EBICglasso) [undirected]
#>   Sample size: 250
#>   Weights: [0.020, 0.325]  |  +13 / -0 edges
#> 
#>   Weight matrix:
#>                Planning Monitoring Evaluating Elaboration Organization Rehearsal
#>   Planning        0.000      0.229      0.130       0.000        0.000     0.000
#>   Monitoring      0.229      0.000      0.300       0.039        0.000     0.000
#>   Evaluating      0.130      0.300      0.000       0.000        0.154     0.000
#>   Elaboration     0.000      0.039      0.000       0.000        0.262     0.158
#>   Organization    0.000      0.000      0.154       0.262        0.000     0.152
#>   Rehearsal       0.000      0.000      0.000       0.158        0.152     0.000
#>   Help_Seeking    0.000      0.000      0.000       0.000        0.000     0.000
#>   Time_Mgmt       0.130      0.000      0.000       0.000        0.000     0.000
#>   Effort_Reg      0.000      0.000      0.093       0.000        0.000     0.000
#>                Help_Seeking Time_Mgmt Effort_Reg
#>   Planning            0.000     0.130      0.000
#>   Monitoring          0.000     0.000      0.000
#>   Evaluating          0.000     0.000      0.093
#>   Elaboration         0.000     0.000      0.000
#>   Organization        0.000     0.000      0.000
#>   Rehearsal           0.000     0.000      0.000
#>   Help_Seeking        0.000     0.065      0.020
#>   Time_Mgmt           0.065     0.000      0.325
#>   Effort_Reg          0.020     0.325      0.000 
#> 
#>   Predictability (R²):
#>   Monitoring    0.173  ████████████████████████████████████████
#>   Evaluating    0.173  ████████████████████████████████████████
#>   Organization  0.137  ████████████████████████████████
#>   Time_Mgmt     0.133  ███████████████████████████████
#>   Effort_Reg    0.122  ████████████████████████████
#>   Planning      0.119  ████████████████████████████
#>   Elaboration   0.115  ██████████████████████████
#>   Rehearsal     0.067  ███████████████
#>   Help_Seeking  0.006  █
#> 
#>   Gamma: 0.50  |  Lambda: 0.1319
```
