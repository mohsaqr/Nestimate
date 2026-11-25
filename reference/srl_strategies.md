# Simulated counts of self-regulated learning strategies

`srl_strategies` is a simulated dataset of how often 250 university
students used nine self-regulated learning (SRL) strategies, one row per
student and one count column per strategy. The strategies form three
clusters: metacognitive (Planning, Monitoring, Evaluating), cognitive
(Elaboration, Organization, Rehearsal) and resource management
(Help_Seeking, Time_Mgmt, Effort_Reg). Strategies of the same cluster
were simulated to correlate moderately and those of different clusters
more weakly; in the bundled sample the correlations range from 0.18 to
0.49 within clusters and from -0.13 to 0.31 across them.

## Usage

``` r
srl_strategies
```

## Format

A data frame with 250 rows and 9 columns, one row per student and one
column per strategy. Every column is numeric (double) and holds a
whole-number count of how often that student used the strategy; the
observed values range from 0 to 37, with no missing values.

## Examples

``` r
net <- build_network(srl_strategies, method = "glasso",
                     params = list(gamma = 0.5))
net
#> Partial Correlation Network (EBICglasso) [undirected]
#>   Sample size: 250
#>   Weights: [0.020, 0.325]  |  +13 / -0 edges
#> 
#>   Strongest edges (10 of 13):
#>    from         to           weight
#>    Time_Mgmt    Effort_Reg   0.325 
#>    Monitoring   Evaluating   0.300 
#>    Elaboration  Organization 0.262 
#>    Planning     Monitoring   0.229 
#>    Elaboration  Rehearsal    0.158 
#>    Evaluating   Organization 0.154 
#>    Organization Rehearsal    0.152 
#>    Planning     Time_Mgmt    0.130 
#>    Planning     Evaluating   0.130 
#>    Evaluating   Effort_Reg   0.093  
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix.
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
