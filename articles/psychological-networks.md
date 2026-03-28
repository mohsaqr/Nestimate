# Psychological Network Analysis with Nestimate

Nestimate supports psychological network analysis (PNA) through partial
correlation and graphical lasso estimation. This vignette uses the
`srl_strategies` dataset — frequency counts of 9 self-regulated learning
strategies for 250 students — to estimate, regularize, and bootstrap a
psychological network.

## Data

The 9 strategies fall into three clusters: metacognitive (Planning,
Monitoring, Evaluating), cognitive (Elaboration, Organization,
Rehearsal), and resource management (Help_Seeking, Time_Mgmt,
Effort_Reg).

``` r
library(Nestimate)
data(srl_strategies)
head(srl_strategies)
#>   Planning Monitoring Evaluating Elaboration Organization Rehearsal
#> 1       13         15         13           5            3        13
#> 2       14         11         12          10           25        19
#> 3       24         20         22          14            3        10
#> 4       19         18         15          17           27        13
#> 5       17         21         15           8            5        12
#> 6        4          6          5          26           24        25
#>   Help_Seeking Time_Mgmt Effort_Reg
#> 1           27        12         21
#> 2           15        17         24
#> 3            6        12         29
#> 4           13        11         20
#> 5            6         7          8
#> 6           14        11         17
```

## Correlation network

The simplest approach estimates pairwise Pearson correlations.

``` r
net_cor <- build_network(srl_strategies, method = "cor")
net_cor
#> Correlation Network [undirected]
#>   Sample size: 250
#>   Weights: [-0.130, 0.485]  |  +26 / -10 edges
#> 
#>   Weight matrix:
#>                Planning Monitoring Evaluating Elaboration Organization Rehearsal
#>   Planning        0.000      0.423      0.358      -0.096       -0.083    -0.019
#>   Monitoring      0.423      0.000      0.485       0.195        0.028     0.132
#>   Evaluating      0.358      0.485      0.000       0.077        0.313     0.076
#>   Elaboration    -0.096      0.195      0.077       0.000        0.432     0.341
#>   Organization   -0.083      0.028      0.313       0.432        0.000     0.339
#>   Rehearsal      -0.019      0.132      0.076       0.341        0.339     0.000
#>   Help_Seeking   -0.108     -0.116      0.023       0.008        0.123    -0.130
#>   Time_Mgmt       0.285      0.015      0.079      -0.033        0.085    -0.106
#>   Effort_Reg     -0.010     -0.008      0.250       0.050        0.135     0.029
#>                Help_Seeking Time_Mgmt Effort_Reg
#>   Planning           -0.108     0.285     -0.010
#>   Monitoring         -0.116     0.015     -0.008
#>   Evaluating          0.023     0.079      0.250
#>   Elaboration         0.008    -0.033      0.050
#>   Organization        0.123     0.085      0.135
#>   Rehearsal          -0.130    -0.106      0.029
#>   Help_Seeking        0.000     0.209      0.176
#>   Time_Mgmt           0.209     0.000      0.467
#>   Effort_Reg          0.176     0.467      0.000
```

## Partial correlation network

Partial correlations control for all other variables, revealing direct
associations.

``` r
net_pcor <- build_network(srl_strategies, method = "pcor")
net_pcor
#> Partial Correlation Network (unregularised) [undirected]
#>   Sample size: 250
#>   Weights: [-0.235, 0.502]  |  +21 / -15 edges
#> 
#>   Weight matrix:
#>                Planning Monitoring Evaluating Elaboration Organization Rehearsal
#>   Planning        0.000      0.268      0.283      -0.103       -0.146     0.046
#>   Monitoring      0.268      0.000      0.432       0.274       -0.213     0.095
#>   Evaluating      0.283      0.432      0.000      -0.156        0.406    -0.098
#>   Elaboration    -0.103      0.274     -0.156       0.000        0.380     0.181
#>   Organization   -0.146     -0.213      0.406       0.380        0.000     0.274
#>   Rehearsal       0.046      0.095     -0.098       0.181        0.274     0.000
#>   Help_Seeking   -0.121     -0.054      0.039       0.004        0.102    -0.144
#>   Time_Mgmt       0.397     -0.007     -0.207      -0.031        0.158    -0.137
#>   Effort_Reg     -0.235     -0.099      0.330       0.049       -0.092     0.096
#>                Help_Seeking Time_Mgmt Effort_Reg
#>   Planning           -0.121     0.397     -0.235
#>   Monitoring         -0.054    -0.007     -0.099
#>   Evaluating          0.039    -0.207      0.330
#>   Elaboration         0.004    -0.031      0.049
#>   Organization        0.102     0.158     -0.092
#>   Rehearsal          -0.144    -0.137      0.096
#>   Help_Seeking        0.000     0.161      0.051
#>   Time_Mgmt           0.161     0.000      0.502
#>   Effort_Reg          0.051     0.502      0.000
```

## Regularized network (EBICglasso)

The graphical lasso applies L1 regularization to the precision matrix,
producing a sparse network. The `gamma` parameter controls sparsity via
EBIC model selection — higher values yield sparser networks.

``` r
net_glasso <- build_network(srl_strategies, method = "glasso",
                            params = list(gamma = 0.5))
net_glasso
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
```

## Predictability

Node predictability measures how well each node is predicted by its
neighbors (R-squared from the network structure).

``` r
pred <- predictability(net_glasso)
round(pred, 3)
#>     Planning   Monitoring   Evaluating  Elaboration Organization    Rehearsal 
#>        0.251        0.316        0.332        0.241        0.279        0.161 
#> Help_Seeking    Time_Mgmt   Effort_Reg 
#>        0.051        0.274        0.252
```

## Bootstrap inference

Non-parametric bootstrap assesses edge stability and significance.

``` r
boot <- boot_glasso(net_glasso, iter = 1000,
                    centrality = c("strength", "expected_influence"),
                    seed = 42)
```

### Edge significance

``` r
summary(boot, type = "edges")
#>                            edge     weight      ci_lower   ci_upper inclusion
#> 36      Time_Mgmt -- Effort_Reg 0.32491515  0.2402787627 0.49995578     1.000
#> 3      Monitoring -- Evaluating 0.30049766  0.2135739099 0.44420273     1.000
#> 10  Elaboration -- Organization 0.26232453  0.1759094880 0.40494836     1.000
#> 1        Planning -- Monitoring 0.22902145  0.1438835150 0.35852942     1.000
#> 14     Elaboration -- Rehearsal 0.15760204  0.0554857192 0.28256604     0.998
#> 9    Evaluating -- Organization 0.15443135  0.0724680604 0.38310566     0.998
#> 15    Organization -- Rehearsal 0.15232840  0.0636729072 0.29191756     0.998
#> 22        Planning -- Time_Mgmt 0.13034652  0.0414121518 0.37922397     0.989
#> 2        Planning -- Evaluating 0.12953651  0.0535907183 0.29342854     0.994
#> 31     Evaluating -- Effort_Reg 0.09309602  0.0000000000 0.30523777     0.953
#> 28    Help_Seeking -- Time_Mgmt 0.06512062  0.0000000000 0.22045877     0.883
#> 5     Monitoring -- Elaboration 0.03881843  0.0000000000 0.27775412     0.842
#> 35   Help_Seeking -- Effort_Reg 0.01956134  0.0000000000 0.16452337     0.639
#> 4       Planning -- Elaboration 0.00000000 -0.1765238132 0.00000000     0.569
#> 6     Evaluating -- Elaboration 0.00000000 -0.1500953615 0.00000000     0.188
#> 7      Planning -- Organization 0.00000000 -0.1700927401 0.00000000     0.597
#> 8    Monitoring -- Organization 0.00000000 -0.2260417137 0.00000000     0.358
#> 11        Planning -- Rehearsal 0.00000000 -0.0422651527 0.03811323     0.119
#> 12      Monitoring -- Rehearsal 0.00000000  0.0000000000 0.14113867     0.465
#> 13      Evaluating -- Rehearsal 0.00000000 -0.0797648695 0.02565436     0.141
#> 16     Planning -- Help_Seeking 0.00000000 -0.1689276000 0.00000000     0.548
#> 17   Monitoring -- Help_Seeking 0.00000000 -0.1325894747 0.00000000     0.525
#> 18   Evaluating -- Help_Seeking 0.00000000  0.0000000000 0.07943598     0.109
#> 19  Elaboration -- Help_Seeking 0.00000000 -0.0438698021 0.04286525     0.116
#> 20 Organization -- Help_Seeking 0.00000000  0.0000000000 0.17705659     0.614
#> 21    Rehearsal -- Help_Seeking 0.00000000 -0.2120800796 0.00000000     0.644
#> 23      Monitoring -- Time_Mgmt 0.00000000 -0.0965496487 0.00000000     0.225
#> 24      Evaluating -- Time_Mgmt 0.00000000 -0.1745304366 0.00000000     0.207
#> 25     Elaboration -- Time_Mgmt 0.00000000 -0.0708145578 0.00000000     0.146
#> 26    Organization -- Time_Mgmt 0.00000000  0.0000000000 0.14618494     0.306
#> 27       Rehearsal -- Time_Mgmt 0.00000000 -0.1601119497 0.00000000     0.527
#> 29       Planning -- Effort_Reg 0.00000000 -0.2377412472 0.00000000     0.503
#> 30     Monitoring -- Effort_Reg 0.00000000 -0.1316710059 0.00000000     0.308
#> 32    Elaboration -- Effort_Reg 0.00000000 -0.0001367775 0.07300666     0.164
#> 33   Organization -- Effort_Reg 0.00000000  0.0000000000 0.09380279     0.336
#> 34      Rehearsal -- Effort_Reg 0.00000000  0.0000000000 0.09946240     0.190
```

### Centrality stability

``` r
summary(boot, type = "centrality")
#> $strength
#>           node      value   ci_lower  ci_upper
#> 1     Planning 0.48890448 0.37253704 1.4995496
#> 2   Monitoring 0.56833754 0.44124222 1.4631744
#> 3   Evaluating 0.67756154 0.53015851 1.6425503
#> 4  Elaboration 0.45874500 0.33341348 1.1808440
#> 5 Organization 0.56908428 0.44195673 1.5442505
#> 6    Rehearsal 0.30993044 0.22306635 0.9836351
#> 7 Help_Seeking 0.08468196 0.01772335 0.8334546
#> 8    Time_Mgmt 0.52038228 0.38168367 1.4554888
#> 9   Effort_Reg 0.43757251 0.32045449 1.2261947
#> 
#> $expected_influence
#>           node      value    ci_lower  ci_upper
#> 1     Planning 0.48890448  0.20691831 0.6142270
#> 2   Monitoring 0.56833754  0.39864460 0.8304885
#> 3   Evaluating 0.67756154  0.53015851 1.1390458
#> 4  Elaboration 0.45874500  0.28614495 0.7207834
#> 5 Organization 0.56908428  0.42522959 0.9582501
#> 6    Rehearsal 0.30993044  0.09868951 0.4621349
#> 7 Help_Seeking 0.08468196 -0.13814207 0.2614040
#> 8    Time_Mgmt 0.52038228  0.36555268 0.8701167
#> 9   Effort_Reg 0.43757251  0.30571131 0.7375021
```
