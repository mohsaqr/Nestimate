# Nested sequences and the permutation test

Do high and low achievers regulate their group work differently, or do
their teams? In `group_regulation_long`, 2,000 students work in 200
teams of 10, and the achievement level is given to the team, not to the
student. A comparison of High and Low is therefore a comparison of
teams. The students are nested in teams. The permutation test treats
each student’s sequence as an independent unit, so a nesting effect,
teammates behaving alike, would enter the comparison. This tutorial runs
the comparison, checks for a nesting effect, and shows how to report it.

``` r

achievers <- build_network(group_regulation_long, method = "relative",
                           actor = "Actor", action = "Action",
                           time = "Time", group = "Achiever")
achievers
```

    ## Group Networks (2 groups, group_col: Achiever)
    ## 
    ##   Group  Nodes  Edges  Weights
    ##   High   9      76     [0.001, 0.576]
    ##   Low    9      75     [0.000, 0.462]

[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
builds one transition network per achievement level from the students’
action sequences. Both networks have the same nine regulation actions
and 75 to 76 non-zero transitions.

## The ordinary test

``` r

permutation(achievers, iter = 1000, seed = 1)
```

    ## Grouped Permutation Test
    ## Groups: High vs Low 
    ## 
    ## -- High vs Low --
    ## Permutation Test: Transition Network (relative probabilities) [directed]
    ##   Iterations: 1000  |  Alpha: 0.05
    ##   Nodes: 9  |  Edges tested: 78  |  Significant: 42
    ##   Global test (networks differ overall?): M = 2.612 (p = 0.000999)  |  S = 0.210 (p = 0.000999)

[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
pools the 2,000 sequences, reassigns them to High and Low 1,000 times,
and compares the observed differences with the reassigned ones. Two
kinds of result are reported. The edge tests ask, for each transition,
whether the High minus Low difference is larger than the reassignments
produce: 42 of 78 transitions are. The global test asks whether the two
networks differ at all, with two statistics from the Network Comparison
Test (van Borkulo et al., 2023): M, the sum of all absolute edge
differences, and S, the largest one. Both have p = 0.000999, the
smallest value 1,000 reassignments can give: no reassignment produced a
difference as large as the observed one.

The test assumes that any student could have been in either group. In
this design only whole teams could, so the next step accounts for the
nesting.

## Reassigning whole teams

``` r

permutation(achievers, iter = 1000, seed = 1, actor = "Group")
```

    ## Grouped Permutation Test
    ## Groups: High vs Low 
    ## 
    ## -- High vs Low --
    ## Permutation Test: Transition Network (relative probabilities) [directed]
    ##   Iterations: 1000  |  Alpha: 0.05  |  Actor: Group (200 actors)
    ##   Nodes: 9  |  Edges tested: 78  |  Significant: 42
    ##   Global test (networks differ overall?): M = 2.612 (p = 0.000999)  |  S = 0.210 (p = 0.000999)
    ##   Nesting in Group: ICC = -0.002 [95% CI -0.006, 0.002]  |  between design
    ##   Design effect (1 = nesting does not matter): edges 1.03  |  global 1.20

`actor = "Group"` tells the test that the sequences are nested in teams.
Each reassignment now moves all ten students of a team together, which
is how the achievement levels were actually given. The observed
differences are the same; only the reference distribution changes.

The edge tests and the global test give the same result as with
student-level reassignment: 42 significant transitions, p = 0.000999 for
M and S. The output also reports the intraclass correlation, ICC =
−0.002 (95% CI −0.006 to 0.002), and the design effect, 1.03 for the
transitions and 1.20 for M. Both are defined and interpreted in the next
section.

A low ICC means that sequences from the same team are about as different
from each other as sequences from different teams: the observations are
close to independent, and an analysis that ignores the nesting gives
nearly the same results as one that accounts for it. Accounting for
nesting becomes necessary when the ICC is not negligible, and more so
when units are large, because the design effect grows with both:
approximately 1 + (m − 1) × ICC for units of size m (Kish, 1965). A
commonly cited rule of thumb is that design effects below 2 do not lead
to seriously misleading inferences (Muthén & Satorra, 1995). Here the
ICC is −0.002 and the design effect 1.03, so accounting for the nesting
was not necessary, and the results with and without it were almost the
same.

## The nesting diagnostics

The intraclass correlation (ICC) is the proportion of the total variance
that lies between units rather than within them (Shrout & Fleiss, 1979):

``` math
\text{ICC} = \frac{\sigma^2_{\text{between}}}{\sigma^2_{\text{between}} + \sigma^2_{\text{within}}}.
```

An ICC of 0 means that members of the same unit are no more alike than
members of different units; an ICC of 1 means that all members of a unit
are identical and all variation lies between units. A value in between
is the proportion of variation due to the unit. In
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
the ICC is computed for each transition: for every student, the share of
their transitions that is, for example, `plan` to `monitor`; the one-way
ANOVA ICC of that share across teams is computed within High and within
Low. The transition ICCs are averaged with weights equal to transition
frequency. The 95% interval comes from the jackknife, leaving out one
team at a time (Efron & Tibshirani, 1993).

The design effect is the ratio of the variance of an estimate under the
nested design to its variance had the units been sampled independently
(Kish, 1965). Here it is the variance of the reassigned differences when
whole teams move, divided by the variance when single students move. A
value of 1 means both reassignments give the same chance variation.

[`permutation_diagnostics()`](https://pak.dynasite.org/Nestimate/reference/permutation_diagnostics.md)
runs the ordinary and the team-level test on the same data and reports
both with these two quantities, one row per compared pair:

``` r

permutation_diagnostics(achievers, actor = "Group", iter = 1000, seed = 1)
```

    ##          pair n_sequences n_actors  design          icc icc_ci_lower
    ## 1 High vs Low        2000      200 between -0.001697397 -0.005641197
    ##   icc_ci_upper deff_edges deff_global p_global_sequence p_global_actor
    ## 1  0.002246402   1.026921    1.200797       0.000999001    0.000999001
    ##   sig_edges_sequence sig_edges_actor edges_changed min_p_actor
    ## 1                 42              42             2 0.000999001

The comparison involves 2,000 sequences in 200 teams. The design is
`between`: every team belongs to one achievement level (`within` would
mean units in both groups, as in repeated measures; `mixed`, some of
each). The ICC is −0.002 with a 95% interval from −0.006 to 0.002, which
indicates little evidence of a nesting effect; a value slightly below
zero is sampling variation around zero. The design effect is 1.03 for
the transitions (median over edges) and 1.20 for the global statistic M.
The global test has p = 0.000999 under both reassignments, and 42
transitions are significant under both. Two transitions change their
conclusion. The smallest p-value the team-level test can produce,
`min_p_actor`, is 0.000999: with 200 teams the number of iterations sets
the limit, not the number of teams.

With `level = "edges"` the same comparison is reported per transition:

``` r

head(permutation_diagnostics(achievers, actor = "Group", iter = 1000,
                             level = "edges", seed = 1), 6)
```

    ##          pair  from         to         diff          icc null_sd_sequence
    ## 1 High vs Low adapt   cohesion -0.014762566  0.010118088       0.03863244
    ## 2 High vs Low adapt  consensus  0.055773975 -0.003796296       0.04560409
    ## 3 High vs Low adapt coregulate -0.029891304  0.002498496       0.01294589
    ## 4 High vs Low adapt    discuss -0.032473790 -0.009052931       0.02104505
    ## 5 High vs Low adapt    emotion  0.030430928 -0.007376103       0.02877059
    ## 6 High vs Low adapt    monitor -0.006957293 -0.003562251       0.01591287
    ##   null_sd_actor      deff p_sequence    p_actor changed
    ## 1    0.04129895 1.1428088 0.68831169 0.72927073   FALSE
    ## 2    0.04255012 0.8705504 0.21378621 0.18981019   FALSE
    ## 3    0.01383698 1.1424014 0.01498501 0.02397602   FALSE
    ## 4    0.01902641 0.8173606 0.11588412 0.09190809   FALSE
    ## 5    0.02727781 0.8989211 0.32267732 0.26373626   FALSE
    ## 6    0.01553230 0.9527393 0.66633367 0.66733267   FALSE

`diff` is the High minus Low transition probability and is the same
under both tests. `icc` is the ICC of that transition.
`null_sd_sequence` and `null_sd_actor` are the spread of the reassigned
differences under student-level and team-level reassignment, and `deff`
is the squared ratio of the two. `p_sequence` and `p_actor` are the two
edge p-values, and `changed` marks a transition that is significant at
0.05 under one test but not the other. In the six rows shown, the ICCs
lie between −0.009 and 0.010 and the design effects between 0.82 and
1.14; the two spreads are close, so the two p-values are close. `adapt`
to `coregulate` has p = 0.015 and p = 0.024 and is significant under
both.

Over all 78 transitions, two change at this seed (`coregulate` to
`adapt`, `emotion` to `plan`), with p-values between 0.047 and 0.072
under both tests. With seeds 2 to 4 a different set of two or three
transitions changes, all with p-values between 0.04 and 0.07. These
changes come from the finite number of reassignments, not from the
nesting. The ICC indicates little evidence of a nesting effect, and the
High–Low difference does not depend on the nesting of students in teams.

## Reporting

> Because the data were nested (students in teams; with repeated
> measurements, sessions in students), the two groups were compared with
> a permutation test that accounted for the nesting: in each of 1,000
> permutations, whole teams rather than single sequences were reassigned
> between the groups, so that sequences from the same team stayed
> together (Good, 2005). Units that belonged to one group were
> reassigned as a whole, and units that contributed to both groups, as
> in a repeated-measures design, had their sequences reassigned within
> the unit. The nesting effect was quantified with the intraclass
> correlation (ICC; Shrout & Fleiss, 1979) and the design effect, the
> ratio of the variance under the nested design to the variance under
> independent sampling (Kish, 1965). The ICC indicated little evidence
> of a nesting effect (ICC = −0.002, 95% CI −0.006 to 0.002; design
> effect 1.03 across 200 teams). The networks of high and low achievers
> differed overall (M = 2.61, p \< .001), and 42 of 78 transitions
> differed at p \< .05.

## When to use which

- [`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md):
  one sequence per student, with no grouping above the student.
- `permutation(actor = )`: several sequences per student, or students in
  teams, classes or other units. It handles units in one group, units in
  both groups, and mixtures, and reports the ICC and design effect.
- `paired = TRUE`: exactly one sequence per student in each of two
  conditions.
- [`permutation_diagnostics()`](https://pak.dynasite.org/Nestimate/reference/permutation_diagnostics.md):
  to show next to each other what the ordinary and the team-level test
  conclude.
- `actor` applies to transition networks (`relative`, `frequency`,
  `co_occurrence`). Correlation networks do not keep the unit
  identifiers after estimation.

## References

Kish, L. (1965). *Survey Sampling*. Wiley.

Muthén, B. O., & Satorra, A. (1995). Complex sample data in structural
equation modeling. *Sociological Methodology*, 25, 267–316.

Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in
assessing rater reliability. *Psychological Bulletin*, 86(2), 420–428.

van Borkulo, C. D., van Bork, R., Boschloo, L., Kossakowski, J. J., Tio,
P., Schoevers, R. A., Borsboom, D., & Waldorp, L. J. (2023). Comparing
network structures on three aspects: A permutation test. *Psychological
Methods*, 28(6), 1273–1285.

Efron, B., & Tibshirani, R. J. (1993). *An Introduction to the
Bootstrap*. Chapman & Hall.

Good, P. (2005). *Permutation, Parametric, and Bootstrap Tests of
Hypotheses* (3rd ed.). Springer.
