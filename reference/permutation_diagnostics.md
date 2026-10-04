# Does Nesting Bias a Permutation Test?

Shows what treating nested sequences as independent would cost. The
comparison is run twice on the same data: once with the ordinary
[`permutation`](https://saqr.me/Nestimate/reference/permutation.md)
test, which shuffles single sequences, and once with `actor`, which
shuffles whole actors (the persons whose sessions they are, or the teams
of students). The result places the two side by side, with the ICC and
design effect that explain any difference between them. See the sections
*Nested data and actor* and *ICC and design effect* of
[`permutation`](https://saqr.me/Nestimate/reference/permutation.md) for
what these quantities mean.

## Usage

``` r
permutation_diagnostics(
  x,
  y = NULL,
  actor,
  iter = 1000L,
  alpha = 0.05,
  level = c("overall", "edges"),
  seed = NULL
)
```

## Arguments

- x:

  A `netobject_group` (every pair of groups is diagnosed) or a
  `netobject` (then `y` is required). Transition methods only
  (`"relative"`, `"frequency"`, `"co_occurrence"`).

- y:

  A `netobject` to compare with `x`, or `NULL`.

- actor:

  Character. Column identifying the actor each sequence belongs to, as
  in
  [`permutation`](https://saqr.me/Nestimate/reference/permutation.md).

- iter:

  Integer. Permutation iterations for each of the two tests. Default
  1000.

- alpha:

  Numeric. Significance level. Default 0.05.

- level:

  Character. `"overall"` (default): one row per group pair. `"edges"`:
  one row per edge per pair.

- seed:

  Integer or NULL. RNG seed; both tests use the same seed.

## Value

A `data.frame`.

With `level = "overall"`, one row per compared pair:

- pair:

  `"<x> vs <y>"`.

- n_sequences, n_actors:

  Sequences and distinct actors in the pair.

- design:

  `"between"` (every actor in one group), `"within"` (every actor in
  both groups) or `"mixed"`.

- icc, icc_ci_lower, icc_ci_upper:

  How alike the sequences of one actor are, with a 95% interval; as
  printed by
  [`permutation`](https://saqr.me/Nestimate/reference/permutation.md).
  `NA` interval with fewer than 3 actors.

- deff_edges:

  Median over edges of the design effect, the ratio of the actor-level
  to the sequence-level null variance of the edge difference (Kish,
  1965).

- deff_global:

  The same ratio for the global `M` statistic.

- p_global_sequence, p_global_actor:

  Permutation p-values of `M` when sequences or whole actors are
  reassigned.

- sig_edges_sequence, sig_edges_actor:

  Edges with `p < alpha` under each test.

- edges_changed:

  Edges significant under one test but not the other.

- min_p_actor:

  Smallest p-value an exact actor-level test can produce,
  `max(1 / arrangements, 1 / (iter + 1))`.

With `level = "edges"`, one row per edge present in either network:
`pair`, `from`, `to`, `diff`, `icc` (per-edge ANOVA ICC, not
bias-corrected; `NA` where the share does not vary), `null_sd_sequence`,
`null_sd_actor`, `deff` (`NaN` where the edge difference never varies
under either null), `p_sequence`, `p_actor`, `changed`.

The ICC and design effects are those of the actor-level run (see the
`clustering` element of
[`permutation`](https://saqr.me/Nestimate/reference/permutation.md));
the p-values and significance counts compare it with a separate ordinary
run. Errors with class `nestimate_actor_unsupported` for association
networks and `nestimate_actor_missing` when `actor` is not a column of
the networks' metadata or sequence data.

## How to read the result

- `deff_edges`, `deff_global`:

  The design effect (Kish, 1965): actor-level over sequence-level null
  variance. 1 means the two shuffles give the same chance variation;
  above 1 the actor-level one varies more, below 1 less.

- `edges_changed`:

  Edges significant under one test but not the other. Edges with
  p-values close to `alpha` can flip from Monte Carlo error alone;
  increase `iter` before reading much into one or two.

- `min_p_actor` above `alpha`:

  Too few actors: the actor-level test cannot reject anything.

## References

Efron, B., & Tibshirani, R. J. (1993). *An Introduction to the
Bootstrap*. Chapman & Hall. (jackknife, ch. 11)

Kish, L. (1965). *Survey Sampling*. Wiley. (design effect)

Anderson, M. J., & ter Braak, C. J. F. (2003). Permutation tests for
multi-factorial analysis of variance. *Journal of Statistical
Computation and Simulation*, 73(2), 85-113.

## See also

[`permutation`](https://saqr.me/Nestimate/reference/permutation.md)

## Examples

``` r
# \donttest{
# Students are nested in teams; Achiever is a team-level label
net <- build_network(group_regulation_long, method = "relative",
                     actor = "Actor", action = "Action", time = "Time",
                     group = "Achiever")
permutation_diagnostics(net, actor = "Group", iter = 200, seed = 1)
#>          pair n_sequences n_actors  design          icc icc_ci_lower
#> 1 High vs Low        2000      200 between -0.001697397 -0.005641197
#>   icc_ci_upper deff_edges deff_global p_global_sequence p_global_actor
#> 1  0.002246402   1.036832    1.285395       0.004975124    0.004975124
#>   sig_edges_sequence sig_edges_actor edges_changed min_p_actor
#> 1                 41              42             3 0.004975124
head(permutation_diagnostics(net, actor = "Group", iter = 200,
                             level = "edges", seed = 1))
#>          pair  from         to         diff          icc null_sd_sequence
#> 1 High vs Low adapt   cohesion -0.014762566  0.010118088       0.03649945
#> 2 High vs Low adapt  consensus  0.055773975 -0.003796296       0.04677298
#> 3 High vs Low adapt coregulate -0.029891304  0.002498496       0.01246761
#> 4 High vs Low adapt    discuss -0.032473790 -0.009052931       0.02085483
#> 5 High vs Low adapt    emotion  0.030430928 -0.007376103       0.03113475
#> 6 High vs Low adapt    monitor -0.006957293 -0.003562251       0.01628239
#>   null_sd_actor      deff  p_sequence    p_actor changed
#> 1    0.03996603 1.1989732 0.666666667 0.71641791   FALSE
#> 2    0.04215318 0.8122143 0.233830846 0.17910448   FALSE
#> 3    0.01370900 1.2090521 0.009950249 0.01990050   FALSE
#> 4    0.01948459 0.8729092 0.119402985 0.09452736   FALSE
#> 5    0.02814810 0.8173486 0.338308458 0.25373134   FALSE
#> 6    0.01559177 0.9169695 0.671641791 0.66666667   FALSE
# }
```
