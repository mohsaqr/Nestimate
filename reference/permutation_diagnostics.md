# Does Nesting Bias a Permutation Test?

`permutation_diagnostics()` shows what treating nested sequences as
independent costs a network comparison. It runs
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
twice on the same networks and seed: once shuffling single sequences,
once with `actor`, shuffling whole actors (the persons whose sessions
the sequences are, or the teams of students). It returns the two results
side by side, one row per compared pair or per edge, with the intraclass
correlation (ICC) and design effect that explain any difference between
them. It takes two transition networks, or every pair of groups of a
`netobject_group`.

A design effect near 1 and few changed edges say that the nesting
changes the comparison little. Where the two tests differ, the
actor-level one is the test that follows the nesting, as
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
sets out for nested data. The sections *Nested data and actor* and *ICC
and design effect* of
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
define these quantities, and *How to read the result* below says what
each column shows.

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

  The networks to compare: a `netobject_group`, whose every pair of
  groups is diagnosed, or a `netobject`, compared with `y`. An `mcml`,
  whose layers hold different states, raises an error of class
  `nestimate_bad_input`. Only the transition methods `"relative"`,
  `"frequency"` and `"co_occurrence"` are supported; any other method,
  association networks and `wtna` among them, raises an error of class
  `nestimate_actor_unsupported`.

- y:

  The `netobject` to compare with when `x` is a `netobject`. Default
  `NULL`, as it must stay when `x` is a group.

- actor:

  The column identifying the actor each sequence belongs to, as in
  [`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md),
  looked up in the networks' metadata or sequence data. A column not
  found, or one with missing values, raises an error of class
  `nestimate_actor_missing`. Required.

- iter:

  The number of permutations for each of the two tests, at least 2.
  Default `1000`.

- alpha:

  The significance level that counts an edge as significant, between 0
  and 1. The p-values are not adjusted for multiplicity. Default `0.05`.

- level:

  `"overall"` (default) returns one row per compared pair; `"edges"` one
  row per edge per pair.

- seed:

  A seed for the shuffles, used by both tests, so the result repeats
  exactly. Default `NULL`.

## Value

A data frame. With `level = "overall"`, one row per compared pair:

- pair:

  `"<x> vs <y>"`, the names of the two groups (`"x vs y"` for two
  networks).

- n_sequences, n_actors:

  The sequences and distinct actors in the pair.

- design:

  `"between"` (every actor in one group), `"within"` (every actor in
  both groups) or `"mixed"`.

- icc, icc_ci_lower, icc_ci_upper:

  How alike the sequences of one actor are: the jackknife bias-corrected
  ICC, as printed by
  [`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md),
  with a 95% interval; with fewer than 3 actors the ICC is the
  uncorrected estimate and the interval is `NA`; the ICC is `NA` when no
  transition share varies.

- deff_edges:

  The median over edges of the design effect, the ratio of the
  actor-level to the sequence-level null variance of the edge difference
  (Kish, 1965).

- deff_global:

  The same ratio for the global `M` statistic.

- p_global_sequence, p_global_actor:

  The permutation p-values of `M` when sequences or whole actors are
  reassigned.

- sig_edges_sequence, sig_edges_actor:

  The number of edges with `p < alpha` under each test.

- edges_changed:

  The number of edges significant under one test but not the other.

- min_p_actor:

  The smallest p-value the actor-level test can reach: one over the
  number of distinct arrangements of the actors, or `1 / (iter + 1)`,
  whichever is larger.

With `level = "edges"`, one row per edge present in either network, with
the columns `pair`, `from`, `to`, `diff` (the observed difference,
`x - y`), `icc` (the edge's ANOVA ICC, not bias-corrected; `NA` where
the share does not vary), `null_sd_sequence` and `null_sd_actor` (the SD
of the shuffled differences under each test), `deff` (`NaN` where the
edge difference never varies under either null), `p_sequence`, `p_actor`
and `changed` (whether the edge is significant under one test only).

The ICC and design effects are those of the actor-level run, as
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
prints them when given `actor`; the p-values and significance counts
compare it with a separate ordinary run.

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

[`permutation`](https://pak.dynasite.org/Nestimate/reference/permutation.md)

## Examples

``` r
# Students are nested in teams; Achiever is a team-level label.
# iter = 50 keeps the example fast; a real analysis uses 1000 or more.
net <- build_network(group_regulation_long, method = "relative",
                     actor = "Actor", action = "Action", time = "Time",
                     group = "Achiever")
permutation_diagnostics(net, actor = "Group", iter = 50, seed = 1)
#>          pair n_sequences n_actors  design          icc icc_ci_lower
#> 1 High vs Low        2000      200 between -0.001697397 -0.005641197
#>   icc_ci_upper deff_edges deff_global p_global_sequence p_global_actor
#> 1  0.002246402   1.069592    1.186787        0.01960784     0.01960784
#>   sig_edges_sequence sig_edges_actor edges_changed min_p_actor
#> 1                 38              42             4  0.01960784
head(permutation_diagnostics(net, actor = "Group", iter = 50,
                             level = "edges", seed = 1))
#>          pair  from         to         diff          icc null_sd_sequence
#> 1 High vs Low adapt   cohesion -0.014762566  0.010118088       0.03366031
#> 2 High vs Low adapt  consensus  0.055773975 -0.003796296       0.04581149
#> 3 High vs Low adapt coregulate -0.029891304  0.002498496       0.01144616
#> 4 High vs Low adapt    discuss -0.032473790 -0.009052931       0.02089267
#> 5 High vs Low adapt    emotion  0.030430928 -0.007376103       0.02759746
#> 6 High vs Low adapt    monitor -0.006957293 -0.003562251       0.01333857
#>   null_sd_actor      deff p_sequence    p_actor changed
#> 1    0.04169614 1.5344598 0.64705882 0.68627451   FALSE
#> 2    0.04532698 0.9789597 0.21568627 0.27450980   FALSE
#> 3    0.01121761 0.9604642 0.01960784 0.01960784   FALSE
#> 4    0.01396427 0.4467334 0.15686275 0.03921569    TRUE
#> 5    0.02870211 1.0816567 0.39215686 0.23529412   FALSE
#> 6    0.01470229 1.2149298 0.74509804 0.60784314   FALSE
```
