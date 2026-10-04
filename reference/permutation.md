# Permutation Test for Network Comparison

Tests whether two networks estimated by
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md)
differ more than chance would produce. The sequences (or rows) of both
networks are pooled, the group labels are shuffled `iter` times, both
networks are re-estimated on every shuffle, and the observed differences
are compared with the shuffled ones. Works with every built-in method
and with registered custom estimators.

## Usage

``` r
permutation(
  x,
  y = NULL,
  iter = 1000L,
  alpha = 0.05,
  paired = FALSE,
  adjust = "none",
  measures = NULL,
  nlambda = 50L,
  seed = NULL,
  actor = NULL
)
```

## Arguments

- x:

  A `netobject` (from
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md))
  or a
  [`net_edge_betweenness`](https://saqr.me/Nestimate/reference/net_edge_betweenness.md)
  object.

- y:

  A `netobject` (from
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md))
  or a
  [`net_edge_betweenness`](https://saqr.me/Nestimate/reference/net_edge_betweenness.md)
  object. Must use the same method and have the same nodes as `x`.
  Default `NULL`: when `x` is a `netobject_group` (or an `mcml`) and `y`
  is left `NULL`, every pair of groups is tested and the result is a
  `net_permutation_group` named `"<group i> vs <group j>"`.

- iter:

  Integer. Number of permutation iterations (default: 1000).

- alpha:

  Numeric. Significance level (default: 0.05).

- paired:

  Logical. If `TRUE`, permute within pairs (requires equal number of
  observations in `x` and `y`). Default: FALSE.

- adjust:

  Character. p-value adjustment method passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) (default:
  `"none"`). Common choices: `"holm"`, `"BH"`, `"bonferroni"`.

- measures:

  Character vector of centrality measures to permutation-test in
  addition to the edges, or `"all"` for every built-in measure. Default
  `NULL` (edges only). When supplied, the result gains a `$centralities`
  block matching the layout of `tna::permutation_test(measures = )`: per
  state and measure it reports the observed difference, an effect size
  (difference / SD of the permutation null), and a permutation p-value,
  all using the same permuted networks as the edge test. Not supported
  for `net_edge_betweenness` inputs.

- nlambda:

  Integer. Number of lambda values for the EBIC-glasso regularisation
  path (only used when `method = "glasso"`). Higher values give finer
  lambda resolution at the cost of speed. Default: 50.

- seed:

  Integer or NULL. RNG seed for reproducibility.

- actor:

  Character or NULL. Name of the column identifying the actor each
  sequence belongs to: the person whose sessions they are, or the team
  of a student. Looked up in the network's `$metadata` (e.g.
  `"student_id"` when sessions are nested in students, `"Group"` for
  students nested in teams) or in its wide sequence data. When supplied,
  whole actors are permuted and the ICC and design effect are reported;
  see the sections *Nested data and actor* and *ICC and design effect*.
  Supported for transition methods (`"relative"`, `"frequency"`,
  `"co_occurrence"`); cannot be combined with `paired = TRUE`, which is
  the special case of one actor per pair. Default `NULL`: sequences are
  permuted individually.

## Value

An object of class `"net_permutation"` containing:

- x:

  The first `netobject`.

- y:

  The second `netobject`.

- diff:

  Observed difference matrix (`x - y`).

- diff_sig:

  Observed difference where `p < alpha`, else 0.

- p_values:

  P-value matrix (adjusted if `adjust != "none"`).

- effect_size:

  Effect size matrix (observed diff / SD of permutation diffs).

- summary:

  Long-format data frame, one row per edge present in either network
  (undirected networks keep one row per unordered pair), with columns
  `from`, `to`, `weight_x`, `weight_y`, `diff`, `effect_size`,
  `p_value`, `sig`.

- global:

  Data frame of the two NCT-style global statistics, one row each:
  `statistic` (`"M"`, the sum of absolute edge differences, and `"S"`,
  the largest absolute edge difference), `observed`, and `p_value` from
  the same permutation null as the edge test. Absent on the
  edge-betweenness path.

- method:

  The network estimation method.

- source_method:

  For edge-betweenness tests, the source network method.

- iter:

  Number of permutation iterations.

- alpha:

  Significance level used.

- paired:

  Whether paired permutation was used.

- adjust:

  p-value adjustment method used.

- actor:

  The `actor` column name, or `NULL`.

- n_actors:

  Number of distinct actors, or `NULL`.

- null_sd:

  Matrix of the SD of each edge difference over the permutation null
  (the effect-size denominator).

- null_sd_m:

  SD of the global `M` statistic over the permutation null. Absent on
  the edge-betweenness path.

- clustering:

  Present only with `actor`. One-row data frame: `n_sequences`,
  `n_actors`, `design` (`"between"`, `"within"`, `"mixed"`), `icc` with
  `icc_ci_lower`/`icc_ci_upper` (how alike sequences of one actor are;
  see
  [`permutation_diagnostics`](https://saqr.me/Nestimate/reference/permutation_diagnostics.md)),
  `deff_edges` (median over edges) and `deff_global` (for `M`): the
  actor-level over the sequence-level null variance, drawn in the same
  run (the design effect; Kish, 1965). `min_p` is the smallest
  attainable p-value.

- clustering_edges:

  Present only with `actor`. One row per edge of `summary`: `from`,
  `to`, `icc`, `null_sd_actor`, `null_sd_sequence`, `deff`.

- centralities:

  Present only when `measures` is supplied. A list with `stats` (one row
  per state-by-measure: `state`, `centrality`, `diff_true`,
  `effect_size`, `p_value`), `diffs_true` (wide observed differences),
  and `diffs_sig` (observed differences where `p < alpha`, else 0).

Grouped input returns a `"net_permutation_group"` (a named list of
`net_permutation` results): one element per matching group name when
both `x` and `y` are `netobject_group`s, or one per group pair when `y`
is `NULL`. Two `wtna_mixed` inputs return a `"wtna_perm_mixed"` with
`$transition` and `$cooccurrence` results.

## What is tested

Two kinds of question are answered from the same shuffles.

- Edge tests:

  One test per edge: is the difference in this edge's weight, `x - y`,
  larger than the shuffles produce? Reported in
  [`summary()`](https://rdrr.io/r/base/summary.html) with an effect size
  (observed difference divided by the SD of the shuffled differences)
  and a p-value
  `(number of shuffles at least as extreme + 1) / (iter + 1)`. With many
  edges, some fall below `alpha` by chance; use `adjust` to correct for
  that.

- Global test:

  One test for the whole network: do the two networks differ at all? Two
  statistics, as in the Network Comparison Test (van Borkulo et al.,
  2023): **M**, the sum of the absolute edge differences (the total
  amount of difference), and **S**, the largest absolute edge
  difference. Being a single test, it needs no multiplicity correction.
  It is shown by [`print()`](https://rdrr.io/r/base/print.html).

The smallest attainable p-value is `1 / (iter + 1)`; with the default
`iter = 1000` it is 0.000999, meaning no shuffle came close.

## Nested data and `actor`

Shuffling treats every sequence as an exchangeable unit. When several
sequences come from the same actor (sessions of one person, students of
one team), `actor` names the column identifying that actor. The shuffle
then respects it (Good, 2005; Anderson & ter Braak, 2003): a person
whose sequences are all in one group moves to the other group as a
whole; a person with sequences in both groups has their labels shuffled
among their own sequences only. Mixed designs combine the two. The
observed differences do not change; only the p-values and effect sizes
do. With few persons there are few distinct ways to shuffle them, and a
warning (class `nestimate_few_actors`) is raised when no p-value could
fall below `alpha`.

`actor` is available for transition networks (`"relative"`,
`"frequency"`, `"co_occurrence"`). Association networks (`"cor"`,
`"pcor"`, `"glasso"`, ...) do not keep row identifiers after estimation
and raise `nestimate_actor_unsupported`.

## ICC and design effect

With `actor`, [`print()`](https://rdrr.io/r/base/print.html) also
reports:

- ICC:

  The intraclass correlation, the proportion of the total variance that
  lies between actors (Shrout & Fleiss, 1979). An ICC close to 0
  indicates little evidence of a nesting effect. Computed as the one-way
  ANOVA ICC of each sequence's transition shares within each group,
  averaged over edges weighted by edge frequency, jackknife
  bias-corrected, with a 95% interval from the leave-one-actor-out
  jackknife (Efron & Tibshirani, 1993).

- Design effect:

  The ratio of the variance of an estimate under the clustered design to
  its variance had the units been sampled independently (Kish, 1965).
  Here: the variance of the shuffled differences when whole actors are
  moved, divided by the variance when single sequences are moved, both
  drawn in the same run. Reported as the median over edges and for the
  global statistic M. For equal numbers of sequences per actor m, Kish
  gives the approximation `1 + (m - 1) * ICC`.

## Reading the printed output


    Permutation Test: Transition Network (relative probabilities) [directed]
      Iterations: 1000  |  Alpha: 0.05  |  Actor: Group (200 actors)
      Nodes: 9  |  Edges tested: 78  |  Significant: 42
      Global test (networks differ overall?): M = 2.612 (p = 0.000999)  |  ...
      Nesting in Group: ICC = -0.002 [95% CI -0.006, 0.002]  |  between design
      Design effect (1 = nesting does not matter): edges 1.03  |  global 1.20

Line 2: settings, and the actor column with its number of actors. Line
3: edges present in either network and how many differ at `alpha`. Line
4: the global test. Lines 5-6, only with `actor`: the ICC with its
interval, whether actors sit in one group (`between`), in both
(`within`) or either (`mixed`), and the design effects.

## Other inputs

For transition methods, per-sequence count matrices are computed once
and each shuffle only re-sums them, which keeps large `iter` fast. For
association methods the estimator is re-run on every shuffle. If a
transition network rests on a single sequence, a warning (class
`nestimate_single_sequence`) says it cannot be validated by resampling.

`permutation()` also accepts two
[`net_edge_betweenness`](https://saqr.me/Nestimate/reference/net_edge_betweenness.md)
objects. It then permutes the source networks, recomputes edge
betweenness for each shuffle, and tests the edge-betweenness
differences. Both objects must come from the same source method and use
the same `invert` setting.

## References

Anderson, M. J., & ter Braak, C. J. F. (2003). Permutation tests for
multi-factorial analysis of variance. *Journal of Statistical
Computation and Simulation*, 73(2), 85-113.

Efron, B., & Tibshirani, R. J. (1993). *An Introduction to the
Bootstrap*. Chapman & Hall.

Good, P. (2005). *Permutation, Parametric, and Bootstrap Tests of
Hypotheses* (3rd ed.). Springer.

Kish, L. (1965). *Survey Sampling*. Wiley.

Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in
assessing rater reliability. *Psychological Bulletin*, 86(2), 420-428.

van Borkulo, C. D., van Bork, R., Boschloo, L., Kossakowski, J. J., Tio,
P., Schoevers, R. A., Borsboom, D., & Waldorp, L. J. (2023). Comparing
network structures on three aspects: A permutation test. *Psychological
Methods*, 28(6), 1273-1285.

## See also

[`permutation_diagnostics`](https://saqr.me/Nestimate/reference/permutation_diagnostics.md)
to compare the actor-level and ordinary tests side by side;
[`bayes_compare`](https://saqr.me/Nestimate/reference/bayes_compare.md)
for the Bayesian complement: instead of "is this difference more extreme
than chance?" it answers "how probable is a difference, and how large?";
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md),
[`bootstrap_network`](https://saqr.me/Nestimate/reference/bootstrap_network.md),
[`print.net_permutation`](https://saqr.me/Nestimate/reference/print.net_permutation.md),
[`summary.net_permutation`](https://saqr.me/Nestimate/reference/summary.net_permutation.md)

## Examples

``` r
s1 <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
s2 <- data.frame(V1 = c("A","C","B"), V2 = c("C","B","A"))
n1 <- build_network(s1, method = "relative")
n2 <- build_network(s2, method = "relative")
perm <- permutation(n1, n2, iter = 10)
# \donttest{
set.seed(1)
d1 <- data.frame(V1 = sample(LETTERS[1:4], 20, TRUE),
                 V2 = sample(LETTERS[1:4], 20, TRUE),
                 V3 = sample(LETTERS[1:4], 20, TRUE))
d2 <- data.frame(V1 = sample(LETTERS[1:4], 20, TRUE),
                 V2 = sample(LETTERS[1:4], 20, TRUE),
                 V3 = sample(LETTERS[1:4], 20, TRUE))
net1 <- build_network(d1, method = "relative")
net2 <- build_network(d2, method = "relative")
perm <- permutation(net1, net2, iter = 100, seed = 42)
print(perm)
#> Permutation Test: Transition Network (relative probabilities) [directed]
#>   Iterations: 100  |  Alpha: 0.05
#>   Nodes: 4  |  Edges tested: 16  |  Significant: 3
#>   Global test (networks differ overall?): M = 3.858 (p = 0.0297)  |  S = 0.567 (p = 0.109)
summary(perm)
#>    from to   weight_x  weight_y        diff effect_size    p_value   sig
#> 1     A  A 0.33333333 0.2222222  0.11111111   0.6504300 0.49504950 FALSE
#> 2     A  B 0.06666667 0.3333333 -0.26666667  -1.9871473 0.04950495  TRUE
#> 3     A  C 0.20000000 0.2222222 -0.02222222  -0.1329339 0.93069307 FALSE
#> 4     A  D 0.40000000 0.2222222  0.17777778   1.0449465 0.25742574 FALSE
#> 5     B  A 0.30769231 0.3333333 -0.02564103  -0.1106148 0.98019802 FALSE
#> 6     B  B 0.38461538 0.0000000  0.38461538   1.9215279 0.03960396  TRUE
#> 7     B  C 0.23076923 0.2500000 -0.01923077  -0.1103038 0.99009901 FALSE
#> 8     B  D 0.07692308 0.4166667 -0.33974359  -1.9910325 0.07920792 FALSE
#> 9     C  A 0.33333333 0.1111111  0.22222222   1.0147102 0.38613861 FALSE
#> 10    C  B 0.66666667 0.3333333  0.33333333   1.3787271 0.24752475 FALSE
#> 11    C  C 0.00000000 0.2222222 -0.22222222  -1.3964532 0.26732673 FALSE
#> 12    C  D 0.00000000 0.3333333 -0.33333333  -1.8265261 0.12871287 FALSE
#> 13    D  A 0.00000000 0.5000000 -0.50000000  -1.8205956 0.10891089 FALSE
#> 14    D  B 0.00000000 0.2000000 -0.20000000  -0.9828935 0.45544554 FALSE
#> 15    D  C 0.33333333 0.2000000  0.13333333   0.5436513 0.56435644 FALSE
#> 16    D  D 0.66666667 0.1000000  0.56666667   2.5117967 0.01980198  TRUE

# Students are nested in teams, and Achiever is a team-level label:
# permute whole teams, not single students
net <- build_network(group_regulation_long, method = "relative",
                     actor = "Actor", action = "Action", time = "Time",
                     group = "Achiever")
permutation(net, iter = 100, actor = "Group", seed = 1)
#> Grouped Permutation Test
#> Groups: High vs Low 
#> 
#> -- High vs Low --
#> Permutation Test: Transition Network (relative probabilities) [directed]
#>   Iterations: 100  |  Alpha: 0.05  |  Actor: Group (200 actors)
#>   Nodes: 9  |  Edges tested: 78  |  Significant: 43
#>   Global test (networks differ overall?): M = 2.612 (p = 0.0099)  |  S = 0.210 (p = 0.0099)
#>   Nesting in Group: ICC = -0.002 [95% CI -0.006, 0.002]  |  between design
#>   Design effect (1 = nesting does not matter): edges 1.03  |  global 1.26
# }
```
