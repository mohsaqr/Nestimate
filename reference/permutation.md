# Permutation Test for Network Comparison

`permutation()` tests whether two networks differ more than chance would
produce. The sequences (or, for association networks, the rows) behind
both networks are pooled, the group labels are shuffled `iter` times,
both networks are re-estimated on every shuffle, and each observed
difference is compared with the shuffled ones. It tests every edge, the
network as a whole and, with `measures`, the nodes' centralities.

It compares two `netobject`s from
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
estimated with the same method on the same nodes, for every built-in
method and for registered estimators whose result keeps the data
(`cleaned_data`, see
[`register_estimator()`](https://pak.dynasite.org/Nestimate/reference/register_estimator.md)).
A `netobject_group` (or an `mcml`, converted with
[`as_tna()`](https://pak.dynasite.org/Nestimate/reference/as_tna.md)) is
tested pair by pair: every pair of its groups when `y` is `NULL`, or the
groups of the same name in a second `netobject_group`. Two `wtna_mixed`
results of
[`wtna()`](https://pak.dynasite.org/Nestimate/reference/wtna.md) are
tested on their transition and co-occurrence networks, and two
[`net_edge_betweenness()`](https://pak.dynasite.org/Nestimate/reference/net_edge_betweenness.md)
results on their edge betweenness.

A small p-value says that the observed difference is larger than the
shuffles produce; it does not say how large or how probable the
difference is, which
[`bayes_compare()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
estimates. With many edges, some fall below `alpha` by chance: set
`adjust` to correct for that. The shuffle treats every sequence as
exchangeable; when sequences are nested in persons or teams, set
`actor`, and see what the nesting changes with
[`permutation_diagnostics()`](https://pak.dynasite.org/Nestimate/reference/permutation_diagnostics.md).

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
  nlambda = NULL,
  seed = NULL,
  actor = NULL
)

# S3 method for class 'net_permutation'
print(x, ...)

# S3 method for class 'net_permutation'
summary(object, ...)

# S3 method for class 'net_permutation_group'
print(x, ...)

# S3 method for class 'net_permutation_group'
summary(object, ...)

# S3 method for class 'wtna_perm_mixed'
print(x, ...)

# S3 method for class 'wtna_perm_mixed'
summary(object, ...)
```

## Arguments

- x:

  The first network: a `netobject` from
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  (a bare `cograph_network` is converted to one), a `netobject_group` or
  `mcml`, a `wtna_mixed`, or a
  [`net_edge_betweenness()`](https://pak.dynasite.org/Nestimate/reference/net_edge_betweenness.md)
  result. For the [`print()`](https://rdrr.io/r/base/print.html)
  methods: a `net_permutation`, `net_permutation_group` or
  `wtna_perm_mixed`.

- y:

  The second network, of the same kind as `x`, estimated with the same
  method on the same nodes (their order may differ); other methods or
  nodes raise an error. Default `NULL`: when `x` is a `netobject_group`
  or an `mcml`, every pair of its groups is tested. When both are
  `netobject_group`s, the groups present in both under the same name are
  tested, and having none in common raises an error.

- iter:

  The number of shuffles, at least 2. Default `1000`. The smallest
  attainable p-value is `1 / (iter + 1)`.

- alpha:

  The significance level, strictly between 0 and 1, that marks an edge
  or centrality as significant (`sig`, `diff_sig`) and that the
  `nestimate_few_actors` warning is checked against. It does not change
  the p-values. Default `0.05`.

- paired:

  Whether the two networks hold paired observations (the i-th sequence
  or row of `x` with the i-th of `y`, such as the same people measured
  twice). Each shuffle then swaps the two members of a pair at random
  instead of shuffling all labels, so `x` and `y` need the same number
  of sequences (after `concat`) or rows; otherwise an error is raised.
  It cannot be combined with `actor`, of which it is the special case of
  one actor per pair. Default `FALSE`.

- adjust:

  The multiplicity correction for the p-values, a method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"none"`
  (default), `"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"`, `"BH"`,
  `"BY"` or `"fdr"`. The correction runs over the tested edges: the
  cells non-zero in either network (a self-loop included), each pair
  counted once in an undirected network, whose matrix stays symmetric.
  Cells zero in both networks keep their raw p-value, 1. The centrality
  p-values are corrected separately, over all states and measures. Any
  other value raises an error of class `nestimate_bad_input` before
  anything runs.

- measures:

  Centrality measures to test as well as the edges: any of
  `"OutStrength"`, `"InStrength"`, `"ClosenessIn"`, `"ClosenessOut"`,
  `"Closeness"`, `"Betweenness"`, `"BetweennessRSP"`, `"Diffusion"` and
  `"Clustering"` (`"InCloseness"` and `"OutCloseness"` are accepted for
  `"ClosenessIn"` and `"ClosenessOut"`), or `"all"` for the nine.
  Default `NULL` (edges only). Every state and measure gets the observed
  difference, an effect size and a p-value from the same shuffled
  networks as the edge test. The centralities are computed as by
  [`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
  with its defaults: self-loops left out, heavier edges as shorter
  steps, and `Diffusion` range-normalised to 0–1 in each network. An
  unknown name raises an error, and so does any `measures` with
  [`net_edge_betweenness()`](https://pak.dynasite.org/Nestimate/reference/net_edge_betweenness.md)
  input.

- nlambda:

  Ignored, kept so existing calls run. Default `NULL`. The shuffled
  glasso networks use the `nlambda` the networks were built with (see
  *Details*); a value other than `NULL` gives a warning of class
  `nestimate_ignored_argument`.

- seed:

  An integer that makes the shuffles reproducible. Default `NULL`. The
  session's own random stream is left as it was; with grouped input each
  pair, or each shared group, is tested from the same seed.

- actor:

  The column identifying the actor each sequence belongs to: the person
  whose sessions they are, or the team of a student. It is looked up in
  the network's metadata (such as `"student_id"` when sessions are
  nested in students, or `"Group"` for students nested in teams) and
  then in its wide sequence data. Whole actors are then shuffled, and
  the ICC and design effect are reported; see the sections *Nested data
  and actor* and *ICC and design effect*. Default `NULL`: sequences are
  shuffled one by one. A value that is not one column name, or `actor`
  with `paired = TRUE`, raises an error of class `nestimate_bad_actor`;
  a column that is absent or has missing values,
  `nestimate_actor_missing`; a column whose values do not give one id
  per sequence (or, with `concat`, that joins rows of different actors),
  `nestimate_actor_misaligned`; a method other than the transition
  methods, `nestimate_actor_unsupported`.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods: not used.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) methods: a
  `net_permutation`, `net_permutation_group` or `wtna_perm_mixed`.

## Value

A `net_permutation` object.
[`print()`](https://rdrr.io/r/base/print.html) shows the settings, the
number of edges tested and significant, the global test and, with
`actor`, the ICC and design effects;
[`summary()`](https://rdrr.io/r/base/summary.html) returns the edge
table, one row per edge present in either network;
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
draws the difference network. Grouped input gives a
`net_permutation_group`, a named list of `net_permutation` results (one
per group pair, named `"<group i> vs <group j>"`, when `y` is `NULL`;
one per shared group name otherwise), whose
[`print()`](https://rdrr.io/r/base/print.html) prints each and whose
[`summary()`](https://rdrr.io/r/base/summary.html) stacks the edge
tables with a leading `group` column naming the pair or group. Two
`wtna_mixed` inputs give a `wtna_perm_mixed`, the `net_permutation`
results of the `transition` and `cooccurrence` networks, whose
[`print()`](https://rdrr.io/r/base/print.html) prints both and whose
[`summary()`](https://rdrr.io/r/base/summary.html) returns a list of
their two edge tables. The
[`print()`](https://rdrr.io/r/base/print.html) methods return their
input invisibly. The `net_permutation`'s fields:

- x, y:

  The two networks compared.

- diff:

  The observed difference matrix, `x - y`.

- diff_sig:

  The observed difference where `p < alpha`, else 0.

- p_values:

  The p-value matrix (adjusted when `adjust` is not `"none"`).

- effect_size:

  The effect-size matrix: the observed difference divided by the SD of
  the shuffled differences.

- summary:

  The edge table [`summary()`](https://rdrr.io/r/base/summary.html)
  returns: one row per edge present in either network (an undirected
  network has one row per unordered pair), with columns `from`, `to`,
  `weight_x`, `weight_y`, `diff`, `effect_size`, `p_value` and `sig`.

- global:

  The global test, one row per statistic: `statistic` (`"M"`, the sum of
  absolute differences, and `"S"`, the largest absolute difference),
  `observed`, and `p_value` from the same shuffles as the edge test.
  Absent for edge-betweenness input.

- method:

  The network estimation method (`"edge_betweenness"` for
  edge-betweenness input).

- source_method:

  For edge-betweenness input: the source networks' method.

- edge_betweenness:

  For edge-betweenness input: the `invert` setting and the source method
  used.

- iter, alpha, paired, adjust:

  The settings used.

- actor:

  The `actor` column name, or `NULL`.

- n_actors:

  The number of distinct actors, or `NULL`.

- null_sd:

  The SD of each edge's shuffled differences (the effect-size
  denominator), as a matrix.

- null_sd_m:

  The SD of the global statistic M over the shuffles. Absent for
  edge-betweenness input.

- clustering:

  Only with `actor`: a one-row data frame with `n_sequences`,
  `n_actors`, `design` (`"between"`, `"within"` or `"mixed"`), `icc`
  with `icc_ci_lower` and `icc_ci_upper` (how alike the sequences of one
  actor are), `deff_edges` (the median design effect over edges),
  `deff_global` (the design effect for M) and `min_p`, the smallest
  attainable p-value.

- clustering_edges:

  Only with `actor`: one row per edge of the edge table, with `from`,
  `to`, `icc`, `null_sd_actor`, `null_sd_sequence` and `deff`.

- centralities:

  Only with `measures`: a list of `stats` (one row per state and
  measure, with `state`, `centrality`, `diff_true`, `effect_size` and
  `p_value`), `diffs_true` (the observed differences, one row per state
  and one column per measure) and `diffs_sig` (the same where
  `p < alpha`, else 0).

## Details

For transition methods, per-sequence count matrices are computed once
and each shuffle only re-sums them, which keeps a large `iter` fast; the
networks' `scaling` and `threshold` are applied to every shuffled
network as to the observed one. For the other methods the estimator is
re-run on every shuffle (`cor` and `pcor` through a lighter direct
path). For `glasso`, every shuffled network is estimated as
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
estimated the observed ones, with `x`'s `params` (each sample's own
lambda path, `nlambda`, `lambda.min.ratio`, `gamma`,
`penalize.diagonal`, `refit`), so the observed difference and its null
come from one estimator (van Borkulo et al. 2022). Two glasso networks
built with different settings raise an error of class
`nestimate_bad_input`. A network that rests on one sequence, or one-hot
data with no actor or session, has nothing to shuffle and raises an
error of class `nestimate_single_sequence`.

Wide rows combined with `params = list(concat = k)` are shuffled as the
merged sequences the estimator used, keeping their internal joins. A
merged sequence must lie within one actor when `actor` is set (class
`nestimate_actor_misaligned` otherwise), and paired tests need the same
number of merged sequences. Both networks must use the same `weighted`
setting (class `nestimate_bad_input` otherwise).

A refit that fails or produces non-finite weights makes the permutation
statistic undefined for that shuffle. The test then raises an error of
class `nestimate_permutation_failed` that reports the number of failed
shuffles, and returns no p-values: failed shuffles are neither counted
as non-exceedances nor dropped from the null.

For
[`net_edge_betweenness()`](https://pak.dynasite.org/Nestimate/reference/net_edge_betweenness.md)
input, the source networks are shuffled and edge betweenness is
recomputed on each shuffle. Both objects must come from the same source
method and use the same `invert` setting.

The centrality block follows the layout of
`tna::permutation_test(measures = )`: a centrality whose shuffled
differences do not vary gets an effect size of `NaN` or `Inf`, as there,
rather than the 0 the edge table uses.

## What is tested

Two kinds of question are answered from the same shuffles.

- Edge tests:

  One test per edge: is the difference in this edge's weight, `x - y`,
  larger in absolute value than the shuffles produce?
  [`summary()`](https://rdrr.io/r/base/summary.html) reports each with
  an effect size (the observed difference divided by the SD of the
  shuffled differences, 0 when that SD is 0) and a p-value,
  `(number of shuffles at least as extreme + 1) / (iter + 1)`.

- Global test:

  One test for the whole network: do the two networks differ at all? Two
  statistics: **M**, the sum of the absolute differences over all cells
  of the weight matrix (the total amount of difference; an undirected
  edge counts twice), and **S**, the largest absolute difference. Being
  a single test, it needs no multiplicity correction.
  [`print()`](https://rdrr.io/r/base/print.html) shows it.

The letters follow this function's historical convention. In
[`nct()`](https://pak.dynasite.org/Nestimate/reference/nct.md), M is
instead the largest absolute edge difference and S the difference in
global strength; the sum of absolute edge differences here measures
total edge difference, not a difference in global strength. The smallest
attainable p-value is `1 / (iter + 1)`; with the default `iter = 1000`
it is 0.000999, meaning that no shuffle came close.

## Nested data and `actor`

Shuffling treats every sequence as an exchangeable unit. When several
sequences come from the same actor (sessions of one person, students of
one team), `actor` names the column identifying that actor, and the
shuffle respects it (Good, 2005; Anderson & ter Braak, 2003): an actor
whose sequences are all in one group moves to the other group as a
whole, and an actor with sequences in both groups has their labels
shuffled among their own sequences only. Mixed designs combine the two.
The observed differences do not change; only the p-values and effect
sizes do. With few actors there are few distinct ways to shuffle them,
and a warning of class `nestimate_few_actors` is raised when no p-value
could fall below `alpha`.

`actor` is available for transition networks built from sequences
(`"relative"`, `"frequency"`, `"co_occurrence"`). On one-hot data the
actor column is held per row rather than per unit, so the call raises an
error of class `nestimate_actor_misaligned`. The other methods,
association networks (`"cor"`, `"pcor"`, `"glasso"`, ...) among them, do
not keep row identifiers after estimation and raise an error of class
`nestimate_actor_unsupported`.

## ICC and design effect

With `actor`, [`print()`](https://rdrr.io/r/base/print.html) also
reports:

- ICC:

  The intraclass correlation, the proportion of the total variance that
  lies between actors (Shrout & Fleiss, 1979). An ICC close to 0
  indicates little evidence of a nesting effect. It is the one-way ANOVA
  ICC of each sequence's transition shares within each group, averaged
  over edges weighted by edge frequency and jackknife bias-corrected,
  with a 95% interval from the leave-one-actor-out jackknife (Efron &
  Tibshirani, 1993).

- Design effect:

  The ratio of the variance of an estimate under the clustered design to
  its variance had the units been sampled independently (Kish, 1965).
  Here it is the variance of the shuffled differences when whole actors
  are moved, divided by the variance when single sequences are moved,
  both drawn in the same run, reported as the median over edges and for
  the global statistic M. For equal numbers of sequences per actor m,
  Kish gives the approximation `1 + (m - 1) * ICC`.

## Reading the printed output


    Permutation Test: Transition Network (relative probabilities) [directed]
      Iterations: 1000  |  Alpha: 0.05  |  Actor: Group (200 actors)
      Nodes: 9  |  Edges tested: 78  |  Significant: 42
      Global test (networks differ overall?): M = 2.612 (p = 0.000999)  |  ...
      Nesting in Group: ICC = -0.002 [95% CI -0.006, 0.002]  |  between design
      Design effect (1 = nesting does not matter): edges 1.03  |  global 1.20

Line 2: the settings, and the actor column with its number of actors.
Line 3: the edges present in either network and how many differ at
`alpha`. Line 4: the global test. Lines 5 and 6, only with `actor`: the
ICC with its interval, whether actors sit in one group (`between`), in
both (`within`) or either (`mixed`), and the design effects. With
`actor`, when no p-value can fall below `alpha` (too few actors, or too
small an `iter`), a last line says so.

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

[`permutation_diagnostics()`](https://pak.dynasite.org/Nestimate/reference/permutation_diagnostics.md)
to compare the actor-level and ordinary tests side by side;
[`bayes_compare()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
for the Bayesian complement, which instead of "is this difference more
extreme than chance?" answers "how probable is a difference, and how
large?";
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)

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
