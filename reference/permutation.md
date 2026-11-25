# Test whether two networks differ by permutation

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

  The first network:

  - a `netobject` from
    [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
    or a bare `cograph_network`;

  - a `netobject_group`, an `mcml` or a `wtna_mixed`;

  - a
    [`net_edge_betweenness()`](https://pak.dynasite.org/Nestimate/reference/net_edge_betweenness.md)
    result;

  - a `tna` model or `group_tna`, rebuilt as a network from its own
    sequences (see
    [`as_netobject()`](https://pak.dynasite.org/Nestimate/reference/as_netobject.md)).
    Raises `nestimate_actor_missing` when `actor` is given and the model
    keeps no metadata; models built with
    [`tna::prepare_data()`](https://sonsoles.me/tna/reference/prepare_data.html)
    keep it.

  For the [`print()`](https://rdrr.io/r/base/print.html) methods: a
  `permutation()` result.

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
  twice; for one-hot data, the i-th actor in the order the actors first
  appear). Each shuffle then swaps the two members of a pair at random
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
  and actor* and *ICC and design effect*. Several column names,
  outermost first, are nested levels (such as `c("Group", "Student")`):
  labels are shuffled at the level where they are assigned, within the
  units above it. Default `NULL`: sequences are shuffled one by one.
  Raises `nestimate_actor_missing` when the column is absent or
  incomplete, and `nestimate_actor_misaligned` when it does not give one
  actor per sequence.

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

  The `actor` column name (or names, outermost first), or `NULL`.

- n_actors:

  The number of distinct actors (with several levels, at the innermost
  level), or `NULL`.

- actor_levels:

  Only with several `actor` levels: one row per level, with `level` (the
  column), `n_actors` (the units at that level), `n_pure` and
  `n_crossed` (the units the shuffle met at that level, inside crossed
  units of the level above, that sit in one group or in both) and
  `design` (`"between"`, `"within"`, `"mixed"`, or `"not opened"` when
  no unit above was crossed).

- null_sd:

  The SD of each edge's shuffled differences (the effect-size
  denominator), as a matrix.

- null_sd_m:

  The SD of the global statistic M over the shuffles. Absent for
  edge-betweenness input.

- clustering:

  Only with `actor`: a one-row data frame with `n_sequences`,
  `n_actors`, `design` (`"between"`, `"within"` or `"mixed"`; with
  several levels, one `"<level>: <design> (<n> actors)"` entry per
  level), `icc` (at the innermost level) with `icc_ci_lower` and
  `icc_ci_upper` (how alike the sequences of one actor are),
  `deff_edges` (the median design effect over edges), `deff_global` (the
  design effect for M) and `min_p`, the smallest p-value the exact test
  over the actors' arrangements can give: one over their number, or
  `1 / (iter + 1)` when larger. The reported p-values are Monte Carlo
  estimates, so with few arrangements one can fall below it by chance.

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
path). A network built with `method = "latency"` (alone or in a group)
raises an error of class `nestimate_latency_resampling`;
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
is the verb that resamples latency networks. For `glasso`, every
shuffled network is estimated as
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
estimated the observed ones, with `x`'s `params` (each sample's own
lambda path, `nlambda`, `lambda.min.ratio`, `gamma`,
`penalize.diagonal`, `refit`), so the observed difference and its null
come from one estimator (van Borkulo et al. 2022). Two glasso networks
built with different settings raise an error of class
`nestimate_bad_input`. A network that rests on one sequence, or one-hot
data with no actor or session, has nothing to shuffle and raises an
error of class `nestimate_single_sequence`. A sequence network with no
observed transition raises an error of class `nestimate_no_transition`.

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

Several columns, outermost first (`actor = c("Group", "Student")`,
students within teams), name nested levels; the unit at each level is
the combination of its column with the ones before it. A level that adds
no units below the one before it (the levels given innermost first)
raises an error of class `nestimate_bad_actor`. A column whose ids
repeat across the units of the next column (a class code reused in every
team) cannot be told apart from a correctly nested one, so the order of
the columns is the caller's statement of the design. The labels are then
shuffled at the level where they are assigned, within the units above it
(Anderson & ter Braak, 2003; Good, 2005): teams whose sequences are all
in one group move as wholes among such teams, exactly as with one
column; a team with sequences in both groups is opened, and inside it
the students all in one group exchange labels with each other (never
with students of another team), while a student in both groups is opened
in turn, down to the last column, inside which single sequences exchange
labels. Set the levels whenever the groups are assigned below the top
level (students of one team in different conditions): `actor = "Group"`
alone would then shuffle single sequences within each team and treat the
sessions of one student as independent. The number of distinct
arrangements, the `nestimate_few_actors` check and `min_p` follow the
whole scheme. One column, or `NULL`, runs exactly the scheme described
above.

When to set `actor`: when the groups are made of different actors (high-
versus low-achieving students, each student in one group), set it.
Shuffling single sequences then treats the sequences of one actor as
independent, and the test is too liberal. In one data set of sessions
nested in students, two groups of different students with no true
difference had 11.7% of edges significant at `alpha = 0.05` when
sessions were shuffled (23% of replicates had at least one significant
edge after Benjamini-Hochberg correction), against 4.6% (0.5%) when
whole students were shuffled. When the groups vary within the same
actors (pre- and post-test sessions of the same students), `actor`
shuffles the labels within each actor; the plain test gave similar size
in that simulation (5.1% against 5.0%). This does not establish validity
for other within-actor designs: the permitted shuffles must respect the
assignment mechanism or exchangeability under the null. Whether nesting
matters is decided by the design effect, not by the ICC alone: in that
data set the transitions of one student correlated at about 0.01, yet
with about 59 transitions per student (size-weighted mean) the Rao-Scott
design effect of the transition table was about 1.5
([`mosaic_test()`](https://pak.dynasite.org/Nestimate/reference/mosaic_test.md)).
The printed ICC is a different quantity, the correlation between the
sequences of one actor (0.05 in that data set, with a printed design
effect of 1.24 for the bootstrap), so Kish's approximation with it
counts sequences per actor, not transitions.

`actor` works for the transition networks built from sequences:
`"relative"`, `"frequency"` and `"co_occurrence"`. The other methods
keep no row identifiers after estimation, and with them `actor` raises
`nestimate_actor_unsupported`; one-hot data, which hold the actor per
row, raise `nestimate_actor_misaligned`.

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
  Tibshirani, 1993). With several `actor` levels it is computed at the
  innermost level (the finest actor, such as the student within a team),
  and [`print()`](https://rdrr.io/r/base/print.html) says so.

- Design effect:

  The ratio of the variance of an estimate under the clustered design to
  its variance had the units been sampled independently (Kish, 1965).
  Here it is the variance of the shuffled differences when whole actors
  are moved, divided by the variance when single sequences are moved
  (with several levels: under the whole nested scheme), both drawn in
  the same run, reported as the median over edges and for the global
  statistic M. For equal numbers of sequences per actor m, Kish gives
  the approximation `1 + (m - 1) * ICC`.

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
`actor`, when `min_p` is above `alpha` (too few actors for the exact
test, or too small an `iter`), a last line says so. With several `actor`
levels, line 2 names the levels (`Actor: Group > Student (nested)`), and
line 5 is replaced by two: the design at every level, with its number of
actors
(`Design: Group: mixed (200 actors), Student: between (2000 actors)`),
and the ICC at the innermost level.

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
perm <- permutation(net1, net2, iter = 20, seed = 42)
print(perm)
#> Permutation Test: Transition Network (relative probabilities) [directed]
#>   Iterations: 20  |  Alpha: 0.05
#>   Nodes: 4  |  Edges tested: 16  |  Significant: 2
#>   Global test (networks differ overall?): M = 3.858 (p = 0.0952)  |  S = 0.567 (p = 0.143)
summary(perm)
#>    from to   weight_x  weight_y        diff effect_size    p_value   sig
#> 1     A  A 0.33333333 0.2222222  0.11111111   0.7072043 0.52380952 FALSE
#> 2     A  B 0.06666667 0.3333333 -0.26666667  -2.0679738 0.04761905  TRUE
#> 3     A  C 0.20000000 0.2222222 -0.02222222  -0.1421124 0.90476190 FALSE
#> 4     A  D 0.40000000 0.2222222  0.17777778   1.1120847 0.23809524 FALSE
#> 5     B  A 0.30769231 0.3333333 -0.02564103  -0.1091872 1.00000000 FALSE
#> 6     B  B 0.38461538 0.0000000  0.38461538   2.0383530 0.09523810 FALSE
#> 7     B  C 0.23076923 0.2500000 -0.01923077  -0.1012939 1.00000000 FALSE
#> 8     B  D 0.07692308 0.4166667 -0.33974359  -1.7822963 0.14285714 FALSE
#> 9     C  A 0.33333333 0.1111111  0.22222222   0.8892224 0.42857143 FALSE
#> 10    C  B 0.66666667 0.3333333  0.33333333   1.2204017 0.28571429 FALSE
#> 11    C  C 0.00000000 0.2222222 -0.22222222  -1.3692569 0.23809524 FALSE
#> 12    C  D 0.00000000 0.3333333 -0.33333333  -1.6800434 0.19047619 FALSE
#> 13    D  A 0.00000000 0.5000000 -0.50000000  -1.6730600 0.14285714 FALSE
#> 14    D  B 0.00000000 0.2000000 -0.20000000  -0.8021641 0.52380952 FALSE
#> 15    D  C 0.33333333 0.2000000  0.13333333   0.4775517 0.66666667 FALSE
#> 16    D  D 0.66666667 0.1000000  0.56666667   2.8652993 0.04761905  TRUE

# Students are nested in teams, and Achiever is a team-level label:
# permute whole teams, not single students
net <- build_network(group_regulation_long, method = "relative",
                     actor = "Actor", action = "Action", time = "Time",
                     group = "Achiever")
permutation(net, iter = 20, actor = "Group", seed = 1)
#> Grouped Permutation Test
#> Groups: High vs Low 
#> 
#> -- High vs Low --
#> Permutation Test: Transition Network (relative probabilities) [directed]
#>   Iterations: 20  |  Alpha: 0.05  |  Actor: Group (200 actors)
#>   Nodes: 9  |  Edges tested: 78  |  Significant: 41
#>   Global test (networks differ overall?): M = 2.612 (p = 0.0476)  |  S = 0.210 (p = 0.0476)
#>   Nesting in Group: ICC = -0.002 [95% CI -0.006, 0.002]  |  between design
#>   Design effect (1 = nesting does not matter): edges 1.08  |  global 0.59
# }
```
