# Nestimate 1.0.1

## Breaking changes since CRAN 0.9.24

* `plot_mosaic()` is no longer exported (it masked `tna::plot_mosaic()`);
  use `mosaic_plot()`, whose data-frame method takes `x_var`, `y_var` and
  `weight` (see the 0.9.29 notes).
* `sequence_plot()` gained `periods`, `fill`, `actor`, `action`, `time` and
  `order` after `group`, so arguments after `group` must be passed by name.
* `hypergraph_centrality()` returns only the clique centrality by default;
  `$Z` and `$H` are computed only when requested and only for uniform
  hypergraphs.

## New features

* `stability_landscape()` draws a transition network as a ball-in-cup
  landscape in which each state's cup depth is the log of its expected stay,
  `-log(1 - P_ii)` (stay = e^depth), with intervals from resampling whole
  sequences, the share of time in each state, its exit probabilities and,
  for three states, the net circulation. It takes a network, a grouped
  network (for example `build_network()` of `build_clusters()`) or a
  sequence data frame or matrix; `states =` sets the order of the cups.
  `print()`, `summary()` (`table = "cups"`, `"exits"`, `"circulation"`) and
  `plot()`.
* `landscape_surface()` draws the cup depths as a 3D surface over a
  per-sequence stability measure from `sequence_metrics()` (self-loop
  tendency, entropy, complexity, ...), estimated from the nearest sequences
  at each value of the measure, with the most stable end in front.
  `summary()` compares the two ends of each measure; `plot()` draws one
  `persp()` surface per network and measure.
* `landscape_trend()` traces the cup depths along the program phase
  (sliding windows of time steps) or along a numeric per-sequence
  covariate such as a grade (nearest-sequence slices), one line per state
  with intervals from resampling sequences.

## Bug fixes

* `build_network()`: wide data with number-coded states and a text column
  (an id) built the network from the text column alone, an all-zero network
  with no warning. A numbered run of numeric position columns is now the
  sequence, and the other columns are metadata.
* The sequence verbs (`extract_sequences()`, `sequence_patterns()`,
  `sequence_metrics()`, `pattern_outcome()`) read the data of a correlation,
  glasso, Ising, MGM or window network as states; they now raise
  `nestimate_bad_input`. A network grouped by two columns (`group =
  c("Achiever", "Course")`) no longer fails; its groups are labelled in one
  column, `Achiever-Course`.
* `sequence_metrics()`: a row with no observed state has `NA`, not 0, for
  `valid_proportion` and `cyclic_feedback_strength`.
* `sequence_patterns(pattern =)` accepts the spaced notation
  `"A -> * -> B"`.
* `nct()`: an unnamed matrix compared with a named one is matched by
  position again, as documented.
* `build_hon()`, `build_mogen()` and `path_counts()`: a missing value inside
  a sequence ends a run of observed states. `build_hon()` stopped with an
  error, `build_mogen()` counted `NA` as a state, and `path_counts()` joined
  the two sides of the gap.
* `association_rules()` reads factor columns of a wide data frame by their
  labels, not their integer codes.
* `build_simplicial(type = "pathway")` leaves out HYPA paths that never
  occur: since 0.9.26 HYPA scores every possible path, and an unobserved
  one flagged "under" became a simplex of states never seen together.
* `persistent_homology()`: the rows of `$persistence` are numbered without
  gaps.
* `labels =` on number-coded states: every state column is relabelled
  together, so `bootstrap_network()` no longer drops the cells of columns
  that held no relabelled state.
* `landscape_surface()` and `plot(landscape_trend(), type = "surface")` draw
  a flat surface instead of failing; a constant measure or covariate raises
  `nestimate_bad_input`; slices resting on fewer than five effective
  sequences warn (`nestimate_thin_slice`); input with no observed
  transition raises `nestimate_landscape_unestimable`; a two-state
  occupancy landscape reports its single ridge once.
* `plot()` of a `passage_time()` result drew the matrix transposed: the
  tile in row A, column B showed the passage time from B to A. Rows are now
  the starting states and columns the targets, as documented, for a single
  network and for each panel of a grouped one.
* `pattern_outcome()`: an `adjust` column named like one of the model's own
  terms (`.p1`, `.y`, ...) was replaced by that term, so the requested
  adjustment was silently lost; an `adjust` name with a space failed to
  parse. Every model column now has an internal name. `outcome` defaults to
  `"last_obs"`, as in codyna and as the 0.9.29 notes stated.
* `sequence_patterns()`: a gapped pattern (`A->*->B`) no longer spans a
  missing position, matching the template search; a repeated `len` or
  `gap` (`len = c(2, 2)`) no longer counts each occurrence twice; a grouped
  search that finds no pattern returns zero rows with all its columns and
  no warning; a template ending in `->` is refused.
* `build_network()` and `extract_sequences()`: an all-missing position
  column stored as logical `NA` no longer stops the numbered position names
  from forming a run, which had turned a character covariate into states.
* `plot()` of a grouped `markov_stability()` result handles more than nine
  groups: colours cycle and shapes shift, so every group keeps its own
  colour-shape pair.
* `extract_sequences()` and the sequence verbs raise `nestimate_bad_input`
  for a model built from a matrix, which holds no sequences.

## Documentation

* `sequence_metrics()`: the emergent-state rule for a return to the first
  state is a spell of at least 3 that is more than half the opening spell
  (codyna's rule); the help page said "over twice".

# Nestimate 0.9.29

## Sequence patterns

* `extract_sequences()` returns the sequences behind any input as a wide
  data frame, one row per sequence: prepared data, a long event log
  (`actor`/`action`/`time`), a wide data frame or character matrix, a
  `netobject`, a grouped `netobject_group` (the group column is added) or a
  `tna` model. `meta = TRUE` adds each sequence's covariates.
* `sequence_patterns()` discovers ordered patterns: contiguous n-grams,
  gapped patterns (`A->*->B`) and repeated runs. For each it reports the
  occurrences, the share of patterns of that length, the number and share
  of sequences that contain it (support), and its lift. The counts, support
  and lift agree with `codyna::discover_patterns()`. With `group` (a column,
  a vector, or the grouping of a `netobject_group`) each pattern also gets
  its support per group and a test of presence by group: for two groups the
  odds ratio, 95% CI and p-value of Fisher's exact test, for more the
  chi-square and Cramer's V; p-values are BH-adjusted. `top` returns the
  first rows.
* `sequence_patterns(pattern =)` searches for a template such as
  `"discuss->*->consensus"` (each `*` is one position of any state) and
  returns every concrete instance; `group = "last_obs"` groups sequences by
  their final state. Both agree with `codyna::discover_patterns()`.
* `sequence_metrics()` computes 23 per-sequence indices (validity, spells,
  entropy, Simpson diversity, transitions, initial-state persistence and
  decay, dominant and emergent states, cyclic feedback, integrative
  potential, complexity), one row per sequence. They agree with
  `codyna::sequence_indices()` on 20 indices; three are corrected:
  `dominant_max_spell` is per sequence (codyna repeats the last sequence's
  value), `longitudinal_entropy` is normalised by the number of states
  (codyna counts `NA` as a state), and `complexity_index` uses that
  entropy once (codyna normalises it twice) and is 0 for a sequence of one
  spell, as in `TraMineR::seqici()` (codyna: `NA`). Empty sequences give
  `NA`.
  A missing position is not an observation: the opening spell starts at the
  first observed state, transitions are pairs of adjacent observed states
  (`self_loop_tendency + transition_rate = 1`), and `integrative_potential`
  weights observed positions only, as `TraMineR::seqintegr()` does. The
  indices shared with TraMineR agree with it on sequences without internal
  gaps.
* `extract_sequences()` and the sequence verbs read TraMineR's void (`%`)
  and missing (`*`) markers as `NA`, as `build_network()` does; an
  `stslist` such as tna's `engagement` no longer gains a `%` state.
* `pattern_outcome()` reproduces `codyna::analyze_outcome()`, including its
  default `outcome = "last_obs"`: the `n`
  patterns ranked highest by `priority` (default: codyna's chi-square)
  enter one logistic regression of a two-valued outcome on their presence,
  with an optional random intercept (`group`). Coefficients agree with
  codyna's fit. It returns a tidy table, one row per pattern, with support
  per outcome value, the selection statistic, estimate, Wald CI, p-value and
  odds ratio. It fits the selected pattern even where codyna 1.0.0 fits a
  data column instead (a lone pattern of one length loses its name, and a
  data column named `V<k>` takes its place).

## Validation tests aligned with the current tna release

* `plot_mosaic()` is no longer exported: it masked `tna::plot_mosaic()`
  whenever both packages were attached, so tna code (and apps such as
  Dynalytics) calling `plot_mosaic(model)` failed with `argument "data" is
  missing`. Its data-frame mosaic is now a method of `mosaic_plot()`:
  `plot_mosaic(df, x = "group", y = "state", weight = "count")` becomes
  `mosaic_plot(df, x_var = "group", y_var = "state", weight = "count")`
  (`x_label`/`y_label` become `xlab`/`ylab`, as in the other methods). A
  network, table or matrix given a data-frame-only argument (`x_var`,
  `y_var`, `weight`, `fill`, `colors`, `show_labels`, and `label_size`
  except with `style = "flat"`, which uses it for its tile labels) raises
  `nestimate_bad_input` instead of ignoring it, and so does an unknown
  argument of the data-frame method.
* `casedrop_reliability()` compares only the edges present in the original
  network. A subsample cannot contain a transition that is absent from the
  full data, so the empty cells agreed in every draw and inflated the
  correlation: on a sparse 8-state chain the CS-coefficient was 0.9 (the
  top of the grid) with a mean correlation of 0.993 after dropping 90% of
  the sequences. `n_edges` now counts the present edges, and an undefined
  CS-coefficient (edge weights that do not vary) is `NA` instead of 0. Such
  a result now has the usual summary columns with `NA` statistics, so
  `summary()` (including the grouped summary) and `plot()` work: the plot
  states that the CS-coefficient is undefined, and the grouped plot names the
  networks it leaves out. This case is common now, for example a fixed
  cycle, where every present edge has probability 1. The documentation no
  longer claims that case-dropping is valid regardless of within-actor
  correlation.
* `permutation(measures = "Diffusion")` (and the node table of
  `compare_networks()`) compares Diffusion after min-max scaling to 0-1, as
  `net_centrality()` reports it and as tna does from 1.4.1.
* The range normalization behind `normalize = TRUE` and Diffusion returns
  `NaN` for a numerically constant vector instead of stretching rounding
  error to 0-1, as tna's does.
* `markov_order_test()` splits sequences at missing values, so no window
  spans an `NA`; an internal `NA` previously became a context named "NA".
  Results on data without internal missing values are unchanged (identical
  under the same seed). Its documentation no longer calls the
  within-context permutation exact: the windows of a sequence overlap, so
  the permutation distribution approximates the null (Besag & Mondal, 2013).

## Bug fixes

* `plot()` works on `passage_time()` and `markov_stability()` results for a
  grouped network: one heatmap panel per group on a shared colour scale,
  and a per-metric dot plot with one point per group. Both used to fall
  through to `plot.default()` and fail. The single-network stability plot
  uses one fill; its per-state palette left a ninth state without colour.
* `build_mmm(covariate_effect = "posthoc")`: `print()` and `summary()`
  labelled the covariates "integrated into EM -- influences cluster
  membership" whatever the setting. The fit records `covariate_effect`,
  and a post-hoc analysis is now labelled as one.
* The `print()` methods of `entropy_trajectory()` and `entropy_bayes()` no
  longer tell the user to reach into the object with `$`; they point to
  `summary()` and `plot()`.
* `build_network()` no longer reads a character covariate of wide sequence
  data as states. An outcome or group label such as `High`/`Low` used to
  become nodes of the network and vanish from `$metadata`. When the
  position columns are named as a run (`T1`, `T2`, ...; `V1`, ...; `1`,
  `2`, ...), every other character or factor column is now metadata, like
  numeric columns already were, even when its values coincide with state
  names. The network equals the one built with `metadata_cols`. Without
  such a run, or with two runs of equal size, nothing changes.
* Association `permutation()` refuses inference with
  `nestimate_permutation_failed` when any permuted network cannot be
  estimated or has non-finite weights, reporting the number of failed draws.
  Failed draws no longer silently reduce p-values or invalidate null moments.
* Transition resampling honors `weighted` and `concat`. Concatenated rows
  are resampled as the merged sequences used by the estimator, including
  their joining transitions. With `actor`, a merged sequence must stay
  within one actor; paired tests require matching merged-sequence counts.
* Relabelling numeric sequence states, start/end markers and onehot columns
  keeps stored data, alphabets and column selectors consistent for resampling.
  Wide sequence positions retain their column names; association measurements
  retain their numeric values.
* Higher-order builders and `sequence_compare()` preserve literal numeric
  states in Nestimate models and decode positional codes only for tna input
  or explicitly encoded converted data. Character digit alphabets such as
  `"1"` through `"12"` are decoded correctly.
* Grouped wide sequences use the same resolved state columns as their child
  estimates. Numeric covariates do not enter the shared alphabet, and factor
  columns keep their labels.
* MGM's numeric `level` is validated after excluding group, actor and
  metadata columns. MGM actor identifiers are excluded from estimation.
* Custom variable estimators retain their numeric cleaned data when column
  names identify the modeled nodes, allowing valid downstream resampling.
* `nct()` rejects `Inf` and `-Inf` instead of silently isolating the affected
  variable; its supported missing-data and constant-variable handling is
  unchanged.
* `plot_mosaic(network, x_label =, y_label =)` applies the labels. Explicit
  primitive-only styling arguments raise a classed error directing the
  caller to `mosaic_plot()`.
* `association_rules()` reads a character matrix of states (such as
  `trajectories`) as wide sequences, like a character data frame. It used to
  treat every matrix as binary: with missing cells it failed with
  "subscript out of bounds", and without them it silently mined the column
  names. Binary input (matrix or data frame, 0/1 or logical) with a missing
  cell or a value other than 0/1 now raises `nestimate_bad_input`; logical
  data-frame columns are read as binary items instead of the strings
  `"TRUE"`/`"FALSE"`.

## Documentation

* MGM's missing-data policy is mean/modal single imputation, as documented
  by psychnets for its `"pairwise"` option; it is not pairwise deletion.
* Permutation M/S labels describe this function's historical total/max edge
  difference statistics and distinguish them from NCT's structure/strength
  convention. Comparison tables no longer call total edge difference
  "global strength". The NCT tutorial and local display checks use current
  notation and capitalization.

# Nestimate 0.9.28

## Higher-order methods: corrections from a review against the source papers

These change results. Each is tested against an independent value (a
hand-derived oracle or a published example), and each new test fails on
the previous code.

* `hypergraph_centrality()` computes the clique-motif centrality (CEC)
  directly as the Perron vector of each connected component. The power
  iteration it replaces oscillated when the clique expansion was bipartite
  (any triangle-free hypergraph), and returned a different wrong vector
  depending on whether `max_iter` was odd or even, without a warning. A
  star of pairs now gives the exact vector for every `max_iter`.
* `hypergraph_centrality()` now defaults to `type = "clique"`. `"Z"` and
  `"H"` are defined by Benson (2019) for uniform hypergraphs only: with
  mixed hyperedge sizes the equation mixes polynomial degrees and its
  solution depends on the normalisation inside the iteration. Requesting
  them on a non-uniform hypergraph raises `nestimate_nonuniform_hypergraph`.
  On uniform hypergraphs they are unchanged and reproduce Benson's Table 2.
  An unconverged Z/H iteration warns with `nestimate_no_converge`.
* `q_analysis()` counts, at each level q, only the maximal simplices of
  dimension at least q (Atkin's Q-vector). It previously counted every
  maximal simplex at every level, so lower-dimensional simplices appeared as
  extra components at high q.
* `persistent_homology()` builds the filtration one dimension above
  `max_dim` (the Ripser convention), so classes of dimension `max_dim` can
  die. Previously every top-dimension cycle was reported as essential; with
  the default `max_dim = 3`, dense 7-state transition networks showed
  spurious never-dying bars in every case tried.
* `build_simplicial()` (clique and VR) records when `max_dim` cut a larger
  clique. `betti_numbers()` then warns (`nestimate_truncated_complex`) that
  the top Betti number is an upper bound, and `persistent_homology()` warns
  the same way when given such a complex, whose top-dimension classes
  cannot die.
* `build_mogen()` and `markov_order_test()` gain `dof = c("paths",
  "observed")`. The default `"paths"` counts each layer's parameters as
  Scholtes (2017, Eq. 7) and pathpy do: every transition of the paths of
  length k that the observed first-order topology allows. The previous rule,
  now `dof = "observed"`, counted only observed transitions; it can fall as
  the order rises and under-penalised high orders on sparse data, so AIC,
  BIC and the likelihood-ratio test favoured higher orders. Selected orders
  can change.
* `build_mogen(criterion = "lrt")` selects the largest order whose test
  rejects (Scholtes 2017) instead of stopping at the first non-rejection,
  and returns every p-value as `lrt_p`.
* `build_hon(method = "hon", min_freq >= 2)`: a source whose transitions
  are all below `min_freq` is no longer an extension candidate nor a
  starting rule, as in the reference pyHON. Such sources could previously
  add HON nodes that pyHON does not create. `min_freq = 1` and HON+ are
  unaffected.
* The Z and H iterations of `hypergraph_centrality()` normalise the tensor
  update before adding the stabilising shift, and count as converged only
  when the eigen-equation residual is small. On large uniform hyperedges the
  unnormalised shift swamped the update, so the iteration stopped at its
  start and returned the uniform vector without a warning (a 3-petal
  sunflower of 20-member hyperedges gave a hub/petal ratio of 1 instead of
  Benson's sqrt(3)).
* Path degrees of freedom are summed as non-negative per-state terms, so
  they stay exact up to the double range; beyond it `build_mogen()` and
  `markov_order_test()` stop with `nestimate_dof_overflow` instead of
  carrying `Inf`/`NaN` into AIC and BIC.
* Truncation detection reads the simplices already built (does a top
  simplex extend by one vertex?) instead of an exact maximum-clique search,
  which took seconds on dense graphs even at `max_dim = 0`; persistence on a
  distance matrix no longer computes it at all. The persistence warning for
  a truncated complex rechecks within `max_scale`.
* `build_hypa()` classifies a path as over-represented when
  `p_over <= alpha`, i.e. HYPA >= 1 - alpha as in LaRock et al. (2020); it
  was `< alpha`. This only changes ties at exactly `alpha`.

## Documentation

* `build_mogen()` is documented as the multi-order model of Scholtes (2017),
  pathpy's `MultiOrderModel`; it is not the multi-order generative model of
  Gote et al. (2023), which adds start and end states.
* `build_honem()` states what it computes (decayed matrix powers of the HON
  transition matrix over all HON nodes) and how that differs from Saebi et
  al. (2020), whose order-indexed matrices over first-order nodes are not
  fully specified in the paper.

# Nestimate 0.9.27

## New features

* `sequence_plot()` gains `periods =` and `fill =`. For a network built with
  `actor =` and `session = "day"` (or any period column),
  `sequence_plot(net, periods = "day")` draws one row per actor and splits
  the time axis into equal-width blocks, one per period, with ticks at the
  period boundaries and each label under its block. The number of steps may
  differ by actor and by period: `fill = "stretch"` (default) spreads an
  actor's steps over the whole block (one step fills it, three take a third
  each); `fill = "pad"` gives every block as many fixed slots as the busiest
  actor-period and leaves the rest empty, so a slot is the same step
  position for everyone. Factor levels set the period order and an unused
  level shows as an empty block. Works for `type = "heatmap"` and
  `"index"`.
* `sequence_plot()` accepts long data directly: `actor =`, `action =`,
  `time =` and `order =` name the columns, as in `build_network()`, e.g.
  `sequence_plot(log, actor = "student", action = "action", periods = "day")`.
  With `periods`, each actor-period is one sequence whatever the time gaps.

## Breaking changes

* `nct()`: the M and S letters were the wrong way round. In van Borkulo et
  al. (2022) and `NetworkComparisonTest`, M is the network-structure test
  (largest absolute edge difference) and S the global-strength test; `nct()`
  stored and printed them swapped. The tests are now `$structure` (M) and
  `$strength` (S, with `by_group`, each network's global strength); reading
  `$M` / `$S` raises a `nestimate_renamed` error naming the field to use.
* `summary.net_nct()` now returns the M and S tests (columns `test`,
  `statistic`, `observed`, `p_value`, `significant`) instead of the per-edge
  table with `m_stat` / `s_stat` attributes (those attributes had the M/S swap
  too). The per-edge table is `as.data.frame()`.

## Bug fixes

* `plot_mosaic()` given a network (netobject, group, mcml, htna) hands it
  to `mosaic_plot()`; it failed with "is.data.frame(data) is not TRUE".
  Data-frame calls are unchanged.

* `cluster_network(cluster_by = "mmm")` passes `covariate_effect` (and any
  other `build_mmm()` argument) to `build_mmm()`; it failed with "Unknown
  argument(s) in cluster_network(): covariate_effect".

* `build_mcml()` decodes a `tna` object's integer-coded sequences (it
  failed with "Unknown nodes in clusters").
* `summary.mcml()` works without within-cluster layers
  (`compute_within = FALSE`, where it crashed), and `within_total`
  includes self-transitions. For raw counts it equals the macro diagonal;
  separately normalized micro and macro layers can have different totals.
* `entropy_trajectory()` windows run through time. Transitions were pooled
  actor by actor, so a window could start in one actor's late events and
  end in another's early ones.
* Markov-chain measures (`transition_entropy()`, `entropy_network()`,
  `passage_time()`, `markov_stability()`, `chain_structure()`) refuse
  correlation-type networks and negative weights (`nestimate_not_markov`)
  instead of row-normalising them into a meaningless chain. Co-occurrence
  networks (a random walk) are still accepted.
* `path_dependence()` no longer drops states whose names contain "NA"
  (e.g. "FINAL", "ANALYSIS"); missing cells are tested directly.
* `sequence_compare()` on a network whose states are numbers (1, 3, 5)
  reads them as states; they were decoded as positions, relabelling one
  state and dropping another.
* `build_mmm()`: BIC/AIC/ICL count covariate-dependent mixing by design
  column (a 4-level factor is 3 columns); and a sequence whose first state
  is missing adds no initial-state term instead of being scored as
  starting in the first state.
* `set_state_colors()` on a `netobject_group` stores the palette as an
  attribute; it was added as a list element, after which every plot,
  print and summary of the group failed.
* `sequence_plot()` heatmaps pin each state to its colour; when the drawn
  part held only some states (e.g. after `trim`), cells were coloured as
  other states. Panels with a single sequence (one row, a singleton
  cluster, one actor with `periods`) draw instead of erroring.

* `build_network()` refuses a `compare_mmm()` table with a classed error
  (`nestimate_bad_input`) giving the call to run, with the lowest-BIC k:
  `build_network(build_mmm(data, k = 2))`. The
  table is a data frame, so it used to be read as sequence data and
  "built" into a network whose states were its log-likelihoods and BICs.

* `permutation()` on cor, pcor and glasso networks pools the two datasets
  by column name. It pooled them by position, so the same variables in a
  different column order were mixed in the permutation null (one edge's
  p-value moved from 0.47 to 0.71 on identical data).
* `boot_glasso()` edge and centrality difference tests count a tie (both
  values equal, e.g. two absent edges) in both tails:
  p = 2 * min(P(d >= 0), P(d <= 0)), capped at 1. Ties used to count as
  evidence of a difference; with one edge present, 29 of 45 edge pairs
  were "significantly different" (now 9, the present edge against each
  absent one).
* Betweenness, edge betweenness and the shortest-path counts behind them
  use one relative tolerance for "same length". A tie that floating-point
  addition breaks (1/10 + 1/15 vs 1/6) was counted as one path but
  credited as a full path, so node b got 1 where igraph gives 0.5.
* Ising networks can be resampled: `bootstrap_network()`,
  `network_reliability()`, `centrality_stability()` and
  `casedrop_reliability()` returned NaN or CS = 0 for every Ising network,
  because the estimator rejected the stored data matrix and every draw
  failed silently. The resamplers now also warn
  (`nestimate_resample_failed`) with the number of draws whose network
  could not be estimated, instead of dropping them without a word.
* `build_hon()`, `build_hypa()`, `build_mogen()`, `path_counts()` and
  `state_frequencies()` read data frames column by column, so a factor
  column next to character columns keeps its labels instead of becoming
  states "1", "2", ...

* `build_network()` input handling:
  * `actor =` now names the ID column in wide and association data, as it
    already did in long data. Before, the ID values became nodes (wide
    sequences) or the ID column a variable (cor, pcor, glasso, Ising, which
    errored on a non-binary ID).
  * A numeric column in character sequence data (age, score) goes to
    `$metadata` instead of becoming states.
  * cor, pcor and glasso convert data to state frequencies only when it
    has at least two categorical columns and fewer than two numeric ones
    (sequences, possibly with one numeric covariate, which is left out).
    One character column among numeric variables (e.g. gender) used to
    turn the whole network into a network of that column's values; the
    numeric variables are now the nodes. Declared `metadata_cols` are
    dropped before this check.
  * `labels =` now also relabels `$data`, `$initial`, `$frequency_matrix`
    and `$node_groups`, so `bootstrap_network()` and the other resamplers
    no longer return all-zero matrices for a relabelled network, and
    `predictability()` no longer reports `NA` RMSE.
  * `prepare()`'s `unix_time_unit`, `is_unix_time` and `custom_format`
    passed to `build_network()` reach `prepare()` (they were silently
    ignored), and `prepare()` warns (`nestimate_time_unit`) when numeric
    times read in the given unit fall after the year 2500, the sign of
    millisecond timestamps read as seconds.
  * `level = "both"` keeps `actor` and `labels`; a grouped call keeps
    `labels`; an `NA` in the `group` column is an error
    (`nestimate_bad_input`) instead of a hidden all-`NA` group.
* Start/end markers (`start = TRUE`, `end = TRUE`) are now honoured by the
  precomputed resampling path behind `bootstrap_network()`,
  `permutation()`, `centrality_stability()`, `network_reliability()` and
  `casedrop_reliability()`. Before, the Start row and End column were zero
  in every replicate (bootstrap mean, SD and CI of 0; every Start edge
  "significant" in a permutation test).

* `build_mmm()` / `cluster_mmm()` with `covariates` and
  `covariate_effect = "em"`: the covariate-dependent mixing step is now an
  exact M-step, `nnet::multinom()` fitted to the posterior probabilities.
  The former hand-written Newton step had no step control; with `k >= 3`
  it diverged (coefficients around 1e7 to 1e9, an empty cluster, a
  log-likelihood below the fit without covariates) while reporting
  `converged = TRUE`. The likelihood now rises at every EM iteration, and
  `converged` also requires the last mixing step to have converged.
  Covariate fits change; `$covariates$beta` columns are now named.

* `nct()` compared the two groups' variables by position: a second group whose
  columns arrived in another order was tested against the wrong variables with
  no error. Columns that differ in name or order now raise
  `nestimate_bad_input`, as psychnets' `net_compare()` does.
* `nct()` failed with the uninformative `nearPD` error "infinite or missing
  values in 'x'" on any missing value, on a variable that is constant in a
  group, and on a permuted split in which a rare binary variable becomes
  constant. New `na_method = c("pairwise", "listwise")` (default
  `"pairwise"`: pairwise correlations, median co-observation count as the
  sample size; listwise drops incomplete rows, whole pairs when `paired =
  TRUE`). A constant variable enters the network as an isolated node (the
  psychnets convention), with a `nestimate_constant_variable` warning when it
  is constant in an observed group; the statistic is then defined on every
  permuted split, so the test stays exact. A pair sharing fewer than two
  observed rows raises `nestimate_degenerate_sample`.
* `nct()` drops non-numeric columns (an id or group label) instead of failing
  in `cor()`, validates `iter` as a single whole number (2.7 was silently
  truncated), and raises `nestimate_bad_input` for malformed arguments.
* `summary.net_nct()` labelled every edge `edge_1`, `edge_2`, ... with `to =
  NA` when the data had no column names; edges are now labelled `V1`, `V2`,
  ... `$E$edge_names` is always a `from`/`to` data frame (it was an
  `expand.grid()` result with `Var1`/`Var2`, or `NULL`).
* Printing an `nct()` result shows everything: both groups by name, M and S
  with each group's global strength, and every edge (largest difference
  first) with both groups' weights (`print(x, n = 10)` shortens the list) --
  the same layout as psychnets' `net_compare()`. New
  `as.data.frame.net_nct()`: one row per edge (`from`, `to`, `weight_1`,
  `weight_2`, `difference`, `p_value`, `significant`).
* `nct()` treats a network weight of |w| <= 1e-10 as no edge. psychnets
  before 0.6.0 can return an absent edge as a ~1e-17 residue, and two empty
  networks then tested as different (M = 5.6e-17 and a p-value below 1; 10 of
  40 random n = 60, p = 12 null datasets were affected with psychnets 0.5.2).
* On complete data every `nct()` number is unchanged (frozen in
  `tests/testthat/fixtures/nct-baseline.rds`).

# Nestimate 0.9.26

## Breaking changes

* `build_hypa()` now defaults to `p_adjust = "none"`, following LaRock et al.
  (2020), who classify paths against the threshold `alpha` without
  multiple-testing correction. Pass `p_adjust = "BH"` for the previous
  behaviour.

## Bug fixes

* `build_hypa()` used a miscalibrated null model. The propensity matrix was
  `s_out * s_in` restricted to observed paths, normalised over the whole
  graph and never refitted, so expected counts did not preserve node
  strengths: high-strength contexts were systematically under-represented and
  rare contexts over-represented, and memoryless (first-order Markov) data
  produced dozens of anomalies. Xi is now built on every possible De Bruijn
  edge (an observed k-gram continued by an observed first-order transition)
  and refitted until expected in- and out-strengths match the observed ones,
  following LaRock et al. (2020) and the reference `tlarock/hypa` code.

* `build_hypa()` fits the propensity matrix by continuous iterative
  proportional fitting instead of a rounded, keep-best sweep. Expected counts
  now preserve every De Bruijn node's in- and out-strength to numerical
  precision (order 2 equals the first-order Markov expectation exactly), and
  `converged` is reported only when the tolerance is met; otherwise a warning
  of class `nestimate_hypa_no_converge` is raised. Stalled fits were
  previously reported as converged. Expected counts can differ from the
  original Python implementation, which rounds inside every sweep, by a few
  thousandths.

* `build_hypa()` validates its input: `min_count` must be a single whole
  number, `alpha` a single number in (0, 0.5) and `order` free of duplicates
  (duplicates previously doubled every row and anomaly). Internal missing
  states now split a trajectory instead of becoming a state called "NA", and
  labels containing the internal separators `\x01`/`\x02` raise an error of
  class `nestimate_hypa_bad_label`.

* `build_hypa()` now uses the over-representation rule of LaRock et al.
  (2020): a path is over-represented when its HYPA score `P(X <= f)` exceeds
  `1 - alpha`, so `p_over` is the strict tail `P(X > f) = 1 - p_under`
  (previously the inclusive `P(X >= f)`).

* `build_hypa()` can now flag paths that never or rarely occur as
  under-represented. `min_count` gates over-representation on the observed
  count and under-representation on the expected count; previously a path
  observed fewer than `min_count` times could never be a deficit.
  Possible but unobserved paths are now scored (`observed = 0`), so
  `$scores` and `n_edges` contain more rows. At order 2 the expected count
  reduces to the first-order Markov expectation.

* `print.q_analysis()` reported the fragmentation point from the highest q
  down, because `q_vector` is stored from the highest q; it now reports the
  lowest q with more than one component (e.g. q = 2, not q = 4, for a clique
  complex whose maximal simplices meet only in triangles).

## Improvements

* `verify_simplicial()` accepts a network object (`netobject`, `tna`,
  `cograph_network`, `net_hon`) as well as a matrix, as `build_simplicial()`
  does.

# Nestimate 0.9.24

## Documentation

* Fixes to the 0.9.23 help-page merge, found in review:
  * Pages are named after their function again: `build_mcml` and
    `plot_state_frequencies` had been named after a `print()` method, and
    their usage now lists the function before its methods.
  * The descriptions of 51 methods (what a plot draws, what a print shows)
    are restored, as a "Methods" section on each function's page.
  * Each page's `...` entry states what each method does with extra
    arguments: methods that reject them say so, `plot.net_bayes()` names
    `cograph::splot()`, and the rest say they are ignored.
  * `compare_networks`: the `...` and `labels` entries no longer repeat;
    `net_edge_betweenness` and `net_pruning_details` no longer list a
    non-existent class.
  * Every S3 method has a documented return value.

## Behaviour notes

* The `mgm` estimator rejects `scale = FALSE` with a classed error
  (`nestimate_mgm_unscaled`); this started in 0.9.23 and was not listed
  there. Missing values use mean/modal single imputation under psychnets'
  named `"pairwise"` policy, now passed explicitly.

## Tests

* Two more reference-implementation tests (a hand-written linear solve for
  `passage_time()`, a manual `multinom` re-computation) and unused reference
  helpers moved to the local suite; section headers that claimed
  equivalence tests are renamed. New tests: mgm on data with missing values,
  and the `scale = FALSE` error.

# Nestimate 0.9.23

## Internal changes

* The `mgm` and `ising` estimators delegate to `psychnets::mgm_fit()` and
  `psychnets::ising_fit()` (glmnet engine), joining `cor`, `pcor` and `glasso`.
  Results are unchanged: identical to the former in-package code (max absolute
  difference 0 for Ising over 180 configurations, 5.6e-17 for MGM over 48).
  Nestimate keeps type detection, input validation and the result fields. The
  Ising result no longer carries the unused `asymm_weights` and
  `lambda_selected` fields. The unused internal moderated-MGM code is removed.

## Documentation

* Each function has one help page. The 181 `print()`, `summary()`, `plot()`,
  `format()` and `as.data.frame()` methods are documented on the page of the
  function that creates their object (`?nct` covers `print.net_nct()` and
  `summary.net_nct()`), as is `compare_model()`'s `netobject_group` method.
  286 help pages become 139; `?print.net_nct` still opens the right page.
* `compare_model()` and `as_netdifference()` have examples.

## Bug fixes

* `build_network(method = "mgm")` no longer fails on data with missing values
  ("number of observations in y not equal to the number of rows of x").

## Tests

* Equivalence tests against other packages (tna, igraph, mgm, stringdist,
  `glm()`) and naive-loop reference implementations are no longer shipped.
* CRAN runs a core of 15 test files (data preparation, estimation registry,
  extraction, Markov, link prediction, Bayesian verbs, pruning, and the htna
  contract); every other file is skipped on CRAN and runs locally and in CI.
  The CRAN test run drops from about 300 s to 22 s.
* Slow examples run with fewer iterations (`boot_glasso()`, `nct()`,
  `permutation_diagnostics()`, `build_gimme()`) and `?sequence_plot` draws
  fewer figures.

# Nestimate 0.9.21

## Bug fixes

* `summary()` of a `predict_links()` result reports `NA` scores for a method
  with no predictions (every possible link already exists). It returned
  `NaN` / `-Inf` / `Inf` with one warning per column.

## Documentation

* Every help page now has example code that runs under `R CMD check`. 31 pages
  had all of it inside `\donttest{}` or `\dontrun{}`; fast examples are
  unwrapped, and slow ones (`nct()`, `permutation_diagnostics()`,
  `entropy_bayes()`) run live with fewer iterations.
* The `predict_links()` family of examples uses fixed data that leaves links
  to predict, and the `loading_stability()` / `as_networks()` examples use items
  with real two-factor structure, so neither emits warnings.
* pkgdown: the `cograph-tutorial-nestimate` article was a 47.8 MB page (300 dpi
  figures inlined as base64); it is now about 70 KB with figures as separate
  files. The README no longer lists outdated imports.

# Nestimate 0.9.20

## Documentation

* The documentation site is now `https://pak.dynasite.org/Nestimate/`
  (`DESCRIPTION` URL, `_pkgdown.yml`, README links). `saqr.me/Nestimate` was a
  stale copy. The README lists the `permutation-nested-data` and
  `compare-networks` articles.

# Nestimate 0.9.18

## Documentation

* The pkgdown site builds again. Four exported topics (`state_colors()`,
  `set_state_colors()`, `composites()`, `item_loadings()`) were missing from the
  reference index, which had failed every site build since 0.9.10. The
  `compare-networks` and `permutation-nested-data` articles now publish.

# Nestimate 0.9.17

## New features

* `compare_networks()` gains `actor =`, passed to the permutation backend:
  whole actors are reassigned between the networks, the printed header names
  the actor column, and `global_differences()` gains the rows `ICC`,
  `Design effect (edges)` and `Design effect (M)` (category `"Nesting"`).
  Requires `test = "permutation"` (error class
  `nestimate_compare_actor_needs_permutation`).

# Nestimate 0.9.16

## Breaking changes (development versions only)

* The nesting argument added in 0.9.13 (`permutation()`,
  `permutation_diagnostics()`) and 0.9.15 (`bootstrap_network()`) is renamed
  from `block` to `actor`, the column identifying whose sequences they are:
  the same vocabulary as `build_network(actor = )`. Results carry `actor`
  and `n_actors`; diagnostic columns use `_sequence` (sequences reassigned)
  and `_actor` (actors reassigned), e.g. `p_global_sequence` /
  `p_global_actor`; conditions are `nestimate_bad_actor`,
  `nestimate_actor_missing`, `nestimate_actor_misaligned`,
  `nestimate_actor_unsupported` and `nestimate_few_actors`. The printout
  reads `Actor: Group (200 actors)` and `Nesting in Group: ICC = ...`.
  Neither name was ever on CRAN.

## Bug fixes

* `build_network()`: columns named in `metadata_cols` (or left out of
  `state_cols`) are no longer read as sequence positions in wide data. They
  were moved to `$metadata` only after estimation, so their values became
  states and created spurious transitions, in single and grouped networks.
  For `cor`, `pcor`, `glasso`, `ising` and `mgm`, declared metadata columns
  are no longer entered as variables.

* An ICC that cannot be estimated (one block, or no variation) is reported
  as not estimable instead of `NaN [NA, NA]` by `permutation()`,
  `bootstrap_network()` and `permutation_diagnostics()`.

# Nestimate 0.9.15

## New features

* `bootstrap_network()` gains `block =`, the column identifying the unit
  that sequences are nested in. Whole units are resampled with replacement,
  keeping their sequences together (cluster bootstrap, top level only;
  Davison & Hinkley, 1997; Field & Welsh, 2007). The result reports the ICC
  (Shrout & Fleiss, 1979) and the design effect (Kish, 1965) for the edges,
  printed under the bootstrap summary and returned as `clustering` and
  `clustering_edges`. With one sequence per unit it reproduces the ordinary
  bootstrap exactly. Transition networks only.

# Nestimate 0.9.14

## Documentation

* New pkgdown article "Comparing two networks with compare_networks()":
  high vs low achievers in `group_regulation_long`, every summary table and
  all eight plot views with a short reading of each.

# Nestimate 0.9.13

## New features

* `permutation()` gains `block =`, the column identifying the unit that
  sequences are nested in (sessions in students, students in teams). Whole
  units are reassigned between the groups; units with sequences in both
  groups, as in repeated-measures designs, have their sequences reassigned
  within the unit; mixed designs combine the two (Good, 2005; Anderson &
  ter Braak, 2003). Supported for transition networks (`relative`,
  `frequency`, `co_occurrence`); association networks raise
  `nestimate_block_unsupported`. A warning of class `nestimate_few_blocks`
  is raised when no p-value could fall below `alpha`.

* With `block`, `permutation()` reports the nesting effect: the intraclass
  correlation (Shrout & Fleiss, 1979) with a jackknife 95% interval, and the
  design effect (Kish, 1965), the blocked over the unblocked null variance,
  for the edges and the global statistic M.

* `permutation_diagnostics(x, block =)` runs the ordinary and the blocked
  test on the same data and returns them side by side as a tidy data frame,
  one row per pair or, with `level = "edges"`, per edge.

## Improvements

* `print()` of a `net_permutation` shows the global test (M and S, with
  p-values) and, with `block`, the ICC and design effects. A grouped result
  prints every pair in full.

* The `permutation()` help page is reorganised into sections: what is
  tested, nested data and `block`, ICC and design effect, reading the
  printed output.

* Building a transition network from a single long sequence now raises a
  message instead of a warning; running `bootstrap_network()` or
  `permutation()` on such a network still warns. Both carry the class
  `nestimate_single_sequence`.

## Bug fixes

* `\%` in markdown roxygen silently cut the rest of a line from the help
  pages of `certainty()`, `build_gimme()`, `build_mcml()` and
  `sequence_plot()`; the text is restored.

* Missing space in the `permutation()` print header.

# Nestimate 0.9.12

## New verbs

* `set_state_colors(x, colors)` attaches a palette to a `netobject`,
  `netobject_group`, `mcml` or `htna`, and every figure drawn from that object
  then uses it -- `sequence_plot()`, `distribution_plot()`,
  `plot_state_frequencies()` **and** `cograph::splot()`. The cograph half goes
  through the documented `meta$splot` producer contract (`defaults$node_fill`,
  stamped in node order), so `splot(net)` needs no arguments and renders
  byte-identically to `splot(net, node_fill = ...)`. A `state_colors(x) <- `
  replacement form is also available.

* `state_colors(x)` reads the resolved palette back as a tidy data.frame: one
  row per key, with `state`, `color` and `source` (`"set"` or `"default"`). It
  is what the figures actually draw, not just what was passed in.

## Improvements

* States a named palette does not mention are now dealt the Okabe-Ito colours
  the palette has **not** used, one each, instead of taking their positional
  slot. Pinning `plan = "#0072B2"` no longer leaves another state defaulting to
  that same blue. Affects the named form only; an unnamed `state_colors` and
  the all-default palette are unchanged.

* An attached palette resolves **once**, in the object's own key order, so
  every figure drawn from it colours a state identically. Previously only the
  colours you set were carried and each figure dealt the remaining defaults in
  its own sort order -- `plot_state_frequencies()` sorts states by frequency
  and `sequence_plot()` alphabetically, so one state could come out amber in
  one figure and black in the other.

* In-tile labels in `plot_state_frequencies()` take whichever ink (dark or
  white) has the higher WCAG contrast ratio against the tile, instead of one
  fixed grey. A percentage on a dark tile -- black, navy, dark wine -- is now
  readable, whether the colour came from the default palette or from a palette
  you set.

# Nestimate 0.9.11

## Fixes

* A named `state_colors` no longer has to match the figure exactly. Names the
  plot does not draw are dropped with a `message()` naming them, so one
  project-wide palette -- states from a wider coding scheme, clusters a given
  call merged away -- can be handed to every plot and each takes the keys that
  apply to it. 0.9.10 raised an error instead, which made a shared palette
  unusable.

# Nestimate 0.9.10

## Improvements

* `state_colors` accepts a **named** vector everywhere `sequence_plot()`,
  `distribution_plot()` and `plot_state_frequencies()` draw, and the names are
  a lookup rather than a positional list: only the keys you name are
  overridden, every other key keeps its default, and the vector may be shorter
  than the number of states. An unnamed vector is still positional, unchanged.

* For an `mcml`, that lookup covers the whole multichannel figure, not only
  the states. A cluster name colours its `Summary` band, its channel strip and
  its faded band in the other panels; a group merged by `combine =` is named
  by its label (the list name, or `"A + B"`); `rest_label` is a key too. So
  `sequence_plot(fit, combine = list(Task = c("Cognitive", "Regulation")),
  state_colors = c(Task = "#0072B2"))` recolours the combined cluster and
  leaves the rest of the palette alone. Previously the cluster and combined
  keys were drawn from a fixed internal palette with no way to set them.

* A `state_colors` name that matches no key in the figure is now an error
  naming the available keys, instead of being silently ignored.

# Nestimate 0.9.9

## New verbs

* `build_mcml()` gains `combine =` and `expand =`. On new input they change
  the partition before estimation, in any `clusters` form and for every input
  type: `build_mcml(data, clusters = cl, combine = c("A", "B"))` equals a build
  with `A` and `B` merged in `cl` (the merged cluster lists its states in
  cluster-name order). On an existing `mcml` they re-partition it:
  `build_mcml(mc, combine =)` merges clusters, `build_mcml(mc, expand =)`
  splits clusters into one cluster per state, and `build_mcml(mc, clusters =)`
  applies a new partition. The
  model (macro network, within-cluster networks, sequences) is re-estimated
  from the sequences the `mcml` carries, with its original `type`, `method`
  and `directed`. With the partition unchanged the result equals the input;
  `expand = "all"` reproduces the node-level transition network of the same
  `type`. Raises `nestimate_mcml_no_sequences` for an `mcml` built from a
  matrix or an edge list (which keeps only within-cluster edges), and errors
  when sequence-shaping arguments (`trim`, `exclude`, `end`, `labels`,
  `actor`, ...) are passed with a re-partition, since the carried sequences
  already reflect them.
* `session_ids()` names the session behind every sequence of a network built
  from long data, and of a `build_mmm()` or `build_clusters()` fit on such a
  network. It returns one row per sequence in model order: `sequence`, the
  `actor` and `session` columns under their own names, `session_label`, and,
  for a fit, `cluster` (plus `posterior` for a mixture). The fit's assignments
  can then be joined to the input by those columns instead of parsing the
  `"actor | session"` label, which breaks when an id contains `" | "`. Raises
  `nestimate_no_session_ids` for wide-data input or for fits made before this
  version, and `nestimate_session_ids_misaligned` when the metadata and the
  sequences differ in number.
* `item_loadings()` returns the tidy item-diagnostic table of a
  `build_mcml_pc()` fit (node, cluster, loading, weight, sign, max_cross,
  cross_cluster, misfit); `misfit = TRUE`/`FALSE` filters it.
* `composites()` returns the per-respondent cluster scores of a
  `build_mcml_pc()` fit: one row per input row (input order and row names,
  `NA` where all of a cluster's items are missing) and one column per
  cluster. Raises `nestimate_no_composites` for the descriptive aggregations
  (`"average"`, `"escoufier"`, `"cancor"`), which form no score.

## Changes

* `sequence_plot()` on an `mcml` gains `combine =`: named clusters are merged
  into one channel (a character vector for one group, a list for several;
  list names label the merged channels, default `"A + B"`). The merged group
  acts as one cluster across the figure (one panel, one Summary key, one
  faded band) and can itself be opened with `expand =`.
* `sequence_plot()` on an `mcml` gains `rest = c("clusters", "pooled",
  "none")`: how a cluster's panel shows the time spent in other clusters (one
  faded band per cluster, one pooled grey band, or blank, leaving
  only the panel's own states). Applies to the carpet and distribution views.
  `rest_label =` (default `"Other states"`) sets the legend text: the pooled
  band takes it as is, per-cluster bands read `"Social (Other states)"`. This
  replaces the former `"(elsewhere)"` wording.
* `sequence_plot()` on an `mcml` now honours `na =` in the distribution view:
  `na = FALSE` drops the `NA` (ended) band and shows each time point as shares
  of the sequences still running, as `distribution_plot()` already did.
* `macro_network()` accepts an `mcml_pc` fit and returns its cluster-level
  network; `expand =` on an `mcml_pc` raises `nestimate_no_expand`, and
  `method =` or `...` error (the estimator is set in `build_mcml_pc()`).
* `sequence_plot()` errors when `combine`, `expand`, `rest` or `rest_label`
  is passed for input that is not an `mcml`, instead of ignoring them.
* `print.mcml_pc()` and its build-time warnings name `item_loadings()`
  instead of pointing at `$loadings`.

* `prepare()` (and so `build_network()` on long data) keeps the `session`
  column(s) in the per-sequence metadata under their own names, and returns
  the metadata explicitly in sequence row order (it was assembled with
  `merge(sort = FALSE)`, whose order is unspecified).
* `build_mmm()` and `build_clusters()` keep the input network's `$metadata`
  (restricted to the fitted rows when `build_mmm()` drops sequences with
  missing covariates).

## Fixes

* `sequence_plot()` on an `mcml` with `type = "heatmap"`/`"index"` and
  `expand =` drew the other-cluster wash as blank cells: the wash was keyed by
  the Summary keys, which are states once a cluster is expanded. It is now
  keyed by cluster.
* `sequence_plot()` on an `mcml` with `type = "distribution"` and `expand =`
  no longer fails with "subscript out of bounds" in the default
  (`normalize = FALSE`) view. The Summary band now opens the expanded
  cluster into its states, and the other panels draw it as one faded
  `"<cluster> (<rest_label>)"` band.
* `sequence_plot()` on an `mcml` with a single channel (one cluster, or every
  cluster merged by `combine`) no longer fails in the carpet view with
  "replacement has 1 row, data has 0".
* `sequence_plot()` on an `mcml`: with `expand =`, each cluster's faded band
  in the other panels now has its own colour (expanded clusters all shared
  one).

# Nestimate 0.9.5

## New verbs

* `macro_network()` returns the macro (cluster-level) network of an `mcml`
  with one or more named clusters expanded back into their member states and
  every other cluster left collapsed — a network at mixed resolution. It
  re-counts from the recoded sequence data rather than splitting the k x k
  aggregate, which cannot be disaggregated; a matrix-derived `mcml` raises
  `nestimate_no_expand_source`. `$node_groups` maps each expanded state back
  to its parent cluster, so the result plots grouped.
* `extract_pathways()` cuts a long event log into pathways and returns one
  row per pathway. Three cuts via `type =`: `"unit"` (one pathway per group),
  `"segments"` (one per contiguous run), `"anchored"` (spans around an anchor
  event). `resolve =` appends a resolution label as the closing state.
* `outcome_model()` regresses a unit-level outcome on sequence or network
  predictors — pattern indicators, `simplicial_features()` output, or any
  numeric covariate — with family auto-detection, an optional `lme4` random
  intercept, `select = "split"` hold-out selection, BH-corrected p-values,
  confidence intervals and odds ratios. `effects_table()` is the tidy
  accessor; `summary()` returns the same table.
* `simplicial_features()` returns topological summaries of one or many
  networks as a tidy data.frame, one row per network per threshold, ready to
  use as regression predictors. Accepts a `netobject`, `netobject_group`,
  `mcml`, a square weight matrix, or a named list of any of these.

## Extended

* `build_mcml()` gains the `exclude`, `trim`, `end` and `end_by` sequence
  arguments, applied in that fixed order.
* `as_tna()` is now a generic. `as_tna.mcml(expand =)` delegates to
  `macro_network()` for the macro layer only; per-cluster layers are
  untouched. An `as_tna.default()` covers everything else.
* `sequence_plot()` and its `mcml` method gain
  `panel = c("both", "summary", "channels")`.

## Fixes

* `frequencies(format = "frequency")` returned zero rows for any long-format
  data whose id column was not integer-valued. The sequence key was coerced
  with `as.integer()` before the merge, turning every character id into `NA`;
  the merge then matched nothing and failed silently. Sequence keys are now
  aligned in their character form.
* `sequence_plot()` on an `mcml` errored whenever a singleton cluster was
  named after its own state: the shared fill scale built its levels with
  `factor(levels = c(states, clusters))`, and duplicate levels are an error.
  Levels and values are now de-duplicated.
* `as.data.frame()` methods for `net_hypergraph_transduction` and
  `net_hypergraph_cluster` placed `what` in the generic's `row.names` slot,
  so `as.data.frame(x, "scores")` bound `"scores"` to `row.names`. The
  generic's arguments now come first and `row.names` is honoured
  (`R CMD check` "S3 generic/method consistency").

# Nestimate 0.9.4

* pkgdown reference index repaired: the five topics added in 0.9.3 were
  unindexed and `build_reference_index()` failed. The S3 Methods section now
  wildcards `as.data.frame.` alongside `print.`, `summary.` and `plot.`, so a
  future result class with the house-standard accessor cannot repeat it.

# Nestimate 0.9.3

## New verb: `compare_networks()`

* Compares two or more networks (`netobject`, `netobject_group`, `tna`,
  `group_tna`, `mcml`, `cograph_network`, matrices, or a list of them) in one
  object. `summary()` returns the tidy one-row-per-pair overview, as every
  other Nestimate `summary()` does; the full tables are their own verbs:
  `edge_differences()`, `node_differences()`, `global_differences()`,
  `network_metrics()`, each with a `pair` column. All pairs by default, or
  every network against one `reference =`.
* Ratios are guarded (`NA`, never `Inf`/`NaN`); every cell is kept, absent
  edges appear at weight 0; `higher` (colour) and `status` (shape) columns
  drive the plots.
* `plot(x)` draws one view per call, selected with `type =` (the package-wide
  argument name; `what =` is accepted as an alias): `"networks"` (default;
  each network drawn once via `cograph::splot()`), `"difference"` (the signed
  difference network of each pair), `"edges"` (ranked dumbbell), `"nodes"`,
  `"global"`, `"heatmap"`, `"scatter"`, and `"inference"` (a forest of the
  edge differences on the difference scale, with credible intervals when the
  Bayesian backend ran, a filled marker for a significant edge and the
  p-value in a right-hand column; it raises
  `nestimate_compare_no_test` when `test = "none"`). One sign-to-colour
  contract (`#4A6FE3` = first/reference higher, `#D33F6A` = second higher)
  backed by line type and marker shape; the heatmap scale follows the data.
* `plot(x, combined = FALSE)` splits a multi-pair view into a named list of
  single-pair plots (one per pair) instead of facetting them into one figure,
  matching the `combined` argument of `plot.net_reliability()`,
  `plot.simplicial_complex()` and friends. The base-graphics views draw one
  panel per page.
* Optional inference on the same tables: `test = "permutation"`, `"bayes"`
  (with an optional `rope`), `"bootstrap"`, combinable. Non-significant
  results are faded, never deleted.
* `permutation()` gains a `$global` data frame with NCT-style `M` (sum of
  absolute differences) and `S` (largest absolute difference) statistics and
  their permutation p-values, computed from the same null.
* `compare_model()` is unchanged and will be soft-deprecated once
  `compare_networks()` has been through one release.

## Hypergraph suite

* `build_hypergraph()` promotes a network's k-cliques (k >= 3) to hyperedges,
  following Burgio, Matamalas, Gomez and Arenas (2020); underlying pairwise
  edges are retained. `clique_expansion()` projects a hypergraph back to a
  pairwise network in one `tcrossprod()` call, closing the
  event data -> `bipartite_groups()` -> hypergraph -> network cycle.
* `hypergraph_measures()` returns the structural-statistics suite (Lee, Choe
  and Shin 2024); an empty hypergraph returns trivial zeros rather than
  erroring. `hypergraph_centrality()` computes Benson's (2019) three
  eigenvector centralities.
* `hypergraph_laplacian()` computes the normalized Laplacian, either the
  Zhou, Huang and Scholkopf (2006) form on the binary incidence pattern or
  the Hayashi, Aksoy, Park and Park (2020) weighted form. Built on it:
  `hypergraph_cluster()` (spectral clustering) and
  `hypergraph_transduction()` (semi-supervised label propagation), each with
  print, summary, plot and `as.data.frame()` methods.

# Nestimate 0.9.1

tna-parity release: sequence-side gaps against `tna::build_model()` and
`tna::prepare_data()` closed.

## New estimators

* `build_network(method = "ngram")` (params `n_gram`, default 2),
  `"gap"` (params `max_gap`, default 1) and `"reverse"` (params `weighted`)
  mirror the tna types `"n-gram"`, `"gap"` and `"reverse"`. Weights are
  numerically identical to tna on shared inputs (tested). All three accept
  wide or long input, `start`/`end` boundary markers, `group =` dispatch and
  `scaling = "normalize"` for row probabilities.
* New aliases: `"n-gram"` / `"n_gram"` -> `"ngram"`, `"co-occurrence"` ->
  `"co_occurrence"`.

## Timezone-safe timestamp parsing

* `prepare()` and `build_network()` gain `timezone = "UTC"` (Olson name).
  Naive timestamps are interpreted in that zone; ISO-8601 `Z`/`UTC`/`GMT`
  markers and numeric offsets (`+0200`, `+02:00`) are converted from their
  offset. Parsing no longer depends on the machine's local time zone.
* Fixes: `...Z` timestamps parsed to `NA` (the marker was stripped before
  matching), offsets were silently ignored, and wall-clock readings inside a
  DST gap could turn into spurious session boundaries. Unparsable values now
  raise an error instead of becoming session breaks.

# Nestimate 0.9.0

Delegation release: Nestimate stops owning psychometric-network math and
delegates it to its two clean-room home packages, `psychnets` (cross-sectional)
and `idiographic` (temporal). No public API changes: every signature, default,
return shape, and — verified against frozen pre-delegation baselines — every
number is preserved.

## Delegated to psychnets (new hard dependency, >= 0.5.2)

* The `cor`, `pcor`, and `glasso` estimators keep Nestimate's input layer,
  validation, and return contract; the math now runs in
  `psychnets::cor_network()` / `pcor_network()` / `ebic_glasso()`. Verified
  field-complete against the frozen baseline: `cor` exact, `pcor` within
  5.6e-17, `glasso` exact on weights, precision, selected lambda, and the
  EBIC path (including `penalize.diagonal` and `refit` branches).
* `nct()` delegates its inner EBIC-glasso solve; the NCT-specific `nearPD`
  symmetrization stays local. Seeded runs: networks and p-values exact.
* `boot_glasso()` and `permutation()` (glasso branch) delegate the
  per-resample solve via `psychnets::ebic_glasso(lambda_path = )`, keeping
  the fixed-path-across-resamples semantics. Seeded `permutation()` is
  byte-identical; `boot_glasso()` is identical except wall-clock timing.
* The internal pure-R glasso kernel (`glasso_pure.R`) and its EBIC helpers
  are deleted; psychnets owns that math now.
* The `ising` and `mgm` estimators remain local for now (their delegation
  needs psychnets to expose per-node lambda selection and a scale
  passthrough) and are unchanged.

## Delegated to idiographic (new hard dependency, >= 0.3.4)

* `build_mlvar()` keeps its signature, S3 methods, and return object; the
  lmer estimation pipeline now runs in `idiographic::fit_mlvar()`. Verified
  identical at tolerance 0 (object, class, print, summary) across
  lag/standardize/day/beep configurations, and against `mlVAR::mlVAR()` at
  8.8e-16 over 954 checks.
* When the between-subjects network is not estimable (a random-intercept SD
  of zero), `build_mlvar()` now warns before returning the zero matrix; it
  previously returned it silently. Numbers are unchanged.
* `build_gimme()` delegates its whole search to `idiographic::fit_gimme()`.
  **This is the one delegation that changes results**: idiographic's search
  reproduces the upstream `gimme` package (>= 10.0) exactly, which
  Nestimate's own search did not — the individual-level path search in
  particular under-detected person-specific paths relative to upstream.
  Group-level results are typically unchanged; individual-level path sets
  can grow. The signature and the `net_gimme` field contract are unchanged
  (the object gains idiographic's netobject fields and now renders directly
  with cograph); `print`/`summary`/`plot` dispatch to idiographic's
  methods. Treat pre-0.9.0 `build_gimme()` individual-level results as
  superseded.

# Nestimate 0.8.5

Bug-fix release. 0.8.4 was published on r-universe but never reached CRAN.

## Plotting in non-UTF-8 locales

* `plot()` on `net_entropy_bayes` and `net_sequence_comparison` results no
  longer fails where the graphics device cannot represent the arrow glyph.
  Edge labels and the comparison subtitle previously used `U+2192`/`U+2190`,
  which a non-UTF-8 locale cannot convert; the graphics engine raised an error
  (`conversion failure ... in 'mbcsToSbcs'`) rather than substituting a
  character, so the plot could not be drawn at all. Both now use the ASCII
  arrow notation (`A -> B`) already used by the higher-order network verbs.
  Console output is unaffected.

* `plot()` on `net_entropy_bayes` results no longer emits a deprecation warning
  under ggplot2 4.0, which removed `geom_errorbarh()`.

## Rendering entropy networks with cograph

* `entropy_network()` now renders with its intended transition-network styling
  on cograph 2.4.4, the current CRAN version. The object states its full style
  through cograph's `meta$splot` producer contract instead of relying on
  cograph recognising the `"entropy"` method name, so no cograph update is
  required. Output is identical under cograph 2.4.4 and 2.4.5.

## Documentation

* The `print()` and `plot()` examples for the MMM clustering attribute now
  reach it through `build_network()`, matching the 0.8.4 change that made
  `cluster_mmm()` return the fitted `net_mmm` object.

# Nestimate 0.8.4

## Transition matrix entropy

* New entropy suite around the entropy rate of a transition matrix
  (Shannon 1948; Cover & Thomas 2006, ch. 4; the transition entropy of
  Krejtz et al. 2015 and the real-time mobile transition matrix entropy of
  Krejtz et al. 2025):
  - `transition_entropy()` — entropy rate H (stationary distribution from
    the eigendecomposition at lambda = 1), stationary entropy, redundancy,
    per-state branching entropies, all with the normalized (scale-free)
    variants; print/summary/plot and `netobject_group` dispatch.
  - `entropy_network()` — the exact edge-level decomposition of H: each
    edge carries its term pi_i P_ij log(1/P_ij); weights sum to H. Displays
    the summands of the entropy-rate equation — no new quantity is
    estimated. `scaling = "share"` relabels edges as percentages of H
    (sum = 100). Additional weights: `"surprisal"` (optionally
    `scaling = "chance"`) and `"production"` (irreversibility). The result
    is a regular `netobject` that inherits its source network's styling and
    ships an entropy house style via cograph's `meta$splot` producer
    contract, which states the styling outright rather than relying on
    cograph recognising the method name (cograph >= 2.4.4).
  - `entropy_trajectory()` — sliding-window entropy over the transition
    stream (the windowed design of Krejtz et al. 2025): tidy per-window
    table, per-group trajectories, loess-trend plot; the per-window
    estimator weights rows by observed occupancy (robust to non-ergodic
    window fragments).
  - `entropy_bayes()` — Dirichlet-posterior estimation: credible intervals
    for H, per-state entropies, and per-edge contributions; edges flagged
    credible when their share of H credibly exceeds `min_share`; `$model`
    holds the pruned stable entropy network.
* New vignette: *Transition Matrix Entropy*.

## Mixed Markov clustering contract

* `cluster_mmm()` now returns the fitted `net_mmm` clustering object, retaining
  assignments, posterior probabilities, mixing proportions, fit criteria, and
  fitted component models. Network materialization remains the responsibility
  of `build_network(fit)` or the one-step
  `cluster_network(..., cluster_by = "mmm")` workflow.
* `as_htna()` gains a `net_mmm` method. An MMM fit created from an HTNA model is
  materialized into an `htna_group` without rerunning the MMM fit, while the
  original actor partition and clustering diagnostics are preserved.

# Nestimate 0.8.3

## HTNA expansion

* `as_htna()` still rebuilds one full node-level network from the original
  source, preserving every between-cluster transition, and now completes the
  canonical HTNA contract. Its result inherits from `htna`, `netobject`, and
  `cograph_network`; stores character actor labels in `$node_groups$group`, a
  factor in `$nodes$groups`, and actor order in `$actor_levels`; and retains
  actor-order metadata on `$node_groups` for lossless partition round trips,
  plus the legacy `$nodes$cluster` and `"cluster_members"` metadata.

# Nestimate 0.8.2

## Session grouping

* `prepare()` now identifies sessions from the observed combinations of the
  `actor` and `session` columns instead of `base::interaction()`. Three
  defects are fixed:

  - **Integer overflow.** `interaction()` codes a combination over the
    marginal level space, which exceeds `.Machine$integer.max` once both
    columns pass 46,341 distinct values. The resulting `NA`s were pasted into
    the literal string `"NA"`, merging unrelated events into one pseudo-session.
    On 50,000 actor-session pairs this silently discarded 14% of sessions and
    manufactured transitions that no input sequence contained, including
    self-loops on terminal states.
  - **Separator collisions.** Identifiers were pasted with `" | "` before being
    used as a grouping and metadata merge key, so `("a | b", "c")` and
    `("a", "b | c")` collapsed into a single session. Grouping now keys on the
    original columns; the readable label is display-only and is exposed as
    `.session_label` in `meta_data`.
  - **Missing identifiers.** Missing values in a grouping column silently
    changed the session count. They now raise an error naming the columns.

  Group numbering reproduces `interaction()`'s ordering, so prepared row order
  and finite same-seed bootstrap results are unchanged.

* `time_threshold = FALSE` switches session-interval splitting off, so each
  actor (or actor-session) forms a single sequence regardless of gap length.
  Accepted by `prepare()`, `build_network()` and `build_mcml()`.

## Testing

* Added a randomized grouping sweep over 1,000 seeded datasets covering
  actor/session column counts, identifier vocabularies including separator and
  UTF-8 cases, time on and off, tied timestamps, explicit order columns,
  shuffled input, and time-gap structures that split, do not split, or sit
  exactly on the threshold. Verified against ground truth and against
  `interaction()` for row-order compatibility.

# Nestimate 0.8.1

## HTNA interoperability

* Distance clustering and mixed-Markov clustering now preserve HTNA inputs.
  `build_clusters()` and `cluster_mmm()` carry the node-to-actor partition into
  network materialization, while `cluster_network()` returns an `htna_group`
  directly. Every child remains an `htna` object with
  `$node_groups`, `$nodes$groups`, and `$actor_levels`; clustering assignments,
  posterior probabilities, fit diagnostics, and other outer attributes remain
  attached.

* Added extensive randomized equivalence coverage across distance clustering,
  mixed-Markov clustering, all transition-network estimators, actor-absent
  clusters, and HTNA-versus-plain input paths.

# Nestimate 0.8.0

## Documentation

* The Bayesian verbs get their own reference section, placed directly after
  Network Estimation: `certainty()`, `bayes_compare()`, `subtract_networks()`
  and `as_netdifference()`. They were previously buried in a fourteen-entry
  "Bootstrap & Inference" list. `bootstrap_network()` now points at
  `certainty()` as its closed-form counterpart, and `permutation()` points at
  `bayes_compare()` as its Bayesian complement, so each pair is reachable from
  either side.

* `frequencies()` is no longer marked `\keyword{internal}`. The topic page and
  the exported function share a roxygen topic name, so the keyword from the
  topic block leaked onto the function's own help page even though the function
  is exported (and called by the \pkg{htna} package). `cluster_data()` keeps its
  internal keyword: it is a deprecated alias for `build_clusters()` and is meant
  to stay out of the index.

* Dropped the `utils` help page, which documented no exported object. The
  `@importFrom` directives it carried are retained.

* `audit_codex/` is no longer tracked; it holds generated audit artifacts.

## Dependencies

* `Suggests: cograph (>= 2.4.4)`. The netdifference verbs added in 0.7.8 need
  cograph 2.4.x: CRAN's cograph 2.3.6 contains no `netdifference` support, so
  `cograph::plot_difference()` does not exist there and `cograph::splot()` on a
  `netdifference` falls through to the plain `netobject` renderer and silently
  draws an unsigned network. Nestimate must not be submitted to CRAN before
  cograph 2.4.4 is available there.

# Nestimate 0.7.8

## New features

* `subtract_networks()` / `as_netdifference()` — verbs for the difference
  between two networks. `subtract_networks(x, y)` returns the edge-wise
  difference as a `netdifference` object; `as_netdifference()` promotes an
  existing comparison result to the same class — a `bayes_compare()` result,
  or a `netdifference`, which passes through; anything else errors — so a
  difference computed by any route prints the same way. Adds
  `print.netdifference`.

* `bayes_compare()` accepts two `net_edge_betweenness()` objects (source
  method `"relative"` only). Edge betweenness is recomputed on every
  posterior draw, giving the Bayesian analogue of `permutation()`'s
  edge-betweenness dispatch, with posterior mean betweenness matrices and
  the plug-in `observed_diff`.

* `permutation()` gains a `measures` argument for centrality permutation
  tests, matching the `tna` package's dispatch.

## Enhancements

* `bayes_compare()`'s probability-of-direction column is renamed `pd` ->
  `p_difference` in the `summary()` frame, and the result now carries class
  `c("net_bayes", "netdifference", "net_permutation")` so it dispatches to
  the difference verbs as well as the permutation ones.

* Non-ASCII characters normalized across R sources and man pages.

## Bug fixes

* `centrality_stability()` no longer errors with
  "missing value where TRUE/FALSE needed" when a requested measure is
  undefined on the network (e.g. `Diffusion` is `NaN` on a small cyclic
  net): `sd()` returned `NA`, which poisoned `if (!any(keep))`. Such
  measures now drop like zero-variance ones.

* `centrality_stability()`'s default `measures` is restored to
  `c("InStrength", "OutStrength", "Betweenness")`. 0.7.7 had swapped
  `OutStrength` for `Diffusion`, which broke the \pkg{htna} package: it calls
  `centrality_stability()` with no `measures` and compares the result
  against its own explicit trio. `centrality()` / `net_centrality()` keep
  the `Diffusion` default; only `centrality_stability()` reverts.

* `Suggests: cograph` relaxed from `(>= 2.4.4)` to `(>= 2.3.6)`, the version
  available on CRAN. `Additional_repositories` removed — every declared
  dependency now resolves from CRAN.

# Nestimate 0.7.7

## Enhancements

* `plot()` on a `net_centrality_group` gains `type = "delta"`, showing the
  between-group difference per measure, and now supports three or more
  groups. Zero-valued edges can be blanked with `drop_zero = TRUE`.

# Nestimate 0.7.6

## New features

* `plot()` on a `net_edge_betweenness()` result (`plot.net_edge_betweenness`).

## Enhancements

* `plot()` on centrality results gains alternative views:
  `type = c("bar", "line", "heatmap")` for a single `net_centrality`, and
  `type = c("bar", "line", "delta")` for a `net_centrality_group`. Count-like
  measures get integer axis labels.

# Nestimate 0.7.5

Version bump only; no user-visible changes.

# Nestimate 0.7.4

## New features

* `as_htna()` — builds a grouped node-level network from data and a
  clustering, keeping every node (unlike `cluster_summary()`, which
  collapses to a cluster-level macro summary). Intended for
  `cograph::plot_htna()`.

## Enhancements

* Centrality gains the `tna`-parity measures. `net_centrality(x, measures =
  "all")` now returns `OutStrength`, `InStrength`, `ClosenessIn`,
  `ClosenessOut`, `Closeness`, `Betweenness`, `BetweennessRSP`, `Diffusion`
  and `Clustering` — previously only the strengths, `Closeness` and
  `Betweenness`. Adds `plot.net_centrality` and `plot.net_centrality_group`.

* `sequence_plot()` and the MCML plots gain layout refinements.

# Nestimate 0.7.3

## New features

* `as_netobject()` / `validate_netobject()` — the boundary layer between
  \pkg{psychnet} (which owns the psychometric-network math and emits a lean
  `cograph_network`) and Nestimate (which owns the canonical `netobject`
  schema). `as_netobject()` promotes a `psychnet` result or a bare
  `cograph_network` to the dual-class `c("netobject", "cograph_network")` so it
  dispatches to every Nestimate verb, parking psychnet-specific fields (including
  the GLASSO KKT certificate) under `$meta$psychnet`; `netobject`s pass through
  unchanged. `validate_netobject()` enforces the shared structural contract so
  schema drift on either side fails loudly. `psychnet` is not a declared
  dependency — Nestimate never calls it; the converter works by S3 dispatch on
  whatever `psychnet` object the caller supplies.

* `certainty()` — analytic Bayesian counterpart of `bootstrap_network()` for
  transition networks. Models each state's outgoing transitions as a
  Dirichlet-Multinomial process (Jeffreys prior) and returns posterior mean,
  sd, credible interval and a stability decision per edge in closed form (no
  resampling). Returns the exact `net_bootstrap` object layout and carries
  class `c("net_certainty", "net_bootstrap")`, so it is a drop-in: every
  `net_bootstrap` method works on it. Completes the assessment trio
  certainty / stability (`bootstrap_network`) / reliability (`reliability`).

## Enhancements

* `sequence_plot()` gains a multichannel view for `mcml` objects built from
  sequences. `sequence_plot(fit)` draws one carpet panel per cluster channel
  plus a macro `Summary` panel — each channel's own states solid, the other
  clusters a faded wash, finished cells white, rows aligned by the macro
  sequence. `sequence_plot(fit, type = "distribution")` stacks the prevalence
  (own states + faded other clusters + an explicit `NA` band, to 100%), and
  `normalize = TRUE` gives a TraMineR-style `seqdplot` where each time point
  sums to 1. ggplot-based and dependency-free; returns a `ggplot` object.

* `bayes_compare()` results are now 100% compatible with the `permutation()`
  format: the object carries class `c("net_bayes", "net_permutation")` with all
  `net_permutation` slots (`diff_sig`, `p_values`, `effect_size`, `iter`,
  `alpha`, `paired`, `adjust`), and its `summary` is a superset of
  `summary.net_permutation` (`from, to, weight_x, weight_y, diff, effect_size,
  p_value, sig` plus the Bayesian extras `count_x, count_y, ci_lower, ci_upper,
  ci_width, pd`). A `bayes_compare()` result is now a drop-in wherever a
  `net_permutation` is consumed.

# Nestimate 0.7.2

## New features

* `bayes_compare()` — Bayesian Dirichlet-Multinomial comparison of two
  transition networks, a complement to `permutation()`. Models each source
  state's outgoing transitions as a Dirichlet-Multinomial process (Jeffreys
  prior) and returns, per edge, a posterior mean difference, a credible
  interval, the probability of direction (`pd`) and its two-sided
  p-equivalent. Adds `print`/`summary`/`plot` methods and `netobject_group`
  dispatch (all-pairwise or matched). Method source: Johnston & Jendoubi
  (2026), *How Delivery Mode Reshapes Resource Engagement: A Bayesian
  Differential Network Analysis*, TNA Workshop 2026.

# Nestimate 0.7.1

## New features

* `as_networks()` — promote a `build_mcml_pc()` result into a
  `netobject_group` (the psychometric-network counterpart of `as_tna()`).
  Singleton clusters with no within-network are dropped with a warning;
  an existing `netobject_group` passes through unchanged.

## Documentation

* Vignettes and articles now call package verbs directly instead of
  hand-assembled base-R subsetting rituals: `markov_order_test()` reads
  sequences straight from a fitted network (`markov_order_test(net)`);
  HYPA anomaly tables use `summary(hypa, order_by = "ratio")`; higher-order
  pathways use `pathways(hon, top = )`; grouped-clustering inspection uses
  `cluster_diagnostics()`.

# Nestimate 0.7.0

## New features (experimental)

* `build_mcml_pc()` — MCML aggregation for psychometric networks
  (cor / pcor / EBICglasso). Five aggregation methods with explicitly
  different statuses: `"average"` (descriptive block-mean; works without
  raw data), `"composite"` (cluster scores re-estimated with the chosen
  estimator — a genuine cluster-level network), `"loadings"` (composites
  weighted by mean within-cluster connection strength — Nestimate's own
  weighting, not an EGA reimplementation), `"rv"` (Escoufier's RV matrix
  correlation between blocks), and `"canonical"` (first canonical
  correlation — the upper bound for composite methods). Within-cluster
  networks re-estimated by default (`within = "reestimate"`) since a
  pcor submatrix is not the subsystem's pcor network. Item diagnostics
  in `$loadings`: signed loadings (reverse-keyed items detected via the
  leading eigenvector of the within-block matrix and flipped in
  composites), cross-cluster strengths, and a `misfit` flag when an item
  is more connected to another cluster than its own (warned). Composites
  tolerate missing data (row-wise renormalized weighted means);
  Composite item weights are selectable via `weighting` — ten built-in
  schemes spanning three views of the cluster: the scale as scored
  (`"equal"`, `"item_total"`), the network's view (`"strength"`,
  `"eigen"`, `"closeness"`, `"betweenness"`, `"expected_influence"`,
  `"specificity"` — the misfit margin as a weighting, zeroing items that
  belong as much to another cluster), and the latent-variable view
  (`"pca"`, `"factor"`); plus fully custom weighting via a named numeric
  vector or a `function(W_block, data_block, nodes)`.
  `aggregation = "loadings"` is the alias for composite + strength.
  The `"factor"` weighting exposes its extraction method via
  `fa_method`: `"ml"` (factanal), `"paf"` (iterated principal axis),
  `"minres"` (ULS), or `"cfa"` (one-factor lavaan model; with
  `cor_method = "polychoric"` the categorical DWLS factor model) — all
  operating on the `cor_method`-consistent correlation structure.
  Reverse-keyed handling works under every sign-carrying scheme
  (item-total correlations are computed on eigen-sign-pre-oriented
  columns so a reversed member cannot contaminate small clusters).
  `cor_method = "polychoric"` (via lavaan) supports ordinal items;
  `id_col` drops identifier columns so
  `convert_sequence_format(format = "frequency")` actor-profiles feed
  the function directly (the within-person co-occurrence view of event
  data).
  Returns class `mcml_pc` (macro + within netobjects, all undirected)
  with print/summary/plot; the composite/loadings macro is a full
  netobject, so `bootstrap_network()`, `vertex_bootstrap()`, and
  `vertex_compare()` apply to it directly. `cograph::plot_mcml()`
  (>= 2.3.8) renders the two-layer undirected MCML view. Experimental:
  API and formulas may change.
* `loading_stability()` — case-bootstrap stability of the
  `build_mcml_pc()` composite weights (percentile CIs, sign-flip rates),
  with print and forest-style plot.

# Nestimate 0.6.5

## New features

* `vertex_bootstrap()` — Snijders & Borgatti (1999) vertex bootstrap for
  network-level statistics (density, mean weight, strength centralization,
  weighted reciprocity, plus custom `statistic_fn`). Needs only the weight
  matrix, so it works on data-less netobjects (`build_mlvar()`
  constituents, `as_tna(mcml)` elements, plain matrices) where
  `bootstrap_network()` cannot run. Returns a tidy one-row-per-statistic
  `net_vertex_bootstrap` with print/summary/plot. Self-loops are preserved
  (diagonal carries the resampled vertex's own self-weight); undirected
  replicates stay symmetric.
* `vertex_compare()` — the Snijders & Borgatti two-network test the
  vertex bootstrap was originally proposed for: z-tests and
  normal-approximation CIs for differences in network-level statistics
  between two networks (netobjects, matrices, or precomputed
  `net_vertex_bootstrap` objects). Tidy `net_vertex_comparison` result
  with print/summary/plot (forest plot of differences).
* `bootstrap_network()` and `vertex_bootstrap()` gain
  `ci_method = c("percentile", "basic")`: basic intervals
  (Davison & Hinkley 1997, eq. 5.6) reflect the percentile bounds around
  the observed estimate, correcting first-order bootstrap bias. Default
  remains `"percentile"`.

# Nestimate 0.6.4

## Bug fixes

* `build_mcml()` (sequence and edge-list paths) now records the
  *effective* directedness in `$meta$directed`: `FALSE` when
  `type = "cooccurrence"`, whose weights are symmetrized, instead of
  echoing the `directed` argument unchanged. Renderers that auto-detect
  directedness (e.g., `cograph::plot_mcml()` with `directed = NULL`) now
  draw co-occurrence MCML objects as undirected networks automatically.

# Nestimate 0.6.3

## New features

* `mosaic_analysis(data, var1, var2)` — two-variable mosaic analysis on a
  `data.frame`: chi-square or Fisher test, Cramer's V (df-adjusted effect
  size) and a flat mosaic plot. Returns class `mosaic_analysis` with a tidy
  one-row-per-cell `$counts`, a one-row `$stats`, and print/summary/plot.
  Distinct from `mosaic_plot()`, which draws from a fitted network object.

## Enhancements

* `mosaic_plot()` gains `style = c("classic", "flat")`. The flat style uses
  variable-width columns, white gutters and in-tile or side labels, sharing
  the classic style's geometry and diverging palette; `values = TRUE` prints
  residuals inside the tiles.

# Nestimate 0.6.2

## Bug fixes

* Corrected stale rank-scaling assertions in the test suite. No user-visible
  change.

# Nestimate 0.6.1

## New features

* `build_network()` and the transition wrappers (`build_tna()`,
  `build_ftna()`, `build_atna()`, `build_cna()`) gain `start` and `end`
  boundary markers: `FALSE` (default), `TRUE` (labels `"Start"` / `"End"`)
  or a custom string. `start` prepends a source state to every sequence;
  `end` places a sink in the single cell after each sequence's last non-`NA`
  state (not absorbing — see `mark_terminal_state()` for that). Honoured by
  the `relative`, `frequency`, `co_occurrence` and `attention` estimators;
  other methods error.

* `build_mmm()` gains `covariate_effect`. `"em"` (default) folds covariates
  into the EM as covariate-dependent mixing, changing the fit; `"posthoc"`
  fits a plain mixture and uses covariates only for the after-fit
  multinomial logit, leaving the clustering bit-identical to a
  no-covariate fit.

# Nestimate 0.6.0

## New features

* `magnitude_difference()` compares the frequency (FTNA) and probability
  (TNA) views of a transition network and quantifies the per-edge
  discrepancy on a common scale, with five metrics, four scalings, and two
  polar `plot()` portraits (stacked and circular).
* Full persistent homology with a Vietoris-Rips filtration
  (`persistent_homology()`, `build_simplicial(type = "vr")`) plus diagram
  tools `bottleneck_distance()` and `persistence_landscape()`.
* Network comparison: `compare_model()` (with `netobject_group` dispatch),
  `summary.netobject()`, `plot.net_comparison()`, and `rename_models()` for
  relabelling grouped network objects.

## Documentation

* The pkgdown reference index now lists every exported function;
  `magnitude_difference()`, `casedrop_reliability()`, `build_hypergraph()`,
  `hypergraph_measures()`, and `cluster_data()` were previously absent.

## Packaging

* Removed the `Remotes:` field; `cograph` and `tna` are available from CRAN,
  so no non-CRAN source pin is needed.

# Nestimate 0.5.1

## Audit follow-up (clustering + MCML)
Followed `codex_docs/audit_clustering` and `codex_docs/audit_mcml`
recommendations across two modules. Eleven of thirteen findings
addressed; two deferred pending design decisions on numeric semantics
(`directed = FALSE` raw-data MCML, MMM first-non-NA initial state).

### Bug fixes
* `cluster_network()` now forwards distance-clustering arguments
  (`na_syms`, `weighted`, `lambda`, `seed`, `q`, `p`, `covariates`)
  to `build_clusters()` instead of silently passing them to
  `build_network()`. The split runs on caller `...` only — netobject
  `build_args` continue to flow only to the `build_network()` step,
  protecting attention-method (`atna`) network history from being
  re-routed to weighted Hamming. (audit_clustering item 1)
* `.auto_detect_clusters()` (used by `build_mcml()` and
  `cluster_summary()`) now requires `node_groups` to carry a node
  identifier column when shaped as a data.frame, or be a named atomic
  vector keyed by node label. Previously, a bare `cluster`-only
  data.frame was read positionally — silently mis-assigning nodes
  whenever `node_groups` rows were in a different order than `x$nodes`.
  (audit_mcml item 1)
* `build_clusters()` now rejects all-missing input early with a clear
  message instead of failing indirectly downstream in pam/hclust.
  (audit_clustering item 4)

### New parameters
* `compare_mmm(return_fits = FALSE)` — when `TRUE`, the fitted
  `net_mmm` models are attached as `attr(result, "fits")` keyed by
  `k`, so users can pick the chosen model without re-running EM.
  Default behaviour unchanged. (audit_clustering item 6)

### Improvements
* `build_clusters()` validation messages now name the offending
  argument (`"'k' must be at least 2 (got k = 1)"`) rather than
  dumping the failing predicate. Top-level type checks switched to
  named-condition `stopifnot()` for the same reason.
  (audit_clustering item 2)

### Documentation
* `summary.mcml()` roxygen corrected — was claiming a printing side
  effect that doesn't exist. (audit_mcml item 5)
* `build_mcml()` `clusters = "<col>"` mode now documents its narrow
  contract: assigns each row's group label to both endpoints, so it
  only makes sense for within-group edge lists. (audit_mcml item 2)
* `build_mcml()` `method` parameter doc now steers raw sequence /
  event-log inputs to `"sum"`, since the function counts observed
  transitions. Other methods are for weighted edge lists or
  pre-existing matrices. (audit_mcml item 4)
* `as_tna.mcml()` "Excluded Clusters" section corrected — drop emits
  a `warning()` (was claimed silent) and only fires for `relative`
  method (was claimed unconditional). (audit_mcml item 6)
* `build_clusters()` `na_syms` doc adds an explicit "Missing-value
  distance rule" subsection: NA becomes a comparable sentinel state,
  not pairwise deletion. (audit_clustering item 3)
* `build_mmm()` adds an "Initial states" section explaining
  first-column-verbatim init and that build_mmm does NOT honor
  build_clusters-style `na_syms` — only actual `NA` cells become NA
  inits. (audit_clustering item 5, doc-only path)

### Tests
* +12 new tests pinning the corrected contracts and the documented
  edge cases (misordered `node_groups` alignment, label propagation
  through `state_distribution()`, `as_tna.mcml()` drop-warning
  fixture, MMM first-column NA behaviour, and the four-way
  `cluster_network()` arg-routing contract). Full sweep:
  1628 / 1628 pass, 0 fail.

# Nestimate 0.5.0

## Bug fixes
* `.extract_edges_from_matrix()` no longer drops the diagonal. Netobjects
  built via `.wrap_netobject()` (and therefore everything from
  `build_network()`, `build_mcml()`, `bootstrap_network()`, `build_mmm()`,
  `wtna()`, `as_tna()`) now have `$edges` containing every non-zero
  matrix entry, including self-loops. Previously `$weights` and `$edges`
  were silently inconsistent on any matrix with a non-zero diagonal,
  causing downstream consumers (e.g. `cograph::centrality()` on an
  MCML macro) to under-count node degree by 2.

## New features
* `plot_state_frequencies()` — native S3 generic for state-frequency
  plots across `netobject`, `netobject_group`, `mcml`, and `htna`.
  Defaults to a marimekko (mosaic) layout where column widths reflect
  per-group totals and segment heights reflect within-group state
  proportions; also supports a colored-bars style and a per-group
  faceted marimekko. Uses the package Okabe-Ito palette throughout.
* `plot_mosaic()` — exported low-level marimekko primitive built on
  `geom_rect()` with cumulative-width / cumulative-height geometry.
  Reusable for any tidy `data.frame(group, state, weight)` input.

# Nestimate 0.4.4

## Bug fixes
* `passage_time()` and `markov_stability()` now raise an explicit error
  naming the dead state when a transition-matrix row sums to zero,
  instead of silently propagating `NaN` through `eigen`/`solve`. Zero
  rows mean the chain is not ergodic; mean first passage times are
  undefined. Shared helper `.mpt_normalize_rows()` factored out of both
  entry points.
* `.prepare_association_input()` no longer hard-rejects non-square
  numeric matrices. For association methods (glasso, pcor, cor) the
  netobject's `$data` slot is a numeric matrix (not a data.frame). Any
  downstream caller that row-subsetted `$data` and re-invoked the
  estimator (`centrality_stability()`, `bootstrap_network()`,
  `reliability()`) was silently producing NULL centralities caught by
  `tryCatch`, which surfaced as an "all centrality measures have zero
  variance" warning or all-`NaN` correlations. The matrix branch now
  recognises non-square input as raw observation data and recursively
  re-enters through the data-frame branch. Square symmetric matrices
  (pre-computed correlation / covariance) still go through the
  symmetric-matrix path with the symmetry check intact.

## New parameters
* `build_network()` gains `state_cols` and `metadata_cols` parameters
  (both default `NULL`). Explicit overrides for the state-vs-metadata
  column classifier, which previously used a "values-in-nodes"
  heuristic that silently misclassifies metadata columns whose values
  coincide with node labels (e.g. a `condition` column with levels
  `"A","B","C"` when nodes are `"A","B","C"`). Validation: error on
  overlap between the two vectors, error on column names not present
  in the input data. Forwarded through the `group = ...` recursive
  dispatch so per-group calls honour the override.

## Removed
* `plot.net_link_prediction()` and `plot.mcml()` removed. Nestimate is
  a computation engine — visualization is the user's concern.
  Previously both methods called `cograph::` directly, violating the
  stated dependency invariant (Nestimate -> cograph direction
  forbidden). Users call `cograph::splot(net)` or
  `cograph::plot_mcml(fit)` directly.

## Documentation
* `wtna()` `@param type` now flags that `type = "relative"` combined
  with `method = "cooccurrence"` produces an asymmetric matrix
  (conditional co-occurrence given row state), not a symmetric
  undirected weight matrix. Use `type = "frequency"` if symmetric
  counts are required.

## Testing infrastructure
* New numerical-equivalence tests (gated by `NESTIMATE_EQUIV_TESTS=true`):
  `test-equiv-permutation.R` (vs. `stats::p.adjust` + hand-coded base-R
  permutation loop), `test-equiv-mlvar.R` (vs. `mlVAR::mlVAR` at
  machine precision), `test-equiv-association-rules.R` (vs.
  `arules::apriori`), `test-equiv-link-prediction.R` (vs. clean-room
  matrix algebra + `igraph::similarity`), `test-equiv-centrality-stability.R`
  (vs. `bootnet::corStability`). Total ~162k per-value comparisons;
  all within machine precision except centrality-stability which uses
  a documented drop-grid tolerance because bootnet uses `igraph`
  path-based centrality and Nestimate uses Floyd-Warshall.
* New HON-family equivalence tests under `local_testing_and_equivalence/`
  validating HON, HONEM, HYPA, MOGen, and hypergraph against pathpy
  2.2.0 (via reticulate), `BiasedUrn`, `RSpectra`, and `HyperG`. Not
  shipped in the R-package `tests/` directory; added to `.Rbuildignore`.
* Branch-matrix coverage added for the four many-mode APIs (`wtna`,
  `bootstrap_network`, `build_clusters`, `sequence_plot`) — systematic
  cross-product tests over all combinations of mode parameters to
  catch regressions where one branch silently diverges.

# Nestimate 0.4.3

## CRAN resubmission (addresses incoming-check NOTEs on 0.4.2)
* Tarball now ships `build/vignette.rds` (the vignette index). Previous 0.4.2 build used `R CMD build --no-build-vignettes`, which preserved pre-built `inst/doc/*.html` but stripped the index — CRAN flagged "VignetteBuilder field but no prebuilt vignette index."
* `test-gimme.R` now `skip_on_cran()`. GIMME tests fit a lavaan SEM per subject and took ~50s locally (2-3× on Windows), pushing total check time to 11 min on win-devel. Full test suite still runs in CI and local dev.

# Nestimate 0.4.2

## CRAN resubmission
* Full `--as-cran --run-donttest` audit pass.
* Purged stale `.Rcheck/` and `Meta/` build artifacts from working tree; added explicit `^Nestimate\.Rcheck$` and `^\.\.Rcheck$` entries to `.Rbuildignore` as belt-and-suspenders against repeat-submission contamination.

# Nestimate 0.4.1

## CRAN resubmission
* Pre-built vignettes included in `inst/doc/` as required by CRAN.
* Fixed 301-redirect URLs in README.
* Added `skip_on_cran()` to slow test block to keep check time under 10 minutes.

# Nestimate 0.4.0

## New functions
* `build_mlvar()` — multilevel VAR networks from ESM/EMA panel data. Estimates temporal (directed), contemporaneous (undirected), and between-subjects (undirected) networks matching `mlVAR::mlVAR()` at machine precision.
* `build_mmm()` / `compare_mmm()` — mixture of Markov models via EM, with BIC/AIC/ICL model selection and optional covariate regression.
* `cooccurrence()` — standalone co-occurrence network builder supporting 6 input formats and 8 similarity methods.
* `sequence_compare()` — k-gram pattern comparison across groups with optional permutation testing.
* `sequence_plot()` / `distribution_plot()` — base-R sequence index and state distribution plots with clustering integration.
* `build_simplicial()`, `persistent_homology()`, `q_analysis()` — topological analysis of networks via simplicial complexes.
* `nct()` — Network Comparison Test matching `NetworkComparisonTest::NCT()` at machine precision.
* `build_gimme()` — group iterative mean estimation for idiographic networks via lavaan.
* `passage_time()`, `markov_stability()` — Markov chain passage times and stability analysis.
* `predict_links()` / `evaluate_links()` — link prediction with 6 structural similarity methods.
* `association_rules()` — Apriori association rule mining from sequences or binary matrices.
* `predictability()` — node predictability for glasso/pcor/cor networks.
* `build_hon()`, `build_honem()`, `build_hypa()`, `build_mogen()` — higher-order network methods (HON, HONEM, HYPA, MOGen) now `cograph_network`-compatible.

## New datasets
* `human_long`, `ai_long` — canonical long-format human–AI pair programming interaction sequences (10,796 turns, 429 sessions).
* `chatgpt_srl` — ChatGPT-generated SRL scale scores for psychological network analysis.
* `trajectories` — 138-student engagement trajectory matrix (15 timepoints, 3 states).

## API
* `build_clusters()`, `network_reliability()`, `permutation()`, and `prepare()` replace earlier internal names for consistency with the `build_*` naming convention.
* `mgm` estimator added (`method = "mgm"`) for mixed continuous + categorical data via nodewise lasso, matching `mgm::mgm()` at machine precision.

## Bug fixes
* `build_mmm()` no longer crashes on platforms where `parallel::detectCores()` returns `NA` (macOS ARM64 CRAN check failure).
* `gimme` convergence filter now correctly handles all typed `NA` variants (`NA_character_`, `NA_real_`, etc.).
* `NaN` values in numeric metadata aggregation (all-`NA` sessions) normalized to `NA_real_`.
* HYPA p-values corrected; `hypa_score` column renamed to `p_value`.

## CRAN compliance
* `.data` pronoun added to `globalVariables()`.
* `base::.rowSums()` / `base::.colSums()` replaced with `rowSums()` / `colSums()`.
* `dev.new()` guarded by `interactive()` — no side effects under knitr or CI.
* Equivalence test files excluded from the built tarball.

## Performance
* `do.call(rbind, ...)` replaced with `data.table::rbindlist()` in `mcml.R` and `sequence_compare.R`.

# Nestimate 0.3.4

* HYPA: Renamed `hypa_score` column to `p_value` for clarity. Added `$over`, `$under`, `$n_over`, `$n_under` fields to `net_hypa` objects. Scores are now pre-sorted with anomalous paths first.
* HYPA: `summary.net_hypa()` now shows over/under-represented paths separately with a configurable `n` parameter.
* `pathways.netobject()`: New S3 method to extract higher-order pathways directly from a netobject (builds HON or HYPA internally).
* `path_counts()`: Now handles NAs in trajectories by stripping them before k-gram counting.

# Nestimate 0.2.15

* Preparing for publication

# Nestimate 0.2.0

* Reduced hard dependencies from 6 to 4 Imports (ggplot2, glasso, data.table, cluster).
* Removed igraph from Imports — `centrality_stability()` and `boot_glasso()` now accept a `centrality_fn` parameter for external centrality computation.
* Removed tna from Imports — moved to Suggests (only used for input class detection).
* Implemented `graphical_var()` from scratch using coordinate descent lasso + graphical lasso with EBIC model selection, eliminating the graphicalVAR dependency.
* Dropped `ml_graphical_var()` — users should use `mlvar()` for multilevel VAR.
* Removed cograph plot wrappers — `plot.netobject()`, `plot.net_bootstrap()`, `plot.net_permutation()`, `plot.net_hon()`, `plot.net_hypa()` and `as_cograph()` removed. Users call cograph plotting functions directly on netobjects.
* Added `attention` estimator for decay-weighted transition networks.
* Increased test coverage from 84.5% to 96.1% (2780 tests).
* Passes R CMD check with 0 errors, 0 warnings, 0 notes.

# Nestimate 0.1.0

* Initial release. Split from Saqrlab v0.3.0.
* Core estimation via `build_network()` with 8 built-in estimators.
* Bootstrap inference (`bootstrap_network()`), permutation testing (`permutation()`), EBICglasso bootstrap (`boot_glasso()`).
* Higher-order networks: HON, HONEM, HYPA, MOGen.
* GIMME, MCML, multilevel VAR, graphical VAR.
* Temporal network analysis and velocity TNA.
* Dual-class `c("netobject", "cograph_network")` output for cograph compatibility.
