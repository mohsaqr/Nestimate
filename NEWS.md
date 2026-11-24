# Nestimate 1.1.0

Changes since the CRAN release 0.9.24.

## Breaking: results with a new shape

These results are now tidy data frames, one row per observation, with a
`group` column where there are groups. `format = "wide"` returns the old
one-column-per-group table; `summary()` returns every row as a plain data
frame, and `print()` shows the first `n` rows.

* `sequence_metrics(group = )`: one row per index and group, with the
  difference, Hedges' g and its interval, eta squared and p-values.
* `sequence_patterns(group = )`: one row per pattern and group, with
  `group_support`.
* `net_centrality()` of a group of networks: one data frame with `group`
  and `state` columns; it returned a list of data frames.
* `summary()` of a grouped network: one row per metric and group.
* `summary()` of `passage_time()`, `transition_entropy()` and
  `chain_structure()` returns one data frame; `table =` gives the second
  table (`"pairs"`, `"chain"`, `"absorption"`). They returned a list.
  Grouped results stack the table with a `group` column.
* `summary()` of `markov_stability()` returns the table and prints nothing.
* `vertex_compare()` names its columns `observed_x` and `observed_y`.
* `summary(nct(), table = "edges")` returns the edge table; `nct()` stores
  the M and S tests as `$structure` and `$strength`.
* `print()` of a network lists its strongest edges instead of the weight
  matrix.

## Rerun an analysis made with 0.9.24 if it uses

* data whose column roles came from their names, or long-log times not in a
  whole format with a four-digit year. Declare `action`, `actor`, `time`,
  `state_cols`, `metadata_cols` or `cols`; one-hot data needs
  `format = "onehot"` or `codes`.
* wide data with metadata, `actor`, `time` or `order` columns, which were
  read as states.
* `permutation()` or `compare_networks(test = "permutation")` on networks
  built with different estimator settings (say Pearson and Spearman): the
  difference was reported as a difference between the groups. Such pairs
  are now refused.
* `boot_glasso()`, `build_mcml_pc()` or `loading_stability()`: they now use
  the network's own estimator settings.
* `casedrop_reliability()` or `centrality_stability()`: failed subsets count
  as below `threshold`, and `casedrop_reliability()` compares only the
  existing edges of a transition network (`target = "all"` gives the
  earlier values).
* `build_mogen()` or `markov_order_test()`: `dof = "paths"` is the default
  (`dof = "observed"` gives the earlier values).
* `build_hon(method = "hon")`, which now follows pyHON, or `build_hypa()`,
  whose null is strength-preserving (the earlier null flagged anomalies in
  first-order data).
* `nct()`, whose `$structure` and `$strength` held each other's values in
  0.9.24, or strengths from `summary()`, whose in and out spreads were
  swapped.
* `build_mlvar()` temporal edges, which now run from predictor to outcome.
* `sequence_compare()` effect sizes, which depended on `min_freq`.
* `predict_links()` ties, `evaluate_links()`, `scaling = "normalize"` on an
  undirected network (now directed when the result is asymmetric) or the
  kernel width of `landscape_surface()` and `landscape_trend()`.
* corrected estimates in centrality (`BetweennessRSP`, betweenness ties,
  min-max normalisation), `net_prune()`, `build_mmm()`, `build_clusters()`,
  `build_mcml()`, `passage_time()`, `chain_structure()`,
  `markov_stability()`, `entropy_bayes()`, `hypergraph_centrality()`,
  `q_analysis()`, `persistent_homology()`, `betti_numbers()`,
  `bottleneck_distance()`, `pathways()`, `mogen_transitions()`,
  `extract_initial_probs()` and `mosaic_analysis()`.

## New features

* `build_network(method = "latency")`: a reply network weighted by reply
  speed (Templeton et al. 2022).
* `tna` and `group_tna` models are accepted wherever a network or sequences
  are.
* `stability_landscape()`, `landscape_surface()` and `landscape_trend()`
  draw a transition network as a ball-in-cup landscape.
* `extract_sequences()`, `sequence_patterns()`, `sequence_metrics()` and
  `pattern_outcome()`: the sequences, patterns and indices behind any
  input. `sequence_metrics()` compares every index between groups, within
  actors when an actor appears in several groups.
* `mark_phases()` splits each actor's events into `n` phases, for
  `build_network(group = "phase")`.
* `mosaic_test()`, and `nested`, `actor` and `report` in `mosaic_plot()` and
  `mosaic_analysis()`, for counts nested in actors.
* `node_clusters()`, `label_predictions()` and `label_scores()`.
* `actor` drops whole actors in `centrality_stability()` and
  `casedrop_reliability()`, and takes nested levels in `permutation()`.
* `entropy_trajectory()` reads sequences as well as an event log, and its
  `plot()` draws a bootstrap confidence band around each trend.
* `sequence_plot()` gains `periods` (weeks or months), `period_unit` and
  `fill`; it and `distribution_plot()` gain `cols` and a two-part `trim`.
* `extract_edges()` and `extract_initial_probs()` read a group of networks;
  `summary()` of `magnitude_difference()` returns its edges.
* `compare_networks()` stars the tested differences in
  `plot(type = "difference")`.
* `nct()` gains `na_method` and accepts missing values and constant
  variables.
* `build_network(params = list(window_size = ))` counts sequences between
  blocks of positions, as `tna::tna()` does (it was ignored).
* The `build_*()` shortcuts take `group`, `scaling`, `threshold` and
  `params` as named arguments; `build_mcml()` gains `from`, `to`, `weight`
  and `cols`; `cooccurrence()` and `association_rules()` gain
  `format = "binary"`.
* New dataset `daily_actions`.

## Other changes

* `print()` of `sequence_compare()` shows one row per pattern and group,
  rounded, with a `group` column; `summary(format = "long")` returns that
  table for every pattern. The result and `summary()` keep their shape.
* One colour rule for states: with no colours set, every figure draws the
  states in Okabe-Ito colours dealt in alphabetical order, and objects
  derived from a network keep its palette. The clusters of an mcml and the
  groups of an `entropy_trajectory()` plot take the colours the states
  leave. `plot_state_frequencies(colors = )` is now
  `state_colors =` (`colors` still works, with a warning).
* Plots that used ggplot's default hues, viridis or green and red use
  Okabe-Ito and the package's blue-red diverging colours, with shapes or
  linetypes where colour alone told groups apart.
* Sequence verbs refuse one-hot and association networks; resampling verbs
  refuse one-sequence data and networks built with `level`, and all but
  `bootstrap_network()` refuse latency networks; Markov measures refuse
  correlation networks.
* `entropy_trajectory()` refuses a `time` column with missing values.
* `prepare_onehot()` and `plot_mosaic()` are no longer exported.
* `sequence_plot()` arguments after `group` must be named.
* `seed =` no longer moves the session's random stream.

## Bug fixes

* `plot()` of `entropy_trajectory()` dropped a group with few windows; it is
  now drawn as its windows and named in a warning.
* `plot()` of `bayes_compare()` on a group of networks failed.
* `build_network()` stored averaged metadata as a one-dimensional array, so
  `landscape_trend(along = )` failed.
* Colours set with `set_state_colors()` were lost by `build_clusters()`,
  `build_mmm()`, `cluster_network()` and other derived objects, and were
  matched against the wrong name column.
* `subset()` of `compare_networks()` tables printed at two decimals.
* One-state data failed in `bootstrap_network()`, `path_dependence()` and
  `predict_links()`.
* `labels =` did not relabel stored data.
* `sequence_plot()` and `distribution_plot()` refused a `data.table` or
  tibble.
* `cluster_network()` did not pass its arguments.
* `build_hon(method = "hon+")` missed some higher-order rules.
* The pyramid plot of `sequence_compare()` cut off its counts.

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
  there. Missing values are handled pairwise, now passed explicitly.

## Tests

* Two more reference-implementation tests (a hand-written linear solve for
  `passage_time()`, a manual `multinom` re-computation) and unused reference
  helpers moved to the local suite; section headers that claimed
  equivalence tests are renamed. New tests: mgm on data with missing values,
  and the `scale = FALSE` error.

# Nestimate 0.3.0

* First release on CRAN.
