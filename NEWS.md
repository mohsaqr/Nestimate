# Nestimate 1.2.0

This version collects everything since 0.9.24, the last release on CRAN.
Some results are refined and some functions now ask you to say what your
columns are; see "Refined results" and "What works differently".

## New

* **Latency networks.** `build_network(method = "latency")` builds a reply
  network in which faster replies carry more weight (Templeton et al. 2022).
* **Stability landscapes.** `stability_landscape()` draws each state as a
  cup whose depth shows how long people stay in it. `landscape_trend()`
  follows the depths over time or along a covariate, and
  `landscape_surface()` draws them along a measure of stability.
* **Sequence patterns and indices.** `sequence_patterns()` finds the
  patterns in sequences and tests them between groups. `sequence_metrics()`
  computes 23 indices for each sequence (entropy, turbulence, complexity,
  spells and others) and compares them between groups. `pattern_outcome()`
  relates patterns to an outcome, and `extract_sequences()` returns the
  sequences behind any input.
* **Phases.** `mark_phases()` splits each person's events into phases, so
  `build_network(group = "phase")` gives one network per phase.
* **Nested data.** `mosaic_test()` tests counts nested in people or teams,
  and `mosaic_plot()` and `mosaic_analysis()` take `nested` and `actor`.
  `permutation()`, `centrality_stability()` and `casedrop_reliability()`
  take `actor` to resample whole people.
* **Entropy over time.** `entropy_trajectory()` reads sequences as well as
  event logs, and its plot shows a bootstrap band around each trend.
* **Sequence plots by period.** `sequence_plot()` can show actions per day,
  week or month (`periods`), and both `sequence_plot()` and
  `distribution_plot()` can drop rare columns (`trim`).
* **tna models** can be passed wherever a network or sequences are
  accepted.
* **New dataset** `daily_actions`: the learning actions of 30 students over
  seven days (simulated).
* Smaller additions: `node_clusters()`, `label_predictions()` and
  `label_scores()`; `nct()` handles missing values and constant variables;
  `extract_edges()` and `extract_initial_probs()` read a group of networks;
  the `build_*()` shortcuts take `group`, `scaling`, `threshold` and
  `params`; `build_mcml()` reads edge lists; `cooccurrence()` and
  `association_rules()` read 0/1 tables (`format = "binary"`);
  `compare_networks()` marks significant differences in its difference
  plot; `window_size` in `build_network(params = )` now works as in tna.

## Refined results

These functions give refined results since 0.9.24.

* `nct()` reported the structure test as the strength test and the other
  way round, and `summary()` swapped the in- and out-strength spreads.
* `build_mlvar()` drew temporal edges in the wrong direction (from outcome
  to predictor).
* `sequence_compare()` effect sizes changed with `min_freq`.
* `build_hypa()` flagged anomalies in first-order data, which has none;
  it now uses a strength-preserving null. `build_hon(method = "hon")` now
  matches pyHON, and `method = "hon+"` missed some rules.
* `permutation()` and `compare_networks(test = "permutation")` accepted two
  networks estimated in different ways (say Pearson and Spearman) and
  reported the difference as a difference between the groups. They now
  refuse.
* `boot_glasso()`, `build_mcml_pc()` and `loading_stability()` ignored the
  network's own estimator settings.
* Numbers were corrected in centrality (`BetweennessRSP`, betweenness ties,
  min-max normalisation), `net_prune()`, `build_mmm()`, `build_clusters()`,
  `build_mcml()`, `passage_time()`, `chain_structure()`,
  `markov_stability()`, `entropy_bayes()`, `hypergraph_centrality()`,
  `q_analysis()`, `persistent_homology()`, `betti_numbers()`,
  `bottleneck_distance()`, `pathways()`, `mogen_transitions()`,
  `extract_initial_probs()`, `mosaic_analysis()`, `predict_links()` (ties),
  `evaluate_links()`, and the smoothing of `landscape_surface()` and
  `landscape_trend()`.

## What works differently

* **Say which column is which.** Nestimate no longer guesses a column's
  role from its name. Give `action`, `actor` and `time` for an event log,
  or `state_cols` for wide data, and `format = "onehot"` for 0/1 code
  columns. Extra columns in wide data (metadata, `actor`, `time`, `order`)
  are no longer read as states. Times must include a four-digit year.
* **Results are tables.** These used to return lists or wide tables and now
  return one data frame, one row per item, with a `group` column when
  groups are compared: `sequence_metrics(group = )`,
  `sequence_patterns(group = )`, `net_centrality()` on several networks,
  and `summary()` of a grouped network, `passage_time()`,
  `transition_entropy()`, `chain_structure()`, `markov_stability()` and
  `nct()`. `format = "wide"` gives the old layout; `table = ` picks the
  second table where there is one. Printing shows the first rows,
  `summary()` all of them.
* **Printing a network** lists its strongest edges instead of the whole
  weight matrix.
* **New defaults.** `build_mogen()` and `markov_order_test()` use
  `dof = "paths"` (`dof = "observed"` gives the old values).
  `casedrop_reliability()` compares only the edges a transition network
  has (`target = "all"` gives the old values), and it and
  `centrality_stability()` count a subset that fails as below the
  threshold.
* **Renamed.** `vertex_compare()` columns are `observed_x` and
  `observed_y`; `plot_state_frequencies(colors = )` is `state_colors =`
  (the old name still works, with a warning).
* **Refused inputs.** Functions that read sequences refuse one-hot
  networks and association networks (correlation, partial correlation,
  glasso, Ising, mixed graphical models); resampling functions refuse a
  single sequence and networks built with `level`; latency networks can
  only be bootstrapped; Markov measures refuse correlation networks;
  `entropy_trajectory()` refuses missing times.
* `scaling = "normalize"` on an undirected network returns a directed
  network when the rows no longer match the columns.
* Arguments of `sequence_plot()` after `group` must be named.
* `prepare_onehot()` and `plot_mosaic()` are no longer exported.
* `seed =` no longer changes your session's random numbers.

## Documentation

* "Comparing sequences between groups" is now an article on the
  package website instead of a vignette installed with the package;
  its address, PDF and DOI are unchanged.
* Examples use 20 resamples (200 posterior draws) so they run quickly;
  a real analysis uses 1,000 or more.

## Colours

Every figure colours states the same way: with no colours set, states get
colour-blind-safe Okabe-Ito colours in alphabetical order, and anything
derived from a network keeps its colours. Plots that used rainbow, viridis
or red-green now use these colours, with shapes or line types where colour
alone separated groups.

## Other fixes

* `plot()` of `entropy_trajectory()` dropped groups with few windows.
* `plot()` of `bayes_compare()` failed on a group of networks.
* `landscape_trend(along = )` failed on averaged metadata.
* Colours set with `set_state_colors()` were lost after clustering.
* `subset()` of `compare_networks()` tables printed only two decimals.
* Data with a single state failed in `bootstrap_network()`,
  `path_dependence()` and `predict_links()`.
* `labels =` did not relabel the stored data.
* `sequence_plot()` and `distribution_plot()` refused a data.table or
  tibble.
* `cluster_network()` ignored its arguments.
* The pyramid plot of `sequence_compare()` cut off its counts, and its
  printout now shows one row per pattern and group.

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
