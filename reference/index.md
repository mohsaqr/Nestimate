# Package index

## Network Estimation

Core functions for building networks from data

- [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  [`print(`*`<netobject>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  [`print(`*`<netobject_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  [`print(`*`<netobject_ml>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  [`summary(`*`<netobject>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  [`summary(`*`<netobject_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  [`print(`*`<summary.netobject>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  [`` `[`( ``*`<summary.netobject_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  [`print(`*`<summary.netobject_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  : Estimate a network
- [`build_tna()`](https://pak.dynasite.org/Nestimate/reference/build_tna.md)
  : Estimate a transition network analysis (TNA) model
- [`build_ftna()`](https://pak.dynasite.org/Nestimate/reference/build_ftna.md)
  : Estimate a frequency transition network analysis (FTNA) model
- [`build_atna()`](https://pak.dynasite.org/Nestimate/reference/build_atna.md)
  : Estimate an attention transition network analysis (ATNA) model
- [`build_cna()`](https://pak.dynasite.org/Nestimate/reference/build_cna.md)
  : Estimate a co-occurrence network analysis (CNA) model
- [`build_ngram()`](https://pak.dynasite.org/Nestimate/reference/build_ngram.md)
  : Estimate an n-gram transition network analysis model
- [`build_gap()`](https://pak.dynasite.org/Nestimate/reference/build_gap.md)
  : Estimate a gap transition network analysis model
- [`build_reverse()`](https://pak.dynasite.org/Nestimate/reference/build_reverse.md)
  : Estimate a reverse transition network analysis model
- [`build_latency()`](https://pak.dynasite.org/Nestimate/reference/build_latency.md)
  : Estimate a latency network analysis model
- [`build_wtna()`](https://pak.dynasite.org/Nestimate/reference/build_wtna.md)
  : Estimate a window transition network analysis (WTNA) model
- [`wtna()`](https://pak.dynasite.org/Nestimate/reference/wtna.md)
  [`print(`*`<wtna_mixed>`*`)`](https://pak.dynasite.org/Nestimate/reference/wtna.md)
  : Build a window-based transition network from one-hot data
- [`build_cor()`](https://pak.dynasite.org/Nestimate/reference/build_cor.md)
  : Estimate a correlation network
- [`build_pcor()`](https://pak.dynasite.org/Nestimate/reference/build_pcor.md)
  : Estimate a partial correlation network
- [`build_glasso()`](https://pak.dynasite.org/Nestimate/reference/build_glasso.md)
  : Estimate an EBIC graphical lasso network
- [`build_ising()`](https://pak.dynasite.org/Nestimate/reference/build_ising.md)
  : Estimate an Ising model
- [`build_mgm()`](https://pak.dynasite.org/Nestimate/reference/build_mgm.md)
  : Estimate a mixed graphical model (MGM)
- [`estimate_network()`](https://pak.dynasite.org/Nestimate/reference/estimate_network.md)
  : Estimate a network (deprecated)
- [`register_estimator()`](https://pak.dynasite.org/Nestimate/reference/register_estimator.md)
  : Register a network estimator
- [`get_estimator()`](https://pak.dynasite.org/Nestimate/reference/get_estimator.md)
  : Retrieve a registered estimator
- [`list_estimators()`](https://pak.dynasite.org/Nestimate/reference/list_estimators.md)
  : List the registered estimators
- [`remove_estimator()`](https://pak.dynasite.org/Nestimate/reference/remove_estimator.md)
  : Remove a registered estimator
- [`cooccurrence()`](https://pak.dynasite.org/Nestimate/reference/cooccurrence.md)
  : Build a co-occurrence network of items in transactions
- [`build_mlvar()`](https://pak.dynasite.org/Nestimate/reference/build_mlvar.md)
  [`print(`*`<net_mlvar>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mlvar.md)
  [`summary(`*`<net_mlvar>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mlvar.md)
  [`plot(`*`<net_mlvar>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mlvar.md)
  : Estimate multilevel vector autoregression networks
- [`build_gimme()`](https://pak.dynasite.org/Nestimate/reference/build_gimme.md)
  [`coefs(`*`<net_gimme>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_gimme.md)
  : Estimate person-specific networks with GIMME

## Bayesian Inference

Dirichlet-Multinomial posterior inference for transition networks.
[`certainty()`](https://pak.dynasite.org/Nestimate/reference/certainty.md)
is the closed-form counterpart of
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md);
[`bayes_compare()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
is the complement of
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md).

- [`certainty()`](https://pak.dynasite.org/Nestimate/reference/certainty.md)
  [`print(`*`<net_certainty>`*`)`](https://pak.dynasite.org/Nestimate/reference/certainty.md)
  : Estimate the Bayesian certainty of transition network edges
- [`bayes_compare()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
  [`print(`*`<net_bayes>`*`)`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
  [`summary(`*`<net_bayes>`*`)`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
  [`plot(`*`<net_bayes>`*`)`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
  [`print(`*`<net_bayes_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
  [`plot(`*`<net_bayes_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
  [`summary(`*`<net_bayes_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
  : Compare two transition networks with a Bayesian Dirichlet model
- [`subtract_networks()`](https://pak.dynasite.org/Nestimate/reference/subtract_networks.md)
  [`print(`*`<netdifference>`*`)`](https://pak.dynasite.org/Nestimate/reference/subtract_networks.md)
  : Subtract one network from another
- [`as_netdifference()`](https://pak.dynasite.org/Nestimate/reference/as_netdifference.md)
  : Convert a network comparison to a difference network

## Higher-Order Networks

Methods for capturing higher-order dependencies

- [`build_hon()`](https://pak.dynasite.org/Nestimate/reference/build_hon.md)
  [`print(`*`<net_hon>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_hon.md)
  [`summary(`*`<net_hon>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_hon.md)
  : Build a higher-order network from sequences
- [`build_honem()`](https://pak.dynasite.org/Nestimate/reference/build_honem.md)
  [`print(`*`<net_honem>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_honem.md)
  [`summary(`*`<net_honem>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_honem.md)
  [`plot(`*`<net_honem>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_honem.md)
  : Embed the nodes of a higher-order network with HONEM
- [`build_hypa()`](https://pak.dynasite.org/Nestimate/reference/build_hypa.md)
  [`print(`*`<net_hypa>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_hypa.md)
  [`summary(`*`<net_hypa>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_hypa.md)
  : Detect over- and under-represented paths with HYPA
- [`build_mogen()`](https://pak.dynasite.org/Nestimate/reference/build_mogen.md)
  [`print(`*`<net_mogen>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mogen.md)
  [`summary(`*`<net_mogen>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mogen.md)
  [`plot(`*`<net_mogen>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mogen.md)
  : Fit a multi-order generative model to sequences
- [`pathways()`](https://pak.dynasite.org/Nestimate/reference/pathways.md)
  : Extract pathways from a higher-order network
- [`mogen_transitions()`](https://pak.dynasite.org/Nestimate/reference/mogen_transitions.md)
  : Extract the transitions of a MOGen model
- [`path_counts()`](https://pak.dynasite.org/Nestimate/reference/path_counts.md)
  : Count the paths of k states in trajectories
- [`bipartite_groups()`](https://pak.dynasite.org/Nestimate/reference/bipartite_groups.md)
  : Build a hypergraph from group memberships
- [`clique_expansion()`](https://pak.dynasite.org/Nestimate/reference/clique_expansion.md)
  : Project a hypergraph onto a pairwise network
- [`hypergraph_centrality()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_centrality.md)
  : Compute the eigenvector centralities of a hypergraph
- [`build_hypergraph()`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md)
  [`print(`*`<net_hypergraph>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md)
  [`summary(`*`<net_hypergraph>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md)
  : Build a hypergraph from a network's cliques
- [`hypergraph_measures()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_measures.md)
  [`print(`*`<hypergraph_measures>`*`)`](https://pak.dynasite.org/Nestimate/reference/hypergraph_measures.md)
  : Compute the structural measures of a hypergraph
- [`hypergraph_cluster()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_cluster.md)
  [`print(`*`<net_hypergraph_cluster>`*`)`](https://pak.dynasite.org/Nestimate/reference/hypergraph_cluster.md)
  [`summary(`*`<net_hypergraph_cluster>`*`)`](https://pak.dynasite.org/Nestimate/reference/hypergraph_cluster.md)
  [`node_clusters()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_cluster.md)
  [`plot(`*`<net_hypergraph_cluster>`*`)`](https://pak.dynasite.org/Nestimate/reference/hypergraph_cluster.md)
  : Cluster the nodes of a hypergraph spectrally
- [`hypergraph_transduction()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_transduction.md)
  [`print(`*`<net_hypergraph_transduction>`*`)`](https://pak.dynasite.org/Nestimate/reference/hypergraph_transduction.md)
  [`summary(`*`<net_hypergraph_transduction>`*`)`](https://pak.dynasite.org/Nestimate/reference/hypergraph_transduction.md)
  [`label_predictions()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_transduction.md)
  [`label_scores()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_transduction.md)
  [`plot(`*`<net_hypergraph_transduction>`*`)`](https://pak.dynasite.org/Nestimate/reference/hypergraph_transduction.md)
  : Classify unlabelled hypergraph nodes by label spreading
- [`hypergraph_laplacian()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_laplacian.md)
  : Compute the normalised Laplacian of a hypergraph

## Markov Analysis

Order, structure, entropy, and stability of Markov chains

- [`chain_structure()`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)
  [`print(`*`<chain_structure>`*`)`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)
  [`plot(`*`<chain_structure>`*`)`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)
  [`summary(`*`<chain_structure>`*`)`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)
  [`print(`*`<chain_structure_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)
  [`summary(`*`<chain_structure_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)
  [`print(`*`<summary_chain_structure>`*`)`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)
  : Classify the states and structure of a Markov chain
- [`markov_order_test()`](https://pak.dynasite.org/Nestimate/reference/markov_order_test.md)
  [`print(`*`<net_markov_order>`*`)`](https://pak.dynasite.org/Nestimate/reference/markov_order_test.md)
  [`print(`*`<net_markov_order_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/markov_order_test.md)
  [`summary(`*`<net_markov_order>`*`)`](https://pak.dynasite.org/Nestimate/reference/markov_order_test.md)
  [`plot(`*`<net_markov_order>`*`)`](https://pak.dynasite.org/Nestimate/reference/markov_order_test.md)
  : Test the Markov order of a sequential process
- [`markov_stability()`](https://pak.dynasite.org/Nestimate/reference/markov_stability.md)
  [`print(`*`<net_markov_stability>`*`)`](https://pak.dynasite.org/Nestimate/reference/markov_stability.md)
  [`print(`*`<net_markov_stability_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/markov_stability.md)
  [`summary(`*`<net_markov_stability>`*`)`](https://pak.dynasite.org/Nestimate/reference/markov_stability.md)
  [`summary(`*`<net_markov_stability_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/markov_stability.md)
  [`plot(`*`<net_markov_stability>`*`)`](https://pak.dynasite.org/Nestimate/reference/markov_stability.md)
  [`plot(`*`<net_markov_stability_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/markov_stability.md)
  : Describe the stability of each state of a Markov chain
- [`passage_time()`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)
  [`print(`*`<net_mpt>`*`)`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)
  [`print(`*`<net_mpt_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)
  [`summary(`*`<net_mpt>`*`)`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)
  [`summary(`*`<net_mpt_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)
  [`plot(`*`<net_mpt>`*`)`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)
  [`plot(`*`<net_mpt_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)
  : Compute the mean first passage times of a Markov chain
- [`path_dependence()`](https://pak.dynasite.org/Nestimate/reference/path_dependence.md)
  [`print(`*`<net_path_dependence>`*`)`](https://pak.dynasite.org/Nestimate/reference/path_dependence.md)
  [`summary(`*`<net_path_dependence>`*`)`](https://pak.dynasite.org/Nestimate/reference/path_dependence.md)
  [`print(`*`<summary.net_path_dependence>`*`)`](https://pak.dynasite.org/Nestimate/reference/path_dependence.md)
  [`plot(`*`<net_path_dependence>`*`)`](https://pak.dynasite.org/Nestimate/reference/path_dependence.md)
  : Measure the path dependence of each context at order k
- [`transition_entropy()`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md)
  [`print(`*`<net_transition_entropy>`*`)`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md)
  [`print(`*`<net_transition_entropy_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md)
  [`summary(`*`<net_transition_entropy_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md)
  [`plot(`*`<net_transition_entropy_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md)
  [`summary(`*`<net_transition_entropy>`*`)`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md)
  [`plot(`*`<net_transition_entropy>`*`)`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md)
  : Compute the transition entropy of a Markov chain
- [`entropy_network()`](https://pak.dynasite.org/Nestimate/reference/entropy_network.md)
  : Decompose the entropy rate into a network of transitions
- [`entropy_trajectory()`](https://pak.dynasite.org/Nestimate/reference/entropy_trajectory.md)
  [`print(`*`<net_entropy_trajectory>`*`)`](https://pak.dynasite.org/Nestimate/reference/entropy_trajectory.md)
  [`summary(`*`<net_entropy_trajectory>`*`)`](https://pak.dynasite.org/Nestimate/reference/entropy_trajectory.md)
  [`plot(`*`<net_entropy_trajectory>`*`)`](https://pak.dynasite.org/Nestimate/reference/entropy_trajectory.md)
  : Track transition entropy over time in sliding windows
- [`entropy_bayes()`](https://pak.dynasite.org/Nestimate/reference/entropy_bayes.md)
  [`print(`*`<net_entropy_bayes>`*`)`](https://pak.dynasite.org/Nestimate/reference/entropy_bayes.md)
  [`print(`*`<net_entropy_bayes_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/entropy_bayes.md)
  [`summary(`*`<net_entropy_bayes_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/entropy_bayes.md)
  [`plot(`*`<net_entropy_bayes_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/entropy_bayes.md)
  [`summary(`*`<net_entropy_bayes>`*`)`](https://pak.dynasite.org/Nestimate/reference/entropy_bayes.md)
  [`plot(`*`<net_entropy_bayes>`*`)`](https://pak.dynasite.org/Nestimate/reference/entropy_bayes.md)
  : Estimate transition entropy with Bayesian credible intervals
- [`stability_landscape()`](https://pak.dynasite.org/Nestimate/reference/stability_landscape.md)
  [`summary(`*`<net_landscape>`*`)`](https://pak.dynasite.org/Nestimate/reference/stability_landscape.md)
  [`print(`*`<net_landscape>`*`)`](https://pak.dynasite.org/Nestimate/reference/stability_landscape.md)
  [`plot(`*`<net_landscape>`*`)`](https://pak.dynasite.org/Nestimate/reference/stability_landscape.md)
  : Estimate the stability landscape of a transition network
- [`landscape_surface()`](https://pak.dynasite.org/Nestimate/reference/landscape_surface.md)
  [`summary(`*`<net_landscape_surface>`*`)`](https://pak.dynasite.org/Nestimate/reference/landscape_surface.md)
  [`print(`*`<net_landscape_surface>`*`)`](https://pak.dynasite.org/Nestimate/reference/landscape_surface.md)
  [`plot(`*`<net_landscape_surface>`*`)`](https://pak.dynasite.org/Nestimate/reference/landscape_surface.md)
  : Trace a stability landscape across a sequence measure
- [`landscape_trend()`](https://pak.dynasite.org/Nestimate/reference/landscape_trend.md)
  [`summary(`*`<net_landscape_trend>`*`)`](https://pak.dynasite.org/Nestimate/reference/landscape_trend.md)
  [`print(`*`<net_landscape_trend>`*`)`](https://pak.dynasite.org/Nestimate/reference/landscape_trend.md)
  [`plot(`*`<net_landscape_trend>`*`)`](https://pak.dynasite.org/Nestimate/reference/landscape_trend.md)
  : Trace a stability landscape along time or a covariate

## Network Pruning

Non-destructive edge pruning and restoration

- [`net_prune()`](https://pak.dynasite.org/Nestimate/reference/net_prune.md)
  : Prune weak or non-significant edges from a network
- [`net_deprune()`](https://pak.dynasite.org/Nestimate/reference/net_deprune.md)
  : Undo the pruning of a network
- [`net_reprune()`](https://pak.dynasite.org/Nestimate/reference/net_reprune.md)
  : Re-apply the pruning of a network
- [`net_pruning_details()`](https://pak.dynasite.org/Nestimate/reference/net_pruning_details.md)
  [`print(`*`<net_pruning_details>`*`)`](https://pak.dynasite.org/Nestimate/reference/net_pruning_details.md)
  : List the edges a pruning removed

## Bootstrap & Inference

Statistical inference for network estimation

- [`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  [`print(`*`<net_bootstrap>`*`)`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  [`summary(`*`<net_bootstrap>`*`)`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  [`print(`*`<net_bootstrap_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  [`summary(`*`<net_bootstrap_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  [`print(`*`<wtna_boot_mixed>`*`)`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  [`summary(`*`<wtna_boot_mixed>`*`)`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  : Bootstrap a network estimate
- [`vertex_bootstrap()`](https://pak.dynasite.org/Nestimate/reference/vertex_bootstrap.md)
  [`print(`*`<net_vertex_bootstrap>`*`)`](https://pak.dynasite.org/Nestimate/reference/vertex_bootstrap.md)
  [`summary(`*`<net_vertex_bootstrap>`*`)`](https://pak.dynasite.org/Nestimate/reference/vertex_bootstrap.md)
  [`plot(`*`<net_vertex_bootstrap>`*`)`](https://pak.dynasite.org/Nestimate/reference/vertex_bootstrap.md)
  : Bootstrap network-level statistics by resampling nodes
- [`vertex_compare()`](https://pak.dynasite.org/Nestimate/reference/vertex_compare.md)
  [`print(`*`<net_vertex_comparison>`*`)`](https://pak.dynasite.org/Nestimate/reference/vertex_compare.md)
  [`summary(`*`<net_vertex_comparison>`*`)`](https://pak.dynasite.org/Nestimate/reference/vertex_compare.md)
  [`plot(`*`<net_vertex_comparison>`*`)`](https://pak.dynasite.org/Nestimate/reference/vertex_compare.md)
  : Test whether network-level statistics differ between two networks
- [`boot_glasso()`](https://pak.dynasite.org/Nestimate/reference/boot_glasso.md)
  [`print(`*`<boot_glasso>`*`)`](https://pak.dynasite.org/Nestimate/reference/boot_glasso.md)
  [`summary(`*`<boot_glasso>`*`)`](https://pak.dynasite.org/Nestimate/reference/boot_glasso.md)
  [`plot(`*`<boot_glasso>`*`)`](https://pak.dynasite.org/Nestimate/reference/boot_glasso.md)
  : Bootstrap a graphical lasso network
- [`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
  [`print(`*`<net_permutation>`*`)`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
  [`summary(`*`<net_permutation>`*`)`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
  [`print(`*`<net_permutation_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
  [`summary(`*`<net_permutation_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
  [`print(`*`<wtna_perm_mixed>`*`)`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
  [`summary(`*`<wtna_perm_mixed>`*`)`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
  : Test whether two networks differ by permutation
- [`permutation_diagnostics()`](https://pak.dynasite.org/Nestimate/reference/permutation_diagnostics.md)
  : Diagnose the effect of nesting on a permutation test
- [`nct()`](https://pak.dynasite.org/Nestimate/reference/nct.md)
  [`print(`*`<net_nct>`*`)`](https://pak.dynasite.org/Nestimate/reference/nct.md)
  [`summary(`*`<net_nct>`*`)`](https://pak.dynasite.org/Nestimate/reference/nct.md)
  [`as.data.frame(`*`<net_nct>`*`)`](https://pak.dynasite.org/Nestimate/reference/nct.md)
  : Compare two networks with the network comparison test
- [`compare_model()`](https://pak.dynasite.org/Nestimate/reference/compare_model.md)
  [`print(`*`<net_comparison>`*`)`](https://pak.dynasite.org/Nestimate/reference/compare_model.md)
  [`plot(`*`<net_comparison>`*`)`](https://pak.dynasite.org/Nestimate/reference/compare_model.md)
  : Compare two networks descriptively
- [`compare_networks()`](https://pak.dynasite.org/Nestimate/reference/compare_networks.md)
  [`print(`*`<net_network_comparison>`*`)`](https://pak.dynasite.org/Nestimate/reference/compare_networks.md)
  [`plot(`*`<net_network_comparison>`*`)`](https://pak.dynasite.org/Nestimate/reference/compare_networks.md)
  : Compare two or more networks
- [`summary(`*`<net_network_comparison>`*`)`](https://pak.dynasite.org/Nestimate/reference/comparison_tables.md)
  [`edge_differences()`](https://pak.dynasite.org/Nestimate/reference/comparison_tables.md)
  [`node_differences()`](https://pak.dynasite.org/Nestimate/reference/comparison_tables.md)
  [`global_differences()`](https://pak.dynasite.org/Nestimate/reference/comparison_tables.md)
  [`network_metrics()`](https://pak.dynasite.org/Nestimate/reference/comparison_tables.md)
  [`print(`*`<net_table>`*`)`](https://pak.dynasite.org/Nestimate/reference/comparison_tables.md)
  [`` `[`( ``*`<net_table>`*`)`](https://pak.dynasite.org/Nestimate/reference/comparison_tables.md)
  : Tabulate the results of a network comparison
- [`rename_models()`](https://pak.dynasite.org/Nestimate/reference/rename_models.md)
  : Rename the networks of a network group
- [`magnitude_difference()`](https://pak.dynasite.org/Nestimate/reference/magnitude_difference.md)
  [`print(`*`<magnitude_difference>`*`)`](https://pak.dynasite.org/Nestimate/reference/magnitude_difference.md)
  [`summary(`*`<magnitude_difference>`*`)`](https://pak.dynasite.org/Nestimate/reference/magnitude_difference.md)
  [`plot(`*`<magnitude_difference>`*`)`](https://pak.dynasite.org/Nestimate/reference/magnitude_difference.md)
  : Measure how far frequency and probability networks disagree

## Reliability & Stability

Assess reliability and stability of network estimates

- [`network_reliability()`](https://pak.dynasite.org/Nestimate/reference/network_reliability.md)
  [`print(`*`<net_reliability>`*`)`](https://pak.dynasite.org/Nestimate/reference/network_reliability.md)
  [`summary(`*`<net_reliability>`*`)`](https://pak.dynasite.org/Nestimate/reference/network_reliability.md)
  [`plot(`*`<net_reliability>`*`)`](https://pak.dynasite.org/Nestimate/reference/network_reliability.md)
  : Estimate the split-half reliability of a network
- [`casedrop_reliability()`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  [`print(`*`<net_casedrop_reliability>`*`)`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  [`summary(`*`<net_casedrop_reliability>`*`)`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  [`print(`*`<net_casedrop_reliability_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  [`summary(`*`<net_casedrop_reliability_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  [`print(`*`<summary.net_casedrop_reliability_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  [`plot(`*`<net_casedrop_reliability>`*`)`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  [`plot(`*`<net_casedrop_reliability_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  : Estimate the case-dropping stability of edge weights
- [`centrality_stability()`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md)
  [`print(`*`<net_stability>`*`)`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md)
  [`print(`*`<net_stability_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md)
  [`summary(`*`<net_stability_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md)
  [`summary(`*`<net_stability>`*`)`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md)
  [`plot(`*`<net_stability>`*`)`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md)
  : Estimate the CS-coefficient of centrality indices
- [`loading_stability()`](https://pak.dynasite.org/Nestimate/reference/loading_stability.md)
  [`print(`*`<pc_loading_stability>`*`)`](https://pak.dynasite.org/Nestimate/reference/loading_stability.md)
  [`summary(`*`<pc_loading_stability>`*`)`](https://pak.dynasite.org/Nestimate/reference/loading_stability.md)
  [`plot(`*`<pc_loading_stability>`*`)`](https://pak.dynasite.org/Nestimate/reference/loading_stability.md)
  : Bootstrap the item weights of a psychometric MCML fit

## Clustering & Grouping

Cluster-based and multilevel network analysis

- [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  [`print(`*`<net_clustering>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  [`summary(`*`<net_clustering>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  [`plot(`*`<net_clustering>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  [`print(`*`<tidy_covariates>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  : Cluster sequences by dissimilarity
- [`cluster_data()`](https://pak.dynasite.org/Nestimate/reference/cluster_data.md)
  : Cluster sequences (deprecated)
- [`cluster_choice()`](https://pak.dynasite.org/Nestimate/reference/cluster_choice.md)
  [`print(`*`<cluster_choice>`*`)`](https://pak.dynasite.org/Nestimate/reference/cluster_choice.md)
  [`summary(`*`<cluster_choice>`*`)`](https://pak.dynasite.org/Nestimate/reference/cluster_choice.md)
  [`plot(`*`<cluster_choice>`*`)`](https://pak.dynasite.org/Nestimate/reference/cluster_choice.md)
  : Compare clusterings across k, dissimilarities and methods
- [`cluster_diagnostics()`](https://pak.dynasite.org/Nestimate/reference/cluster_diagnostics.md)
  [`print(`*`<net_cluster_diagnostics>`*`)`](https://pak.dynasite.org/Nestimate/reference/cluster_diagnostics.md)
  [`plot(`*`<net_cluster_diagnostics>`*`)`](https://pak.dynasite.org/Nestimate/reference/cluster_diagnostics.md)
  [`summary(`*`<net_cluster_diagnostics>`*`)`](https://pak.dynasite.org/Nestimate/reference/cluster_diagnostics.md)
  : Collect the quality measures of a clustering
- [`cluster_summary()`](https://pak.dynasite.org/Nestimate/reference/cluster_summary.md)
  : Aggregate a network's weights to clusters of nodes
- [`cluster_mmm()`](https://pak.dynasite.org/Nestimate/reference/cluster_mmm.md)
  : Cluster sequences with a mixed Markov model
- [`cluster_network()`](https://pak.dynasite.org/Nestimate/reference/cluster_network.md)
  : Cluster sequences and estimate one network per cluster
- [`build_mcml()`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
  [`print(`*`<mcml_layer>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
  [`print(`*`<mcml>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
  [`summary(`*`<mcml>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
  : Build an MCML model of cluster and within-cluster transitions
- [`build_mcml_pc()`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
  [`print(`*`<mcml_pc>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
  [`summary(`*`<mcml_pc>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
  [`plot(`*`<mcml_pc>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
  : Build an MCML model of a psychometric network
- [`composites()`](https://pak.dynasite.org/Nestimate/reference/composites.md)
  : Compute cluster scores from a psychometric MCML fit
- [`item_loadings()`](https://pak.dynasite.org/Nestimate/reference/item_loadings.md)
  : Tabulate the item loadings of a psychometric MCML fit
- [`macro_network()`](https://pak.dynasite.org/Nestimate/reference/macro_network.md)
  : Build the cluster-level network with chosen clusters expanded
- [`build_mmm()`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
  [`print(`*`<net_mmm>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
  [`summary(`*`<net_mmm>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
  [`plot(`*`<net_mmm>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
  [`print(`*`<net_mmm_clustering>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
  [`plot(`*`<net_mmm_clustering>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
  : Fit a mixed Markov model
- [`compare_mmm()`](https://pak.dynasite.org/Nestimate/reference/compare_mmm.md)
  [`print(`*`<mmm_compare>`*`)`](https://pak.dynasite.org/Nestimate/reference/compare_mmm.md)
  [`summary(`*`<mmm_compare>`*`)`](https://pak.dynasite.org/Nestimate/reference/compare_mmm.md)
  [`plot(`*`<mmm_compare>`*`)`](https://pak.dynasite.org/Nestimate/reference/compare_mmm.md)
  : Compare mixed Markov models across numbers of subgroups
- [`session_ids()`](https://pak.dynasite.org/Nestimate/reference/session_ids.md)
  : Identify the session behind each sequence

## Simplicial Complex Analysis

Topological analysis of networks

- [`build_simplicial()`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md)
  [`print(`*`<simplicial_complex>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md)
  [`plot(`*`<simplicial_complex>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md)
  : Build a simplicial complex
- [`simplicial_features()`](https://pak.dynasite.org/Nestimate/reference/simplicial_features.md)
  : Tabulate the topological features of one or many networks
- [`persistent_homology()`](https://pak.dynasite.org/Nestimate/reference/persistent_homology.md)
  [`print(`*`<persistent_homology>`*`)`](https://pak.dynasite.org/Nestimate/reference/persistent_homology.md)
  [`plot(`*`<persistent_homology>`*`)`](https://pak.dynasite.org/Nestimate/reference/persistent_homology.md)
  : Compute the persistent homology of a network
- [`bottleneck_distance()`](https://pak.dynasite.org/Nestimate/reference/bottleneck_distance.md)
  : Compute the bottleneck distance between persistence diagrams
- [`persistence_landscape()`](https://pak.dynasite.org/Nestimate/reference/persistence_landscape.md)
  [`print(`*`<persistence_landscape>`*`)`](https://pak.dynasite.org/Nestimate/reference/persistence_landscape.md)
  [`plot(`*`<persistence_landscape>`*`)`](https://pak.dynasite.org/Nestimate/reference/persistence_landscape.md)
  : Compute the persistence landscape of a persistence diagram
- [`q_analysis()`](https://pak.dynasite.org/Nestimate/reference/q_analysis.md)
  [`print(`*`<q_analysis>`*`)`](https://pak.dynasite.org/Nestimate/reference/q_analysis.md)
  [`plot(`*`<q_analysis>`*`)`](https://pak.dynasite.org/Nestimate/reference/q_analysis.md)
  : Describe the q-connectivity of a simplicial complex
- [`betti_numbers()`](https://pak.dynasite.org/Nestimate/reference/betti_numbers.md)
  : Count the holes of a simplicial complex
- [`euler_characteristic()`](https://pak.dynasite.org/Nestimate/reference/euler_characteristic.md)
  : Compute the Euler characteristic of a simplicial complex
- [`simplicial_degree()`](https://pak.dynasite.org/Nestimate/reference/simplicial_degree.md)
  : Count the simplices that contain each node
- [`verify_simplicial()`](https://pak.dynasite.org/Nestimate/reference/verify_simplicial.md)
  [`print(`*`<net_simplicial_check>`*`)`](https://pak.dynasite.org/Nestimate/reference/verify_simplicial.md)
  : Verify a simplicial complex against igraph

## Data Preparation

Convert and prepare data for network estimation

- [`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md)
  [`print(`*`<nestimate_data>`*`)`](https://pak.dynasite.org/Nestimate/reference/prepare.md)
  : Prepare an event log for network estimation
- [`human_long`](https://pak.dynasite.org/Nestimate/reference/long-data.md)
  [`ai_long`](https://pak.dynasite.org/Nestimate/reference/long-data.md)
  : Coded turns of human-AI vibe coding sessions
- [`prepare_for_tna()`](https://pak.dynasite.org/Nestimate/reference/prepare_for_tna.md)
  : Prepare sequence data for the tna package
- [`action_to_onehot()`](https://pak.dynasite.org/Nestimate/reference/action_to_onehot.md)
  : Convert an action column to one-hot columns
- [`wide_to_long()`](https://pak.dynasite.org/Nestimate/reference/wide_to_long.md)
  : Convert wide sequences to long format
- [`long_to_wide()`](https://pak.dynasite.org/Nestimate/reference/long_to_wide.md)
  : Convert long event data to wide sequences
- [`convert_sequence_format()`](https://pak.dynasite.org/Nestimate/reference/convert_sequence_format.md)
  : Convert sequence data to counts, one-hot codes or transitions
- [`actor_endpoints()`](https://pak.dynasite.org/Nestimate/reference/actor_endpoints.md)
  : Tabulate the first and last states of each sequence
- [`mark_first_state()`](https://pak.dynasite.org/Nestimate/reference/mark_first_state.md)
  : Fill leading missing cells with a start state
- [`mark_terminal_state()`](https://pak.dynasite.org/Nestimate/reference/mark_terminal_state.md)
  : Fill trailing missing cells with an end state
- [`mark_phases()`](https://pak.dynasite.org/Nestimate/reference/mark_phases.md)
  : Split each actor's sequence into phases

## Utilities

Helper functions and extractors

- [`predictability()`](https://pak.dynasite.org/Nestimate/reference/predictability.md)
  : Compute the predictability of each node
- [`frequencies()`](https://pak.dynasite.org/Nestimate/reference/frequencies.md)
  [`summary(`*`<nest_transition_counts>`*`)`](https://pak.dynasite.org/Nestimate/reference/frequencies.md)
  : Count the transitions between states in sequences
- [`state_frequencies()`](https://pak.dynasite.org/Nestimate/reference/state_frequencies.md)
  : Count state frequencies in trajectories
- [`net_aggregate_weights()`](https://pak.dynasite.org/Nestimate/reference/net_aggregate_weights.md)
  : Aggregate edge weights into one value
- [`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
  [`print(`*`<net_centrality>`*`)`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
  [`summary(`*`<net_centrality>`*`)`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
  [`print(`*`<net_centrality_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
  [`summary(`*`<net_centrality_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
  [`` `[`( ``*`<net_centrality_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
  [`` `[`( ``*`<net_centrality>`*`)`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
  [`plot(`*`<net_centrality>`*`)`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
  [`plot(`*`<net_centrality_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
  : Compute centrality measures for a network
- [`net_edge_betweenness()`](https://pak.dynasite.org/Nestimate/reference/net_edge_betweenness.md)
  [`plot(`*`<net_edge_betweenness>`*`)`](https://pak.dynasite.org/Nestimate/reference/net_edge_betweenness.md)
  : Weight a network's edges by their betweenness
- [`coefs()`](https://pak.dynasite.org/Nestimate/reference/coefs.md) :
  Tabulate the temporal coefficients of an mlVAR fit
- [`as_tna()`](https://pak.dynasite.org/Nestimate/reference/as_tna.md) :
  Convert the layers of an MCML model to a group of networks
- [`as_htna()`](https://pak.dynasite.org/Nestimate/reference/as_htna.md)
  : Build a node-level network with clusters as node groups
- [`as_networks()`](https://pak.dynasite.org/Nestimate/reference/as_networks.md)
  : Convert a psychometric MCML fit to a group of networks
- [`as_netobject()`](https://pak.dynasite.org/Nestimate/reference/as_netobject.md)
  : Convert a network from another package to a netobject
- [`validate_netobject()`](https://pak.dynasite.org/Nestimate/reference/validate_netobject.md)
  : Validate the structure of a network object
- [`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md)
  : Extract the edge list of a network
- [`extract_initial_probs()`](https://pak.dynasite.org/Nestimate/reference/extract_initial_probs.md)
  [`print(`*`<nest_initial_probs>`*`)`](https://pak.dynasite.org/Nestimate/reference/extract_initial_probs.md)
  [`summary(`*`<nest_initial_probs>`*`)`](https://pak.dynasite.org/Nestimate/reference/extract_initial_probs.md)
  : Extract the initial-state probabilities of a network
- [`extract_transition_matrix()`](https://pak.dynasite.org/Nestimate/reference/extract_transition_matrix.md)
  [`print(`*`<nest_transition_matrix>`*`)`](https://pak.dynasite.org/Nestimate/reference/extract_transition_matrix.md)
  [`summary(`*`<nest_transition_matrix>`*`)`](https://pak.dynasite.org/Nestimate/reference/extract_transition_matrix.md)
  : Extract the transition matrix of a network
- [`state_colors()`](https://pak.dynasite.org/Nestimate/reference/state_colors.md)
  : Get the state colours of an object
- [`set_state_colors()`](https://pak.dynasite.org/Nestimate/reference/set_state_colors.md)
  [`` `state_colors<-`() ``](https://pak.dynasite.org/Nestimate/reference/set_state_colors.md)
  : Set the state colours of a network object

## Sequence Analysis

Sequence visualization, pattern comparison, and association mining

- [`sequence_plot()`](https://pak.dynasite.org/Nestimate/reference/sequence_plot.md)
  [`print(`*`<mcml_sequence_plot>`*`)`](https://pak.dynasite.org/Nestimate/reference/sequence_plot.md)
  : Draw a sequence index plot, heatmap or distribution plot
- [`distribution_plot()`](https://pak.dynasite.org/Nestimate/reference/distribution_plot.md)
  : Draw a state distribution plot
- [`plot_state_frequencies()`](https://pak.dynasite.org/Nestimate/reference/plot_state_frequencies.md)
  [`print(`*`<nestimate_facet_plot>`*`)`](https://pak.dynasite.org/Nestimate/reference/plot_state_frequencies.md)
  [`print(`*`<nestimate_facet_list>`*`)`](https://pak.dynasite.org/Nestimate/reference/plot_state_frequencies.md)
  [`print(`*`<state_freq>`*`)`](https://pak.dynasite.org/Nestimate/reference/plot_state_frequencies.md)
  [`plot(`*`<state_freq>`*`)`](https://pak.dynasite.org/Nestimate/reference/plot_state_frequencies.md)
  : Plot state frequencies as treemaps or bars
- [`state_distribution()`](https://pak.dynasite.org/Nestimate/reference/state_distribution.md)
  : Tabulate the state distribution of each group
- [`mosaic_plot()`](https://pak.dynasite.org/Nestimate/reference/mosaic_plot.md)
  : Draw a mosaic plot of transition or co-occurrence counts
- [`mosaic_analysis()`](https://pak.dynasite.org/Nestimate/reference/mosaic_analysis.md)
  [`plot(`*`<mosaic_analysis>`*`)`](https://pak.dynasite.org/Nestimate/reference/mosaic_analysis.md)
  [`print(`*`<mosaic_analysis>`*`)`](https://pak.dynasite.org/Nestimate/reference/mosaic_analysis.md)
  [`summary(`*`<mosaic_analysis>`*`)`](https://pak.dynasite.org/Nestimate/reference/mosaic_analysis.md)
  : Test and draw the association of two categorical variables
- [`mosaic_test()`](https://pak.dynasite.org/Nestimate/reference/mosaic_test.md)
  : Test the association of a mosaic with and without clustering
- [`extract_sequences()`](https://pak.dynasite.org/Nestimate/reference/extract_sequences.md)
  : Extract sequences in wide format
- [`sequence_patterns()`](https://pak.dynasite.org/Nestimate/reference/sequence_patterns.md)
  [`print(`*`<net_sequence_patterns>`*`)`](https://pak.dynasite.org/Nestimate/reference/sequence_patterns.md)
  [`summary(`*`<net_sequence_patterns>`*`)`](https://pak.dynasite.org/Nestimate/reference/sequence_patterns.md)
  [`` `[`( ``*`<net_sequence_patterns>`*`)`](https://pak.dynasite.org/Nestimate/reference/sequence_patterns.md)
  : Discover sequential patterns
- [`sequence_metrics()`](https://pak.dynasite.org/Nestimate/reference/sequence_metrics.md)
  [`print(`*`<net_sequence_metrics>`*`)`](https://pak.dynasite.org/Nestimate/reference/sequence_metrics.md)
  [`summary(`*`<net_sequence_metrics>`*`)`](https://pak.dynasite.org/Nestimate/reference/sequence_metrics.md)
  [`` `[`( ``*`<net_sequence_metrics>`*`)`](https://pak.dynasite.org/Nestimate/reference/sequence_metrics.md)
  : Compute descriptive indices for each sequence
- [`sequence_compare()`](https://pak.dynasite.org/Nestimate/reference/sequence_compare.md)
  [`print(`*`<net_sequence_comparison>`*`)`](https://pak.dynasite.org/Nestimate/reference/sequence_compare.md)
  [`summary(`*`<net_sequence_comparison>`*`)`](https://pak.dynasite.org/Nestimate/reference/sequence_compare.md)
  [`plot(`*`<net_sequence_comparison>`*`)`](https://pak.dynasite.org/Nestimate/reference/sequence_compare.md)
  : Compare subsequence patterns between groups
- [`extract_pathways()`](https://pak.dynasite.org/Nestimate/reference/extract_pathways.md)
  : Cut an event log into pathways
- [`association_rules()`](https://pak.dynasite.org/Nestimate/reference/association_rules.md)
  [`print(`*`<net_association_rules>`*`)`](https://pak.dynasite.org/Nestimate/reference/association_rules.md)
  [`summary(`*`<net_association_rules>`*`)`](https://pak.dynasite.org/Nestimate/reference/association_rules.md)
  [`plot(`*`<net_association_rules>`*`)`](https://pak.dynasite.org/Nestimate/reference/association_rules.md)
  : Discover association rules in transactions or sequences

## Link Prediction

Predict and evaluate missing connections

- [`predict_links()`](https://pak.dynasite.org/Nestimate/reference/predict_links.md)
  [`print(`*`<net_link_prediction>`*`)`](https://pak.dynasite.org/Nestimate/reference/predict_links.md)
  [`summary(`*`<net_link_prediction>`*`)`](https://pak.dynasite.org/Nestimate/reference/predict_links.md)
  : Predict missing or future links in a network
- [`evaluate_links()`](https://pak.dynasite.org/Nestimate/reference/evaluate_links.md)
  : Evaluate link predictions against known edges

## Outcome Modelling

Relate network and sequence features to an outcome

- [`pattern_outcome()`](https://pak.dynasite.org/Nestimate/reference/pattern_outcome.md)
  : Model a binary outcome from sequential patterns
- [`outcome_model()`](https://pak.dynasite.org/Nestimate/reference/outcome_model.md)
  [`print(`*`<net_outcome_model>`*`)`](https://pak.dynasite.org/Nestimate/reference/outcome_model.md)
  [`summary(`*`<net_outcome_model>`*`)`](https://pak.dynasite.org/Nestimate/reference/outcome_model.md)
  [`plot(`*`<net_outcome_model>`*`)`](https://pak.dynasite.org/Nestimate/reference/outcome_model.md)
  : Model a unit-level outcome from sequence or network predictors
- [`effects_table()`](https://pak.dynasite.org/Nestimate/reference/effects_table.md)
  : Tabulate the effects of an outcome model

## Data

Example datasets

- [`human_long`](https://pak.dynasite.org/Nestimate/reference/long-data.md)
  [`ai_long`](https://pak.dynasite.org/Nestimate/reference/long-data.md)
  : Coded turns of human-AI vibe coding sessions
- [`daily_actions`](https://pak.dynasite.org/Nestimate/reference/daily_actions.md)
  : Simulated learning actions over seven days
- [`srl_strategies`](https://pak.dynasite.org/Nestimate/reference/srl_strategies.md)
  : Simulated counts of self-regulated learning strategies
- [`learning_activities`](https://pak.dynasite.org/Nestimate/reference/learning_activities.md)
  : Online learning activity indicators
- [`group_regulation_long`](https://pak.dynasite.org/Nestimate/reference/group_regulation_long.md)
  : Group regulation actions in collaborative learning
- [`chatgpt_srl`](https://pak.dynasite.org/Nestimate/reference/chatgpt_srl.md)
  : ChatGPT scores on self-regulated learning scales
- [`trajectories`](https://pak.dynasite.org/Nestimate/reference/trajectories.md)
  : Student engagement trajectories
