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
  [`print(`*`<summary.netobject_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
  : Build a Network
- [`estimate_network()`](https://pak.dynasite.org/Nestimate/reference/estimate_network.md)
  : Estimate a Network (Deprecated)
- [`register_estimator()`](https://pak.dynasite.org/Nestimate/reference/register_estimator.md)
  : Register a Network Estimator
- [`get_estimator()`](https://pak.dynasite.org/Nestimate/reference/get_estimator.md)
  : Retrieve a Registered Estimator
- [`list_estimators()`](https://pak.dynasite.org/Nestimate/reference/list_estimators.md)
  : List All Registered Estimators
- [`remove_estimator()`](https://pak.dynasite.org/Nestimate/reference/remove_estimator.md)
  : Remove a Registered Estimator
- [`build_tna()`](https://pak.dynasite.org/Nestimate/reference/build_tna.md)
  : Build a Transition Network (TNA)
- [`build_atna()`](https://pak.dynasite.org/Nestimate/reference/build_atna.md)
  : Build an Attention-Weighted Transition Network (ATNA)
- [`build_ftna()`](https://pak.dynasite.org/Nestimate/reference/build_ftna.md)
  : Build a Frequency Transition Network (FTNA)
- [`build_cna()`](https://pak.dynasite.org/Nestimate/reference/build_cna.md)
  : Build a Co-occurrence Network (CNA)
- [`build_cor()`](https://pak.dynasite.org/Nestimate/reference/build_cor.md)
  : Build a Correlation Network
- [`build_pcor()`](https://pak.dynasite.org/Nestimate/reference/build_pcor.md)
  : Build a Partial Correlation Network
- [`build_glasso()`](https://pak.dynasite.org/Nestimate/reference/build_glasso.md)
  : Build a Graphical Lasso Network (EBICglasso)
- [`build_ising()`](https://pak.dynasite.org/Nestimate/reference/build_ising.md)
  : Build an Ising Network
- [`wtna()`](https://pak.dynasite.org/Nestimate/reference/wtna.md)
  [`print(`*`<wtna_mixed>`*`)`](https://pak.dynasite.org/Nestimate/reference/wtna.md)
  : Window-based Transition Network Analysis
- [`cooccurrence()`](https://pak.dynasite.org/Nestimate/reference/cooccurrence.md)
  : Build a Co-occurrence Network
- [`build_mlvar()`](https://pak.dynasite.org/Nestimate/reference/build_mlvar.md)
  [`print(`*`<net_mlvar>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mlvar.md)
  [`summary(`*`<net_mlvar>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mlvar.md)
  : Build a Multilevel Vector Autoregression (mlVAR) network
- [`build_gimme()`](https://pak.dynasite.org/Nestimate/reference/build_gimme.md)
  : GIMME: Group Iterative Multiple Model Estimation

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
  : Analytic certainty of network edges (Bayesian Dirichlet-Multinomial)
- [`bayes_compare()`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
  [`print(`*`<net_bayes>`*`)`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
  [`summary(`*`<net_bayes>`*`)`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
  [`plot(`*`<net_bayes>`*`)`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
  [`print(`*`<net_bayes_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
  [`summary(`*`<net_bayes_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/bayes_compare.md)
  : Bayesian Dirichlet-Multinomial comparison of two transition networks
- [`subtract_networks()`](https://pak.dynasite.org/Nestimate/reference/subtract_networks.md)
  [`print(`*`<netdifference>`*`)`](https://pak.dynasite.org/Nestimate/reference/subtract_networks.md)
  : Subtract one network from another
- [`as_netdifference()`](https://pak.dynasite.org/Nestimate/reference/as_netdifference.md)
  : Coerce an inferential comparison to a network difference

## Higher-Order Networks

Methods for capturing higher-order dependencies

- [`build_hon()`](https://pak.dynasite.org/Nestimate/reference/build_hon.md)
  [`print(`*`<net_hon>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_hon.md)
  [`summary(`*`<net_hon>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_hon.md)
  : Build a Higher-Order Network (HON)
- [`build_honem()`](https://pak.dynasite.org/Nestimate/reference/build_honem.md)
  [`print(`*`<net_honem>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_honem.md)
  [`summary(`*`<net_honem>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_honem.md)
  [`plot(`*`<net_honem>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_honem.md)
  : Build HONEM Embeddings for Higher-Order Networks
- [`build_hypa()`](https://pak.dynasite.org/Nestimate/reference/build_hypa.md)
  [`print(`*`<net_hypa>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_hypa.md)
  [`summary(`*`<net_hypa>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_hypa.md)
  : Detect Path Anomalies via HYPA
- [`build_mogen()`](https://pak.dynasite.org/Nestimate/reference/build_mogen.md)
  [`print(`*`<net_mogen>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mogen.md)
  [`summary(`*`<net_mogen>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mogen.md)
  [`plot(`*`<net_mogen>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mogen.md)
  : Build Multi-Order Generative Model (MOGen)
- [`pathways()`](https://pak.dynasite.org/Nestimate/reference/pathways.md)
  : Extract Pathways from Higher-Order Network Objects
- [`mogen_transitions()`](https://pak.dynasite.org/Nestimate/reference/mogen_transitions.md)
  : Extract Transition Table from a MOGen Model
- [`path_counts()`](https://pak.dynasite.org/Nestimate/reference/path_counts.md)
  : Count Path Frequencies in Trajectory Data
- [`bipartite_groups()`](https://pak.dynasite.org/Nestimate/reference/bipartite_groups.md)
  : Hypergraph from bipartite group / event data
- [`clique_expansion()`](https://pak.dynasite.org/Nestimate/reference/clique_expansion.md)
  : Clique expansion of a hypergraph
- [`hypergraph_centrality()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_centrality.md)
  : Hypergraph eigenvector centralities
- [`build_hypergraph()`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md)
  [`print(`*`<net_hypergraph>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md)
  [`summary(`*`<net_hypergraph>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_hypergraph.md)
  : Higher-order hypergraph from a network's clique structure
- [`hypergraph_measures()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_measures.md)
  [`print(`*`<hypergraph_measures>`*`)`](https://pak.dynasite.org/Nestimate/reference/hypergraph_measures.md)
  : Structural measures for a hypergraph
- [`hypergraph_cluster()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_cluster.md)
  [`print(`*`<net_hypergraph_cluster>`*`)`](https://pak.dynasite.org/Nestimate/reference/hypergraph_cluster.md)
  [`summary(`*`<net_hypergraph_cluster>`*`)`](https://pak.dynasite.org/Nestimate/reference/hypergraph_cluster.md)
  [`as.data.frame(`*`<net_hypergraph_cluster>`*`)`](https://pak.dynasite.org/Nestimate/reference/hypergraph_cluster.md)
  [`plot(`*`<net_hypergraph_cluster>`*`)`](https://pak.dynasite.org/Nestimate/reference/hypergraph_cluster.md)
  : Spectral clustering of hypergraph vertices
- [`hypergraph_transduction()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_transduction.md)
  [`print(`*`<net_hypergraph_transduction>`*`)`](https://pak.dynasite.org/Nestimate/reference/hypergraph_transduction.md)
  [`summary(`*`<net_hypergraph_transduction>`*`)`](https://pak.dynasite.org/Nestimate/reference/hypergraph_transduction.md)
  [`as.data.frame(`*`<net_hypergraph_transduction>`*`)`](https://pak.dynasite.org/Nestimate/reference/hypergraph_transduction.md)
  [`plot(`*`<net_hypergraph_transduction>`*`)`](https://pak.dynasite.org/Nestimate/reference/hypergraph_transduction.md)
  : Transductive label spreading on a hypergraph
- [`hypergraph_laplacian()`](https://pak.dynasite.org/Nestimate/reference/hypergraph_laplacian.md)
  : Normalized hypergraph Laplacian

## Markov Analysis

Order, structure, entropy, and stability of Markov chains

- [`chain_structure()`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)
  [`print(`*`<chain_structure>`*`)`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)
  [`plot(`*`<chain_structure>`*`)`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)
  [`summary(`*`<chain_structure>`*`)`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)
  [`print(`*`<chain_structure_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)
  [`summary(`*`<chain_structure_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)
  [`print(`*`<summary_chain_structure>`*`)`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md)
  : Qualitative structure of a discrete-time Markov chain
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
  [`plot(`*`<net_markov_stability>`*`)`](https://pak.dynasite.org/Nestimate/reference/markov_stability.md)
  : Markov Stability Analysis
- [`passage_time()`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)
  [`print(`*`<net_mpt>`*`)`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)
  [`print(`*`<net_mpt_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)
  [`summary(`*`<net_mpt>`*`)`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)
  [`print(`*`<summary.net_mpt>`*`)`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)
  [`plot(`*`<net_mpt>`*`)`](https://pak.dynasite.org/Nestimate/reference/passage_time.md)
  : Mean First Passage Times
- [`path_dependence()`](https://pak.dynasite.org/Nestimate/reference/path_dependence.md)
  [`print(`*`<net_path_dependence>`*`)`](https://pak.dynasite.org/Nestimate/reference/path_dependence.md)
  [`summary(`*`<net_path_dependence>`*`)`](https://pak.dynasite.org/Nestimate/reference/path_dependence.md)
  [`print(`*`<summary.net_path_dependence>`*`)`](https://pak.dynasite.org/Nestimate/reference/path_dependence.md)
  [`plot(`*`<net_path_dependence>`*`)`](https://pak.dynasite.org/Nestimate/reference/path_dependence.md)
  : Per-Context Path Dependence at Order k
- [`transition_entropy()`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md)
  [`print(`*`<net_transition_entropy>`*`)`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md)
  [`print(`*`<net_transition_entropy_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md)
  [`summary(`*`<net_transition_entropy>`*`)`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md)
  [`print(`*`<summary.net_transition_entropy>`*`)`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md)
  [`plot(`*`<net_transition_entropy>`*`)`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md)
  : Transition Entropy of a Markov Chain
- [`entropy_network()`](https://pak.dynasite.org/Nestimate/reference/entropy_network.md)
  : Transition Entropy Network
- [`entropy_trajectory()`](https://pak.dynasite.org/Nestimate/reference/entropy_trajectory.md)
  [`print(`*`<net_entropy_trajectory>`*`)`](https://pak.dynasite.org/Nestimate/reference/entropy_trajectory.md)
  [`summary(`*`<net_entropy_trajectory>`*`)`](https://pak.dynasite.org/Nestimate/reference/entropy_trajectory.md)
  [`plot(`*`<net_entropy_trajectory>`*`)`](https://pak.dynasite.org/Nestimate/reference/entropy_trajectory.md)
  : Sliding-Window Transition Entropy Trajectory
- [`entropy_bayes()`](https://pak.dynasite.org/Nestimate/reference/entropy_bayes.md)
  [`print(`*`<net_entropy_bayes>`*`)`](https://pak.dynasite.org/Nestimate/reference/entropy_bayes.md)
  [`print(`*`<net_entropy_bayes_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/entropy_bayes.md)
  [`summary(`*`<net_entropy_bayes>`*`)`](https://pak.dynasite.org/Nestimate/reference/entropy_bayes.md)
  [`plot(`*`<net_entropy_bayes>`*`)`](https://pak.dynasite.org/Nestimate/reference/entropy_bayes.md)
  : Bayesian Transition Entropy

## Network Pruning

Non-destructive edge pruning and restoration

- [`net_prune()`](https://pak.dynasite.org/Nestimate/reference/net_prune.md)
  : Prune a Network's Edges
- [`net_deprune()`](https://pak.dynasite.org/Nestimate/reference/net_deprune.md)
  : Undo Network Pruning
- [`net_reprune()`](https://pak.dynasite.org/Nestimate/reference/net_reprune.md)
  : Re-apply Network Pruning
- [`net_pruning_details()`](https://pak.dynasite.org/Nestimate/reference/net_pruning_details.md)
  [`print(`*`<net_pruning_details>`*`)`](https://pak.dynasite.org/Nestimate/reference/net_pruning_details.md)
  : Report Network Pruning Details

## Bootstrap & Inference

Statistical inference for network estimation

- [`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  [`print(`*`<net_bootstrap>`*`)`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  [`summary(`*`<net_bootstrap>`*`)`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  [`print(`*`<net_bootstrap_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  [`summary(`*`<net_bootstrap_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  [`print(`*`<wtna_boot_mixed>`*`)`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  [`summary(`*`<wtna_boot_mixed>`*`)`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
  : Bootstrap a Network Estimate

- [`vertex_bootstrap()`](https://pak.dynasite.org/Nestimate/reference/vertex_bootstrap.md)
  [`print(`*`<net_vertex_bootstrap>`*`)`](https://pak.dynasite.org/Nestimate/reference/vertex_bootstrap.md)
  [`summary(`*`<net_vertex_bootstrap>`*`)`](https://pak.dynasite.org/Nestimate/reference/vertex_bootstrap.md)
  [`plot(`*`<net_vertex_bootstrap>`*`)`](https://pak.dynasite.org/Nestimate/reference/vertex_bootstrap.md)
  : Vertex Bootstrap for Network-Level Statistics

- [`vertex_compare()`](https://pak.dynasite.org/Nestimate/reference/vertex_compare.md)
  [`print(`*`<net_vertex_comparison>`*`)`](https://pak.dynasite.org/Nestimate/reference/vertex_compare.md)
  [`summary(`*`<net_vertex_comparison>`*`)`](https://pak.dynasite.org/Nestimate/reference/vertex_compare.md)
  [`plot(`*`<net_vertex_comparison>`*`)`](https://pak.dynasite.org/Nestimate/reference/vertex_compare.md)
  : Compare Network-Level Statistics of Two Networks

- [`boot_glasso()`](https://pak.dynasite.org/Nestimate/reference/boot_glasso.md)
  [`print(`*`<boot_glasso>`*`)`](https://pak.dynasite.org/Nestimate/reference/boot_glasso.md)
  [`summary(`*`<boot_glasso>`*`)`](https://pak.dynasite.org/Nestimate/reference/boot_glasso.md)
  [`plot(`*`<boot_glasso>`*`)`](https://pak.dynasite.org/Nestimate/reference/boot_glasso.md)
  : Bootstrap for Regularized Partial Correlation Networks

- [`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
  [`print(`*`<net_permutation>`*`)`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
  [`summary(`*`<net_permutation>`*`)`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
  [`print(`*`<net_permutation_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
  [`summary(`*`<net_permutation_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
  [`print(`*`<wtna_perm_mixed>`*`)`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
  [`summary(`*`<wtna_perm_mixed>`*`)`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
  : Permutation Test for Network Comparison

- [`permutation_diagnostics()`](https://pak.dynasite.org/Nestimate/reference/permutation_diagnostics.md)
  : Does Nesting Bias a Permutation Test?

- [`nct()`](https://pak.dynasite.org/Nestimate/reference/nct.md)
  [`print(`*`<net_nct>`*`)`](https://pak.dynasite.org/Nestimate/reference/nct.md)
  [`summary(`*`<net_nct>`*`)`](https://pak.dynasite.org/Nestimate/reference/nct.md)
  : Network Comparison Test

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
  : Tables of a network comparison

- [`rename_models()`](https://pak.dynasite.org/Nestimate/reference/rename_models.md)
  :

  Rename the models of a `netobject_group`

- [`magnitude_difference()`](https://pak.dynasite.org/Nestimate/reference/magnitude_difference.md)
  [`print(`*`<magnitude_difference>`*`)`](https://pak.dynasite.org/Nestimate/reference/magnitude_difference.md)
  [`plot(`*`<magnitude_difference>`*`)`](https://pak.dynasite.org/Nestimate/reference/magnitude_difference.md)
  : Magnitude difference between the frequency and probability views

## Reliability & Stability

Assess reliability and stability of network estimates

- [`network_reliability()`](https://pak.dynasite.org/Nestimate/reference/network_reliability.md)
  [`print(`*`<net_reliability>`*`)`](https://pak.dynasite.org/Nestimate/reference/network_reliability.md)
  [`summary(`*`<net_reliability>`*`)`](https://pak.dynasite.org/Nestimate/reference/network_reliability.md)
  [`plot(`*`<net_reliability>`*`)`](https://pak.dynasite.org/Nestimate/reference/network_reliability.md)
  : Split-Half Reliability for Network Estimates
- [`casedrop_reliability()`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  [`print(`*`<net_casedrop_reliability>`*`)`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  [`summary(`*`<net_casedrop_reliability>`*`)`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  [`print(`*`<net_casedrop_reliability_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  [`summary(`*`<net_casedrop_reliability_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  [`print(`*`<summary.net_casedrop_reliability_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  [`plot(`*`<net_casedrop_reliability>`*`)`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  [`plot(`*`<net_casedrop_reliability_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/casedrop_reliability.md)
  : Edge-weight Case-dropping Stability
- [`centrality_stability()`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md)
  [`print(`*`<net_stability>`*`)`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md)
  [`print(`*`<net_stability_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md)
  [`summary(`*`<net_stability_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md)
  [`summary(`*`<net_stability>`*`)`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md)
  [`plot(`*`<net_stability>`*`)`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md)
  : Centrality Stability Coefficient (CS-coefficient)
- [`loading_stability()`](https://pak.dynasite.org/Nestimate/reference/loading_stability.md)
  [`print(`*`<pc_loading_stability>`*`)`](https://pak.dynasite.org/Nestimate/reference/loading_stability.md)
  [`plot(`*`<pc_loading_stability>`*`)`](https://pak.dynasite.org/Nestimate/reference/loading_stability.md)
  : Composite-Weight Stability Under Case Resampling

## Clustering & Grouping

Cluster-based and multilevel network analysis

- [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  [`print(`*`<net_clustering>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  [`summary(`*`<net_clustering>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  [`plot(`*`<net_clustering>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  [`print(`*`<tidy_covariates>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
  : Cluster Sequences by Dissimilarity
- [`cluster_data()`](https://pak.dynasite.org/Nestimate/reference/cluster_data.md)
  : Cluster sequence data (deprecated alias)
- [`cluster_choice()`](https://pak.dynasite.org/Nestimate/reference/cluster_choice.md)
  [`print(`*`<cluster_choice>`*`)`](https://pak.dynasite.org/Nestimate/reference/cluster_choice.md)
  [`summary(`*`<cluster_choice>`*`)`](https://pak.dynasite.org/Nestimate/reference/cluster_choice.md)
  [`plot(`*`<cluster_choice>`*`)`](https://pak.dynasite.org/Nestimate/reference/cluster_choice.md)
  : Cluster Choice – sweep k, dissimilarity and method
- [`cluster_diagnostics()`](https://pak.dynasite.org/Nestimate/reference/cluster_diagnostics.md)
  [`print(`*`<net_cluster_diagnostics>`*`)`](https://pak.dynasite.org/Nestimate/reference/cluster_diagnostics.md)
  [`plot(`*`<net_cluster_diagnostics>`*`)`](https://pak.dynasite.org/Nestimate/reference/cluster_diagnostics.md)
  [`as.data.frame(`*`<net_cluster_diagnostics>`*`)`](https://pak.dynasite.org/Nestimate/reference/cluster_diagnostics.md)
  : Cluster Diagnostics
- [`cluster_summary()`](https://pak.dynasite.org/Nestimate/reference/cluster_summary.md)
  : Cluster Summary Statistics
- [`cluster_mmm()`](https://pak.dynasite.org/Nestimate/reference/cluster_mmm.md)
  : Cluster sequences using Mixed Markov Models
- [`cluster_network()`](https://pak.dynasite.org/Nestimate/reference/cluster_network.md)
  : Cluster data and build per-cluster networks in one step
- [`print(`*`<mcml_layer>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
  [`build_mcml()`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
  [`print(`*`<mcml>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
  [`summary(`*`<mcml>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md)
  : Build MCML from Raw Transition Data
- [`build_mcml_pc()`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
  [`print(`*`<mcml_pc>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
  [`summary(`*`<mcml_pc>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
  [`plot(`*`<mcml_pc>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mcml_pc.md)
  : Multi-Cluster Multi-Level Aggregation for Psychometric Networks
- [`composites()`](https://pak.dynasite.org/Nestimate/reference/composites.md)
  : Cluster Scores From a Psychometric MCML Fit
- [`item_loadings()`](https://pak.dynasite.org/Nestimate/reference/item_loadings.md)
  : Item Diagnostics From a Psychometric MCML Fit
- [`macro_network()`](https://pak.dynasite.org/Nestimate/reference/macro_network.md)
  : Cluster-Level Network, With One Cluster Expanded
- [`build_mmm()`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
  [`print(`*`<net_mmm>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
  [`summary(`*`<net_mmm>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
  [`plot(`*`<net_mmm>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
  [`print(`*`<net_mmm_clustering>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
  [`plot(`*`<net_mmm_clustering>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
  : Fit a Mixed Markov Model
- [`compare_mmm()`](https://pak.dynasite.org/Nestimate/reference/compare_mmm.md)
  [`print(`*`<mmm_compare>`*`)`](https://pak.dynasite.org/Nestimate/reference/compare_mmm.md)
  [`summary(`*`<mmm_compare>`*`)`](https://pak.dynasite.org/Nestimate/reference/compare_mmm.md)
  [`plot(`*`<mmm_compare>`*`)`](https://pak.dynasite.org/Nestimate/reference/compare_mmm.md)
  : Compare MMM fits across different k
- [`session_ids()`](https://pak.dynasite.org/Nestimate/reference/session_ids.md)
  : The session behind each sequence

## Simplicial Complex Analysis

Topological analysis of networks

- [`build_simplicial()`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md)
  [`print(`*`<simplicial_complex>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md)
  [`plot(`*`<simplicial_complex>`*`)`](https://pak.dynasite.org/Nestimate/reference/build_simplicial.md)
  : Build a Simplicial Complex
- [`simplicial_features()`](https://pak.dynasite.org/Nestimate/reference/simplicial_features.md)
  : Tidy Topological Features for One or Many Networks
- [`persistent_homology()`](https://pak.dynasite.org/Nestimate/reference/persistent_homology.md)
  [`print(`*`<persistent_homology>`*`)`](https://pak.dynasite.org/Nestimate/reference/persistent_homology.md)
  [`plot(`*`<persistent_homology>`*`)`](https://pak.dynasite.org/Nestimate/reference/persistent_homology.md)
  : Persistent Homology
- [`bottleneck_distance()`](https://pak.dynasite.org/Nestimate/reference/bottleneck_distance.md)
  : Bottleneck Distance Between Persistence Diagrams
- [`persistence_landscape()`](https://pak.dynasite.org/Nestimate/reference/persistence_landscape.md)
  [`print(`*`<persistence_landscape>`*`)`](https://pak.dynasite.org/Nestimate/reference/persistence_landscape.md)
  [`plot(`*`<persistence_landscape>`*`)`](https://pak.dynasite.org/Nestimate/reference/persistence_landscape.md)
  : Persistence Landscape
- [`q_analysis()`](https://pak.dynasite.org/Nestimate/reference/q_analysis.md)
  [`print(`*`<q_analysis>`*`)`](https://pak.dynasite.org/Nestimate/reference/q_analysis.md)
  [`plot(`*`<q_analysis>`*`)`](https://pak.dynasite.org/Nestimate/reference/q_analysis.md)
  : Q-Analysis
- [`betti_numbers()`](https://pak.dynasite.org/Nestimate/reference/betti_numbers.md)
  : Betti Numbers
- [`euler_characteristic()`](https://pak.dynasite.org/Nestimate/reference/euler_characteristic.md)
  : Euler Characteristic
- [`simplicial_degree()`](https://pak.dynasite.org/Nestimate/reference/simplicial_degree.md)
  : Simplicial Degree
- [`verify_simplicial()`](https://pak.dynasite.org/Nestimate/reference/verify_simplicial.md)
  : Verify Simplicial Complex Against igraph

## Data Preparation

Convert and prepare data for network estimation

- [`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md)
  [`print(`*`<nestimate_data>`*`)`](https://pak.dynasite.org/Nestimate/reference/prepare.md)
  : Prepare Event Log Data for Network Estimation
- [`prepare_for_tna()`](https://pak.dynasite.org/Nestimate/reference/prepare_for_tna.md)
  : Prepare Data for TNA Analysis
- [`action_to_onehot()`](https://pak.dynasite.org/Nestimate/reference/action_to_onehot.md)
  : Convert Action Column to One-Hot Encoding
- [`prepare_onehot()`](https://pak.dynasite.org/Nestimate/reference/prepare_onehot.md)
  : Import One-Hot Encoded Data into Sequence Format
- [`wide_to_long()`](https://pak.dynasite.org/Nestimate/reference/wide_to_long.md)
  : Convert Wide Sequences to Long Format
- [`long_to_wide()`](https://pak.dynasite.org/Nestimate/reference/long_to_wide.md)
  : Convert Long Format to Wide Sequences
- [`convert_sequence_format()`](https://pak.dynasite.org/Nestimate/reference/convert_sequence_format.md)
  : Convert Sequence Data to Different Formats
- [`actor_endpoints()`](https://pak.dynasite.org/Nestimate/reference/actor_endpoints.md)
  : Tidy per-actor endpoint summary of a wide-format sequence dataset
- [`mark_first_state()`](https://pak.dynasite.org/Nestimate/reference/mark_first_state.md)
  : Mark leading-NA cells with an explicit state label
- [`mark_terminal_state()`](https://pak.dynasite.org/Nestimate/reference/mark_terminal_state.md)
  : Mark terminal-NA cells with an explicit state label

## Utilities

Helper functions and extractors

- [`predictability()`](https://pak.dynasite.org/Nestimate/reference/predictability.md)
  : Compute Node Predictability
- [`frequencies()`](https://pak.dynasite.org/Nestimate/reference/frequencies.md)
  [`summary(`*`<nest_transition_counts>`*`)`](https://pak.dynasite.org/Nestimate/reference/frequencies.md)
  : Build a Transition Frequency Matrix
- [`state_frequencies()`](https://pak.dynasite.org/Nestimate/reference/state_frequencies.md)
  : Compute State Frequencies from Trajectory Data
- [`net_aggregate_weights()`](https://pak.dynasite.org/Nestimate/reference/net_aggregate_weights.md)
  : Aggregate Edge Weights
- [`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
  [`plot(`*`<net_centrality>`*`)`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
  [`plot(`*`<net_centrality_group>`*`)`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
  : Compute Centrality Measures for a Network
- [`net_edge_betweenness()`](https://pak.dynasite.org/Nestimate/reference/net_edge_betweenness.md)
  [`plot(`*`<net_edge_betweenness>`*`)`](https://pak.dynasite.org/Nestimate/reference/net_edge_betweenness.md)
  : Edge Betweenness Network
- [`coefs()`](https://pak.dynasite.org/Nestimate/reference/coefs.md) :
  Tidy coefficients from a fitted mlvar model
- [`as_tna()`](https://pak.dynasite.org/Nestimate/reference/as_tna.md) :
  Promote the Layers of an mcml to Networks
- [`as_htna()`](https://pak.dynasite.org/Nestimate/reference/as_htna.md)
  : Build a grouped node-level network (htna) from data and a clustering
- [`as_networks()`](https://pak.dynasite.org/Nestimate/reference/as_networks.md)
  : Promote a psychometric MCML result to a network group
- [`as_netobject()`](https://pak.dynasite.org/Nestimate/reference/as_netobject.md)
  : Coerce a network object to a Nestimate netobject
- [`validate_netobject()`](https://pak.dynasite.org/Nestimate/reference/validate_netobject.md)
  : Validate a netobject / cograph_network against the shared schema
- [`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md)
  : Extract Edge List with Weights
- [`extract_initial_probs()`](https://pak.dynasite.org/Nestimate/reference/extract_initial_probs.md)
  [`summary(`*`<nest_initial_probs>`*`)`](https://pak.dynasite.org/Nestimate/reference/extract_initial_probs.md)
  : Extract Initial Probabilities from Model
- [`extract_transition_matrix()`](https://pak.dynasite.org/Nestimate/reference/extract_transition_matrix.md)
  [`summary(`*`<nest_transition_matrix>`*`)`](https://pak.dynasite.org/Nestimate/reference/extract_transition_matrix.md)
  : Extract Transition Matrix from Model
- [`state_colors()`](https://pak.dynasite.org/Nestimate/reference/state_colors.md)
  : The state colours an object will draw with
- [`set_state_colors()`](https://pak.dynasite.org/Nestimate/reference/set_state_colors.md)
  [`` `state_colors<-`() ``](https://pak.dynasite.org/Nestimate/reference/set_state_colors.md)
  : Set the state colours carried by a network object

## Sequence Analysis

Sequence visualization, pattern comparison, and association mining

- [`sequence_plot()`](https://pak.dynasite.org/Nestimate/reference/sequence_plot.md)
  [`print(`*`<mcml_sequence_plot>`*`)`](https://pak.dynasite.org/Nestimate/reference/sequence_plot.md)
  : Sequence Plot (heatmap, index, or distribution)
- [`distribution_plot()`](https://pak.dynasite.org/Nestimate/reference/distribution_plot.md)
  : State Distribution Plot Over Time
- [`print(`*`<nestimate_facet_plot>`*`)`](https://pak.dynasite.org/Nestimate/reference/plot_state_frequencies.md)
  [`print(`*`<nestimate_facet_list>`*`)`](https://pak.dynasite.org/Nestimate/reference/plot_state_frequencies.md)
  [`plot_state_frequencies()`](https://pak.dynasite.org/Nestimate/reference/plot_state_frequencies.md)
  [`print(`*`<state_freq>`*`)`](https://pak.dynasite.org/Nestimate/reference/plot_state_frequencies.md)
  [`plot(`*`<state_freq>`*`)`](https://pak.dynasite.org/Nestimate/reference/plot_state_frequencies.md)
  [`as.data.frame(`*`<state_freq>`*`)`](https://pak.dynasite.org/Nestimate/reference/plot_state_frequencies.md)
  : Plot State Frequency Distributions
- [`state_distribution()`](https://pak.dynasite.org/Nestimate/reference/state_distribution.md)
  : Per-Class State Distribution as a Tidy Data Frame
- [`plot_mosaic()`](https://pak.dynasite.org/Nestimate/reference/plot_mosaic.md)
  : Draw a Marimekko / Mosaic Plot from a Tidy Data Frame
- [`mosaic_plot()`](https://pak.dynasite.org/Nestimate/reference/mosaic_plot.md)
  : Mosaic Plot of a Network's Transition or Co-occurrence Counts
- [`mosaic_analysis()`](https://pak.dynasite.org/Nestimate/reference/mosaic_analysis.md)
  [`plot(`*`<mosaic_analysis>`*`)`](https://pak.dynasite.org/Nestimate/reference/mosaic_analysis.md)
  [`print(`*`<mosaic_analysis>`*`)`](https://pak.dynasite.org/Nestimate/reference/mosaic_analysis.md)
  [`summary(`*`<mosaic_analysis>`*`)`](https://pak.dynasite.org/Nestimate/reference/mosaic_analysis.md)
  : Two-variable mosaic analysis (chi-square test + flat mosaic)
- [`sequence_compare()`](https://pak.dynasite.org/Nestimate/reference/sequence_compare.md)
  [`print(`*`<net_sequence_comparison>`*`)`](https://pak.dynasite.org/Nestimate/reference/sequence_compare.md)
  [`summary(`*`<net_sequence_comparison>`*`)`](https://pak.dynasite.org/Nestimate/reference/sequence_compare.md)
  [`plot(`*`<net_sequence_comparison>`*`)`](https://pak.dynasite.org/Nestimate/reference/sequence_compare.md)
  : Compare Subsequence Patterns Between Groups
- [`extract_pathways()`](https://pak.dynasite.org/Nestimate/reference/extract_pathways.md)
  : Cut an Event Log into Pathways
- [`association_rules()`](https://pak.dynasite.org/Nestimate/reference/association_rules.md)
  [`print(`*`<net_association_rules>`*`)`](https://pak.dynasite.org/Nestimate/reference/association_rules.md)
  [`summary(`*`<net_association_rules>`*`)`](https://pak.dynasite.org/Nestimate/reference/association_rules.md)
  [`plot(`*`<net_association_rules>`*`)`](https://pak.dynasite.org/Nestimate/reference/association_rules.md)
  : Discover Association Rules from Sequential or Transaction Data

## Link Prediction

Predict and evaluate missing connections

- [`predict_links()`](https://pak.dynasite.org/Nestimate/reference/predict_links.md)
  [`print(`*`<net_link_prediction>`*`)`](https://pak.dynasite.org/Nestimate/reference/predict_links.md)
  [`summary(`*`<net_link_prediction>`*`)`](https://pak.dynasite.org/Nestimate/reference/predict_links.md)
  : Predict Missing or Future Links in a Network
- [`evaluate_links()`](https://pak.dynasite.org/Nestimate/reference/evaluate_links.md)
  : Evaluate Link Predictions Against Known Edges

## Outcome Modelling

Relate network and sequence features to an outcome

- [`outcome_model()`](https://pak.dynasite.org/Nestimate/reference/outcome_model.md)
  [`print(`*`<net_outcome_model>`*`)`](https://pak.dynasite.org/Nestimate/reference/outcome_model.md)
  [`summary(`*`<net_outcome_model>`*`)`](https://pak.dynasite.org/Nestimate/reference/outcome_model.md)
  [`plot(`*`<net_outcome_model>`*`)`](https://pak.dynasite.org/Nestimate/reference/outcome_model.md)
  : Model Unit-Level Outcomes from Sequence or Network Predictors
- [`effects_table()`](https://pak.dynasite.org/Nestimate/reference/effects_table.md)
  : Effect Table of a Fitted Outcome Model

## Data

Example datasets

- [`human_long`](https://pak.dynasite.org/Nestimate/reference/long-data.md)
  [`ai_long`](https://pak.dynasite.org/Nestimate/reference/long-data.md)
  : Human-AI Vibe Coding Interaction Data (Long Format)
- [`srl_strategies`](https://pak.dynasite.org/Nestimate/reference/srl_strategies.md)
  : Self-Regulated Learning Strategy Frequencies
- [`learning_activities`](https://pak.dynasite.org/Nestimate/reference/learning_activities.md)
  : Online Learning Activity Indicators
- [`group_regulation_long`](https://pak.dynasite.org/Nestimate/reference/group_regulation_long.md)
  : Group Regulation in Collaborative Learning (Long Format)
- [`chatgpt_srl`](https://pak.dynasite.org/Nestimate/reference/chatgpt_srl.md)
  : ChatGPT Self-Regulated Learning Scale Scores
- [`trajectories`](https://pak.dynasite.org/Nestimate/reference/trajectories.md)
  : Student Engagement Trajectories
