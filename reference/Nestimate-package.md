# Nestimate: Dynamic, Probabilistic, and Higher-Order Network Analysis

Nestimate estimates networks from sequences, event logs, one-hot codes
and panel data, compares them, and tests how far they can be trusted.
One entry point,
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
estimates transition networks (transition probabilities, frequencies,
co-occurrence and attention-weighted transitions; Saqr et al. 2025,
[doi:10.1145/3706468.3706513](https://doi.org/10.1145/3706468.3706513) )
and psychological networks (correlation, partial correlation, graphical
lasso and Ising models, Saqr, Beck and López-Pernas 2024,
[doi:10.1007/978-3-031-54464-4_19](https://doi.org/10.1007/978-3-031-54464-4_19)
; and mixed graphical models), and returns a `netobject`. The other
verbs, by the question they answer (most read a `netobject`; the
panel-data ones read the data directly):

- inference on a network:
  [`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md),
  [`boot_glasso()`](https://pak.dynasite.org/Nestimate/reference/boot_glasso.md),
  split-half reliability with
  [`network_reliability()`](https://pak.dynasite.org/Nestimate/reference/network_reliability.md)
  and centrality stability with
  [`centrality_stability()`](https://pak.dynasite.org/Nestimate/reference/centrality_stability.md);

- comparison of networks:
  [`compare_networks()`](https://pak.dynasite.org/Nestimate/reference/compare_networks.md),
  and the permutation tests
  [`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
  and [`nct()`](https://pak.dynasite.org/Nestimate/reference/nct.md);

- the structure of a network:
  [`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md),
  [`chain_structure()`](https://pak.dynasite.org/Nestimate/reference/chain_structure.md),
  [`passage_time()`](https://pak.dynasite.org/Nestimate/reference/passage_time.md),
  [`markov_stability()`](https://pak.dynasite.org/Nestimate/reference/markov_stability.md),
  [`stability_landscape()`](https://pak.dynasite.org/Nestimate/reference/stability_landscape.md)
  and
  [`transition_entropy()`](https://pak.dynasite.org/Nestimate/reference/transition_entropy.md);

- subgroups: clustering with
  [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md),
  mixed Markov models with
  [`build_mmm()`](https://pak.dynasite.org/Nestimate/reference/build_mmm.md)
  and multi-cluster multi-layer networks with
  [`build_mcml()`](https://pak.dynasite.org/Nestimate/reference/build_mcml.md);

- higher-order networks:
  [`build_hon()`](https://pak.dynasite.org/Nestimate/reference/build_hon.md),
  higher-order network embedding with
  [`build_honem()`](https://pak.dynasite.org/Nestimate/reference/build_honem.md),
  hyper-path anomalies with
  [`build_hypa()`](https://pak.dynasite.org/Nestimate/reference/build_hypa.md),
  the multi-order generative model
  [`build_mogen()`](https://pak.dynasite.org/Nestimate/reference/build_mogen.md),
  and
  [`markov_order_test()`](https://pak.dynasite.org/Nestimate/reference/markov_order_test.md);

- panel data: multilevel vector autoregression with
  [`build_mlvar()`](https://pak.dynasite.org/Nestimate/reference/build_mlvar.md)
  and person-specific networks with
  [`build_gimme()`](https://pak.dynasite.org/Nestimate/reference/build_gimme.md).

## References

Saqr, M., López-Pernas, S., Törmänen, T., Kaliisa, R., Misiejuk, K., &
Tikka, S. (2025). Transition network analysis: A novel framework for
modeling, visualizing, and identifying the temporal patterns of learners
and learning processes. In *Proceedings of the 15th International
Learning Analytics and Knowledge Conference (LAK25)*. ACM.
[doi:10.1145/3706468.3706513](https://doi.org/10.1145/3706468.3706513)

Saqr, M., Beck, E., & López-Pernas, S. (2024). Psychological networks: A
modern approach to analysis of learning and complex learning processes.
In M. Saqr & S. López-Pernas (Eds.), *Learning Analytics Methods and
Tutorials: A Practical Guide Using R*. Springer.
[doi:10.1007/978-3-031-54464-4_19](https://doi.org/10.1007/978-3-031-54464-4_19)

## See also

Useful links:

- <https://github.com/mohsaqr/Nestimate>

- <https://pak.dynasite.org/Nestimate/>

- Report bugs at <https://github.com/mohsaqr/Nestimate/issues>

## Author

**Maintainer**: Mohammed Saqr <saqr@saqr.me> \[copyright holder\]

Authors:

- Sonsoles López-Pernas <sonsoles.lopez@uef.fi>

- Kamila Misiejuk <kamila.misiejuk@fernuni-hagen.de>
