#' @description
#' Nestimate estimates networks from sequences, event logs, one-hot codes and
#' panel data, compares them, and tests how far they can be trusted. One
#' entry point, [build_network()], estimates transition networks
#' (transition probabilities, frequencies, co-occurrence and
#' attention-weighted transitions; Saqr et al. 2025, \doi{10.1145/3706468.3706513}) and
#' psychological networks (correlation, partial correlation, graphical lasso
#' and Ising models, Saqr, Beck and L\enc{ó}{o}pez-Pernas 2024,
#' \doi{10.1007/978-3-031-54464-4_19}; and mixed graphical models), and
#' returns a `netobject`. The other verbs, by the question they answer (most
#' read a `netobject`; the panel-data ones read the data directly):
#' \itemize{
#'   \item inference on a network: [bootstrap_network()], [boot_glasso()],
#'     split-half reliability with [network_reliability()] and centrality
#'     stability with [centrality_stability()];
#'   \item comparison of networks: [compare_networks()], and the
#'     permutation tests [permutation()] and [nct()];
#'   \item the structure of a network: [net_centrality()],
#'     [chain_structure()], [passage_time()], [markov_stability()],
#'     [stability_landscape()] and [transition_entropy()];
#'   \item subgroups: clustering with [build_clusters()], mixed Markov models
#'     with [build_mmm()] and multi-cluster multi-layer networks with
#'     [build_mcml()];
#'   \item higher-order networks: [build_hon()], higher-order network
#'     embedding with [build_honem()], hyper-path anomalies with
#'     [build_hypa()], the multi-order generative model [build_mogen()], and
#'     [markov_order_test()];
#'   \item panel data: multilevel vector autoregression with [build_mlvar()]
#'     and person-specific networks with [build_gimme()].
#' }
#'
#' @references
#' Saqr, M., L\enc{ó}{o}pez-Pernas, S., T\enc{ö}{o}rm\enc{ä}{a}nen, T.,
#' Kaliisa, R., Misiejuk, K., & Tikka, S. (2025). Transition network
#' analysis: A novel framework for modeling, visualizing, and identifying the
#' temporal patterns of learners and learning processes. In
#' \emph{Proceedings of the 15th International Learning Analytics and
#' Knowledge Conference (LAK25)}. ACM. \doi{10.1145/3706468.3706513}
#'
#' Saqr, M., Beck, E., & L\enc{ó}{o}pez-Pernas, S. (2024). Psychological
#' networks: A modern approach to analysis of learning and complex learning
#' processes. In M. Saqr & S. L\enc{ó}{o}pez-Pernas (Eds.), \emph{Learning
#' Analytics Methods and Tutorials: A Practical Guide Using R}. Springer.
#' \doi{10.1007/978-3-031-54464-4_19}
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import ggplot2
#' @importFrom stats cor var sd median quantile setNames aggregate
#'   approx spline cmdscale coef confint cov fitted lm na.omit
#'   p.adjust pnorm predict pt qnorm qt residuals rnorm runif
#'   rbinom rpois reshape symnum ecdf complete.cases cov2cor
#' @importFrom utils head tail combn modifyList
#' @importFrom data.table data.table setDT rbindlist
## usethis namespace: end
NULL

.onLoad <- function(libname, pkgname) {
  .register_builtin_estimators() # nocov
}

# Global variable declarations for R CMD check
utils::globalVariables(c(
  ".data",
  "from", "to", "weight", "value", "name", "variable",
  "edge", "node", "state", "time", "onset", "terminus",
  "lower", "upper", "significant", "p_value", "ci_lower", "ci_upper",
  "metric", "centrality", "vertex", "group",
  # data.table NSE column names
  ".orig_row", ".seq_grp", ".grp_key",
  # permutation() data.table columns
  "weight_x", "weight_y",
  # mmm ggplot2 aes variables
  "max_posterior", "cluster", "k", "criterion",
  # simplicial plot aes variables
  "threshold", "betti", "dim_label", "birth", "death",
  "persistence", "components", "count", "total",
  # association_rules plot aes variables
  "support", "confidence", "lift"
))
