# ---- Permutation Test for Network Comparison ----

#' Permutation Test for Network Comparison
#'
#' @description
#' Tests whether two networks estimated by \code{\link{build_network}}
#' differ more than chance would produce. The sequences (or rows) of both
#' networks are pooled, the group labels are shuffled \code{iter} times, both
#' networks are re-estimated on every shuffle, and the observed differences
#' are compared with the shuffled ones. Works with every built-in method and
#' with registered custom estimators.
#'
#' @section What is tested:
#' Two kinds of question are answered from the same shuffles.
#' \describe{
#'   \item{Edge tests}{One test per edge: is the difference in this edge's
#'     weight, \code{x - y}, larger than the shuffles produce? Reported in
#'     \code{summary()} with an effect size (observed difference divided by
#'     the SD of the shuffled differences) and a p-value
#'     \code{(number of shuffles at least as extreme + 1) / (iter + 1)}. With
#'     many edges, some fall below \code{alpha} by chance; use \code{adjust}
#'     to correct for that.}
#'   \item{Global test}{One test for the whole network: do the two networks
#'     differ at all? Two statistics, as in the Network Comparison Test
#'     (van Borkulo et al., 2023): \strong{M}, the sum of the absolute edge
#'     differences (the total amount of difference), and \strong{S}, the
#'     largest absolute edge difference. Being a single test, it needs no
#'     multiplicity correction. It is shown by \code{print()}.}
#' }
#' The smallest attainable p-value is \code{1 / (iter + 1)}; with the
#' default \code{iter = 1000} it is 0.000999, meaning no shuffle came close.
#'
#' @section Nested data and \code{actor}:
#' Shuffling treats every sequence as an exchangeable unit. When several
#' sequences come from the same actor (sessions of one person, students of
#' one team), \code{actor} names the column identifying that actor. The
#' shuffle then respects it (Good, 2005; Anderson & ter Braak, 2003): a
#' person whose sequences are all in one group moves to the other group
#' as a whole; a person with sequences in both groups has their labels
#' shuffled among their own sequences only. Mixed designs combine the two.
#' The observed differences do not change; only the p-values and effect
#' sizes do. With few persons there are few distinct ways to shuffle them,
#' and a warning (class \code{nestimate_few_actors}) is raised when no
#' p-value could fall below \code{alpha}.
#'
#' \code{actor} is available for transition networks (\code{"relative"},
#' \code{"frequency"}, \code{"co_occurrence"}). Association networks
#' (\code{"cor"}, \code{"pcor"}, \code{"glasso"}, ...) do not keep row
#' identifiers after estimation and raise \code{nestimate_actor_unsupported}.
#'
#' @section ICC and design effect:
#' With \code{actor}, \code{print()} also reports:
#' \describe{
#'   \item{ICC}{The intraclass correlation, the proportion of the total
#'     variance that lies between actors (Shrout & Fleiss, 1979). An ICC
#'     close to 0 indicates little evidence of a nesting effect. Computed as
#'     the one-way ANOVA ICC of each sequence's transition shares within each
#'     group, averaged over edges weighted by edge frequency, jackknife
#'     bias-corrected, with a 95% interval from the leave-one-actor-out
#'     jackknife (Efron & Tibshirani, 1993).}
#'   \item{Design effect}{The ratio of the variance of an estimate under the
#'     clustered design to its variance had the units been sampled
#'     independently (Kish, 1965). Here: the variance of the shuffled
#'     differences when whole actors are moved, divided by the variance when
#'     single sequences are moved, both drawn in the same run. Reported as
#'     the median over edges and for the global statistic M. For equal
#'     numbers of sequences per actor m, Kish gives the approximation \code{1 + (m - 1) * ICC}.}
#' }
#'
#' @section Reading the printed output:
#' \preformatted{
#' Permutation Test: Transition Network (relative probabilities) [directed]
#'   Iterations: 1000  |  Alpha: 0.05  |  Actor: Group (200 actors)
#'   Nodes: 9  |  Edges tested: 78  |  Significant: 42
#'   Global test (networks differ overall?): M = 2.612 (p = 0.000999)  |  ...
#'   Nesting in Group: ICC = -0.002 [95\% CI -0.006, 0.002]  |  between design
#'   Design effect (1 = nesting does not matter): edges 1.03  |  global 1.20
#' }
#' Line 2: settings, and the actor column with its number of actors.
#' Line 3: edges present in either network and how many differ at
#' \code{alpha}. Line 4: the global test. Lines 5-6, only with
#' \code{actor}: the ICC with its interval, whether actors sit in one group
#' (\code{between}), in both (\code{within}) or either (\code{mixed}), and
#' the design effects.
#'
#' @section Other inputs:
#' For transition methods, per-sequence count matrices are computed once
#' and each shuffle only re-sums them, which keeps large \code{iter} fast.
#' For association methods the estimator is re-run on every shuffle. If a
#' transition network rests on a single sequence, a warning (class
#' \code{nestimate_single_sequence}) says it cannot be validated by
#' resampling.
#'
#' \code{permutation()} also accepts two \code{\link{net_edge_betweenness}}
#' objects. It then permutes the source networks, recomputes edge
#' betweenness for each shuffle, and tests the edge-betweenness
#' differences. Both objects must come from the same source method and use
#' the same \code{invert} setting.
#'
#' @param x A \code{netobject} (from \code{\link{build_network}}) or a
#'   \code{\link{net_edge_betweenness}} object.
#' @param y A \code{netobject} (from \code{\link{build_network}}) or a
#'   \code{\link{net_edge_betweenness}} object.
#'   Must use the same method and have the same nodes as \code{x}.
#'   Default \code{NULL}: when \code{x} is a \code{netobject_group} (or an
#'   \code{mcml}) and \code{y} is left \code{NULL}, every pair of groups is
#'   tested and the result is a \code{net_permutation_group} named
#'   \code{"<group i> vs <group j>"}.
#' @param iter Integer. Number of permutation iterations (default: 1000).
#' @param alpha Numeric. Significance level (default: 0.05).
#' @param paired Logical. If \code{TRUE}, permute within pairs (requires
#'   equal number of observations in \code{x} and \code{y}). Default: FALSE.
#' @param adjust Character. p-value adjustment method passed to
#'   \code{\link[stats]{p.adjust}} (default: \code{"none"}). Common choices:
#'   \code{"holm"}, \code{"BH"}, \code{"bonferroni"}.
#' @param measures Character vector of centrality measures to permutation-test
#'   in addition to the edges, or \code{"all"} for every built-in measure.
#'   Default \code{NULL} (edges only). When supplied, the result gains a
#'   \code{$centralities} block matching the layout of
#'   \code{tna::permutation_test(measures = )}: per state and measure it reports
#'   the observed difference, an effect size (difference / SD of the permutation
#'   null), and a permutation p-value, all using the same permuted networks as
#'   the edge test. Not supported for \code{net_edge_betweenness} inputs.
#' @param nlambda Integer. Number of lambda values for the EBIC-glasso
#'   regularisation path (only used when \code{method = "glasso"}).
#'   Higher values give finer lambda resolution at the cost of speed.
#'   Default: 50.
#' @param seed Integer or NULL. RNG seed for reproducibility.
#' @param actor Character or NULL. Name of the column identifying the actor
#'   each sequence belongs to: the person whose sessions they are, or the
#'   team of a student. Looked up in the network's \code{$metadata} (e.g.
#'   \code{"student_id"} when sessions are nested in students,
#'   \code{"Group"} for students nested in teams) or in its wide sequence
#'   data. When supplied, whole actors are permuted and the ICC and design
#'   effect are reported; see the sections \emph{Nested data and actor} and
#'   \emph{ICC and design effect}. Supported for transition
#'   methods (\code{"relative"}, \code{"frequency"}, \code{"co_occurrence"});
#'   cannot be combined with \code{paired = TRUE}, which is the special case
#'   of one actor per pair. Default \code{NULL}: sequences are permuted
#'   individually.
#'
#' @return An object of class \code{"net_permutation"} containing:
#' \describe{
#'   \item{x}{The first \code{netobject}.}
#'   \item{y}{The second \code{netobject}.}
#'   \item{diff}{Observed difference matrix (\code{x - y}).}
#'   \item{diff_sig}{Observed difference where \code{p < alpha}, else 0.}
#'   \item{p_values}{P-value matrix (adjusted if \code{adjust != "none"}).}
#'   \item{effect_size}{Effect size matrix (observed diff / SD of permutation diffs).}
#'   \item{summary}{Long-format data frame, one row per edge present in
#'     either network (undirected networks keep one row per unordered
#'     pair), with columns \code{from}, \code{to}, \code{weight_x},
#'     \code{weight_y}, \code{diff}, \code{effect_size}, \code{p_value},
#'     \code{sig}.}
#'   \item{global}{Data frame of the two NCT-style global statistics, one
#'     row each: \code{statistic} (\code{"M"}, the sum of absolute edge
#'     differences, and \code{"S"}, the largest absolute edge difference),
#'     \code{observed}, and \code{p_value} from the same permutation null
#'     as the edge test. Absent on the edge-betweenness path.}
#'   \item{method}{The network estimation method.}
#'   \item{source_method}{For edge-betweenness tests, the source network method.}
#'   \item{iter}{Number of permutation iterations.}
#'   \item{alpha}{Significance level used.}
#'   \item{paired}{Whether paired permutation was used.}
#'   \item{adjust}{p-value adjustment method used.}
#'   \item{actor}{The \code{actor} column name, or \code{NULL}.}
#'   \item{n_actors}{Number of distinct actors, or \code{NULL}.}
#'   \item{null_sd}{Matrix of the SD of each edge difference over the
#'     permutation null (the effect-size denominator).}
#'   \item{null_sd_m}{SD of the global \code{M} statistic over the
#'     permutation null. Absent on the edge-betweenness path.}
#'   \item{clustering}{Present only with \code{actor}. One-row data frame:
#'     \code{n_sequences}, \code{n_actors}, \code{design}
#'     (\code{"between"}, \code{"within"}, \code{"mixed"}), \code{icc} with
#'     \code{icc_ci_lower}/\code{icc_ci_upper} (how alike sequences of one
#'     actor are; see \code{\link{permutation_diagnostics}}),
#'     \code{deff_edges} (median over edges) and \code{deff_global} (for
#'     \code{M}): the actor-level over the sequence-level null variance, drawn in the
#'     same run (the design effect; Kish, 1965). \code{min_p} is the
#'     smallest attainable p-value.}
#'   \item{clustering_edges}{Present only with \code{actor}. One row per
#'     edge of \code{summary}: \code{from}, \code{to}, \code{icc},
#'     \code{null_sd_actor}, \code{null_sd_sequence}, \code{deff}.}
#'   \item{centralities}{Present only when \code{measures} is supplied. A list
#'     with \code{stats} (one row per state-by-measure: \code{state},
#'     \code{centrality}, \code{diff_true}, \code{effect_size}, \code{p_value}),
#'     \code{diffs_true} (wide observed differences), and \code{diffs_sig}
#'     (observed differences where \code{p < alpha}, else 0).}
#' }
#' Grouped input returns a \code{"net_permutation_group"} (a named list of
#' \code{net_permutation} results): one element per matching group name
#' when both \code{x} and \code{y} are \code{netobject_group}s, or one per
#' group pair when \code{y} is \code{NULL}. Two \code{wtna_mixed} inputs
#' return a \code{"wtna_perm_mixed"} with \code{$transition} and
#' \code{$cooccurrence} results.
#'
#' @examples
#' s1 <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
#' s2 <- data.frame(V1 = c("A","C","B"), V2 = c("C","B","A"))
#' n1 <- build_network(s1, method = "relative")
#' n2 <- build_network(s2, method = "relative")
#' perm <- permutation(n1, n2, iter = 10)
#' \donttest{
#' set.seed(1)
#' d1 <- data.frame(V1 = sample(LETTERS[1:4], 20, TRUE),
#'                  V2 = sample(LETTERS[1:4], 20, TRUE),
#'                  V3 = sample(LETTERS[1:4], 20, TRUE))
#' d2 <- data.frame(V1 = sample(LETTERS[1:4], 20, TRUE),
#'                  V2 = sample(LETTERS[1:4], 20, TRUE),
#'                  V3 = sample(LETTERS[1:4], 20, TRUE))
#' net1 <- build_network(d1, method = "relative")
#' net2 <- build_network(d2, method = "relative")
#' perm <- permutation(net1, net2, iter = 100, seed = 42)
#' print(perm)
#' summary(perm)
#'
#' # Students are nested in teams, and Achiever is a team-level label:
#' # permute whole teams, not single students
#' net <- build_network(group_regulation_long, method = "relative",
#'                      actor = "Actor", action = "Action", time = "Time",
#'                      group = "Achiever")
#' permutation(net, iter = 100, actor = "Group", seed = 1)
#' }
#'
#' @references
#' Anderson, M. J., & ter Braak, C. J. F. (2003). Permutation tests for
#' multi-factorial analysis of variance. \emph{Journal of Statistical
#' Computation and Simulation}, 73(2), 85-113.
#'
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An Introduction to the
#' Bootstrap}. Chapman & Hall.
#'
#' Good, P. (2005). \emph{Permutation, Parametric, and Bootstrap Tests of
#' Hypotheses} (3rd ed.). Springer.
#'
#' Kish, L. (1965). \emph{Survey Sampling}. Wiley.
#'
#' Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in
#' assessing rater reliability. \emph{Psychological Bulletin}, 86(2),
#' 420-428.
#'
#' van Borkulo, C. D., van Bork, R., Boschloo, L., Kossakowski, J. J.,
#' Tio, P., Schoevers, R. A., Borsboom, D., & Waldorp, L. J. (2023).
#' Comparing network structures on three aspects: A permutation test.
#' \emph{Psychological Methods}, 28(6), 1273-1285.
#'
#' @seealso \code{\link{permutation_diagnostics}} to compare the
#'   actor-level and ordinary tests side by side; \code{\link{bayes_compare}} for the Bayesian complement: instead of
#'   "is this difference more extreme than chance?" it answers "how probable is
#'   a difference, and how large?";
#'   \code{\link{build_network}}, \code{\link{bootstrap_network}},
#'   \code{\link{print.net_permutation}},
#'   \code{\link{summary.net_permutation}}
#'
#' @importFrom stats p.adjust sd
#' @export
permutation <- function(x, y = NULL,
                             iter = 1000L,
                             alpha = 0.05,
                             paired = FALSE,
                             adjust = "none",
                             measures = NULL,
                             nlambda = 50L,
                             seed = NULL,
                             actor = NULL) {

  # ---- wtna_mixed dispatch: permute both components ----
  if (inherits(x, "wtna_mixed") || inherits(y, "wtna_mixed")) {
    if (!inherits(x, "wtna_mixed") || !inherits(y, "wtna_mixed")) {
      stop("Both x and y must be wtna_mixed objects.", call. = FALSE)
    }
    result <- list(
      transition = permutation(
        x$transition, y$transition, iter = iter, alpha = alpha,
        paired = paired, adjust = adjust, measures = measures, nlambda = nlambda,
        seed = seed, actor = actor
      ),
      cooccurrence = permutation(
        x$cooccurrence, y$cooccurrence, iter = iter, alpha = alpha,
        paired = paired, adjust = adjust, measures = measures,
        nlambda = nlambda, seed = seed, actor = actor
      )
    )
    class(result) <- "wtna_perm_mixed"
    return(result)
  }

  # ---- mcml dispatch: convert to netobject_group via as_tna ----
  if (inherits(x, "mcml")) x <- as_tna(x)
  if (inherits(y, "mcml")) y <- as_tna(y)

  # ---- Single netobject_group: all-pairs permutation tests ----
  if (inherits(x, "netobject_group") && is.null(y)) {
    grp_names <- names(x)
    n_grps <- length(grp_names)
    if (n_grps < 2L) {
      stop("Need at least 2 groups for pairwise permutation tests.",
           call. = FALSE)
    }
    pairs <- combn(n_grps, 2L)
    results <- lapply(seq_len(ncol(pairs)), function(k) {
      i <- pairs[1L, k]
      j <- pairs[2L, k]
      permutation(x[[i]], x[[j]], iter = iter, alpha = alpha,
                  paired = paired, adjust = adjust, measures = measures,
                  nlambda = nlambda, seed = seed, actor = actor)
    })
    pair_labels <- vapply(seq_len(ncol(pairs)), function(k) {
      paste(grp_names[pairs[1L, k]], "vs", grp_names[pairs[2L, k]])
    }, character(1))
    names(results) <- pair_labels
    class(results) <- c("net_permutation_group", "list")
    return(results)
  }

  # ---- netobject_group dispatch: permute each matching element ----
  if (inherits(x, "netobject_group") && inherits(y, "netobject_group")) {
    common <- intersect(names(x), names(y))
    if (length(common) == 0L) {
      stop("No matching group names between x and y.", call. = FALSE)
    }
    results <- lapply(common, function(nm) {
      permutation(x[[nm]], y[[nm]], iter = iter, alpha = alpha,
                  paired = paired, adjust = adjust, measures = measures,
                  nlambda = nlambda, seed = seed, actor = actor)
    })
    names(results) <- common
    class(results) <- c("net_permutation_group", "list")
    return(results)
  }

  # ---- Coerce cograph_network inputs ----
  if (inherits(x, "cograph_network")) x <- .as_netobject(x)
  if (inherits(y, "cograph_network")) y <- .as_netobject(y)

  # ---- Edge-betweenness dispatch: permute source networks, compare EB ----
  if (inherits(x, "net_edge_betweenness") ||
      inherits(y, "net_edge_betweenness")) {
    return(.permutation_edge_betweenness(
      x = x, y = y, iter = iter, alpha = alpha, paired = paired,
      adjust = adjust, measures = measures, nlambda = nlambda, seed = seed,
      actor = actor
    ))
  }

  # ---- Input validation ----
  stopifnot(
    inherits(x, "netobject"),
    inherits(y, "netobject"),
    is.numeric(iter), length(iter) == 1, iter >= 2,
    is.numeric(alpha), length(alpha) == 1, alpha > 0, alpha < 1,
    is.logical(paired), length(paired) == 1,
    is.character(adjust), length(adjust) == 1
  )
  iter <- as.integer(iter)
  .check_permutation_block(actor, paired)

  if (is.null(x$data)) {
    stop("'x' does not contain $data. Rebuild with build_network().",
         call. = FALSE)
  }
  if (is.null(y$data)) {
    stop("'y' does not contain $data. Rebuild with build_network().",
         call. = FALSE)
  }

  if (x$method != y$method) {
    stop("Methods must match: x uses '", x$method,
         "', y uses '", y$method, "'.", call. = FALSE)
  }

  if (!setequal(x$nodes$label, y$nodes$label)) {
    stop("Nodes must be the same in both networks.", call. = FALSE)
  }

  # Ensure same node order
  nodes <- x$nodes$label
  if (!identical(x$nodes$label, y$nodes$label)) {
    y$weights <- y$weights[nodes, nodes]
  }

  method <- .resolve_method_alias(x$method)
  directed <- x$directed
  n_nodes <- length(nodes)

  if (paired) {
    if (nrow(x$data) != nrow(y$data)) {
      stop("Paired test requires equal number of observations in x and y.",
           call. = FALSE)
    }
  }

  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1)
    set.seed(seed)
  }

  # ---- Observed difference ----
  obs_diff <- x$weights - y$weights

  # ---- Centrality differences (optional, tna-parity) ----
  if (!is.null(measures)) {
    if (length(measures) == 1L && identical(tolower(measures), "all")) {
      measures <- .centrality_all_measures()
    }
    bad <- setdiff(measures, .centrality_builtin_measures())
    if (length(bad) > 0L) {
      stop("Unknown measures: ", paste(bad, collapse = ", "),
           ". Options: ", paste(.centrality_builtin_measures(), collapse = ", "),
           call. = FALSE)
    }
    obs_cent_diff <- .perm_cent_diff_mat(x$weights, y$weights, nodes,
                                         directed, measures)
  } else {
    obs_cent_diff <- NULL
  }

  # ---- Dispatch permutation ----
  has_data_x <- is.data.frame(x$data) && ncol(x$data) > 0L
  has_data_y <- is.data.frame(y$data) && ncol(y$data) > 0L
  if (method %in% c("relative", "frequency", "co_occurrence")) {
    if (!has_data_x || !has_data_y) {
      stop("Permutation test requires the original data stored in the netobject. ",
           "For wtna/cna networks, use wtna() directly instead of ",
           "build_network(method='cna').", call. = FALSE)
    }
    x$data <- .resampling_transition_data(x$data, x$metadata, x$params)
    y$data <- .resampling_transition_data(y$data, y$metadata, y$params)
    n_seq_x <- .transition_resampling_n_sequences(x$data, x$params)
    n_seq_y <- .transition_resampling_n_sequences(y$data, y$params)
    if ((!is.na(n_seq_x) && n_seq_x <= 1L) ||
        (!is.na(n_seq_y) && n_seq_y <= 1L)) {
      .single_sequence_notice()
    }
    perm_result <- .permutation_transition(
      x = x, y = y, nodes = nodes, method = method,
      iter = iter, paired = paired,
      measures = measures, directed = directed, obs_cent_diff = obs_cent_diff,
      block = actor, alpha = alpha
    )
  } else {
    .stop_block_unsupported(actor, method)
    perm_result <- .permutation_association(
      x = x, y = y, nodes = nodes, method = method,
      iter = iter, paired = paired, nlambda = nlambda,
      measures = measures, directed = directed, obs_cent_diff = obs_cent_diff
    )
  }

  # ---- P-values ----
  # (sum(|perm_diff| >= |obs_diff|) + 1) / (iter + 1)
  obs_flat <- as.vector(obs_diff)
  p_values_flat <- (perm_result$exceed_counts + 1L) / (iter + 1L)

  # Apply multiple comparison correction
  p_values_flat <- p.adjust(p_values_flat, method = adjust)

  p_mat <- matrix(p_values_flat, n_nodes, n_nodes,
                  dimnames = list(nodes, nodes))

  # ---- Effect size ----
  # Cohen's d style: observed_diff / sd(perm_diffs)
  perm_sd <- perm_result$perm_sd
  perm_sd[perm_sd == 0] <- NA_real_
  es_flat <- obs_flat / perm_sd
  es_flat[is.na(es_flat)] <- 0
  es_mat <- matrix(es_flat, n_nodes, n_nodes,
                   dimnames = list(nodes, nodes))

  # ---- Significant diff ----
  sig_mask <- (p_mat < alpha) * 1
  diff_sig <- obs_diff * sig_mask

  # ---- Summary ----
  summary_df <- .build_permutation_summary(
    obs_diff = obs_diff,
    p_mat = p_mat,
    es_mat = es_mat,
    x_matrix = x$weights,
    y_matrix = y$weights,
    nodes = nodes,
    directed = directed,
    alpha = alpha
  )

  # ---- Assemble result ----
  result <- list(
    x           = x,
    y           = y,
    diff        = obs_diff,
    diff_sig    = diff_sig,
    p_values    = p_mat,
    effect_size = es_mat,
    summary     = summary_df,
    method      = method,
    iter        = iter,
    alpha       = alpha,
    paired      = paired,
    adjust      = adjust,
    actor       = actor,
    n_actors    = perm_result$n_blocks,
    global      = .build_permutation_global(perm_result$global, iter),
    null_sd     = matrix(perm_result$perm_sd, n_nodes, n_nodes,
                         dimnames = list(nodes, nodes)),
    null_sd_m   = perm_result$global$m_null_sd
  )
  if (!is.null(perm_result$clustering)) {
    clus <- .build_permutation_clustering(perm_result, summary_df, nodes, iter)
    result$clustering <- clus$overall
    result$clustering_edges <- clus$edges
  }

  # ---- Centrality permutation block (tna-parity) ----
  if (!is.null(measures)) {
    result$centralities <- .build_permutation_centralities(
      obs_cent_diff = obs_cent_diff,
      cent_exceed = perm_result$cent_exceed,
      cent_sum = perm_result$cent_sum,
      cent_sumsq = perm_result$cent_sumsq,
      states = nodes, measures = measures,
      iter = iter, alpha = alpha, adjust = adjust
    )
  }

  class(result) <- "net_permutation"
  result
}

# ---- Global (NCT-style) statistics ----

#' Global permutation statistics: M (sum |diff|) and S (max |diff|)
#'
#' Same permutation null as the edge test; p = (exceed + 1) / (iter + 1).
#' @noRd
.build_permutation_global <- function(g, iter) {
  data.frame(
    statistic = c("M", "S"),
    observed = c(g$m_obs, g$s_obs),
    p_value = (c(g$m_exceed, g$s_exceed) + 1L) / (iter + 1L),
    stringsAsFactors = FALSE
  )
}

# ---- Centrality permutation helpers (tna-parity) ----

#' Centralities of a permuted weight matrix, tna::centralities() settings
#' @noRd
.perm_centralities_mat <- function(mat, nodes, directed, measures) {
  dimnames(mat) <- list(nodes, nodes)
  res <- .compute_centralities(
    mat, nodes, directed, measures,
    loops = FALSE, normalize = FALSE, invert = TRUE,
    normalize_diffusion = FALSE
  )
  vapply(measures, function(m) {
    v <- res[[.centrality_canonical_measure(m)]]
    if (is.null(v)) rep(NA_real_, length(nodes)) else unname(v)
  }, numeric(length(nodes)))
}

#' Observed/permuted centrality difference matrix (states x measures)
#' @noRd
.perm_cent_diff_mat <- function(wx, wy, nodes, directed, measures) {
  cx <- .perm_centralities_mat(wx, nodes, directed, measures)
  cy <- .perm_centralities_mat(wy, nodes, directed, measures)
  m <- cx - cy
  dim(m) <- c(length(nodes), length(measures))
  dimnames(m) <- list(nodes, measures)
  m
}

#' Assemble the $centralities block exactly like tna::permutation_test()
#' @noRd
.build_permutation_centralities <- function(obs_cent_diff, cent_exceed,
                                            cent_sum, cent_sumsq, states,
                                            measures, iter, alpha, adjust) {
  cent_p <- (cent_exceed + 1) / (iter + 1)
  cent_p[] <- stats::p.adjust(as.vector(cent_p), method = adjust)
  cent_mean <- cent_sum / iter
  cent_sd <- sqrt(pmax(cent_sumsq / iter - cent_mean^2, 0))
  effect <- obs_cent_diff / cent_sd          # tna: diff / sd (no guard)
  sig <- obs_cent_diff * (cent_p < alpha)
  state_f <- factor(states, levels = states)
  stats_df <- expand.grid(state = state_f, centrality = measures,
                          KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  stats_df$centrality <- factor(stats_df$centrality, levels = measures)
  stats_df$diff_true <- as.vector(obs_cent_diff)
  stats_df$effect_size <- as.vector(effect)
  stats_df$p_value <- as.vector(cent_p)
  diffs_true <- data.frame(state = state_f, as.data.frame(obs_cent_diff),
                           check.names = FALSE)
  diffs_sig <- data.frame(state = state_f, as.data.frame(sig),
                          check.names = FALSE)
  rownames(diffs_true) <- NULL
  rownames(diffs_sig) <- NULL
  list(stats = stats_df, diffs_true = diffs_true, diffs_sig = diffs_sig)
}


# ---- Edge-betweenness permutation path ----

#' Source-network proxy for a net_edge_betweenness object
#' @noRd
.edge_betweenness_source_net <- function(x) {
  if (!inherits(x, "net_edge_betweenness")) {
    stop("Both x and y must be net_edge_betweenness objects.",
         call. = FALSE)
  }
  source_method <- x$edge_betweenness$source_method
  if (is.null(source_method) || length(source_method) != 1L ||
      is.na(source_method) || !nzchar(source_method)) {
    stop(
      "net_edge_betweenness object does not carry its source method. ",
      "Recreate it with net_edge_betweenness() before permutation().",
      call. = FALSE
    )
  }
  src <- x
  src$method <- source_method
  if (!is.null(x$source_weights)) {
    src$weights <- x$source_weights
  }
  src$params <- x$params %||% list()
  src$scaling <- x$scaling
  src$threshold <- x$threshold %||% 0
  src
}

#' Permutation test for edge-betweenness differences
#' @noRd
.permutation_edge_betweenness <- function(x, y, iter, alpha, paired,
                                          adjust, measures, nlambda, seed,
                                          actor = NULL) {
  if (!inherits(x, "net_edge_betweenness") ||
      !inherits(y, "net_edge_betweenness")) {
    stop("Both x and y must be net_edge_betweenness objects.",
         call. = FALSE)
  }
  if (!is.null(measures)) {
    stop("`measures` is not supported for edge-betweenness permutation tests.",
         call. = FALSE)
  }

  stopifnot(
    is.numeric(iter), length(iter) == 1, iter >= 2,
    is.numeric(alpha), length(alpha) == 1, alpha > 0, alpha < 1,
    is.logical(paired), length(paired) == 1,
    is.character(adjust), length(adjust) == 1
  )
  iter <- as.integer(iter)
  .check_permutation_block(actor, paired)

  if (is.null(x$data)) {
    stop("'x' does not contain $data. Rebuild with build_network().",
         call. = FALSE)
  }
  if (is.null(y$data)) {
    stop("'y' does not contain $data. Rebuild with build_network().",
         call. = FALSE)
  }

  if (!setequal(x$nodes$label, y$nodes$label)) {
    stop("Nodes must be the same in both networks.", call. = FALSE)
  }

  nodes <- x$nodes$label
  if (!identical(x$nodes$label, y$nodes$label)) {
    y$weights <- y$weights[nodes, nodes]
  }

  x_src <- .edge_betweenness_source_net(x)
  y_src <- .edge_betweenness_source_net(y)
  method <- .resolve_method_alias(x_src$method)
  y_method <- .resolve_method_alias(y_src$method)
  if (method != y_method) {
    stop("Source methods must match: x uses '", method,
         "', y uses '", y_method, "'.", call. = FALSE)
  }

  directed <- x$directed
  n_nodes <- length(nodes)

  if (paired && nrow(x$data) != nrow(y$data)) {
    stop("Paired test requires equal number of observations in x and y.",
         call. = FALSE)
  }

  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1)
    set.seed(seed)
  }

  obs_diff <- x$weights - y$weights
  invert <- isTRUE(x$edge_betweenness$invert)
  if (!identical(isTRUE(y$edge_betweenness$invert), invert)) {
    stop("x and y must use the same edge-betweenness `invert` setting.",
         call. = FALSE)
  }
  transform_eb <- function(mat) {
    dimnames(mat) <- list(nodes, nodes)
    .edge_betweenness(mat, invert = invert)
  }

  has_data_x <- is.data.frame(x_src$data) && ncol(x_src$data) > 0L
  has_data_y <- is.data.frame(y_src$data) && ncol(y_src$data) > 0L
  if (method %in% c("relative", "frequency", "co_occurrence")) {
    if (!has_data_x || !has_data_y) {
      stop("Permutation test requires the original data stored in the netobject. ",
           "For wtna/cna networks, use wtna() directly instead of ",
           "build_network(method='cna').", call. = FALSE)
    }
    x_src$data <- .resampling_transition_data(
      x_src$data, x_src$metadata, x_src$params
    )
    y_src$data <- .resampling_transition_data(
      y_src$data, y_src$metadata, y_src$params
    )
    n_seq_x <- .transition_resampling_n_sequences(x_src$data, x_src$params)
    n_seq_y <- .transition_resampling_n_sequences(y_src$data, y_src$params)
    if ((!is.na(n_seq_x) && n_seq_x <= 1L) ||
        (!is.na(n_seq_y) && n_seq_y <= 1L)) {
      .single_sequence_notice()
    }
    perm_result <- .permutation_transition(
      x = x_src, y = y_src, nodes = nodes, method = method,
      iter = iter, paired = paired,
      measures = NULL, directed = directed, obs_cent_diff = NULL,
      transform = transform_eb, obs_diff = obs_diff,
      block = actor, alpha = alpha
    )
  } else {
    .stop_block_unsupported(actor, method)
    perm_result <- .permutation_association(
      x = x_src, y = y_src, nodes = nodes, method = method,
      iter = iter, paired = paired, nlambda = nlambda,
      measures = NULL, directed = directed, obs_cent_diff = NULL,
      transform = transform_eb, obs_diff = obs_diff
    )
  }

  obs_flat <- as.vector(obs_diff)
  p_values_flat <- (perm_result$exceed_counts + 1L) / (iter + 1L)
  p_values_flat <- p.adjust(p_values_flat, method = adjust)
  p_mat <- matrix(p_values_flat, n_nodes, n_nodes,
                  dimnames = list(nodes, nodes))

  perm_sd <- perm_result$perm_sd
  perm_sd[perm_sd == 0] <- NA_real_
  es_flat <- obs_flat / perm_sd
  es_flat[is.na(es_flat)] <- 0
  es_mat <- matrix(es_flat, n_nodes, n_nodes,
                   dimnames = list(nodes, nodes))

  diff_sig <- obs_diff * ((p_mat < alpha) * 1)
  summary_df <- .build_permutation_summary(
    obs_diff = obs_diff,
    p_mat = p_mat,
    es_mat = es_mat,
    x_matrix = x$weights,
    y_matrix = y$weights,
    nodes = nodes,
    directed = directed,
    alpha = alpha
  )

  result <- list(
    x = x,
    y = y,
    diff = obs_diff,
    diff_sig = diff_sig,
    p_values = p_mat,
    effect_size = es_mat,
    summary = summary_df,
    method = "edge_betweenness",
    source_method = method,
    iter = iter,
    alpha = alpha,
    paired = paired,
    adjust = adjust,
    actor = actor,
    n_actors = perm_result$n_blocks,
    null_sd = matrix(perm_result$perm_sd, n_nodes, n_nodes,
                     dimnames = list(nodes, nodes)),
    edge_betweenness = list(invert = invert, source_method = method)
  )
  if (!is.null(perm_result$clustering)) {
    clus <- .build_permutation_clustering(perm_result, summary_df, nodes, iter)
    result$clustering <- clus$overall
    result$clustering_edges <- clus$edges
  }
  class(result) <- "net_permutation"
  result
}


# ---- Transition fast path ----

#' Permutation test for transition networks via pre-computed counts
#' @noRd
.permutation_transition <- function(x, y, nodes, method, iter, paired,
                                    measures = NULL, directed = TRUE,
                                    obs_cent_diff = NULL,
                                    transform = NULL,
                                    obs_diff = NULL,
                                    block = NULL, alpha = 0.05) {
  n_nodes <- length(nodes)
  nbins <- n_nodes * n_nodes
  is_relative <- method == "relative"
  do_cent <- !is.null(measures)
  if (is.null(transform)) transform <- identity
  if (do_cent) {
    cent_abs_true <- abs(obs_cent_diff)
    cent_exceed <- matrix(0L, n_nodes, length(measures))
    cent_sum <- matrix(0, n_nodes, length(measures))
    cent_sumsq <- matrix(0, n_nodes, length(measures))
  }

  # Pre-compute per-sequence counts for both groups
  trans_x <- .precompute_per_sequence(x$data, method, x$params, nodes)
  trans_y <- .precompute_per_sequence(y$data, method, y$params, nodes)

  n_x <- nrow(trans_x)
  n_y <- nrow(trans_y)

  # Pool sequences
  pooled <- rbind(trans_x, trans_y)
  n_total <- n_x + n_y

  # Blocked design: whole persons/teams move together (see ?permutation)
  block_ids <- if (is.null(block)) NULL else
    c(.permutation_block_ids(x, block, n_x, "x"),
      .permutation_block_ids(y, block, n_y, "y"))
  design <- if (is.null(block)) NULL else .block_design(
    block_ids = block_ids, is_x = rep(c(TRUE, FALSE), c(n_x, n_y)),
    block = block, alpha = alpha
  )
  if (!is.null(design)) {
    # unblocked reference null, drawn alongside, gives the design effect
    ref_sum <- numeric(nbins)
    ref_sumsq <- numeric(nbins)
    ref_m_sum <- 0
    ref_m_sumsq <- 0
    split_diff <- function(in_x) {
      mat_x <- .postprocess_counts(colSums(pooled[in_x, , drop = FALSE]),
                                   n_nodes, is_relative, x$scaling, x$threshold)
      mat_y <- .postprocess_counts(colSums(pooled[-in_x, , drop = FALSE]),
                                   n_nodes, is_relative, y$scaling, y$threshold)
      as.vector(transform(mat_x)) - as.vector(transform(mat_y))
    }
  }

  # Observed diff (recomputed from counts for consistency)
  obs_flat <- as.vector(obs_diff %||% (x$weights - y$weights))

  # Running counters
  exceed_counts <- integer(nbins)
  sum_diffs <- numeric(nbins)
  sum_diffs_sq <- numeric(nbins)
  # Global statistics (NCT-style): M = sum |diff|, S = max |diff|
  m_obs <- sum(abs(obs_flat))
  s_obs <- max(abs(obs_flat))
  m_exceed <- 0L
  s_exceed <- 0L
  m_sum <- 0
  m_sumsq <- 0

  for (i in seq_len(iter)) {
    if (paired) {
      # Paired: randomly swap x/y within each pair
      swaps <- sample(c(TRUE, FALSE), n_x, replace = TRUE)
      idx_x <- ifelse(swaps, seq(n_x + 1L, n_total), seq_len(n_x))
      idx_y <- ifelse(swaps, seq_len(n_x), seq(n_x + 1L, n_total))
      counts_x <- colSums(pooled[idx_x, , drop = FALSE])
      counts_y <- colSums(pooled[idx_y, , drop = FALSE])
    } else if (!is.null(design)) {
      # Blocked: reassign whole blocks / shuffle within crossed blocks
      perm_is_x <- .block_permute(design)
      counts_x <- colSums(pooled[perm_is_x, , drop = FALSE])
      counts_y <- colSums(pooled[!perm_is_x, , drop = FALSE])
    } else {
      # Unpaired: shuffle group labels
      idx_x <- sample.int(n_total, n_x)
      counts_x <- colSums(pooled[idx_x, , drop = FALSE])
      counts_y <- colSums(pooled[-idx_x, , drop = FALSE])
    }

    # Post-process each group to get network matrix
    mat_x <- .postprocess_counts(counts_x, n_nodes, is_relative,
                                 x$scaling, x$threshold)
    mat_y <- .postprocess_counts(counts_y, n_nodes, is_relative,
                                 y$scaling, y$threshold)

    mat_x <- transform(mat_x)
    mat_y <- transform(mat_y)
    perm_diff <- as.vector(mat_x) - as.vector(mat_y)

    # Accumulate
    exceed_counts <- exceed_counts + (abs(perm_diff) >= abs(obs_flat))
    sum_diffs <- sum_diffs + perm_diff
    sum_diffs_sq <- sum_diffs_sq + perm_diff^2
    m_exceed <- m_exceed + (sum(abs(perm_diff)) >= m_obs)
    s_exceed <- s_exceed + (max(abs(perm_diff)) >= s_obs)
    m_sum <- m_sum + sum(abs(perm_diff))
    m_sumsq <- m_sumsq + sum(abs(perm_diff))^2

    if (!is.null(design)) {
      ref_diff <- split_diff(sample.int(n_total, n_x))
      ref_sum <- ref_sum + ref_diff
      ref_sumsq <- ref_sumsq + ref_diff^2
      ref_m_sum <- ref_m_sum + sum(abs(ref_diff))
      ref_m_sumsq <- ref_m_sumsq + sum(abs(ref_diff))^2
    }

    if (do_cent) {
      cd <- .perm_cent_diff_mat(mat_x, mat_y, nodes, directed, measures)
      cent_exceed <- cent_exceed + (abs(cd) >= cent_abs_true)
      cent_sum <- cent_sum + cd
      cent_sumsq <- cent_sumsq + cd^2
    }
  }

  # SD of permutation diffs
  perm_mean <- sum_diffs / iter
  perm_sd <- sqrt(pmax(sum_diffs_sq / iter - perm_mean^2, 0))

  out <- list(
    exceed_counts = exceed_counts,
    perm_sd = perm_sd,
    global = list(m_obs = m_obs, s_obs = s_obs,
                  m_exceed = m_exceed, s_exceed = s_exceed,
                  m_null_sd = sqrt(max(m_sumsq / iter - (m_sum / iter)^2, 0))),
    n_blocks = design$n_blocks
  )
  if (!is.null(design)) {
    ref_mean <- ref_sum / iter
    out$clustering <- list(
      design = design,
      n_sequences = n_total,
      icc = .block_icc(pooled, block_ids, ifelse(design$is_x, "x", "y")),
      ref_sd = sqrt(pmax(ref_sumsq / iter - ref_mean^2, 0)),
      ref_m_sd = sqrt(max(ref_m_sumsq / iter - (ref_m_sum / iter)^2, 0))
    )
  }
  if (do_cent) {
    out$cent_exceed <- cent_exceed
    out$cent_sum <- cent_sum
    out$cent_sumsq <- cent_sumsq
  }
  out
}


# ---- Blocked (cluster-level) permutation ----

#' Validate the `block` argument
#' @noRd
.check_permutation_block <- function(block, paired) {
  if (is.null(block)) return(invisible(NULL))
  if (!is.character(block) || length(block) != 1L || is.na(block) ||
      !nzchar(block)) {
    stop(errorCondition("`actor` must be a single column name or NULL.",
                        class = "nestimate_bad_actor", call = NULL))
  }
  if (isTRUE(paired)) {
    stop(errorCondition(
      paste0("`actor` and `paired = TRUE` cannot be combined: a paired ",
             "design is the special case of one actor per pair. Use one."),
      class = "nestimate_bad_actor", call = NULL))
  }
  invisible(NULL)
}

#' Blocked permutation needs per-row identifiers, which association
#' networks do not keep after estimation
#' @noRd
.stop_block_unsupported <- function(block, method) {
  if (is.null(block)) return(invisible(NULL))
  stop(errorCondition(
    sprintf(paste0("`actor` is supported for transition networks ",
                   "(relative, frequency, co_occurrence), not method '%s': ",
                   "its netobject does not keep row identifiers."), method),
    class = "nestimate_actor_unsupported", call = NULL))
}

#' Block id of every sequence row of a netobject
#'
#' Looked up in `$metadata` (aligned with `$data` rows for long-format
#' builds) or in the wide sequence data itself.
#' @noRd
.permutation_block_ids <- function(net, block, n_rows, label) {
  meta <- net$metadata
  ids <- if (is.data.frame(meta) && block %in% names(meta)) {
    meta[[block]]
  } else if (is.data.frame(net$data) && block %in% names(net$data)) {
    net$data[[block]]
  } else {
    stop(errorCondition(
      sprintf("`actor = \"%s\"` is not a column of %s's metadata or sequence data.",
              block, label),
      class = "nestimate_actor_missing", call = NULL))
  }
  if (length(ids) != n_rows) {
    stop(errorCondition(
      sprintf(paste0("`actor = \"%s\"` has %d values but %s has %d sequences; ",
                     "this needs one actor id per sequence."),
              block, length(ids), label, n_rows),
      class = "nestimate_actor_misaligned", call = NULL))
  }
  if (anyNA(ids)) {
    stop(errorCondition(
      sprintf("`actor = \"%s\"` has %d missing values in %s.",
              block, sum(is.na(ids)), label),
      class = "nestimate_actor_missing", call = NULL))
  }
  as.character(ids)
}

#' Precompute the blocked permutation scheme
#'
#' Blocks entirely inside one group ("pure") are reassigned as whole units;
#' blocks with sequences in both groups ("crossed") keep their per-block
#' group counts and shuffle labels internally. Warns when the design admits
#' too few distinct arrangements to reach `alpha`.
#' @noRd
.block_design <- function(block_ids, is_x, block, alpha, warn = TRUE) {
  n_in_block <- tapply(is_x, block_ids, length)
  n_x_in_block <- tapply(is_x, block_ids, sum)
  is_pure <- n_x_in_block == 0L | n_x_in_block == n_in_block
  pure <- names(n_in_block)[is_pure]
  crossed <- names(n_in_block)[!is_pure]

  pure_rows <- which(block_ids %in% pure)
  crossed_rows <- which(block_ids %in% crossed)
  crossed_block <- block_ids[crossed_rows]

  # log number of distinct label arrangements the scheme can produce
  n_pure_x <- sum(n_x_in_block[pure] > 0L)
  log_arrangements <- lchoose(length(pure), n_pure_x) +
    sum(lchoose(n_in_block[crossed], n_x_in_block[crossed]))
  if (warn && log_arrangements < log(1 / alpha)) {
    warning(warningCondition(
      sprintf(paste0("`actor = \"%s\"` allows only %d distinct permutations, ",
                     "so no p-value can fall below %.3g (alpha = %g). ",
                     "More actors are needed for this test to detect anything."),
              block, as.integer(round(exp(log_arrangements))),
              exp(-log_arrangements), alpha),
      class = "nestimate_few_actors", call = NULL))
  }

  list(
    is_x = is_x,
    pure_rows = pure_rows,
    pure_index = match(block_ids[pure_rows], pure),
    pure_label = unname(n_x_in_block[pure] > 0L),
    crossed_rows = crossed_rows,
    crossed_block = crossed_block,
    crossed_sorted = crossed_rows[order(crossed_block)],
    n_blocks = length(n_in_block),
    n_pure = length(pure),
    n_crossed = length(crossed),
    log_arrangements = log_arrangements
  )
}

#' One blocked permutation of the group labels
#' @noRd
.block_permute <- function(design) {
  perm_is_x <- design$is_x
  if (length(design$pure_label) > 1L) {
    shuffled <- design$pure_label[sample.int(length(design$pure_label))]
    perm_is_x[design$pure_rows] <- shuffled[design$pure_index]
  }
  if (length(design$crossed_rows) > 0L) {
    # rows grouped by block in random within-block order receive the
    # block-grouped original labels: a uniform shuffle inside each block
    target <- design$crossed_rows[order(design$crossed_block,
                                        stats::runif(length(design$crossed_rows)))]
    perm_is_x[target] <- design$is_x[design$crossed_sorted]
  }
  perm_is_x
}


#' Tidy clustering diagnostics of a blocked permutation
#'
#' Design effect = blocked null variance / unblocked null variance, per edge
#' and for the global M statistic (Kish, 1965). One-row overall table plus a
#' one-row-per-edge table restricted to the edges in the summary.
#' @noRd
.build_permutation_clustering <- function(perm_result, summary_df, nodes, iter) {
  clus <- perm_result$clustering
  design <- clus$design
  n_nodes <- length(nodes)
  as_mat <- function(v) matrix(v, n_nodes, n_nodes, dimnames = list(nodes, nodes))
  sd_blocked <- as_mat(perm_result$perm_sd)
  sd_unblocked <- as_mat(clus$ref_sd)
  icc_edges <- matrix(clus$icc$per_edge, n_nodes, n_nodes, byrow = TRUE,
                      dimnames = list(nodes, nodes))
  at <- cbind(summary_df$from, summary_df$to)
  edges <- data.frame(
    from = summary_df$from, to = summary_df$to,
    icc = icc_edges[at],
    null_sd_actor = sd_blocked[at], null_sd_sequence = sd_unblocked[at],
    deff = (sd_blocked[at] / sd_unblocked[at])^2,
    stringsAsFactors = FALSE
  )
  overall <- data.frame(
    n_sequences = clus$n_sequences,
    n_actors = design$n_blocks,
    design = if (design$n_crossed == 0L) "between"
             else if (design$n_pure == 0L) "within" else "mixed",
    icc = clus$icc$estimate,
    icc_ci_lower = clus$icc$ci[1L], icc_ci_upper = clus$icc$ci[2L],
    deff_edges = stats::median(edges$deff[is.finite(edges$deff)]),
    deff_global = (perm_result$global$m_null_sd / clus$ref_m_sd)^2,
    min_p = max(exp(-design$log_arrangements), 1 / (iter + 1)),
    stringsAsFactors = FALSE
  )
  list(overall = overall, edges = edges)
}


#' Convert flat count vector to network matrix with post-processing
#' @noRd
.postprocess_counts <- function(counts, n_nodes, is_relative,
                                scaling, threshold) {
  mat <- matrix(counts, n_nodes, n_nodes, byrow = TRUE)
  if (is_relative) {
    rs <- rowSums(mat)
    nz <- rs > 0
    mat[nz, ] <- mat[nz, ] / rs[nz]
  }
  if (!is.null(scaling)) mat <- .apply_scaling(mat, scaling) # nocov start
  if (threshold > 0) mat[abs(mat) < threshold] <- 0 # nocov end
  mat
}


# ---- Association path (optimized) ----

#' Permutation test for association networks
#'
#' Pre-cleans pooled data once, then uses lightweight per-iteration
#' estimation (direct cor/solve/glasso calls) to avoid repeated
#' input validation overhead. For custom/unknown estimators, falls
#' back to full estimator calls.
#' @noRd
.permutation_association <- function(x, y, nodes, method, iter, paired,
                                    nlambda = 50L, measures = NULL,
                                    directed = FALSE, obs_cent_diff = NULL,
                                    transform = NULL,
                                    obs_diff = NULL) {
  n_nodes <- length(nodes)
  nbins <- n_nodes * n_nodes
  do_cent <- !is.null(measures)
  if (is.null(transform)) transform <- identity
  if (do_cent) {
    cent_abs_true <- abs(obs_cent_diff)
    cent_exceed <- matrix(0L, n_nodes, length(measures))
    cent_sum <- matrix(0, n_nodes, length(measures))
    cent_sumsq <- matrix(0, n_nodes, length(measures))
  }

  # $data is already cleaned by the estimator (numeric matrix, no NAs,
  # no zero-variance columns) - just pool directly
  n_x <- nrow(x$data)
  n_y <- nrow(y$data)
  pooled_mat <- rbind(x$data, y$data)
  n_total <- n_x + n_y

  # Extract params
  params_x <- x$params
  cor_method <- .param_get(params_x, "cor_method", "pearson")
  threshold_x <- x$threshold
  threshold_y <- y$threshold
  scaling_x <- x$scaling
  scaling_y <- y$scaling
  obs_flat <- as.vector(obs_diff %||% (x$weights - y$weights))
  m_obs <- sum(abs(obs_flat))
  s_obs <- max(abs(obs_flat))
  m_exceed <- 0L
  s_exceed <- 0L
  m_sum <- 0
  m_sumsq <- 0

  # Select fast path based on method
  use_fast <- method %in% c("cor", "pcor", "glasso")

  # Pre-compute the glasso lambda path once from the pooled correlation;
  # every permutation iteration reuses it via psychnets::ebic_glasso
  if (use_fast && method == "glasso") {
    gamma <- .param_get(params_x, "gamma", 0.5)
    penalize_diag <- .param_get(params_x, "penalize.diagonal", FALSE)

    S_pooled <- cor(pooled_mat, method = cor_method)
    perm_rholist <- .compute_lambda_path(S_pooled, nlambda, 0.01)
    p_glasso <- ncol(pooled_mat)
  }

  # Build the per-iteration estimator function
  if (use_fast) {
    estimate_from_rows <- switch(method,
      cor = function(mat_subset) {
        S <- cor(mat_subset, method = cor_method)
        diag(S) <- 0
        S
      },
      pcor = function(mat_subset) {
        S <- cor(mat_subset, method = cor_method)
        Wi <- tryCatch(solve(S), error = function(e) NULL)
        if (is.null(Wi)) return(NULL) # nocov
        .precision_to_pcor(Wi, threshold = 0)
      },
      glasso = function(mat_subset) {
        S <- cor(mat_subset, method = cor_method)
        n_obs <- nrow(mat_subset)
        # Solver delegated to psychnets; the pooled lambda path is fixed
        # once and reused per iteration, refit = FALSE reproduces the
        # former select-from-path semantics exactly. Iteration-level
        # failures stay NULL by the permutation contract (counted and
        # reported upstream).
        fit <- tryCatch(
          psychnets::ebic_glasso(cor_matrix = S, n = n_obs, gamma = gamma,
                                 lambda_path = perm_rholist,
                                 penalize_diagonal = penalize_diag,
                                 refit = FALSE, native = TRUE),
          error = function(e) NULL
        )
        if (is.null(fit)) return(NULL) # nocov
        .precision_to_pcor(fit$precision, threshold = 0)
      }
    )
  } else {
    # Fallback: full estimator for custom methods
    estimator <- get_estimator(method)
    estimate_from_rows <- function(mat_subset) {
      df <- as.data.frame(mat_subset)
      est <- tryCatch(
        do.call(estimator$fn, c(list(data = df), params_x)),
        error = function(e) NULL
      )
      if (is.null(est)) return(NULL)
      mat <- est$matrix
      if (!identical(rownames(mat), nodes)) {
        common <- intersect(nodes, rownames(mat))
        if (length(common) < n_nodes) return(NULL)
        mat <- mat[nodes, nodes] # nocov
      }
      mat
    }
  }

  # Running counters
  exceed_counts <- integer(nbins)
  sum_diffs <- numeric(nbins)
  sum_diffs_sq <- numeric(nbins)

  for (i in seq_len(iter)) {
    if (paired) {
      swaps <- sample(c(TRUE, FALSE), n_x, replace = TRUE)
      idx_x <- ifelse(swaps, seq(n_x + 1L, n_total), seq_len(n_x))
      idx_y <- ifelse(swaps, seq_len(n_x), seq(n_x + 1L, n_total))
    } else {
      idx_x <- sample.int(n_total, n_x)
      idx_y <- seq_len(n_total)[-idx_x]
    }

    mat_x <- estimate_from_rows(pooled_mat[idx_x, , drop = FALSE])
    mat_y <- estimate_from_rows(pooled_mat[idx_y, , drop = FALSE])

    if (is.null(mat_x) || is.null(mat_y)) next

    # Apply scaling and threshold
    if (!is.null(scaling_x)) mat_x <- .apply_scaling(mat_x, scaling_x) # nocov
    if (threshold_x > 0) mat_x[abs(mat_x) < threshold_x] <- 0
    if (!is.null(scaling_y)) mat_y <- .apply_scaling(mat_y, scaling_y) # nocov
    if (threshold_y > 0) mat_y[abs(mat_y) < threshold_y] <- 0

    mat_x <- transform(mat_x)
    mat_y <- transform(mat_y)
    perm_diff <- as.vector(mat_x) - as.vector(mat_y)

    exceed_counts <- exceed_counts + (abs(perm_diff) >= abs(obs_flat))
    sum_diffs <- sum_diffs + perm_diff
    sum_diffs_sq <- sum_diffs_sq + perm_diff^2
    m_exceed <- m_exceed + (sum(abs(perm_diff)) >= m_obs)
    s_exceed <- s_exceed + (max(abs(perm_diff)) >= s_obs)
    m_sum <- m_sum + sum(abs(perm_diff))
    m_sumsq <- m_sumsq + sum(abs(perm_diff))^2

    if (do_cent) {
      cd <- .perm_cent_diff_mat(mat_x, mat_y, nodes, directed, measures)
      cent_exceed <- cent_exceed + (abs(cd) >= cent_abs_true)
      cent_sum <- cent_sum + cd
      cent_sumsq <- cent_sumsq + cd^2
    }
  }

  perm_mean <- sum_diffs / iter
  perm_sd <- sqrt(pmax(sum_diffs_sq / iter - perm_mean^2, 0))

  out <- list(
    exceed_counts = exceed_counts,
    perm_sd = perm_sd,
    global = list(m_obs = m_obs, s_obs = s_obs,
                  m_exceed = m_exceed, s_exceed = s_exceed,
                  m_null_sd = sqrt(max(m_sumsq / iter - (m_sum / iter)^2, 0)))
  )
  if (do_cent) {
    out$cent_exceed <- cent_exceed
    out$cent_sum <- cent_sum
    out$cent_sumsq <- cent_sumsq
  }
  out
}


# ---- Summary builder ----

#' Build long-format summary data frame from permutation test results
#' @noRd
.build_permutation_summary <- function(obs_diff, p_mat, es_mat,
                                       x_matrix, y_matrix,
                                       nodes, directed, alpha) {
  n <- length(nodes)
  dt <- data.table::data.table(
    from        = rep(nodes, each = n),
    to          = rep(nodes, times = n),
    weight_x    = as.vector(t(x_matrix)),
    weight_y    = as.vector(t(y_matrix)),
    diff        = as.vector(t(obs_diff)),
    effect_size = as.vector(t(es_mat)),
    p_value     = as.vector(t(p_mat)),
    sig         = as.vector(t(p_mat)) < alpha
  )

  # Filter: keep edges present in either network
  if (directed) {
    dt <- dt[weight_x != 0 | weight_y != 0]
  } else {
    dt <- dt[(weight_x != 0 | weight_y != 0) & from <= to]
  }

  as.data.frame(dt)
}


# ---- S3 Methods ----

#' Print Method for net_permutation
#'
#' @param x A \code{net_permutation} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' s1 <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
#' s2 <- data.frame(V1 = c("A","C","B"), V2 = c("C","B","A"))
#' n1 <- build_network(s1, method = "relative")
#' n2 <- build_network(s2, method = "relative")
#' perm <- permutation(n1, n2, iter = 10)
#' print(perm)
#' \donttest{
#' set.seed(1)
#' d1 <- data.frame(V1 = c("A","B","A"), V2 = c("B","C","B"),
#'                  V3 = c("C","A","C"))
#' d2 <- data.frame(V1 = c("C","A","C"), V2 = c("A","B","A"),
#'                  V3 = c("B","C","B"))
#' net1 <- build_network(d1, method = "relative")
#' net2 <- build_network(d2, method = "relative")
#' perm <- permutation(net1, net2, iter = 20, seed = 1)
#' print(perm)
#' }
#'
#' @export
print.net_permutation <- function(x, ...) {
  method_labels <- c(
    relative      = "Transition Network (relative probabilities)",
    frequency     = "Transition Network (frequency counts)",
    co_occurrence = "Co-occurrence Network",
    glasso        = "Partial Correlation Network (EBICglasso)",
    pcor          = "Partial Correlation Network (unregularised)",
    cor           = "Correlation Network",
    attention     = "Attention Network (decay-weighted transitions)",
    wtna          = "Window TNA (transitions)",
    edge_betweenness = "Edge-Betweenness Network"
  )
  label <- if (x$method %in% names(method_labels)) {
    method_labels[[x$method]]
  } else {
    sprintf("Network (method: %s)", x$method)
  }

  dir_label <- if (x$x$directed) " [directed]" else " [undirected]"

  cat("Permutation Test: ", label, dir_label, "\n", sep = "")
  cat(sprintf("  Iterations: %d  |  Alpha: %.2f",
              x$iter, x$alpha))
  if (x$paired) cat("  |  Paired")
  if (!is.null(x$actor)) {
    cat(sprintf("  |  Actor: %s (%d actors)", x$actor, x$n_actors))
  }
  if (x$adjust != "none") cat(sprintf("  |  Adjust: %s", x$adjust))
  cat("\n")

  n_sig <- sum(x$summary$sig)
  n_total <- nrow(x$summary)
  cat(sprintf("  Nodes: %d  |  Edges tested: %d  |  Significant: %d\n",
              x$x$n_nodes, n_total, n_sig))

  if (!is.null(x$global)) {
    g <- x$global
    cat(sprintf(paste0("  Global test (networks differ overall?): ",
                       "M = %s (p = %s)  |  S = %s (p = %s)\n"),
                format(round(g$observed[1], 3), nsmall = 3),
                format(g$p_value[1], digits = 3),
                format(round(g$observed[2], 3), nsmall = 3),
                format(g$p_value[2], digits = 3)))
  }
  if (!is.null(x$clustering)) {
    cl <- x$clustering
    cat(sprintf("  Nesting in %s: %s  |  %s design\n", x$actor,
                .format_icc(cl$icc, cl$icc_ci_lower, cl$icc_ci_upper),
                cl$design))
    cat(sprintf(paste0("  Design effect (1 = nesting does not matter): ",
                       "edges %.2f  |  global %.2f\n"),
                cl$deff_edges, cl$deff_global))
    if (cl$min_p > x$alpha) {
      cat(sprintf("  Note: with %d actors no p-value can fall below %s\n",
                  cl$n_actors, format(cl$min_p, digits = 3)))
    }
  }

  invisible(x)
}


#' Summary Method for net_permutation
#'
#' @param object A \code{net_permutation} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The \code{$summary} data frame: one row per edge present in
#'   either network, with columns \code{from}, \code{to},
#'   \code{weight_x}, \code{weight_y}, \code{diff}, \code{effect_size},
#'   \code{p_value}, \code{sig}.
#'
#' @examples
#' s1 <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
#' s2 <- data.frame(V1 = c("A","C","B"), V2 = c("C","B","A"))
#' n1 <- build_network(s1, method = "relative")
#' n2 <- build_network(s2, method = "relative")
#' perm <- permutation(n1, n2, iter = 10)
#' summary(perm)
#' \donttest{
#' set.seed(1)
#' d1 <- data.frame(V1 = c("A","B","A"), V2 = c("B","C","B"),
#'                  V3 = c("C","A","C"))
#' d2 <- data.frame(V1 = c("C","A","C"), V2 = c("A","B","A"),
#'                  V3 = c("B","C","B"))
#' net1 <- build_network(d1, method = "relative")
#' net2 <- build_network(d2, method = "relative")
#' perm <- permutation(net1, net2, iter = 20, seed = 1)
#' summary(perm)
#' }
#'
#' @export
summary.net_permutation <- function(object, ...) {
  object$summary
}


#' Print Method for net_permutation_group
#'
#' @param x A \code{net_permutation_group} object.
#' @param ... Additional arguments (ignored).
#' @return \code{x} invisibly.
#' @examples
#' s1 <- data.frame(V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
#'   V3 = c("C","A","C","B"), grp = c("X","X","Y","Y"))
#' s2 <- data.frame(V1 = c("C","A","C","B"), V2 = c("A","B","A","C"),
#'   V3 = c("B","C","B","A"), grp = c("X","X","Y","Y"))
#' nets1 <- build_network(s1, method = "relative", group = "grp")
#' nets2 <- build_network(s2, method = "relative", group = "grp")
#' perm  <- permutation(nets1, nets2, iter = 10)
#' print(perm)
#' \donttest{
#' set.seed(1)
#' s1 <- data.frame(V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
#'                  V3 = c("C","A","C","B"), grp = c("X","X","Y","Y"))
#' s2 <- data.frame(V1 = c("C","A","C","B"), V2 = c("A","B","A","C"),
#'                  V3 = c("B","C","B","A"), grp = c("X","X","Y","Y"))
#' nets1 <- build_network(s1, method = "relative", group = "grp")
#' nets2 <- build_network(s2, method = "relative", group = "grp")
#' perm  <- permutation(nets1, nets2, iter = 20, seed = 1)
#' print(perm)
#' }
#' @export
print.net_permutation_group <- function(x, ...) {
  cat("Grouped Permutation Test\n")
  cat("Groups:", paste(names(x), collapse = ", "), "\n")
  lapply(names(x), function(nm) {
    cat("\n-- ", nm, " --\n", sep = "")
    print(x[[nm]], ...)
  })
  invisible(x)
}

#' Summary Method for net_permutation_group
#'
#' Returns a combined summary data frame across all groups.
#'
#' @param object A \code{net_permutation_group} object.
#' @param ... Additional arguments (ignored).
#' @return The per-group summaries stacked into one data frame: the
#'   columns of \code{\link{summary.net_permutation}} prefixed by a
#'   \code{group} column naming the group (or group pair) each row came
#'   from.
#' @examples
#' s1 <- data.frame(V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
#'   V3 = c("C","A","C","B"), grp = c("X","X","Y","Y"))
#' s2 <- data.frame(V1 = c("C","A","C","B"), V2 = c("A","B","A","C"),
#'   V3 = c("B","C","B","A"), grp = c("X","X","Y","Y"))
#' nets1 <- build_network(s1, method = "relative", group = "grp")
#' nets2 <- build_network(s2, method = "relative", group = "grp")
#' perm  <- permutation(nets1, nets2, iter = 10)
#' summary(perm)
#' \donttest{
#' set.seed(1)
#' s1 <- data.frame(V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
#'                  V3 = c("C","A","C","B"), grp = c("X","X","Y","Y"))
#' s2 <- data.frame(V1 = c("C","A","C","B"), V2 = c("A","B","A","C"),
#'                  V3 = c("B","C","B","A"), grp = c("X","X","Y","Y"))
#' nets1 <- build_network(s1, method = "relative", group = "grp")
#' nets2 <- build_network(s2, method = "relative", group = "grp")
#' perm  <- permutation(nets1, nets2, iter = 20, seed = 1)
#' summary(perm)
#' }
#' @export
summary.net_permutation_group <- function(object, ...) {
  do.call(rbind, lapply(names(object), function(nm) {
    df      <- object[[nm]]$summary
    df$group <- nm
    df[c("group", setdiff(names(df), "group"))]
  }))
}


#' Print Method for wtna_perm_mixed
#'
#' @param x A \code{wtna_perm_mixed} object.
#' @param ... Additional arguments (ignored).
#' @return The input object, invisibly.
#' @export
print.wtna_perm_mixed <- function(x, ...) {
  cat("Mixed WTNA Permutation Test (transition + co-occurrence)\n")
  cat("-- Transition (directed) --\n")
  print(x$transition)
  cat("-- Co-occurrence (undirected) --\n")
  print(x$cooccurrence)
  invisible(x)
}


#' Summary Method for wtna_perm_mixed
#'
#' @param object A \code{wtna_perm_mixed} object.
#' @param ... Additional arguments (ignored).
#' @return A list with transition and co-occurrence permutation summaries.
#' @export
summary.wtna_perm_mixed <- function(object, ...) {
  list(
    transition = summary(object$transition),
    cooccurrence = summary(object$cooccurrence)
  )
}
