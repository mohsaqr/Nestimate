# ---- Bayesian Dirichlet-Multinomial Network Comparison ----
#
# A Bayesian complement to permutation() for differential analysis of
# transition networks. Where the permutation test returns a binary
# significance decision and a point estimate, this models the outgoing
# transitions from each state as a Dirichlet-Multinomial process and
# returns, for every edge, a posterior mean difference and a credible
# interval that quantify both the magnitude and the precision of the
# difference directly from the transition counts.
#
# Reference: Johnston, L. & Jendoubi, T. (2026). How Delivery Mode Reshapes
# Resource Engagement: A Bayesian Differential Network Analysis. TNA Workshop
# 2026. Supporting theory: Gelman et al. (2013, Bayesian Data Analysis).
# The Jeffreys prior (alpha = 0.5) is non-informative and standard for
# categorical data (Jeffreys 1946).

#' Bayesian Dirichlet-Multinomial comparison of two transition networks
#'
#' @description
#' `bayes_compare()` compares two transition networks edge by edge and
#' returns, for every transition, the posterior mean difference in
#' transition probability, a credible interval, and the probability that
#' the difference lies on its side of zero. The transitions out of each
#' state are modelled as a multinomial draw with a Dirichlet prior on the
#' transition probabilities, so the posterior for the transitions out of
#' state \eqn{i} is \eqn{\mathrm{Dirichlet}(c_i + \alpha)}{Dirichlet(c_i +
#' alpha)}, where \eqn{c_i} are the observed outgoing counts and
#' \eqn{\alpha}{alpha} is `prior`. Each edge probability is then marginally
#' Beta-distributed: the posterior mean difference has a closed form, and
#' the credible interval comes from Monte Carlo draws. It takes two networks
#' estimated by [build_network()] with method `"relative"` or `"frequency"`
#' (for either, the difference is in transition probabilities estimated from
#' the counts), a group of networks, or two edge-betweenness networks from
#' \code{\link[=net_edge_betweenness]{net_edge_betweenness()}} (see
#' *Details*). An `mcml` object raises an error of class
#' `nestimate_bad_input`: its layers hold different states and carry no
#' transition counts.
#'
#' The comparison complements [permutation()]: the permutation test asks
#' whether a difference is more extreme than chance, and this comparison
#' gives the plausible range of the difference and how precisely the counts
#' estimate it. An edge whose source state has few outgoing transitions gets
#' a wide interval even when its row-normalised probabilities look far
#' apart. An edge is flagged as credibly different when its interval
#' excludes zero and its mean and nearest bound pass `mean_threshold` and
#' `bound_threshold`; the thresholds are conventions, and no correction is
#' made for comparing many edges at once. The reported `p_value` (the
#' Bayesian p-equivalent \eqn{2(1 - p_d)}{2 * (1 - pd)}, with \eqn{p_d}{pd}
#' the probability of direction) summarises posterior mass and is not a
#' frequentist p-value. The model treats every transition as an independent
#' draw, so it does not account for transitions being nested in sequences,
#' which a permutation of sequences with [permutation()] respects.
#'
#' @details
#' For two `net_edge_betweenness` objects (source method `"relative"`
#' only), edge betweenness is a nonlinear function of the whole transition
#' matrix, so Beta marginals do not apply. Instead the full transition
#' matrix is drawn from each network's row-wise Dirichlet posterior and edge
#' betweenness is recomputed on every draw, the Bayesian analogue of
#' [permutation()] on edge-betweenness networks. The result summarises the
#' posterior of \eqn{EB(x) - EB(y)}: `diff` is the posterior mean
#' difference, `prob_x` and `prob_y` the posterior mean betweenness
#' matrices, and `observed_diff` the plug-in difference of the two input
#' networks. Both inputs must use the same `invert` setting, and
#' `mean_threshold` and `bound_threshold` are then on the betweenness
#' scale.
#'
#' @param x The first network: a `netobject`
#'   estimated with `"relative"` or `"frequency"` or one of their aliases, a
#'   `netobject_group`, or a `net_edge_betweenness` object. A bare
#'   `cograph_network` is recounted from the sequences it carries, and used
#'   only when they reproduce its weights; otherwise, and for an `mcml`, an
#'   error of class `nestimate_bad_input` is raised. Both networks must use
#'   the same method and hold the same nodes; any other method raises an
#'   error, because the model needs transition counts. For the `print()` and `plot()` methods:
#'   a `net_bayes`, or for `print()` a `net_bayes_group`.
#' @param y The second network, of the same kind as `x`. Default `NULL`:
#'   with a `netobject_group` as `x`, every pair of its groups is compared.
#'   With two `netobject_group`s, the groups sharing a name are compared,
#'   and an error is raised when none do.
#' @param prior The Dirichlet concentration added to every cell of the
#'   count matrices, a positive number. Default `0.5`, the Jeffreys prior;
#'   `1` gives a uniform (Laplace) prior. Every cell gets it, so a
#'   transition observed in neither network still has a small posterior
#'   probability.
#' @param draws The number of Monte Carlo posterior draws behind the
#'   credible intervals, at least 2. Default `10000`. For edge-betweenness
#'   networks the betweenness is recomputed on every draw, so time grows
#'   with `draws`.
#' @param ci The mass of the equal-tailed credible interval, a number
#'   between 0 and 1. Default `0.95`.
#' @param mean_threshold The absolute posterior mean difference an edge
#'   must exceed to be flagged as credibly different. Default `0.01`, on
#'   the transition-probability scale.
#' @param bound_threshold The absolute value the interval bound nearest
#'   zero must exceed for an edge to be flagged. It guards against
#'   differences that are detectable but negligibly small. Default `0.001`.
#' @param seed An integer that makes the draws reproducible; the session's
#'   own random stream is left as it was.
#'   Default `NULL`. For a group, it is set again before each comparison.
#'
#' @return A `net_bayes` object (class
#'   `c("net_bayes", "netdifference", "net_permutation")`), or, for a
#'   `netobject_group` input, a `net_bayes_group`: a list of `net_bayes`
#'   objects named by comparison (`"A vs B"`), or by group name when two
#'   groups are compared. `print()` shows the method, the prior, draws,
#'   interval and thresholds, and the numbers of nodes, of edges compared
#'   and of edges credibly different (for a group, it lists the
#'   comparisons); [summary()] returns one row per edge observed in
#'   at least one network (for a group, with a `comparison` column); and
#'   `plot()` draws the difference network. A `net_bayes` carries the same
#'   fields as a [permutation()] result, so it is accepted wherever a
#'   `net_permutation` is, and also carries a difference matrix for the
#'   `netdifference` verbs ([as_netdifference()], cograph's difference
#'   plots). Its fields:
#' \describe{
#'   \item{x, y}{The two input networks.}
#'   \item{diff}{Posterior mean difference matrix (\code{prob_x - prob_y});
#'     the analogue of the permutation test's observed difference.}
#'   \item{difference_matrix}{The same matrix as \code{diff}, for cograph's
#'     \code{netdifference} helpers.}
#'   \item{diff_sig}{The difference where \code{sig} is \code{TRUE}, else 0.}
#'   \item{p_values}{The two-sided Bayesian p-equivalent, in the field a
#'     \code{net_permutation} consumer reads as p-values (see
#'     \code{p_bayes}).}
#'   \item{effect_size}{Posterior mean difference over its posterior
#'     standard deviation (0 where that is 0).}
#'   \item{ci_lower, ci_upper}{Credible-interval bound matrices.}
#'   \item{p_difference}{Probability of direction: the share of posterior
#'     mass on the dominant side of zero, in \eqn{[0.5, 1]}.}
#'   \item{p_bayes}{The same matrix as \code{p_values}: the two-sided
#'     Bayesian p-\emph{equivalent}
#'     \eqn{2(1-\mathrm{p\_difference})}{2 * (1 - p_difference)}. It
#'     summarises posterior mass, not a frequentist tail probability, so it
#'     is not a p-value and should not be reported as one.}
#'   \item{prob_x, prob_y}{Posterior mean transition-probability matrices
#'     (posterior mean betweenness for edge-betweenness input).}
#'   \item{sig}{Logical matrix: the interval excludes zero and the mean
#'     and nearest bound exceed their thresholds.}
#'   \item{summary}{The per-edge table that \code{summary()} returns, with
#'     the columns of a \code{net_permutation} summary (\code{from},
#'     \code{to}, \code{weight_x}, \code{weight_y}, \code{diff},
#'     \code{effect_size}, \code{p_value}, \code{sig}) plus
#'     \code{count_x}, \code{count_y}, \code{ci_lower}, \code{ci_upper},
#'     \code{ci_width} and \code{p_difference}.}
#'   \item{method, iter, alpha, paired, adjust}{Permutation-compatible
#'     settings (\code{iter = draws}, \code{alpha = 1 - ci},
#'     \code{paired = FALSE}, \code{adjust = "none"}); \code{method} is
#'     \code{"edge_betweenness"} for edge-betweenness input.}
#'   \item{prior, draws, ci, mean_threshold, bound_threshold}{The Bayesian
#'     settings.}
#'   \item{observed_diff, source_method, edge_betweenness}{For
#'     edge-betweenness input only: the plug-in difference, the source
#'     networks' method and the betweenness settings.}
#' }
#'
#' @examples
#' s1 <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
#' s2 <- data.frame(V1 = c("A","C","B"), V2 = c("C","B","A"))
#' n1 <- build_network(s1, method = "relative")
#' n2 <- build_network(s2, method = "relative")
#' bayes_compare(n1, n2, draws = 500, seed = 1)
#'
#' @references
#' Johnston, L. & Jendoubi, T. (2026). How Delivery Mode Reshapes Resource
#' Engagement: A Bayesian Differential Network Analysis. TNA Workshop 2026.
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., &
#' Rubin, D. B. (2013). \emph{Bayesian Data Analysis} (3rd ed.). CRC Press.
#'
#' Jeffreys, H. (1946). An invariant form for the prior probability in
#' estimation problems. \emph{Proceedings of the Royal Society of London A},
#' 186(1007), 453-461.
#'
#' @seealso \code{\link{permutation}} for the frequentist complement;
#'   \code{\link{certainty}} for single-network posterior edge intervals;
#'   \code{\link{subtract_networks}} and \code{\link{as_netdifference}} for the
#'   difference verbs; \code{\link{build_network}}
#' @importFrom stats rbeta quantile
#' @export
bayes_compare <- function(x, y = NULL,
                          prior = 0.5,
                          draws = 10000L,
                          ci = 0.95,
                          mean_threshold = 0.01,
                          bound_threshold = 0.001,
                          seed = NULL) {

  # ---- mcml: refused. Its layers (macro + one per cluster) hold different
  # states and carry no transition counts, so no comparison can run. ----
  if (inherits(x, "mcml") || inherits(y, "mcml")) {
    .stop_mcml_input("bayes_compare()", "compares two transition networks")
  }

  # ---- Single netobject_group: all-pairs comparisons ----
  if (inherits(x, "netobject_group") && is.null(y)) {
    grp_names <- names(x)
    if (length(grp_names) < 2L) {
      stop(errorCondition(
        "Need at least 2 groups for pairwise Bayesian comparison.",
        class = "nestimate_bad_input", call = NULL))
    }
    pairs <- utils::combn(length(grp_names), 2L)
    results <- lapply(seq_len(ncol(pairs)), function(k) {
      bayes_compare(x[[pairs[1L, k]]], x[[pairs[2L, k]]], prior = prior,
                    draws = draws, ci = ci, mean_threshold = mean_threshold,
                    bound_threshold = bound_threshold, seed = seed)
    })
    names(results) <- vapply(seq_len(ncol(pairs)), function(k) {
      paste(grp_names[pairs[1L, k]], "vs", grp_names[pairs[2L, k]])
    }, character(1))
    class(results) <- c("net_bayes_group", "list")
    return(results)
  }

  # ---- Two netobject_groups: compare each matching element ----
  if (inherits(x, "netobject_group") && inherits(y, "netobject_group")) {
    common <- intersect(names(x), names(y))
    if (length(common) == 0L) {
      stop(errorCondition("No matching group names between x and y.",
        class = "nestimate_bad_input", call = NULL))
    }
    results <- lapply(common, function(nm) {
      bayes_compare(x[[nm]], y[[nm]], prior = prior, draws = draws, ci = ci,
                    mean_threshold = mean_threshold,
                    bound_threshold = bound_threshold, seed = seed)
    })
    names(results) <- common
    class(results) <- c("net_bayes_group", "list")
    return(results)
  }

  # ---- Coerce cograph_network inputs ----
  if (inherits(x, "cograph_network")) x <- .as_netobject(x)
  if (inherits(y, "cograph_network")) y <- .as_netobject(y)

  # ---- Edge-betweenness dispatch: posterior draws of the source networks ----
  if (inherits(x, "net_edge_betweenness") ||
      inherits(y, "net_edge_betweenness")) {
    return(.bayes_compare_edge_betweenness(
      x = x, y = y, prior = prior, draws = draws, ci = ci,
      mean_threshold = mean_threshold, bound_threshold = bound_threshold,
      seed = seed
    ))
  }

  # ---- Input validation ----
  .assert(
    "`x` and `y` must be netobjects (from build_network())" =
      inherits(x, "netobject") && inherits(y, "netobject"),
    "`prior` must be a single positive number" =
      is.numeric(prior) && length(prior) == 1L && prior > 0,
    "`draws` must be a single number >= 2" =
      is.numeric(draws) && length(draws) == 1L && draws >= 2,
    "`ci` must be a single number in (0, 1)" =
      is.numeric(ci) && length(ci) == 1L && ci > 0 && ci < 1,
    "`mean_threshold` must be a single number" =
      is.numeric(mean_threshold) && length(mean_threshold) == 1L,
    "`bound_threshold` must be a single number" =
      is.numeric(bound_threshold) && length(bound_threshold) == 1L
  )
  draws <- as.integer(draws)

  method_x <- .resolve_method_alias(x$method)
  method_y <- .resolve_method_alias(y$method)
  if (method_x != method_y) {
    stop(errorCondition(paste0("Methods must match: x uses '", x$method,
         "', y uses '", y$method, "'."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (!method_x %in% c("relative", "frequency")) {
    stop(errorCondition(
      paste0(
        "bayes_compare() requires a transition-count method ",
        "(\"relative\" or \"frequency\"); got '",
        x$method,
        "'. ",
        "The Dirichlet-Multinomial model needs integer transition counts."),
      class = "nestimate_bad_input", call = NULL))
  }
  x$frequency_matrix <- .network_transition_counts(x, "bayes_compare()")
  y$frequency_matrix <- .network_transition_counts(y, "bayes_compare()")
  if (!setequal(x$nodes$label, y$nodes$label)) {
    stop(errorCondition("Nodes must be the same in both networks.",
      class = "nestimate_bad_input", call = NULL))
  }

  if (!is.null(seed)) {
    .assert("`seed` must be NULL or a single number" =
              is.numeric(seed) && length(seed) == 1L)
    .local_seed(seed)
  }

  # ---- Align count matrices on a common node order ----
  nodes <- x$nodes$label
  count_x <- x$frequency_matrix[nodes, nodes, drop = FALSE]
  count_y <- y$frequency_matrix[nodes, nodes, drop = FALSE]
  k <- length(nodes)

  # ---- Posterior parameters (Dirichlet per source row) ----
  alpha_x <- count_x + prior
  alpha_y <- count_y + prior
  rowA_x <- matrix(rowSums(alpha_x), k, k)   # row-constant concentration sums
  rowA_y <- matrix(rowSums(alpha_y), k, k)

  # Posterior mean transition probabilities (closed form)
  prob_x <- alpha_x / rowA_x
  prob_y <- alpha_y / rowA_y
  diff_mat <- prob_x - prob_y

  # ---- Credible interval via Monte Carlo on Beta marginals ----
  # Each edge (i, j): p ~ Beta(alpha_ij, rowSum_i - alpha_ij)
  res <- .bayes_beta_ci(alpha_x, rowA_x - alpha_x,
                        alpha_y, rowA_y - alpha_y,
                        draws = draws, ci = ci, k = k)
  ci_lower <- res$lower
  ci_upper <- res$upper
  pd       <- res$pd          # probability of direction in [0.5, 1]
  p_bayes  <- 2 * (1 - pd)    # two-sided Bayesian p-equivalent

  # ---- Significance: CI excludes zero, magnitude + nearest bound large ----
  ci_excl <- (ci_lower > 0) | (ci_upper < 0)
  nearest <- pmin(abs(ci_lower), abs(ci_upper))
  sig <- ci_excl &
    (abs(diff_mat) > mean_threshold) &
    (nearest > bound_threshold)

  # ---- Permutation-format fields (so net_bayes is a net_permutation drop-in) ----
  # Cohen's-d-style effect size: posterior mean diff over its posterior sd.
  perm_sd <- res$sd
  perm_sd[perm_sd == 0] <- NA_real_
  effect_size <- diff_mat / perm_sd
  effect_size[is.na(effect_size)] <- 0
  diff_sig <- diff_mat * sig

  # ---- Long-format summary (permutation columns + Bayesian extras) ----
  summary_df <- .build_bayes_summary(
    diff_mat = diff_mat, ci_lower = ci_lower, ci_upper = ci_upper,
    prob_x = prob_x, prob_y = prob_y, count_x = count_x, count_y = count_y,
    effect_size = effect_size, pd = pd, p_bayes = p_bayes,
    sig = sig, nodes = nodes
  )

  result <- list(
    x               = x,
    y               = y,
    diff            = diff_mat,
    difference_matrix = diff_mat,
    diff_sig        = diff_sig,
    p_values        = p_bayes,
    effect_size     = effect_size,
    ci_lower        = ci_lower,
    ci_upper        = ci_upper,
    p_difference    = pd,
    p_bayes         = p_bayes,
    prob_x          = prob_x,
    prob_y          = prob_y,
    sig             = sig,
    summary         = summary_df,
    method          = method_x,
    iter            = draws,
    alpha           = 1 - ci,
    paired          = FALSE,
    adjust          = "none",
    prior           = prior,
    draws           = draws,
    ci              = ci,
    mean_threshold  = mean_threshold,
    bound_threshold = bound_threshold,
    meta            = list(splot = list(renderer = "permutation"))
  )
  class(result) <- c("net_bayes", "netdifference", "net_permutation")
  result
}


#' Bayesian comparison of two edge-betweenness networks
#'
#' Edge betweenness is a nonlinear function of the whole transition matrix,
#' so per-edge Beta marginals do not apply. Instead: draw full transition
#' matrices from each group's row-wise Dirichlet posterior, recompute edge
#' betweenness on every draw (the Bayesian analogue of the permutation
#' dispatch, which recomputes it on every shuffle), and summarise the
#' posterior of \code{EB(x) - EB(y)} per edge.
#' @noRd
.bayes_compare_edge_betweenness <- function(x, y, prior, draws, ci,
                                            mean_threshold, bound_threshold,
                                            seed) {
  if (!inherits(x, "net_edge_betweenness") ||
      !inherits(y, "net_edge_betweenness")) {
    stop(errorCondition("Both x and y must be net_edge_betweenness objects.",
      class = "nestimate_bad_input", call = NULL))
  }
  .assert(
    "`prior` must be a single positive number" =
      is.numeric(prior) && length(prior) == 1L && prior > 0,
    "`draws` must be a single number >= 2" =
      is.numeric(draws) && length(draws) == 1L && draws >= 2,
    "`ci` must be a single number in (0, 1)" =
      is.numeric(ci) && length(ci) == 1L && ci > 0 && ci < 1,
    "`mean_threshold` must be a single number" =
      is.numeric(mean_threshold) && length(mean_threshold) == 1L,
    "`bound_threshold` must be a single number" =
      is.numeric(bound_threshold) && length(bound_threshold) == 1L
  )
  draws <- as.integer(draws)

  x_src <- .edge_betweenness_source_net(x)
  y_src <- .edge_betweenness_source_net(y)
  method <- .resolve_method_alias(x_src$method)
  y_method <- .resolve_method_alias(y_src$method)
  if (method != y_method) {
    stop(errorCondition(paste0("Source methods must match: x uses '", method,
         "', y uses '", y_method, "'."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (method != "relative") {
    stop(errorCondition(
      paste0(
        "Bayesian edge-betweenness comparison requires source method ",
        "\"relative\": the Dirichlet posterior models the transition ",
        "probability rows that the betweenness is computed on; got '",
        method,
        "'."),
      class = "nestimate_bad_input", call = NULL))
  }
  invert <- isTRUE(x$edge_betweenness$invert)
  if (!identical(isTRUE(y$edge_betweenness$invert), invert)) {
    stop(errorCondition(
      "x and y must use the same edge-betweenness `invert` setting.",
      class = "nestimate_bad_input", call = NULL))
  }
  if (is.null(x$frequency_matrix) || is.null(y$frequency_matrix)) {
    stop(errorCondition(
      paste0(
        "Transition counts ($frequency_matrix) not found. ",
        "Recreate with net_edge_betweenness(build_network(...))."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (!setequal(x$nodes$label, y$nodes$label)) {
    stop(errorCondition("Nodes must be the same in both networks.",
      class = "nestimate_bad_input", call = NULL))
  }

  if (!is.null(seed)) {
    .assert("`seed` must be NULL or a single number" =
              is.numeric(seed) && length(seed) == 1L)
    .local_seed(seed)
  }

  nodes <- x$nodes$label
  count_x <- x$frequency_matrix[nodes, nodes, drop = FALSE]
  count_y <- y$frequency_matrix[nodes, nodes, drop = FALSE]
  k <- length(nodes)
  nbins <- k * k
  dn <- list(nodes, nodes)

  alpha_x <- count_x + prior
  alpha_y <- count_y + prior

  # ---- Posterior draws: Dirichlet rows -> edge betweenness per draw ----
  samp <- vapply(seq_len(draws), function(d) {
    px <- .dirichlet_rows(alpha_x)
    py <- .dirichlet_rows(alpha_y)
    ex <- as.vector(.edge_betweenness(px, invert = invert))
    ey <- as.vector(.edge_betweenness(py, invert = invert))
    c(ex - ey, ex, ey)
  }, numeric(3L * nbins))
  diff_draws <- samp[seq_len(nbins), , drop = FALSE]
  eb_x_mean <- matrix(rowMeans(samp[nbins + seq_len(nbins), , drop = FALSE]),
                      k, k, dimnames = dn)
  eb_y_mean <- matrix(rowMeans(samp[2L * nbins + seq_len(nbins), , drop = FALSE]),
                      k, k, dimnames = dn)

  # ---- Posterior summaries (same rules as the transition path) ----
  diff_mat <- matrix(rowMeans(diff_draws), k, k, dimnames = dn)
  probs <- c((1 - ci) / 2, 1 - (1 - ci) / 2)
  qs <- apply(diff_draws, 1L, stats::quantile, probs = probs, names = FALSE)
  ci_lower <- matrix(qs[1L, ], k, k, dimnames = dn)
  ci_upper <- matrix(qs[2L, ], k, k, dimnames = dn)
  prop_pos <- rowMeans(diff_draws > 0)
  pd <- matrix(pmax(prop_pos, 1 - prop_pos), k, k, dimnames = dn)
  p_bayes <- 2 * (1 - pd)
  sd_d <- sqrt(pmax(rowMeans(diff_draws^2) - rowMeans(diff_draws)^2, 0))
  perm_sd <- matrix(sd_d, k, k, dimnames = dn)

  ci_excl <- (ci_lower > 0) | (ci_upper < 0)
  nearest <- pmin(abs(ci_lower), abs(ci_upper))
  sig <- ci_excl &
    (abs(diff_mat) > mean_threshold) &
    (nearest > bound_threshold)

  perm_sd[perm_sd == 0] <- NA_real_
  effect_size <- diff_mat / perm_sd
  effect_size[is.na(effect_size)] <- 0
  diff_sig <- diff_mat * sig

  summary_df <- .build_bayes_summary(
    diff_mat = diff_mat, ci_lower = ci_lower, ci_upper = ci_upper,
    prob_x = eb_x_mean, prob_y = eb_y_mean,
    count_x = count_x, count_y = count_y,
    effect_size = effect_size, pd = pd, p_bayes = p_bayes,
    sig = sig, nodes = nodes
  )

  result <- list(
    x               = x,
    y               = y,
    diff            = diff_mat,
    difference_matrix = diff_mat,
    diff_sig        = diff_sig,
    p_values        = p_bayes,
    effect_size     = effect_size,
    ci_lower        = ci_lower,
    ci_upper        = ci_upper,
    p_difference    = pd,
    p_bayes         = p_bayes,
    prob_x          = eb_x_mean,
    prob_y          = eb_y_mean,
    observed_diff   = x$weights - y$weights[nodes, nodes],
    sig             = sig,
    summary         = summary_df,
    method          = "edge_betweenness",
    source_method   = method,
    edge_betweenness = list(invert = invert, source_method = method),
    iter            = draws,
    alpha           = 1 - ci,
    paired          = FALSE,
    adjust          = "none",
    prior           = prior,
    draws           = draws,
    ci              = ci,
    mean_threshold  = mean_threshold,
    bound_threshold = bound_threshold,
    meta            = list(splot = list(renderer = "permutation"))
  )
  class(result) <- c("net_bayes", "netdifference", "net_permutation")
  result
}


#' Row-wise Dirichlet draw of a transition matrix
#'
#' One sample of the full transition matrix: each row i is drawn from
#' Dirichlet(\code{alpha[i, ]}) via normalised Gamma variates.
#' @noRd
.dirichlet_rows <- function(alpha) {
  g <- matrix(stats::rgamma(length(alpha), shape = alpha),
              nrow(alpha), ncol(alpha), dimnames = dimnames(alpha))
  g / rowSums(g)
}


#' Monte Carlo credible interval for differences of two Beta fields
#'
#' Draws \code{draws} samples per edge from each group's Beta marginal and
#' returns per-edge quantiles of the difference plus the probability of
#' direction. Vectorised: a single \code{rbeta} call fills an
#' (edges x draws) matrix per group.
#' @noRd
.bayes_beta_ci <- function(s1_x, s2_x, s1_y, s2_y, draws, ci, k) {
  nbins <- k * k
  a_x <- as.vector(s1_x); b_x <- as.vector(s2_x)
  a_y <- as.vector(s1_y); b_y <- as.vector(s2_y)

  samp_x <- matrix(stats::rbeta(nbins * draws, rep(a_x, draws), rep(b_x, draws)),
                   nrow = nbins, ncol = draws)
  samp_y <- matrix(stats::rbeta(nbins * draws, rep(a_y, draws), rep(b_y, draws)),
                   nrow = nbins, ncol = draws)
  diff_draws <- samp_x - samp_y

  probs <- c((1 - ci) / 2, 1 - (1 - ci) / 2)
  qs <- apply(diff_draws, 1L, stats::quantile, probs = probs, names = FALSE)

  # Probability of direction: share of posterior mass on the dominant side.
  prop_pos <- rowMeans(diff_draws > 0)
  pd <- pmax(prop_pos, 1 - prop_pos)

  # Posterior sd of the difference (for the permutation-style effect size).
  mean_d <- rowMeans(diff_draws)
  sd_d <- sqrt(pmax(rowMeans(diff_draws^2) - mean_d^2, 0))

  dn <- dimnames(s1_x)
  list(
    lower = matrix(qs[1L, ], k, k, dimnames = dn),
    upper = matrix(qs[2L, ], k, k, dimnames = dn),
    pd    = matrix(pd, k, k, dimnames = dn),
    sd    = matrix(sd_d, k, k, dimnames = dn)
  )
}


#' Build long-format summary for net_bayes
#'
#' Columns are a superset of \code{summary.net_permutation}: it carries
#' \code{from, to, weight_x, weight_y, diff, effect_size, p_value, sig}
#' (the permutation columns, so a net_bayes summary is a drop-in) plus the
#' Bayesian extras \code{count_x, count_y, ci_lower, ci_upper, ci_width,
#' p_difference}.
#' @noRd
.build_bayes_summary <- function(diff_mat, ci_lower, ci_upper,
                                 prob_x, prob_y, count_x, count_y,
                                 effect_size, pd, p_bayes, sig, nodes) {
  n <- length(nodes)
  dt <- data.table::data.table(
    from        = rep(nodes, each = n),
    to          = rep(nodes, times = n),
    weight_x    = as.vector(t(prob_x)),
    weight_y    = as.vector(t(prob_y)),
    count_x     = as.vector(t(count_x)),
    count_y     = as.vector(t(count_y)),
    diff        = as.vector(t(diff_mat)),
    effect_size = as.vector(t(effect_size)),
    ci_lower    = as.vector(t(ci_lower)),
    ci_upper    = as.vector(t(ci_upper)),
    ci_width    = as.vector(t(ci_upper - ci_lower)),
    p_difference = as.vector(t(pd)),
    p_value     = as.vector(t(p_bayes)),
    sig         = as.vector(t(sig))
  )
  # Keep directed edges observed in at least one network
  dt <- dt[count_x > 0 | count_y > 0]
  as.data.frame(dt)
}


# ---- S3 Methods ----

#' @rdname bayes_compare
#' @param ... For `plot()`: further arguments for `cograph::splot()` when
#'   cograph is installed (such as `weight_digits`, the digits the drawn
#'   differences are rounded to, default 2, or `edge_label_digits`, the
#'   digits of the edge labels, default 3); otherwise
#'   ignored. The `print()` and `summary()` methods ignore it.
#' @param object For `summary()`: a `net_bayes` or `net_bayes_group`.
#' @param significant_only For `plot()`: whether to draw only the edges
#'   flagged as credibly different. Default `TRUE`. When no edge is left to
#'   draw, an error suggests `significant_only = FALSE`.
#' @param title For `plot()`: the plot title. Default `NULL`, which gives
#'   `"Bayesian network difference"` (with cograph) or a title naming the
#'   posterior mean difference (without it).
#' @return For the `print()` methods: the input, invisibly. For `summary()`:
#'   the per-edge data frame described under the `summary` field; for a
#'   `net_bayes_group`, the tables of every comparison stacked, with a
#'   `comparison` column first. For `plot()`: with cograph installed, it
#'   draws the difference network with `cograph::splot()`, each edge
#'   labelled with its posterior mean difference and credible interval, and
#'   returns the `cograph_network` invisibly; without cograph, it returns a
#'   `ggplot` chord diagram in which colour gives the direction of the
#'   difference and width its size.
#' @export
print.net_bayes <- function(x, ...) {
  method_labels <- c(
    relative  = "Transition Network (relative probabilities)",
    frequency = "Transition Network (frequency counts)",
    edge_betweenness = "Edge-Betweenness Network"
  )
  label <- if (x$method %in% names(method_labels)) {
    method_labels[[x$method]]
  } else {
    sprintf("Network (method: %s)", x$method)
  }
  cat("Bayesian Dirichlet-Multinomial Comparison:", label, "\n")
  cat(sprintf("  Prior: Dirichlet(%.2f)  |  Draws: %d  |  CI: %d%%\n",
              x$prior, x$draws, round(x$ci * 100)))
  cat(sprintf("  Thresholds: |mean diff| > %.3f, nearest CI bound > %.3f\n",
              x$mean_threshold, x$bound_threshold))
  n_sig <- sum(x$summary$sig)
  n_total <- nrow(x$summary)
  cat(sprintf("  Nodes: %d  |  Edges compared: %d  |  Credibly different: %d\n",
              x$x$n_nodes, n_total, n_sig))
  invisible(x)
}


#' @rdname bayes_compare
#' @export
summary.net_bayes <- function(object, ...) {
  object$summary
}


#' @rdname bayes_compare
#' @export
plot.net_bayes <- function(x, significant_only = TRUE, title = NULL, ...) {
  if (requireNamespace("cograph", quietly = TRUE)) {
    args <- list(...)

    # Per-edge CI vectors must line up with the edge set cograph's
    # net_permutation renderer enumerates: it rounds the display matrix to
    # weight_digits first and, for undirected networks, keeps only the upper
    # triangle - replicate both so the vectors stay aligned.
    weights_display <- if (isTRUE(significant_only)) x$diff_sig else x$diff
    weights_display <- round(weights_display, args[["weight_digits"]] %||% 2)
    directed <- if (!is.null(x$x$directed)) isTRUE(x$x$directed) else TRUE
    edge_pos <- if (directed) {
      which(weights_display != 0, arr.ind = TRUE)
    } else {
      which(weights_display != 0 & upper.tri(weights_display), arr.ind = TRUE)
    }
    if (nrow(edge_pos) == 0L) {
      stop(errorCondition(paste0("No edges to plot",
           if (significant_only) " (none credibly different; try significant_only = FALSE)"
           else "", "."),
        class = "nestimate_bad_input", call = NULL))
    }
    # Exact [[ ]] indexing throughout: `$` on a list partially matches
    # (args$title finds args$title_size when title is absent), which
    # silently skipped the title default here.
    if (is.null(args[["edge_ci_lower"]])) args[["edge_ci_lower"]] <- x$ci_lower[edge_pos]
    if (is.null(args[["edge_ci_upper"]])) args[["edge_ci_upper"]] <- x$ci_upper[edge_pos]
    if (is.null(args[["edge_ci"]])) {
      args[["edge_ci"]] <- pmax(args[["edge_ci_upper"]] - args[["edge_ci_lower"]], 0)
    }
    if (is.null(args[["edge_label_template"]]) && is.null(args[["edge_label_style"]])) {
      args[["edge_label_template"]] <- "{est} {stars}\n[{low}, {up}]"
    }
    if (is.null(args[["edge_label_oneline"]])) args[["edge_label_oneline"]] <- FALSE
    if (is.null(args[["edge_label_digits"]])) args[["edge_label_digits"]] <- 3L
    if (is.null(args[["edge_label_leading_zero"]])) {
      args[["edge_label_leading_zero"]] <- FALSE
    }
    if (is.null(args[["edge_label_stars"]])) args[["edge_label_stars"]] <- TRUE
    if (is.null(args[["edge_label_size"]])) args[["edge_label_size"]] <- 0.36
    if (is.null(args[["edge_ci_alpha"]])) args[["edge_ci_alpha"]] <- 0.10
    if (is.null(args[["edge_ci_scale"]])) args[["edge_ci_scale"]] <- 1.0
    if (is.null(args[["edge_ci_style"]])) args[["edge_ci_style"]] <- 1
    if (is.null(args[["title_size"]])) args[["title_size"]] <- 0.82
    if (is.null(args[["title"]])) {
      args[["title"]] <- title %||% "Bayesian network difference"
    }
    args[["show_nonsig"]] <- !isTRUE(significant_only)

    return(do.call(cograph::splot, c(list(x = x), args)))
  }

  ed <- x$summary
  if (significant_only) ed <- ed[ed$sig, , drop = FALSE]
  if (nrow(ed) == 0L) {
    stop(errorCondition(paste0("No edges to plot",
         if (significant_only) " (none credibly different; try significant_only = FALSE)"
         else "", "."),
      class = "nestimate_bad_input", call = NULL))
  }
  ed$signed <- ed$diff
  ed$value  <- abs(ed$diff)
  if (is.null(title)) {
    title <- "Differential transition network (posterior mean difference)"
  }
  .magdiff_circular(ed, value_col = "signed", abs_col = "value", title = title)
}


#' @rdname bayes_compare
#' @export
print.net_bayes_group <- function(x, ...) {
  cat("Grouped Bayesian Dirichlet-Multinomial Comparison\n")
  cat("Comparisons:", paste(names(x), collapse = ", "), "\n")
  cat("Use summary() for combined edge-level results.\n")
  invisible(x)
}


#' @rdname bayes_compare
#' @export
summary.net_bayes_group <- function(object, ...) {
  do.call(rbind, lapply(names(object), function(nm) {
    df <- object[[nm]]$summary
    df$comparison <- nm
    df[c("comparison", setdiff(names(df), "comparison"))]
  }))
}


# ---- Shared input guards (bayes_compare(), certainty(), predict_links(),
# permutation_diagnostics()) ----

#' Refuse an mcml where one network or a pair is needed
#' @noRd
.stop_mcml_input <- function(caller, what) {
  stop(errorCondition(sprintf(paste(
    "%s %s, and an mcml holds several networks over different states (the",
    "macro network and one per cluster) without transition counts. Build",
    "the networks from the data with build_network(), using `group =` for",
    "groups."), caller, what),
    class = "nestimate_bad_input", call = NULL))
}

#' Transition counts of a relative or frequency network
#'
#' A network from build_network() carries `$frequency_matrix`. One without it
#' (a bare `cograph_network`) is recounted from its sequences by the
#' estimator's own counting (`.te_network_counts()`), and the recount is used
#' only when it reproduces the network's weights: row-normalised counts for
#' "relative", the counts themselves for "frequency".
#' @noRd
.network_transition_counts <- function(x, caller) {
  if (!is.null(x$frequency_matrix)) return(x$frequency_matrix)
  method <- .resolve_method_alias(x$method %||% "")
  has_data <- is.data.frame(x$data) && nrow(x$data) > 0L
  if (!has_data && !identical(method, "frequency")) {
    stop(errorCondition(sprintf(paste(
      "%s needs transition counts, and this network carries neither its",
      "counts nor its sequences. Pass the netobject that build_network()",
      "returned."), caller),
      class = "nestimate_bad_input", call = NULL))
  }
  counts <- .te_network_counts(x)
  states <- x$nodes$label
  counts <- counts[states, states, drop = FALSE]
  expected <- if (identical(method, "relative")) {
    rs <- rowSums(counts)
    counts / ifelse(rs == 0, 1, rs)
  } else {
    counts
  }
  reproduced <- isTRUE(all.equal(
    unname(expected), unname(x$weights[states, states, drop = FALSE]),
    tolerance = sqrt(.Machine$double.eps)))
  if (!reproduced) {
    stop(errorCondition(sprintf(paste(
      "%s needs transition counts, and the sequences this network carries",
      "do not reproduce its weights (a threshold, scaling or other setting",
      "was applied). Pass the netobject that build_network() returned."),
      caller),
      class = "nestimate_bad_input", call = NULL))
  }
  counts
}
