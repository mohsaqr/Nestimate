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
#' Compares two transition networks estimated by \code{\link{build_network}}
#' (method \code{"relative"} or \code{"frequency"}) using a Bayesian
#' Dirichlet-Multinomial model. The outgoing transitions from each source
#' state are modelled as a Multinomial draw with a Dirichlet prior on the
#' transition probabilities. With a Jeffreys prior the posterior for the
#' transitions out of state \eqn{i} is \eqn{\mathrm{Dirichlet}(c_i + \alpha)},
#' where \eqn{c_i} are the observed outgoing counts. Each edge probability is
#' then marginally Beta-distributed, so the posterior mean difference between
#' the two networks is available in closed form and a credible interval is
#' obtained by Monte Carlo.
#'
#' This is a complement to \code{\link{permutation}}: the permutation test
#' answers "is this difference more extreme than chance?"; the Bayesian
#' comparison answers "what is the plausible range of the true difference,
#' and how precisely is it estimated given the counts?". An edge with few
#' outgoing transitions from its source state yields a wide credible
#' interval even when its row-normalised probability looks decisive.
#'
#' @param x A \code{netobject} (from \code{\link{build_network}}), a
#'   \code{netobject_group}, or an \code{mcml} object. Must use a transition
#'   method (\code{"relative"} / \code{"frequency"} and their aliases).
#' @param y A second object of the same kind as \code{x}, or \code{NULL}.
#'   When \code{x} is a \code{netobject_group} and \code{y} is \code{NULL},
#'   all pairwise comparisons among the groups are returned.
#' @param prior Numeric. Dirichlet prior concentration added to every cell
#'   (default \code{0.5}, the Jeffreys prior). Use \code{1} for a uniform
#'   (Laplace) prior.
#' @param draws Integer. Number of Monte Carlo posterior draws used for the
#'   credible intervals (default \code{10000}).
#' @param ci Numeric in (0, 1). Credible interval mass (default \code{0.95}).
#' @param mean_threshold Numeric. An edge is flagged significant only if the
#'   absolute posterior mean difference exceeds this (default \code{0.01}).
#' @param bound_threshold Numeric. An edge is flagged significant only if the
#'   credible-interval bound nearest zero exceeds this in absolute value
#'   (default \code{0.001}). Guards against differences that are detectable
#'   but negligibly small.
#' @param seed Integer or NULL. RNG seed for reproducible credible intervals.
#'
#' @return An object of class \code{c("net_bayes", "net_permutation")}. It carries
#'   the same fields as a \code{\link{permutation}} result, so it is a drop-in
#'   wherever a \code{net_permutation} is consumed, plus Bayesian extras:
#' \describe{
#'   \item{x, y}{The two input \code{netobject}s.}
#'   \item{diff}{Posterior mean difference matrix (\code{prob_x - prob_y});
#'     the analogue of the permutation observed difference.}
#'   \item{diff_sig}{Difference where \code{sig}, else 0.}
#'   \item{p_values}{P-value matrix (the two-sided Bayesian p-equivalent).}
#'   \item{effect_size}{Posterior mean difference over its posterior SD.}
#'   \item{ci_lower, ci_upper}{Credible-interval bound matrices.}
#'   \item{pd}{Probability-of-direction matrix in \eqn{[0.5, 1]}.}
#'   \item{p_bayes}{Alias of \code{p_values} (two-sided Bayesian p, \eqn{2(1-pd)}).}
#'   \item{prob_x, prob_y}{Posterior mean transition-probability matrices.}
#'   \item{sig}{Logical significance matrix (CI excludes zero, mean and
#'     nearest bound exceed their thresholds).}
#'   \item{summary}{Long-format data frame whose columns are a superset of
#'     \code{summary.net_permutation} (\code{from, to, weight_x, weight_y, diff,
#'     effect_size, p_value, sig}) plus \code{count_x, count_y, ci_lower,
#'     ci_upper, ci_width, pd}.}
#'   \item{method, iter, alpha, paired, adjust}{Permutation-compatible settings
#'     (\code{iter = draws}, \code{alpha = 1 - ci}, \code{paired = FALSE},
#'     \code{adjust = "none"}).}
#'   \item{prior, draws, ci, mean_threshold, bound_threshold}{Bayesian settings.}
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
#' @seealso \code{\link{permutation}}, \code{\link{build_network}}
#' @importFrom stats rbeta quantile
#' @export
bayes_compare <- function(x, y = NULL,
                          prior = 0.5,
                          draws = 10000L,
                          ci = 0.95,
                          mean_threshold = 0.01,
                          bound_threshold = 0.001,
                          seed = NULL) {

  # ---- mcml dispatch: convert to netobject_group via as_tna ----
  if (inherits(x, "mcml")) x <- as_tna(x)
  if (inherits(y, "mcml")) y <- as_tna(y)

  # ---- Single netobject_group: all-pairs comparisons ----
  if (inherits(x, "netobject_group") && is.null(y)) {
    grp_names <- names(x)
    if (length(grp_names) < 2L) {
      stop("Need at least 2 groups for pairwise Bayesian comparison.",
           call. = FALSE)
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
      stop("No matching group names between x and y.", call. = FALSE)
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

  # ---- Input validation ----
  stopifnot(
    inherits(x, "netobject"),
    inherits(y, "netobject"),
    is.numeric(prior), length(prior) == 1, prior > 0,
    is.numeric(draws), length(draws) == 1, draws >= 2,
    is.numeric(ci), length(ci) == 1, ci > 0, ci < 1,
    is.numeric(mean_threshold), length(mean_threshold) == 1,
    is.numeric(bound_threshold), length(bound_threshold) == 1
  )
  draws <- as.integer(draws)

  method_x <- .resolve_method_alias(x$method)
  method_y <- .resolve_method_alias(y$method)
  if (method_x != method_y) {
    stop("Methods must match: x uses '", x$method,
         "', y uses '", y$method, "'.", call. = FALSE)
  }
  if (!method_x %in% c("relative", "frequency")) {
    stop("bayes_compare() requires a transition-count method ",
         "(\"relative\" or \"frequency\"); got '", x$method, "'. ",
         "The Dirichlet-Multinomial model needs integer transition counts.",
         call. = FALSE)
  }
  if (is.null(x$frequency_matrix) || is.null(y$frequency_matrix)) {
    stop("Transition counts ($frequency_matrix) not found. ",
         "Rebuild with build_network(method = \"relative\").",
         call. = FALSE)
  }
  if (!setequal(x$nodes$label, y$nodes$label)) {
    stop("Nodes must be the same in both networks.", call. = FALSE)
  }

  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1)
    set.seed(seed)
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
    diff_sig        = diff_sig,
    p_values        = p_bayes,
    effect_size     = effect_size,
    ci_lower        = ci_lower,
    ci_upper        = ci_upper,
    pd              = pd,
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
    bound_threshold = bound_threshold
  )
  class(result) <- c("net_bayes", "net_permutation")
  result
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
#' Bayesian extras \code{count_x, count_y, ci_lower, ci_upper, ci_width, pd}.
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
    pd          = as.vector(t(pd)),
    p_value     = as.vector(t(p_bayes)),
    sig         = as.vector(t(sig))
  )
  # Keep directed edges observed in at least one network
  dt <- dt[count_x > 0 | count_y > 0]
  as.data.frame(dt)
}


# ---- S3 Methods ----

#' Print method for net_bayes
#'
#' @param x A \code{net_bayes} object.
#' @param ... Additional arguments (ignored).
#' @return The input object, invisibly.
#' @examples
#' s1 <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
#' s2 <- data.frame(V1 = c("A","C","B"), V2 = c("C","B","A"))
#' b <- bayes_compare(build_network(s1, method = "relative"),
#'                    build_network(s2, method = "relative"),
#'                    draws = 500, seed = 1)
#' print(b)
#' @export
print.net_bayes <- function(x, ...) {
  method_labels <- c(
    relative  = "Transition Network (relative probabilities)",
    frequency = "Transition Network (frequency counts)"
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


#' Summary method for net_bayes
#'
#' @param object A \code{net_bayes} object.
#' @param ... Additional arguments (ignored).
#' @return A data frame with edge-level posterior differences and intervals.
#' @examples
#' s1 <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
#' s2 <- data.frame(V1 = c("A","C","B"), V2 = c("C","B","A"))
#' b <- bayes_compare(build_network(s1, method = "relative"),
#'                    build_network(s2, method = "relative"),
#'                    draws = 500, seed = 1)
#' summary(b)
#' @export
summary.net_bayes <- function(object, ...) {
  object$summary
}


#' Plot method for net_bayes
#'
#' Draws a differential transition network as a directed chord diagram.
#' Edge colour encodes the signed posterior mean difference (\code{x}
#' stronger vs \code{y} stronger) and edge width its magnitude.
#'
#' @param x A \code{net_bayes} object.
#' @param significant_only Logical. Show only credibly-different edges
#'   (default \code{TRUE}).
#' @param title Optional plot title.
#' @param ... Additional arguments (ignored).
#' @return A \code{ggplot} object.
#' @examples
#' \donttest{
#' s1 <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
#' s2 <- data.frame(V1 = c("A","C","B"), V2 = c("C","B","A"))
#' b <- bayes_compare(build_network(s1, method = "relative"),
#'                    build_network(s2, method = "relative"),
#'                    draws = 500, seed = 1)
#' plot(b, significant_only = FALSE)
#' }
#' @export
plot.net_bayes <- function(x, significant_only = TRUE, title = NULL, ...) {
  ed <- x$summary
  if (significant_only) ed <- ed[ed$sig, , drop = FALSE]
  if (nrow(ed) == 0L) {
    stop("No edges to plot",
         if (significant_only) " (none credibly different; try significant_only = FALSE)"
         else "", ".", call. = FALSE)
  }
  ed$signed <- ed$diff
  ed$value  <- abs(ed$diff)
  if (is.null(title)) {
    title <- "Differential transition network (posterior mean difference)"
  }
  .magdiff_circular(ed, value_col = "signed", abs_col = "value", title = title)
}


#' Print method for net_bayes_group
#'
#' @param x A \code{net_bayes_group} object.
#' @param ... Additional arguments (ignored).
#' @return The input object, invisibly.
#' @examples
#' s <- data.frame(V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
#'                 grp = c("X","X","Y","Y"))
#' nets <- build_network(s, method = "relative", group = "grp")
#' print(bayes_compare(nets, draws = 200, seed = 1))
#' @export
print.net_bayes_group <- function(x, ...) {
  cat("Grouped Bayesian Dirichlet-Multinomial Comparison\n")
  cat("Comparisons:", paste(names(x), collapse = ", "), "\n")
  cat("Use summary() for combined edge-level results.\n")
  invisible(x)
}


#' Summary method for net_bayes_group
#'
#' @param object A \code{net_bayes_group} object.
#' @param ... Additional arguments (ignored).
#' @return A combined data frame with a \code{comparison} column.
#' @examples
#' s <- data.frame(V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
#'                 grp = c("X","X","Y","Y"))
#' nets <- build_network(s, method = "relative", group = "grp")
#' summary(bayes_compare(nets, draws = 200, seed = 1))
#' @export
summary.net_bayes_group <- function(object, ...) {
  do.call(rbind, lapply(names(object), function(nm) {
    df <- object[[nm]]$summary
    df$comparison <- nm
    df[c("comparison", setdiff(names(df), "comparison"))]
  }))
}
