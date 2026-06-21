# ---- Analytic Dirichlet certainty (Bayesian counterpart of bootstrap_network) ----
#
# certainty() is the analytic sibling of bootstrap_network(). Where the
# bootstrap quantifies edge uncertainty by resampling sequences, certainty()
# derives it in closed form: the outgoing transitions from each state follow a
# Dirichlet-Multinomial posterior (Jeffreys prior), so each edge probability is
# marginally Beta-distributed. Posterior mean, sd, credible interval and the
# stability decision all have closed forms -- no resampling.
#
# It returns the EXACT net_bootstrap object layout (same slots, same summary
# columns, same pruned model) and carries class c("net_certainty",
# "net_bootstrap") so every net_bootstrap method works on it as a drop-in.
#
# Caveat (documented): the Dirichlet posterior treats transitions as
# independent, so on strongly heterogeneous data (latent classes) it reports
# more certainty than the sequence bootstrap. Use bootstrap_network() when the
# population is a mixture and sequences are long.

#' Analytic certainty of network edges (Bayesian Dirichlet-Multinomial)
#'
#' @description
#' Closed-form alternative to \code{\link{bootstrap_network}} for transition
#' networks. Models the outgoing transitions from each state as a
#' Dirichlet-Multinomial process: with a Jeffreys prior the posterior for state
#' \eqn{i} is \eqn{\mathrm{Dirichlet}(c_i + \mathrm{prior})}, so each edge is
#' marginally Beta and its posterior mean, standard deviation, credible interval
#' and stability decision are available analytically. No resampling, so it runs
#' in microseconds.
#'
#' The return value has the same structure as \code{\link{bootstrap_network}}
#' (same slots and summary columns) and carries class \code{c("net_certainty",
#' "net_bootstrap")}, so \code{summary()} and any code that consumes a
#' \code{net_bootstrap} object work unchanged.
#'
#' @details
#' Certainty (this function), stability (\code{\link{bootstrap_network}}) and
#' reliability (\code{\link[=network_reliability]{reliability}}) answer different questions about an
#' edge: how precisely it is pinned down by the observed counts, whether it
#' survives resampling the sequences, and whether it is consistent across
#' split-halves. Certainty and stability agree on homogeneous data; certainty is
#' over-confident when the data are a mixture of latent classes, because it
#' treats transitions clustered within a sequence as independent.
#'
#' @param x A \code{netobject} from \code{\link{build_network}} using a
#'   transition-probability method (\code{"relative"} / \code{"tna"}), or a
#'   \code{netobject_group}.
#' @param prior Numeric. Dirichlet prior concentration added to every cell
#'   (default \code{0.5}, the Jeffreys prior).
#' @param ci_level Numeric in (0,1). Tail level for credible intervals and the
#'   stability decision (default \code{0.05}, i.e. a 95\% interval). Named to
#'   match \code{bootstrap_network()}.
#' @param inference Character. \code{"stability"} (default) tests whether the
#'   posterior keeps the edge within a multiplicative \code{consistency_range}
#'   of its weight; \code{"threshold"} tests whether the edge exceeds
#'   \code{edge_threshold}.
#' @param consistency_range Numeric vector of length 2. Multiplicative bounds
#'   for stability inference (default \code{c(0.75, 1.25)}).
#' @param edge_threshold Numeric or NULL. Fixed threshold for
#'   \code{inference = "threshold"}. If NULL, defaults to the 10th percentile of
#'   non-zero edge weights.
#'
#' @return An object of class \code{c("net_certainty", "net_bootstrap")} with the
#'   same fields as \code{\link{bootstrap_network}}: \code{original}, \code{mean},
#'   \code{sd}, \code{p_values}, \code{significant}, \code{ci_lower},
#'   \code{ci_upper}, \code{cr_lower}, \code{cr_upper}, \code{summary},
#'   \code{model}, \code{method}, \code{params}, \code{ci_level},
#'   \code{inference}, \code{consistency_range}, \code{edge_threshold}, plus
#'   \code{prior} and \code{iter = NA} (no iterations).
#'
#' @examples
#' seqs <- data.frame(V1 = c("A","B","A","C","B"), V2 = c("B","C","B","A","C"),
#'                    V3 = c("C","A","C","B","A"))
#' net <- build_network(seqs, method = "relative")
#' cert <- certainty(net)
#' cert
#' summary(cert)
#'
#' @seealso \code{\link{bootstrap_network}}, \code{\link{bayes_compare}},
#'   \code{\link{network_reliability}}
#' @references
#' Johnston, L. & Jendoubi, T. (2026). How Delivery Mode Reshapes Resource
#' Engagement: A Bayesian Differential Network Analysis. TNA Workshop 2026.
#'
#' @importFrom stats qbeta pbeta quantile
#' @export
certainty <- function(x,
                      prior = 0.5,
                      ci_level = 0.05,
                      inference = c("stability", "threshold"),
                      consistency_range = c(0.75, 1.25),
                      edge_threshold = NULL) {
  inference <- match.arg(inference)

  # ---- netobject_group dispatch ----
  if (inherits(x, "netobject_group")) {
    results <- lapply(x, function(net) {
      certainty(net, prior = prior, ci_level = ci_level, inference = inference,
                consistency_range = consistency_range,
                edge_threshold = edge_threshold)
    })
    class(results) <- c("net_certainty_group", "net_bootstrap_group", "list")
    return(results)
  }

  if (inherits(x, "cograph_network")) x <- .as_netobject(x)
  stopifnot(
    inherits(x, "netobject"),
    is.numeric(prior), length(prior) == 1, prior > 0,
    is.numeric(ci_level), length(ci_level) == 1, ci_level > 0, ci_level < 1,
    is.numeric(consistency_range), length(consistency_range) == 2,
    all(consistency_range > 0)
  )

  method <- .resolve_method_alias(x$method)
  if (method != "relative") {
    stop("certainty() needs a transition-probability network ",
         "(method = \"relative\"); got '", x$method, "'. The ",
         "Dirichlet-Multinomial models per-state outgoing probabilities.",
         call. = FALSE)
  }
  if (is.null(x$frequency_matrix)) {
    stop("Transition counts ($frequency_matrix) not found. ",
         "Rebuild with build_network(method = \"relative\").", call. = FALSE)
  }
  if (!is.null(x$scaling) && !identical(x$scaling, "none")) {
    stop("certainty() needs an unscaled transition network; got scaling = '",
         x$scaling, "'. The Dirichlet posterior is on the row-probability ",
         "scale, so the consistency band would compare it against transformed ",
         "weights. Rebuild with build_network(method = \"relative\") and no ",
         "scaling.", call. = FALSE)
  }

  states <- x$nodes$label
  n <- length(states)
  counts <- x$frequency_matrix[states, states, drop = FALSE]
  W <- x$weights[states, states, drop = FALSE]
  dn <- list(states, states)

  # ---- Beta posterior per edge (Dirichlet rows) ----
  a <- counts + prior
  rowA <- matrix(rowSums(a), n, n, dimnames = dn)   # row-constant concentration
  b <- rowA - a
  mean_mat <- a / rowA
  sd_mat <- sqrt(a * b / (rowA^2 * (rowA + 1)))
  ci_lower <- matrix(qbeta(ci_level / 2, a, b), n, n, dimnames = dn)
  ci_upper <- matrix(qbeta(1 - ci_level / 2, a, b), n, n, dimnames = dn)

  # ---- Stability / threshold decision (analytic p = posterior mass) ----
  if (inference == "stability") {
    cr_lo <- pmin(W * consistency_range[1], W * consistency_range[2])
    cr_hi <- pmax(W * consistency_range[1], W * consistency_range[2])
    p_out <- pbeta(cr_lo, a, b) + (1 - pbeta(cr_hi, a, b))   # mass outside band
    p_values <- matrix(p_out, n, n, dimnames = dn)
    cr_lower <- matrix(cr_lo, n, n, dimnames = dn)
    cr_upper <- matrix(cr_hi, n, n, dimnames = dn)
  } else {
    if (is.null(edge_threshold)) {
      nz <- W[W != 0]
      edge_threshold <- if (length(nz)) unname(quantile(nz, 0.10)) else 0
    }
    p_values <- matrix(pbeta(edge_threshold, a, b), n, n, dimnames = dn)  # mass below
    cr_lower <- matrix(0, n, n, dimnames = dn)
    cr_upper <- matrix(0, n, n, dimnames = dn)
  }

  sig_mask <- (p_values < ci_level) & (W != 0)
  significant <- W * sig_mask

  stats <- list(mean = mean_mat, sd = sd_mat, p_values = p_values,
                ci_lower = ci_lower, ci_upper = ci_upper,
                cr_lower = cr_lower, cr_upper = cr_upper)

  # ---- Summary + pruned model (reuse the bootstrap builders) ----
  summary_df <- .build_bootstrap_summary(
    stats = stats, original_matrix = W, states = states,
    directed = x$directed, ci_level = ci_level, inference = inference
  )
  pruned_edges <- .extract_edges_from_matrix(significant, directed = x$directed)
  model <- list(
    weights = significant, nodes = x$nodes, edges = pruned_edges,
    directed = x$directed, method = x$method, params = x$params,
    scaling = x$scaling, threshold = x$threshold, n_nodes = n,
    n_edges = nrow(pruned_edges), level = ci_level,
    meta = list(source = "nestimate", layout = NULL,
                tna = list(method = x$method)),
    node_groups = NULL
  )
  class(model) <- c("netobject", "cograph_network")

  result <- list(
    original          = x,
    mean              = mean_mat,
    sd                = sd_mat,
    p_values          = p_values,
    significant       = significant,
    ci_lower          = ci_lower,
    ci_upper          = ci_upper,
    cr_lower          = cr_lower,
    cr_upper          = cr_upper,
    summary           = summary_df,
    model             = model,
    method            = x$method,
    params            = x$params,
    iter              = NA_integer_,
    ci_level          = ci_level,
    inference         = inference,
    consistency_range = consistency_range,
    edge_threshold    = edge_threshold,
    ci_method         = "analytic",
    prior             = prior
  )
  class(result) <- c("net_certainty", "net_bootstrap")
  result
}


#' Print Method for net_certainty
#'
#' @param x A \code{net_certainty} object.
#' @param ... Additional arguments (ignored).
#' @return The input object, invisibly.
#' @examples
#' seqs <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
#' print(certainty(build_network(seqs, method = "relative")))
#' @export
print.net_certainty <- function(x, ...) {
  method_labels <- c(
    relative = "Transition Network (relative)",
    tna      = "Transition Network (relative)"
  )
  lbl <- method_labels[.resolve_method_alias(x$method)]
  label <- if (is.na(lbl)) sprintf("Network (%s)", x$method) else unname(lbl)
  dir_label <- if (x$original$directed) "directed" else "undirected"
  ci_pct <- sprintf("%d%%", round((1 - x$ci_level) * 100))
  n_orig <- x$original$n_edges
  n_sig <- x$model$n_edges

  s <- x$summary
  if (!is.null(s) && nrow(s) > 0L) {
    sig_s <- s[s$sig, , drop = FALSE]
    if (nrow(sig_s) > 0L) {
      sig_s <- sig_s[order(abs(sig_s$mean), decreasing = TRUE), , drop = FALSE]
      top <- head(sig_s, 5L)
      cat("  Edge                   Mean     95% CI          p\n")
      cat("  -----------------------------------------------\n")
      lbl2  <- sprintf("%-20s", paste0(top$from, " \u2192 ", top$to))
      stars <- ifelse(top$p_value < 0.001, "***",
                      ifelse(top$p_value < 0.01, "** ", "*  "))
      cat(sprintf("  %s  %6.3f  [%5.3f, %5.3f]  %s\n",
                  lbl2, top$mean, top$ci_lower, top$ci_upper, stars),
          sep = "")
      if (nrow(sig_s) > 5L)
        cat(sprintf("  ... and %d more certain edges\n", nrow(sig_s) - 5L))
      cat("\n")
    }
  }

  cat(sprintf("Certainty (Dirichlet)  [%s | %s]\n", label, dir_label))
  cat(sprintf("  Prior      : Dirichlet(%.2f)  |  Nodes : %d\n",
              x$prior, x$original$n_nodes))
  cat(sprintf("  Edges      : %d certain / %d total\n", n_sig, n_orig))
  cat(sprintf("  CI         : %s  |  Inference: %s", ci_pct, x$inference))
  if (x$inference == "stability") {
    cat(sprintf("  |  CR [%.2f, %.2f]",
                x$consistency_range[1L], x$consistency_range[2L]))
  } else {
    cat(sprintf("  |  Threshold: %g", x$edge_threshold))
  }
  cat("\n")
  invisible(x)
}
