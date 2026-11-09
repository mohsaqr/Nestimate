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
#' `certainty()` gives each edge of a transition-probability network a
#' posterior mean, standard deviation, credible interval and a decision on
#' whether the edge is certain, all in closed form, and returns them in the
#' layout of [bootstrap_network()]. It models the outgoing transitions of
#' each state as a Dirichlet-Multinomial process: with a Dirichlet prior of
#' concentration `prior`, the posterior for state \eqn{i} is
#' \eqn{\mathrm{Dirichlet}(c_i + \mathrm{prior})}, where \eqn{c_i} are its
#' transition counts, so each edge probability is marginally Beta. It reads a
#' `netobject` built with `method = "relative"`, or a `netobject_group` of
#' them.
#'
#' With `inference = "stability"`, an edge is certain when the posterior
#' mass outside a band around its weight (`consistency_range`) is below
#' `ci_level`; with `"threshold"`, when the mass below `edge_threshold` is.
#' Certainty says how precisely the observed counts pin an edge down, not
#' whether the edge survives resampling the sequences ([bootstrap_network()])
#' or holds across split halves ([network_reliability()]). Certainty and
#' bootstrap stability agree on homogeneous data. Certainty treats every
#' transition as independent, so it is over-confident when the sequences come
#' from a mixture of latent classes, the more so the longer they are; use
#' [bootstrap_network()] there.
#'
#' @details
#' Nothing is resampled: the posterior quantities are Beta means, variances
#' and quantiles (`stats::qbeta()`, `stats::pbeta()`), so a call takes
#' milliseconds.
#'
#' @param x The network: a `netobject` built with [build_network()] using
#'   `method = "relative"` (or its aliases `"tna"`, `"transition"`), without
#'   `scaling`, or a `netobject_group` of such networks, each handled on its
#'   own. A bare `cograph_network` is recounted from the sequences it
#'   carries, and used only when they reproduce its weights; otherwise, and
#'   for an `mcml`, an error of class `nestimate_bad_input` is raised.
#'   Another method or a scaled network raises an error.
#'   For the `print()` method: a `net_certainty` object.
#' @param prior The Dirichlet concentration added to every cell of the count
#'   matrix, zero cells and self-loops included; a positive number. Default
#'   `0.5`, the Jeffreys prior.
#' @param ci_level The tail level, in (0, 1), of the equal-tailed credible
#'   intervals and of the decision: an edge is certain when its posterior
#'   mass (the `p_value` column) is below it. Default `0.05`, a 95%
#'   interval. Named to match [bootstrap_network()].
#' @param inference The decision rule. `"stability"` (default) takes the
#'   posterior mass outside the band `consistency_range` times the edge's
#'   weight; `"threshold"` takes the posterior mass below `edge_threshold`.
#' @param consistency_range For `inference = "stability"`: the two
#'   multiplicative bounds of the band around each weight, positive numbers
#'   in either order. Default `c(0.75, 1.25)`.
#' @param edge_threshold For `inference = "threshold"`: the weight an edge
#'   must exceed. Default `NULL`: the 10th percentile of the non-zero
#'   weights. Ignored for `"stability"`.
#'
#' @return For a `netobject`: a `net_certainty` object (class
#'   `c("net_certainty", "net_bootstrap")`), so the `net_bootstrap` methods
#'   work on it. `print()` lists up to five certain edges with the largest
#'   absolute posterior means (with their intervals), then the prior, the
#'   count of certain edges and the settings;
#'   [summary()] returns one row per non-zero edge with columns `from`,
#'   `to`, `weight`, `mean`, `sd`, `p_value`, `sig`, `ci_lower` and
#'   `ci_upper`, plus `cr_lower` and `cr_upper` (the band) for
#'   `inference = "stability"`; `cograph::splot()` draws it. The object has
#'   the fields of a [bootstrap_network()] result: `original` (the input
#'   network), `mean`, `sd`, `p_values`, `significant` (the weights of the
#'   certain edges, zero elsewhere), `ci_lower`, `ci_upper`, `cr_lower`,
#'   `cr_upper` (matrices), `summary`, `model` (the network of certain edges
#'   only), `method`, `params`, `ci_level`, `inference`,
#'   `consistency_range` and `edge_threshold`, plus `prior`,
#'   `ci_method = "analytic"` and `iter = NA` (no iterations).
#'
#'   For a `netobject_group`: a named list of those objects, one per
#'   network, of class `c("net_certainty_group", "net_bootstrap_group",
#'   "list")`, printed and summarised by the `net_bootstrap_group` methods
#'   ([summary()] stacks the per-network tables with a `group` column).
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

  if (inherits(x, "mcml")) {
    .stop_mcml_input("certainty()", "reads one transition network")
  }
  if (inherits(x, "cograph_network")) x <- .as_netobject(x)
  .assert(
    inherits(x, "netobject"),
    is.numeric(prior), length(prior) == 1, prior > 0,
    is.numeric(ci_level), length(ci_level) == 1, ci_level > 0, ci_level < 1,
    is.numeric(consistency_range), length(consistency_range) == 2,
    all(consistency_range > 0)
  )

  method <- .resolve_method_alias(x$method)
  if (method != "relative") {
    stop(errorCondition(
      paste0(
        "certainty() needs a transition-probability network ",
        "(method = \"relative\"); got '",
        x$method,
        "'. The ",
        "Dirichlet-Multinomial models per-state outgoing probabilities."),
      class = "nestimate_bad_input", call = NULL))
  }
  x$frequency_matrix <- .network_transition_counts(x, "certainty()")
  if (!is.null(x$scaling) && !identical(x$scaling, "none")) {
    stop(errorCondition(
      paste0(
        "certainty() needs an unscaled transition network; got scaling = '",
        x$scaling,
        "'. The Dirichlet posterior is on the row-probability ",
        "scale, so the consistency band would compare it against transformed ",
        "weights. Rebuild with build_network(method = \"relative\") and no ",
        "scaling."),
      class = "nestimate_bad_input", call = NULL))
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


#' @rdname certainty
#' @param ... For the `print()` method: ignored.
#' @return The `print()` method returns its input invisibly.
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
