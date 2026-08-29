# ==============================================================================
# plot.net_network_comparison(): one view per call.
#
# Views: "networks" (each network once), "difference" (the signed difference
# per pair), "edges", "nodes", "global", "heatmap", "scatter".
#
# Colour contract (every panel, every backend):
#   .cn_col_a  "#4A6FE3"  first / reference network higher
#   .cn_col_b  "#D33F6A"  second network higher
#   .cn_col_eq "grey65"   no difference
# Sign is never carried by colour alone: solid vs dashed segments, circle vs
# square endpoint markers, and the printed signed value repeat it. Evidence
# (when `test` was run) is carried by opacity and a starred label, never by
# deleting an edge. One view per plot() call: no composite dashboard.
# `combined = FALSE` splits a multi-pair view into one plot per pair rather
# than cramming them into facets.
# Network panels are base graphics via cograph::splot() (built-in circular
# drawer only when cograph is absent).
# ==============================================================================

.cn_col_a <- "#4A6FE3"
.cn_col_b <- "#D33F6A"
.cn_col_eq <- "grey65"

.cn_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = ggplot2::margin(t = 0, b = 0),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# Panel grid for a base-graphics page of `n` panels: as square as possible.
.cn_grid <- function(n) {
  nr <- max(1L, floor(sqrt(n)))
  c(nr, ceiling(n / nr))
}

# Colour + linetype scales for the `higher` channel. With a single pair the
# legend shows the real network names; with several pairs the strip carries
# the names and the legend describes the role.
.cn_role_scales <- function(x, pair_names) {
  one_pair <- length(pair_names) == 1L
  if (one_pair) {
    pr <- x$pairs[x$pairs$pair == pair_names, ]
    labs <- c(paste(pr$network_a, "higher"), paste(pr$network_b, "higher"),
              "equal")
    vals <- c(.cn_col_a, .cn_col_b, .cn_col_eq)
    names(vals) <- c(pr$network_a, pr$network_b, "equal")
    lt <- c("solid", "22", "dotted")
    names(lt) <- names(vals)
    key <- names(vals)
  } else {
    key <- c("a", "b", "equal")
    labs <- if (is.null(x$reference)) {
      c("first network of pair higher", "second network of pair higher", "equal")
    } else {
      c(sprintf("%s (reference) higher", x$reference),
        "comparison network higher", "equal")
    }
    vals <- c(a = .cn_col_a, b = .cn_col_b, equal = .cn_col_eq)
    lt <- c(a = "solid", b = "22", equal = "dotted")
  }
  list(
    key = key,
    colour = ggplot2::scale_colour_manual(values = vals, labels = labs,
                                          breaks = key, drop = FALSE,
                                          name = NULL),
    fill = ggplot2::scale_fill_manual(values = vals, labels = labs,
                                      breaks = key, drop = FALSE, name = NULL),
    linetype = ggplot2::scale_linetype_manual(values = lt, labels = labs,
                                              breaks = key, drop = FALSE,
                                              name = NULL)
  )
}

# Map the `higher` factor to the role key used by .cn_role_scales().
.cn_role <- function(d, one_pair) {
  h <- as.character(d$higher)
  if (one_pair) return(factor(h, levels = unique(c(d$network_a[1L], d$network_b[1L], "equal"))))
  factor(ifelse(h == d$network_a, "a", ifelse(h == d$network_b, "b", "equal")),
         levels = c("a", "b", "equal"))
}

# Strip label: the direction of the difference, "early - mid". The pair name
# ("early vs mid") already names the comparison, so repeating it in the strip
# only costs width.
.cn_strip <- function(d) {
  lab <- sprintf("%s - %s", d$network_a, d$network_b)
  factor(lab, levels = unique(lab))
}

.cn_resolve_pairs <- function(x, pair) .cn_match_pairs(x, pair)

# Networks touched by a set of pairs, in the order they were supplied.
.cn_networks_in <- function(x, pairs) {
  pr <- x$pairs[x$pairs$pair %in% pairs, , drop = FALSE]
  intersect(names(x$networks), unique(c(pr$network_a, pr$network_b)))
}

.cn_has_evidence <- function(d) "sig" %in% names(d) && !all(is.na(d$sig))

.cn_evidence_caption <- function(x) {
  if (identical(x$test, "none")) return(NULL)
  sprintf("Evidence: %s, iter = %d, alpha = %s%s. Non-significant results are faded, never removed.",
          paste(x$test, collapse = " + "), x$iter, format(x$alpha),
          if (x$adjust != "none") paste0(", adjust = ", x$adjust) else "")
}

# Facet strips are built from a pair-prefixed factor ("<pair>\r<label>") so
# each panel keeps its own row order; the prefix is stripped at render time.
# Assigning de-duplicated `levels()` instead would MERGE the levels and force
# every panel onto the first panel's order.
.cn_strip_prefix <- function(l) sub("^.*\r", "", l)


#' @rdname compare_networks
#' @param type For `plot()`: one view per call. `"networks"` (default) draws
#'   each network once with `cograph::splot()`; `"difference"` draws the
#'   signed difference network of each pair; `"edges"` is a ranked dumbbell
#'   of the largest edge differences; `"nodes"` the same for centralities;
#'   `"global"` the 22 comparison metrics; `"heatmap"` the signed difference
#'   matrix; `"scatter"` weight against weight; `"inference"` is a forest of
#'   the edge differences on the difference scale, with credible intervals
#'   when the Bayesian backend ran and the p-value printed per edge. It needs
#'   `test != "none"` and raises `nestimate_compare_no_test` otherwise.
#' @param pair Optional selection of comparisons to draw; default all. Any of:
#'   the pair name(s) as printed (`"A vs B"`); the two network names
#'   (`c("A", "B")`, either order); one network name (`"A"`, every pair it
#'   takes part in); or index/indices into the pair table (`1`, `c(1, 3)`).
#'   For `type = "networks"` this selects the networks taking part in the
#'   chosen pairs.
#' @param combined When `TRUE` (default), a multi-pair view is one figure
#'   (facets for the `ggplot` views, one base-graphics page for `"networks"`
#'   and `"difference"`). When `FALSE`, the view is split: the `ggplot`
#'   views return a named list of single-pair plots, one per pair, and the
#'   base-graphics views draw one panel per page.
#' @param top_n Number of edges shown in the edge view (largest absolute
#'   differences first). Default `20`.
#' @param measure Optional centrality measure(s) to restrict the node view.
#' @param labels Logical; print the signed difference on the edge and node
#'   views. Default `TRUE`. The heatmap always shows its values.
#' @param digits Decimals in printed values. Default `2`.
#' @param what Deprecated alias for `type`, kept so existing calls keep
#'   working.
#' @param ... Passed to `cograph::splot()` for the network views (e.g.
#'   `layout`, `node_size`, `minimum`); ignored by the other views.
#' @return `plot()` returns a `ggplot` for `type = "edges"`, `"nodes"`,
#'   `"global"`, `"heatmap"`, `"scatter"` and `"inference"`, or a named list of such plots
#'   (one per pair) when `combined = FALSE`. `type = "networks"` and
#'   `type = "difference"` draw in base graphics -- with `cograph::splot()`
#'   when cograph is installed, otherwise with a built-in circular drawer --
#'   and return `NULL` invisibly.
#' @export
plot.net_network_comparison <- function(x,
                                        type = c("networks", "difference",
                                                 "edges", "nodes", "global",
                                                 "heatmap", "scatter",
                                                 "inference"),
                                        pair = NULL,
                                        combined = TRUE,
                                        top_n = 20L,
                                        measure = NULL,
                                        labels = TRUE,
                                        digits = 2L,
                                        what = NULL,
                                        ...) {
  if (!is.null(what)) type <- what
  type <- match.arg(type, c("networks", "difference", "edges", "nodes",
                            "global", "heatmap", "scatter", "inference"))
  stopifnot(
    "`combined` must be TRUE or FALSE" = isTRUE(combined) || isFALSE(combined),
    "`top_n` must be a single positive integer" =
      is.numeric(top_n) && length(top_n) == 1L && top_n >= 1,
    "`labels` must be TRUE or FALSE" = isTRUE(labels) || isFALSE(labels)
  )
  digits <- .cn_check_digits(digits)
  pairs <- .cn_resolve_pairs(x, pair)

  if (type %in% c("networks", "difference")) {
    return(.cn_plot_base_view(x, type, pairs, combined, ...))
  }
  if (!combined) {
    out <- lapply(pairs, function(p)
      .cn_plot_ggplot_view(x, type, p, top_n, measure, labels, digits))
    names(out) <- pairs
    return(invisible(out))
  }
  .cn_plot_ggplot_view(x, type, pairs, top_n, measure, labels, digits)
}

.cn_plot_ggplot_view <- function(x, type, pairs, top_n, measure, labels,
                                 digits) {
  switch(type,
    edges = .cn_plot_edges(x, pairs, top_n, labels, digits),
    nodes = .cn_plot_nodes(x, pairs, measure, labels, digits),
    global = .cn_plot_global(x, pairs, digits = digits),
    heatmap = .cn_plot_heatmap(x, pairs, digits),
    scatter = .cn_plot_scatter(x, pairs, digits),
    inference = .cn_plot_inference(x, pairs, top_n, labels, digits)
  )
}


# ---- Dumbbell core (edges and nodes share it) --------------------------------

# d: rows with columns label, weight_a, weight_b, diff, higher, network_a,
# network_b, pair, optional sig, bayes_ci_lower/upper.
.cn_dumbbell <- function(x, d, one_pair, x_label, title, subtitle, labels,
                         digits = 2L) {
  sc <- .cn_role_scales(x, unique(d$pair))
  d$role <- .cn_role(d, one_pair)
  d$strip <- .cn_strip(d)
  evidence <- .cn_has_evidence(d)
  d$alpha_lvl <- if (evidence) ifelse(is.na(d$sig) | !d$sig, 0.35, 1) else 1
  net_a <- d$network_a[1L]
  net_b <- d$network_b[1L]
  d$role_a <- if (one_pair) net_a else "first of pair"
  d$role_b <- if (one_pair) net_b else "second of pair"

  p <- ggplot2::ggplot(d, ggplot2::aes(y = .data$label)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$weight_b, xend = .data$weight_a,
                   yend = .data$label, colour = .data$role,
                   linetype = .data$role, alpha = .data$alpha_lvl),
      linewidth = 1.6, lineend = "round"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$weight_b, shape = .data$role_b,
                   alpha = .data$alpha_lvl),
      size = 2.6, fill = "white", stroke = 0.9
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$weight_a, shape = .data$role_a,
                   alpha = .data$alpha_lvl),
      size = 2.6, fill = "black", stroke = 0.9
    ) +
    sc$colour + sc$linetype +
    ggplot2::scale_shape_manual(
      values = stats::setNames(c(19, 22), c(d$role_a[1L], d$role_b[1L])),
      breaks = c(d$role_a[1L], d$role_b[1L]), name = "Network"
    ) +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_y_discrete(labels = .cn_strip_prefix) +
    ggplot2::labs(x = x_label, y = NULL, title = title,
                  subtitle = if (evidence && labels) {
                    paste0(subtitle, if (is.null(subtitle)) "" else "; ",
                           "a starred value is significant")
                  } else subtitle,
                  caption = .cn_evidence_caption(x)) +
    .cn_theme()

  if ("bayes_ci_lower" %in% names(d)) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$weight_b + .data$bayes_ci_lower,
                   xmax = .data$weight_b + .data$bayes_ci_upper,
                   alpha = .data$alpha_lvl),
      width = 0.25, linewidth = 0.4, colour = "grey30"
    )
  }
  if (labels) {
    d$lab <- paste0(.cn_fmt(d$diff, digits, signed = TRUE),
                    if (evidence) ifelse(!is.na(d$sig) & d$sig, "*", "") else "")
    p <- p + ggplot2::geom_text(
      data = d,
      ggplot2::aes(x = pmax(.data$weight_a, .data$weight_b),
                   label = .data$lab),
      hjust = -0.25, size = 3.5
    )
  }
  # Leave room on the right for the printed difference so it is never clipped.
  # The label sits at the largest of the two weights, so the headroom has to
  # cover a whole label ("+0.10", "-0.26*") drawn from the panel's rightmost
  # point -- 0.22 still cut the last character on narrow facets.
  p <- p + ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(mult = if (labels) c(0.06, 0.34) else c(0.05, 0.05))
  )
  if (!one_pair) {
    p <- p + ggplot2::facet_wrap(~ strip, scales = "free_y")
  }
  p
}

.cn_plot_edges <- function(x, pairs, top_n, labels, digits = 2L) {
  # Rows are ordered by |diff| within pair and the level is prefixed with the
  # pair so each facet keeps its own order (the prefix is stripped when the
  # axis is drawn, see .cn_strip_prefix()).
  d <- .cn_top_edges(.cn_pair_subset(x$edges, x, pairs), pairs, top_n)
  one_pair <- length(pairs) == 1L
  pr <- x$pairs[x$pairs$pair == pairs[1L], ]
  n_shown <- min(top_n, max(table(d$pair)))
  subtitle <- if (one_pair) {
    sprintf("Top %d transitions by |%s - %s|", n_shown, pr$network_a, pr$network_b)
  } else {
    sprintf("Top %d transitions per pair by |difference|", n_shown)
  }
  .cn_dumbbell(x, d, one_pair,
               x_label = sprintf("Edge weight (scaling = %s)", x$scaling),
               title = "Largest edge differences", subtitle = subtitle,
               labels = labels, digits = digits)
}

.cn_plot_nodes <- function(x, pairs, measure, labels, digits = 2L) {
  if (is.null(x$nodes)) {
    .cn_stop("No node table: `compare_networks()` was called with no `measures`.",
             "no_nodes")
  }
  d <- .cn_pair_subset(x$nodes, x, pairs)
  if (!is.null(measure)) {
    unknown <- setdiff(measure, x$measures)
    if (length(unknown) > 0L) {
      .cn_stop(paste0("Unknown measure(s): ", paste(unknown, collapse = ", "),
                      ". Available: ", paste(x$measures, collapse = ", "), "."),
               "unknown_measure")
    }
    d <- d[d$measure %in% measure, , drop = FALSE]
  }
  # Node order: mean |diff| across measures and pairs, so every panel of the
  # figure reads on the same vertical scale.
  ord <- stats::aggregate(abs_diff ~ node, data = d, FUN = mean)
  ord <- ord[order(ord$abs_diff, ord$node), ]
  d$label <- factor(d$node, levels = ord$node)
  d$measure <- factor(d$measure, levels = intersect(x$measures, unique(d$measure)))
  d$weight_a <- d$value_a
  d$weight_b <- d$value_b
  one_pair <- length(pairs) == 1L
  p <- .cn_dumbbell(x, d, one_pair, x_label = "Centrality",
                    title = "Node centrality differences",
                    subtitle = if (one_pair) {
                      pr <- x$pairs[x$pairs$pair == pairs, ]
                      sprintf("%s vs %s", pr$network_a, pr$network_b)
                    } else NULL,
                    labels = labels, digits = digits)
  if (one_pair) {
    p + ggplot2::facet_wrap(~ measure, scales = "free_x")
  } else {
    p + ggplot2::facet_grid(strip ~ measure, scales = "free") +
      ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0))
  }
}


# ---- Inference: forest of edge differences -----------------------------------

# Rank the rows of `tab` by |abs_diff| within pair, keep the top `top_n`, and
# build the pair-prefixed row factor the facets need. Shared by the edge
# dumbbell and the inference forest.
.cn_top_edges <- function(tab, pairs, top_n) {
  d <- do.call(rbind, lapply(split(tab, factor(tab$pair, levels = pairs)),
                             function(dd) {
    dd <- dd[order(-dd$abs_diff, dd$from, dd$to), , drop = FALSE]
    utils::head(dd, top_n)
  }))
  d <- d[order(match(d$pair, pairs), d$abs_diff), , drop = FALSE]
  keyed <- paste(d$pair, paste(d$from, "->", d$to), sep = "\r")
  d$label <- factor(keyed, levels = unique(keyed))
  d
}

# The difference scale, not the weight scale: `diff` with its credible
# interval (Bayesian backend) against a zero reference line, so an interval
# that crosses zero is read directly off the axis. A filled marker is a
# significant edge, an open one is not; the p-value is printed per row.
.cn_plot_inference <- function(x, pairs, top_n, labels, digits = 2L) {
  if (identical(x$test, "none")) {
    .cn_stop(paste0("No inference to plot: `compare_networks()` was called ",
                    "with test = \"none\". Re-run it with, for example, ",
                    "test = \"permutation\"."),
             "no_test")
  }
  d <- .cn_top_edges(.cn_pair_subset(x$edges, x, pairs), pairs, top_n)
  one_pair <- length(pairs) == 1L
  sc <- .cn_role_scales(x, pairs)
  d$role <- .cn_role(d, one_pair)
  d$strip <- .cn_strip(d)
  evidence <- .cn_has_evidence(d)
  d$alpha_lvl <- if (evidence) ifelse(is.na(d$sig) | !d$sig, 0.4, 1) else 1
  # Filled marker = significant, open marker = not. Colour still carries the
  # sign, so the fill only ever adds evidence on top of it.
  h <- as.character(d$higher)
  d$point_col <- ifelse(h == d$network_a, .cn_col_a,
                        ifelse(h == d$network_b, .cn_col_b, .cn_col_eq))
  d$point_fill <- if (evidence) {
    ifelse(is.na(d$sig) | !d$sig, "white", d$point_col)
  } else d$point_col

  has_ci <- all(c("bayes_ci_lower", "bayes_ci_upper") %in% names(d))
  # Name the printed quantity after the column it comes from: `bayes_p` is
  # the two-sided Bayesian p-equivalent 2 * (1 - pd), not a frequentist
  # p-value, and printing it as "p" beside a credible interval reads as one.
  p_col <- if ("perm_p" %in% names(d)) "perm_p" else
    if ("bayes_p" %in% names(d)) "bayes_p" else NULL
  p_lab <- if (identical(p_col, "perm_p")) "p" else "p_bayes"

  x_label <- if (one_pair) {
    pr <- x$pairs[x$pairs$pair == pairs, ]
    sprintf("Difference (%s - %s)", pr$network_a, pr$network_b)
  } else "Difference (first - second network of pair)"

  # A bar drawn from zero, not a dot floating beside the line: magnitude is
  # read as length, and the direction of the bar repeats the sign.
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$diff, y = .data$label)) +
    ggplot2::geom_vline(xintercept = 0, colour = .cn_col_eq,
                        linetype = "dotted") +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = .data$diff, yend = .data$label,
                   colour = .data$role, linetype = .data$role,
                   alpha = .data$alpha_lvl),
      linewidth = 1.6, lineend = "round"
    ) + sc$linetype
  if (has_ci) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$bayes_ci_lower, xmax = .data$bayes_ci_upper,
                   alpha = .data$alpha_lvl),
      width = 0.22, linewidth = 0.5, colour = "grey25"
    )
  }
  p <- p +
    ggplot2::geom_point(
      ggplot2::aes(colour = .data$role, fill = .data$point_fill,
                   alpha = .data$alpha_lvl),
      shape = 21, size = 3, stroke = 1.1
    ) +
    sc$colour +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_y_discrete(labels = .cn_strip_prefix) +
    # scale_fill_identity() draws no key, which would leave the legend showing
    # hollow markers -- the plot's own symbol for "not significant".
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(
      shape = 21, fill = c(.cn_col_a, .cn_col_b, .cn_col_eq), alpha = 1))) +
    ggplot2::labs(x = x_label, y = NULL,
                  title = "Edge differences with evidence",
                  subtitle = if (has_ci) {
                    sprintf("Bar = difference from zero, whisker = %d%% credible interval; filled marker = significant",
                            round(100 * (1 - x$alpha)))
                  } else "Bar = difference from zero; filled marker = significant",
                  caption = .cn_evidence_caption(x)) +
    .cn_theme()

  if (labels && !is.null(p_col)) {
    d$lab <- ifelse(is.na(d[[p_col]]), NA_character_,
                    sprintf("%s = %.3f", p_lab, d[[p_col]]))
    # One right-hand column of p-values, as a forest plot reads. Anchoring
    # each label to its own interval instead would pile every negative row
    # onto the zero line.
    p <- p + ggplot2::geom_text(
      data = d[!is.na(d$lab), , drop = FALSE],
      ggplot2::aes(y = .data$label, label = .data$lab),
      x = Inf, hjust = 1.05, size = 3.2, colour = "grey25",
      inherit.aes = FALSE
    ) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.30)))
  }
  if (!one_pair) p <- p + ggplot2::facet_wrap(~ strip, scales = "free_y")
  p
}


# ---- Global metrics ------------------------------------------------------------

# "Identical networks" reads as 1 for correlation/similarity metrics and as 0
# for deviation/dissimilarity metrics; that value anchors the dotted reference
# line each category is read against.
.cn_global_ref <- function(category) {
  ifelse(category %in% c("Correlations", "Similarities",
                         "Pattern Similarities"), 1, 0)
}

.cn_headline_keys <- function() {
  c("pearson", "spearman", "cosine", "jaccard", "mean_abs_diff", "rms_diff",
    "bray_curtis", "sign_agreement", "perm_M", "perm_S",
    "boot_density", "boot_mean_weight", "boot_centralization",
    "boot_reciprocity")
}

.cn_plot_global <- function(x, pairs, headline = FALSE, digits = 2L) {
  d <- .cn_pair_subset(x$global, x, pairs)
  if (headline) d <- d[d$key %in% .cn_headline_keys(), , drop = FALSE]
  d$metric <- factor(d$metric, levels = rev(unique(x$global$metric)))
  d$category <- factor(d$category, levels = unique(x$global$category))
  one_pair <- length(pairs) == 1L
  # Reference lines: 1 for correlations/similarities, 0 elsewhere.
  ref <- data.frame(category = levels(d$category), stringsAsFactors = FALSE)
  ref$at <- .cn_global_ref(ref$category)
  ref$category <- factor(ref$category, levels = levels(d$category))
  ref <- ref[ref$category %in% d$category, , drop = FALSE]
  d$has_ci <- "boot_ci_lower" %in% names(d) & !is.na(d$boot_ci_lower %||% NA)
  d$ptext <- NA_character_
  if ("perm_p" %in% names(d)) d$ptext <- ifelse(is.na(d$perm_p), d$ptext, sprintf("p = %.3f", d$perm_p))
  if ("boot_p" %in% names(d)) d$ptext <- ifelse(is.na(d$boot_p), d$ptext, sprintf("p = %.3f", d$boot_p))
  d$vtext <- .cn_fmt(d$value, digits)

  # Two columns, not five: with free scales each category carries its own
  # metric names, and five side-by-side panels clip both the strip titles and
  # the axis labels.
  n_cat <- nlevels(droplevels(d$category))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$value, y = .data$metric)) +
    ggplot2::geom_vline(data = ref, ggplot2::aes(xintercept = .data$at),
                        colour = .cn_col_eq, linetype = "dotted") +
    ggplot2::facet_wrap(~ category, scales = "free", ncol = min(2L, n_cat)) +
    ggplot2::labs(x = NULL, y = NULL, title = "Global comparison metrics",
                  subtitle = if (one_pair) {
                    pr <- x$pairs[x$pairs$pair == pairs, ]
                    sprintf("%s vs %s", pr$network_a, pr$network_b)
                  } else NULL,
                  caption = .cn_evidence_caption(x)) +
    .cn_theme()
  if ("boot_ci_lower" %in% names(d)) {
    p <- p + ggplot2::geom_errorbar(
      data = d[d$has_ci, , drop = FALSE],
      ggplot2::aes(xmin = .data$boot_ci_lower, xmax = .data$boot_ci_upper),
      width = 0.3, colour = "grey40"
    )
  }
  if (one_pair) {
    p <- p +
      ggplot2::geom_point(size = 2.8, colour = "black") +
      ggplot2::geom_text(ggplot2::aes(label = .data$vtext),
                         hjust = -0.35, size = 3.5)
    if (any(!is.na(d$ptext))) {
      p <- p + ggplot2::geom_text(
        data = d[!is.na(d$ptext), , drop = FALSE],
        ggplot2::aes(label = .data$ptext), hjust = 1.3, size = 3,
        colour = "grey30"
      )
    }
    p <- p + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.25, 0.45)))
  } else {
    pal <- rep(.okabe_ito, length.out = length(pairs))
    shp <- rep(c(16, 17, 15, 18, 8, 7, 9, 10, 11), length.out = length(pairs))
    d$pair <- factor(d$pair, levels = pairs)
    p <- p +
      ggplot2::geom_point(
        data = d,
        ggplot2::aes(colour = .data$pair, shape = .data$pair),
        size = 2.8, position = ggplot2::position_dodge(width = 0.5)
      ) +
      ggplot2::scale_colour_manual(values = pal, name = "Pair") +
      ggplot2::scale_shape_manual(values = shp, name = "Pair") +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.1)))
  }
  p
}


# ---- Heatmap and scatter ----------------------------------------------------------

.cn_plot_heatmap <- function(x, pairs, digits = 2L) {
  d <- .cn_pair_subset(x$edges, x, pairs)
  nodes <- rownames(x$matrices[[1L]])
  d$from <- factor(d$from, levels = rev(nodes))
  d$to <- factor(d$to, levels = nodes)
  d$strip <- .cn_strip(d)
  lim <- max(abs(d$diff), .cn_tol)
  one_pair <- length(pairs) == 1L
  fill_name <- if (one_pair) {
    pr <- x$pairs[x$pairs$pair == pairs, ]
    sprintf("%s - %s", pr$network_a, pr$network_b)
  } else if (is.null(x$reference)) "first - second" else
    sprintf("%s - other", x$reference)
  evidence <- .cn_has_evidence(d)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$to, y = .data$from)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$diff), colour = "white",
                       linewidth = 0.4) +
    ggplot2::scale_fill_gradient2(
      low = .cn_col_b, mid = "white", high = .cn_col_a, midpoint = 0,
      limits = c(-lim, lim), name = fill_name, n.breaks = 5L,
      # A wide, short bar: the default is too narrow for signed labels, which
      # then run into each other.
      guide = ggplot2::guide_colourbar(
        barwidth = ggplot2::unit(10, "lines"),
        barheight = ggplot2::unit(0.7, "lines"),
        title.position = "top"
      )
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = "to", y = "from", title = "Signed edge differences",
                  caption = .cn_evidence_caption(x)) +
    .cn_theme()
  if (evidence) {
    ns <- d[is.na(d$sig) | !d$sig, , drop = FALSE]
    p <- p + ggplot2::geom_tile(data = ns, fill = NA, colour = "grey45",
                                linetype = "dotted", linewidth = 0.5)
  }
  d$lab <- .cn_fmt(d$diff, digits, signed = TRUE)
  p <- p + ggplot2::geom_text(data = d, ggplot2::aes(label = .data$lab), size = 3.3)
  if (!one_pair) p <- p + ggplot2::facet_wrap(~ strip)
  p
}

.cn_plot_scatter <- function(x, pairs, digits = 2L) {
  d <- .cn_pair_subset(x$edges, x, pairs)
  d$strip <- .cn_strip(d)
  one_pair <- length(pairs) == 1L
  g <- .cn_pair_subset(x$global, x, pairs)
  stat_txt <- vapply(pairs, function(pr) {
    gg <- g[g$pair == pr, ]
    sprintf("r = %s, rho = %s, tau = %s",
            .cn_fmt(gg$value[gg$key == "pearson"], digits),
            .cn_fmt(gg$value[gg$key == "spearman"], digits),
            .cn_fmt(gg$value[gg$key == "kendall"], digits))
  }, character(1))
  ann <- data.frame(pair = pairs, txt = stat_txt, stringsAsFactors = FALSE)
  ann <- merge(ann, unique(d[c("pair", "strip")]), by = "pair")
  pr1 <- x$pairs[x$pairs$pair == pairs[1L], ]
  xlab <- if (one_pair) pr1$network_a else if (is.null(x$reference)) "First network of pair" else sprintf("%s (reference)", x$reference)
  ylab <- if (one_pair) pr1$network_b else "Second network of pair"
  status_cols <- stats::setNames(.okabe_ito[c(5, 1, 6, 8)],
                                 c("both", "only_a", "only_b", "neither"))
  status_shapes <- c(both = 16, only_a = 17, only_b = 15, neither = 1)
  status_labs <- c(both = "present in both",
                   only_a = sprintf("only in %s", if (one_pair) pr1$network_a else "first"),
                   only_b = sprintf("only in %s", if (one_pair) pr1$network_b else "second"),
                   neither = "absent in both")
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$weight_a, y = .data$weight_b)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = .cn_col_eq) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                         colour = .cn_col_a, linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$status, shape = .data$status),
                        size = 2.4, alpha = 0.85) +
    ggplot2::geom_text(data = ann, ggplot2::aes(label = .data$txt),
                       x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 3.5,
                       inherit.aes = FALSE) +
    ggplot2::scale_colour_manual(values = status_cols, labels = status_labs,
                                 drop = FALSE, name = "Edge") +
    ggplot2::scale_shape_manual(values = status_shapes, labels = status_labs,
                                drop = FALSE, name = "Edge") +
    ggplot2::labs(x = xlab, y = ylab, title = "Edge weights, network against network") +
    .cn_theme()
  if (!one_pair) p <- p + ggplot2::facet_wrap(~ strip)
  p
}


# ---- Networks: base graphics, cograph::splot() when installed ----------------

# Draw one network (or a signed difference) in base graphics on a circle.
# Used only when cograph is not installed; cograph::splot() is preferred.
.cn_draw_network_base <- function(W, title, signed, directed, ...) {
  n <- nrow(W)
  nodes <- rownames(W)
  ang <- pi / 2 - 2 * pi * (seq_len(n) - 1L) / n
  px <- cos(ang)
  py <- sin(ang)
  graphics::plot.new()
  graphics::plot.window(xlim = c(-1.45, 1.45), ylim = c(-1.45, 1.45), asp = 1)
  graphics::title(main = title, font.main = 1, cex.main = 1.1)
  idx <- which(abs(W) > .cn_tol, arr.ind = TRUE)
  if (!directed) idx <- idx[idx[, 1L] <= idx[, 2L], , drop = FALSE]
  if (nrow(idx) > 0L) {
    w <- W[idx]
    lwd <- 0.6 + 4 * abs(w) / max(abs(w))
    col <- if (signed) ifelse(w > 0, .cn_col_a, .cn_col_b) else "grey40"
    lty <- if (signed) ifelse(w > 0, 1L, 2L) else rep(1L, length(w))
    loop <- idx[, 1L] == idx[, 2L]
    if (any(!loop)) {
      i <- idx[!loop, 1L]; j <- idx[!loop, 2L]
      dx <- px[j] - px[i]; dy <- py[j] - py[i]
      len <- sqrt(dx^2 + dy^2)
      shrink <- 0.16
      graphics::arrows(px[i] + dx / len * shrink, py[i] + dy / len * shrink,
                       px[j] - dx / len * shrink, py[j] - dy / len * shrink,
                       length = if (directed) 0.08 else 0, angle = 20,
                       lwd = lwd[!loop], col = col[!loop], lty = lty[!loop])
    }
    if (any(loop)) {
      i <- idx[loop, 1L]
      graphics::symbols(px[i] * 1.27, py[i] * 1.27, circles = rep(0.09, sum(loop)),
                        inches = FALSE, add = TRUE, fg = col[loop],
                        lwd = lwd[loop], lty = lty[loop])
    }
  }
  graphics::points(px, py, pch = 21, bg = "white", col = "black", cex = 4.2)
  graphics::text(px, py, labels = nodes, cex = 0.8)
  invisible(NULL)
}

# One panel: a network as estimated.
.cn_draw_one_network <- function(x, nm, ...) {
  if (requireNamespace("cograph", quietly = TRUE)) {
    cograph::splot(x$networks[[nm]], title = nm, ...)
  } else {
    .cn_draw_network_base(x$matrices[[nm]], nm, signed = FALSE,
                          directed = x$directed[[nm]])
  }
  invisible(NULL)
}

# One panel: the signed difference of a pair, a - b.
.cn_draw_one_difference <- function(x, pr_name, ...) {
  pr <- x$pairs[x$pairs$pair == pr_name, ]
  ttl <- sprintf("%s - %s", pr$network_a, pr$network_b)
  if (requireNamespace("cograph", quietly = TRUE)) {
    cograph::splot(x$differences[[pr_name]], title = ttl,
                   edge_positive_color = .cn_col_a,
                   edge_negative_color = .cn_col_b, ...)
  } else {
    .cn_draw_network_base(x$matrices[[pr$network_a]] - x$matrices[[pr$network_b]],
                          ttl, signed = TRUE,
                          directed = x$directed[[pr$network_a]] ||
                            x$directed[[pr$network_b]])
  }
  invisible(NULL)
}

# "networks" draws each network once; "difference" draws one signed
# difference per pair. Drawing the a / b / (a - b) triptych per pair instead
# would repeat every network once per pair it takes part in (a four-network
# comparison would draw 18 panels for 4 + 6 distinct ones).
.cn_plot_base_view <- function(x, type, pairs, combined, ...) {
  items <- if (identical(type, "networks")) .cn_networks_in(x, pairs) else pairs
  draw <- if (identical(type, "networks")) .cn_draw_one_network else
    .cn_draw_one_difference
  grid <- if (combined) .cn_grid(length(items)) else c(1L, 1L)
  old <- graphics::par(mfrow = grid, mar = c(1, 1, 2.5, 1))
  on.exit(graphics::par(old), add = TRUE)
  invisible(lapply(items, function(it) draw(x, it, ...)))
  invisible(NULL)
}
