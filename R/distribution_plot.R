# ==============================================================================
# distribution_plot() - state-distribution plot over time (stacked area / bar).
# Companion to sequence_plot(). Optional grouping + ncol/nrow facet grid.
# ==============================================================================

#' State Distribution Plot Over Time
#'
#' Draws how state proportions (or counts) evolve across time points. For
#' each time column, tabulates how many sequences are in each state and
#' renders the result as a stacked area (default) or stacked bar chart.
#' Accepts the same inputs as \code{\link{sequence_plot}}.
#'
#' @param x Wide-format sequence data. Accepts the same inputs as
#'   \code{\link{sequence_plot}}: \code{data.frame}, \code{matrix},
#'   \code{netobject}, \code{net_clustering}, \code{netobject_group},
#'   \code{net_mmm}, or \code{tna}. When clustering info is available,
#'   one panel is drawn per cluster.
#' @param group Optional grouping vector (length \code{nrow(x)}) producing
#'   one panel per group. Ignored if \code{x} is a \code{net_clustering}.
#' @param scale \code{"proportion"} (default) divides each column by its
#'   total so bands fill 0..1. \code{"count"} keeps raw counts.
#' @param geom \code{"area"} (default) draws stacked polygons;
#'   \code{"bar"} draws stacked bars.
#' @param na If \code{TRUE} (default), \code{NA} cells are shown as an
#'   extra band coloured \code{na_color}.
#' @param trim Optional time-axis truncation, to stop a few long
#'   sequences from stretching the plot. \code{NULL} (default) keeps the
#'   full width. A fraction in \code{(0, 1)} drops everything past that
#'   quantile of sequence lengths (e.g. \code{trim = 0.95}); a value
#'   \code{>= 1} is an absolute cut (\code{trim = 50} keeps the first 50
#'   time points). See \code{\link{sequence_plot}}.
#' @param trim_clusterwise Grouped plots only, fractional \code{trim}
#'   only. \code{FALSE} (default) uses one pooled cutoff for every panel
#'   so the time axes stay aligned; \code{TRUE} crops each group to its
#'   own length quantile (panels can differ in width). See
#'   \code{\link{sequence_plot}}.
#' @param state_colors Vector of colours, one per state. Defaults to
#'   Okabe-Ito.
#' @param na_color Colour for the \code{NA} band.
#' @param frame If \code{TRUE} (default), draw a box around each panel.
#' @param width,height Optional device dimensions. See
#'   \code{\link{sequence_plot}}.
#' @param main Plot title.
#' @param show_n Append \code{"(n = N)"} (per-group when grouped) to the
#'   title.
#' @param time_label X-axis label.
#' @param xlab Alias for \code{time_label}.
#' @param y_label Y-axis label. Defaults to \code{"Proportion"} or
#'   \code{"Count"} based on \code{scale}.
#' @param ylab Alias for \code{y_label}.
#' @param tick Show every Nth x-axis label. \code{NULL} = auto.
#' @param ncol,nrow Facet grid dimensions. \code{NULL} = auto:
#'   \code{ncol = ceiling(sqrt(G))}, \code{nrow = ceiling(G / ncol)}.
#'   Ignored when \code{combined = FALSE}.
#' @param combined When \code{TRUE} (default), groups are arranged on one
#'   figure via \code{graphics::layout()}. When \code{FALSE}, each group
#'   is drawn on its own page (one full-size figure per group, with its
#'   own legend). Useful when you want each group at full size in knitr
#'   (\code{fig.show = "asis"}) or to save each as a separate file.
#'   Single-group calls (\code{G == 1}) ignore this argument.
#' @param legend Legend position: \code{"right"} (default),
#'   \code{"bottom"}, or \code{"none"}.
#' @param legend_size Legend text size. \code{NULL} (default) auto-scales
#'   from device width (clamped to \code{[0.65, 1.2]}).
#' @param legend_title Optional legend title.
#' @param legend_ncol Number of legend columns.
#' @param legend_border Swatch border colour.
#' @param legend_bty \code{"n"} (borderless) or \code{"o"} (boxed).
#'
#' @return Invisibly, a list with \code{counts}, \code{proportions},
#'   \code{levels}, \code{palette}, and \code{groups}.
#' @seealso \code{\link{sequence_plot}}, \code{\link{build_clusters}}
#' @examples
#' \donttest{
#' distribution_plot(as.data.frame(trajectories))
#' }
#' @export
distribution_plot <- function(x,
                              group          = NULL,
                              scale          = c("proportion", "count"),
                              geom           = c("area", "bar"),
                              na             = TRUE,
                              trim           = NULL,
                              trim_clusterwise = FALSE,
                              state_colors   = NULL,
                              na_color       = "grey90",
                              frame          = FALSE,
                              width          = NULL,
                              height         = NULL,
                              main           = NULL,
                              show_n         = TRUE,
                              time_label     = "Time",
                              xlab           = NULL,
                              y_label        = NULL,
                              ylab           = NULL,
                              tick           = NULL,
                              ncol           = NULL,
                              nrow           = NULL,
                              combined       = TRUE,
                              legend         = c("right", "bottom", "none"),
                              legend_size    = NULL,
                              legend_title   = NULL,
                              legend_ncol    = NULL,
                              legend_border  = NA,
                              legend_bty     = "n") {

  scale  <- match.arg(scale)
  geom   <- match.arg(geom)
  legend <- match.arg(legend)
  if (!is.null(xlab)) time_label <- xlab
  if (!is.null(ylab)) y_label    <- ylab

  stopifnot(
    is.logical(na),       length(na)       == 1L,
    is.logical(frame),    length(frame)    == 1L,
    is.logical(combined), length(combined) == 1L,
    is.logical(trim_clusterwise), length(trim_clusterwise) == 1L
  )

  if (interactive() && (!is.null(width) || !is.null(height))) {
    grDevices::dev.new(width  = if (!is.null(width))  width  else 7,
                       height = if (!is.null(height)) height else 7,
                       noRStudioGD = TRUE)
  }

  # After any dev.new(), pick the device-sized adaptive legend_size.
  legend_size <- .auto_legend_size(legend_size)
  stopifnot(is.numeric(legend_size), length(legend_size) == 1L,
            legend_size > 0)

  # Extract data and group from various input types
  extracted <- .extract_seqplot_input(x, group)
  x <- extracted$data
  group <- extracted$group

  if (!is.null(group)) {
    stopifnot(length(group) == nrow(x))
    group <- as.factor(group)
  }

  enc        <- .encode_states(x)
  full_codes <- enc$codes
  levels_all <- enc$levels
  col_names  <- colnames(full_codes)
  K_core     <- length(levels_all)
  K          <- K_core + (if (na) 1L else 0L)
  # One pooled cut keeps every panel the same width (aligned axes);
  # trim_clusterwise crops each group to its own length quantile.
  global_cut <- .trim_cut(full_codes, trim)

  palette       <- .state_palette(state_colors, K_core)
  full_palette  <- if (na) c(palette, na_color) else palette
  legend_labels <- if (na) c(levels_all, "NA") else levels_all

  # Factor-rep so grouping dispatch works identically with or without group.
  groups <- if (is.null(group)) factor(rep("all", nrow(full_codes))) else group
  group_levels <- levels(groups)
  G <- length(group_levels)

  count_list <- lapply(group_levels, function(g) {
    sub <- full_codes[groups == g, , drop = FALSE]
    cut <- if (isTRUE(trim_clusterwise)) .trim_cut(sub, trim) else global_cut
    sub <- sub[, seq_len(cut), drop = FALSE]
    tab <- vapply(seq_len(cut), function(t) {
      col <- sub[, t]
      out <- tabulate(col, nbins = K_core)
      if (na) out <- c(out, sum(is.na(col)))
      out
    }, numeric(K))
    rownames(tab) <- legend_labels
    tab
  })
  names(count_list) <- group_levels
  group_sizes <- vapply(group_levels,
                        function(g) sum(groups == g), integer(1))

  prop_list <- lapply(count_list, function(m) {
    cs <- colSums(m)
    cs[cs == 0] <- 1
    sweep(m, 2, cs, "/")
  })
  plot_mats <- if (scale == "proportion") prop_list else count_list
  y_max <- if (scale == "proportion") 1 else max(vapply(plot_mats, max, 0))
  if (is.null(y_label)) {
    y_label <- if (scale == "proportion") "Proportion" else "Count"
  }

  if (is.null(ncol) && is.null(nrow)) {
    ncol <- as.integer(ceiling(sqrt(G)))
    nrow <- as.integer(ceiling(G / ncol))
  } else if (is.null(ncol)) {
    ncol <- as.integer(ceiling(G / nrow))
  } else if (is.null(nrow)) {
    nrow <- as.integer(ceiling(G / ncol))
  }
  stopifnot(ncol * nrow >= G)

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  oma <- .legend_oma_size(legend_labels, legend, legend_size,
                          legend_ncol, legend_title)
  graphics::par(oma = c(if (legend == "bottom") oma[["oma_b"]] else 0.3,
                        0.3, 0.3,
                        if (legend == "right")  oma[["oma_r"]] else 0.3))

  use_layout <- combined && nrow * ncol > 1L
  if (use_layout) {
    layout_mat <- matrix(c(seq_len(G), integer(nrow * ncol - G)),
                         nrow = nrow, ncol = ncol, byrow = TRUE)
    graphics::layout(layout_mat)
  }

  mar_bottom <- if (!is.null(time_label) && nzchar(time_label)) 3 else 1
  mar_top    <- if (!is.null(main) || isTRUE(show_n)) 2.5 else 1.5

  invisible(lapply(seq_len(G), function(g_idx) {
    g   <- group_levels[g_idx]
    mat <- plot_mats[[g_idx]]
    nc  <- ncol(mat)
    ticks <- .tick_positions(nc, tick, col_names[seq_len(nc)])
    # When combined = FALSE, restore oma per page so each group's legend
    # has room (graphics::layout() preserves oma for the whole page; with
    # no layout, every plot.new() starts a fresh page that resets oma).
    if (!combined && g_idx > 1L) {
      graphics::par(oma = c(if (legend == "bottom") oma[["oma_b"]] else 0.3,
                            0.3, 0.3,
                            if (legend == "right")  oma[["oma_r"]] else 0.3))
    }
    graphics::par(mar = c(mar_bottom, 4, mar_top, 1))
    graphics::plot.new()
    graphics::plot.window(xlim = c(0.5, nc + 0.5),
                          ylim = c(0, y_max), xaxs = "i", yaxs = "i")

    cum <- rbind(0, apply(mat, 2, cumsum))
    if (geom == "area") {
      x_poly <- c(0.5, seq_len(nc), nc + 0.5)
      invisible(lapply(seq_len(K), function(k) {
        y_bot <- c(cum[k,     1L], cum[k,     ], cum[k,     nc])
        y_top <- c(cum[k + 1, 1L], cum[k + 1, ], cum[k + 1, nc])
        graphics::polygon(c(x_poly, rev(x_poly)),
                          c(y_bot, rev(y_top)),
                          col = full_palette[k], border = NA)
      }))
    } else {
      ks <- rep(seq_len(K),  times = nc)
      ts <- rep(seq_len(nc), each  = K)
      graphics::rect(xleft   = ts - 0.5,
                     ybottom = cum[cbind(ks,     ts)],
                     xright  = ts + 0.5,
                     ytop    = cum[cbind(ks + 1, ts)],
                     col     = full_palette[ks],
                     border  = NA)
    }
    if (isTRUE(frame)) graphics::box()
    graphics::axis(1, at = ticks$at, labels = ticks$labels,
                   mgp = c(1, 0.5, 0),
                   lwd = if (isTRUE(frame)) 1 else 0, lwd.ticks = 1)
    graphics::axis(2, las = 1, mgp = c(1, 0.5, 0),
                   lwd = if (isTRUE(frame)) 1 else 0, lwd.ticks = 1)
    if (!is.null(time_label) && nzchar(time_label)) {
      graphics::mtext(time_label, side = 1, line = 1.7)
    }
    graphics::mtext(y_label, side = 2, line = 2.5)

    panel_title <- if (G > 1L) {
      base <- sprintf("Cluster %s", g)
      if (isTRUE(show_n)) sprintf("%s  (n = %d)", base,
                                  group_sizes[g_idx]) else base
    } else if (!is.null(main)) {
      if (isTRUE(show_n)) sprintf("%s  (n = %d)", main,
                                  group_sizes[1L]) else main
    } else if (isTRUE(show_n)) sprintf("n = %d", group_sizes[1L]) else NULL
    if (!is.null(panel_title)) {
      graphics::mtext(panel_title, side = 3, line = 0.5, font = 2,
                      cex = if (G > 1L) 0.9 else 1)
    }
    # In combined = FALSE mode, draw a legend on each page (layout() is not
    # used, so the post-loop draw would only land on the last page).
    if (!combined && legend != "none") {
      .draw_legend_in_oma(legend_labels, full_palette, legend, legend_size,
                          legend_ncol, legend_title, legend_border, legend_bty)
    }
  }))
  if (combined && G > 1L && !is.null(main)) {
    graphics::mtext(main, side = 3, line = 1.5, font = 2,
                    outer = TRUE, cex = 1)
  }

  if (combined && legend != "none") {
    .draw_legend_in_oma(legend_labels, full_palette, legend, legend_size,
                        legend_ncol, legend_title, legend_border, legend_bty)
  }

  invisible(list(counts      = count_list,
                 proportions = prop_list,
                 levels      = legend_labels,
                 palette     = full_palette,
                 groups      = group_levels))
}
