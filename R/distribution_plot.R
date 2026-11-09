# ==============================================================================
# distribution_plot() - state-distribution plot over time (stacked area / bar).
# Companion to sequence_plot(). Optional grouping + ncol/nrow facet grid.
# ==============================================================================

#' State Distribution Plot Over Time
#'
#' @description
#' `distribution_plot()` draws how the states are distributed at each time
#' point: for every position column it counts the sequences in each state
#' and stacks the shares (or counts) as areas (default) or bars. It takes the
#' same inputs as [sequence_plot()] (wide sequences, networks built from
#' sequences, clusterings, mixture fits and `tna` models), draws one panel
#' per group or cluster, and is the figure `sequence_plot(type =
#' "distribution")` draws.
#'
#' The plot shows the composition of the sample at each time point, not how
#' individual sequences move: the same bands can come from sequences that
#' stay in one state or from sequences that change state at every step. Look
#' at single sequences with [sequence_plot()] and at the transitions with
#' [build_tna()]. Missing cells form their own band by default (`na`), so
#' sequences that end early stay visible instead of inflating the other
#' shares.
#'
#' @param x The sequences: a wide data frame or matrix (one row per
#'   sequence, one column per time point), a `netobject` built from
#'   sequences, a `netobject_group`, a `net_clustering`, a `net_mmm` or a
#'   `tna` model. A clustering, a mixture fit or a `netobject_group` gives
#'   one panel per cluster or group. A one-hot or association network, and
#'   a frame or matrix holding only 0 and 1 (unless `cols` declares a
#'   frame's columns as positions), raise an error of class
#'   `nestimate_bad_input`.
#' @param group A vector of one group label per sequence (length
#'   `nrow(x)`), drawing one panel per group in the order of its factor
#'   levels (sorted for a character or numeric vector). `NULL` (default)
#'   uses the cluster assignments carried by a `net_clustering`, `net_mmm`
#'   or `netobject_group`; a `group` given overrides them. Panels are
#'   titled `"Cluster <label>"`, whatever the grouping.
#' @param scale What the y-axis shows: `"proportion"` (default) divides
#'   each time point by its total, so the bands fill 0 to 1; `"count"`
#'   keeps the counts, and the axis reaches the tallest stack.
#' @param geom How the bands are drawn: `"area"` (default) as stacked
#'   polygons, `"bar"` as stacked bars.
#' @param na Whether missing cells are drawn as an extra band coloured
#'   `na_color`. Default `TRUE`. With `FALSE` they are left out, and the
#'   proportions are shares of the observed states at each time point.
#' @param trim Cuts the time axis (one number) or the time axis and the
#'   rows (two numbers); `NA` on either leaves that axis alone. Default
#'   `NULL` (no cut). The first number: a fraction in (0, 1) keeps the
#'   busiest time points that together hold that share of the observed
#'   states, wherever they sit, and the kept ones close up under their own
#'   labels (`trim = 0.95`); a value of 1 or more keeps the first that many
#'   time points (`trim = 50`). The second number ends each sequence earlier
#'   before the proportions are taken: a count keeps its first that many
#'   observed states, a fraction its first share of them (`trim = c(0.9,
#'   0.9)`); no sequence is removed.
#' @param trim_clusterwise For grouped plots with a fractional `trim`:
#'   `FALSE` (default) chooses one set of time points over all sequences, so
#'   every panel has the same time axis; `TRUE` chooses them within each
#'   group, so panels can differ in width.
#' @param state_colors Colours for the state bands: an unnamed vector, one
#'   colour per state in sorted state order, or a named vector such as
#'   `c(plan = "#0072B2")`, where only the named states change and the rest
#'   keep the Okabe-Ito defaults. Names this plot does not draw are dropped
#'   with a message, so one palette can serve several figures. `NULL`
#'   (default) uses the palette stored on `x` with [set_state_colors()], or
#'   the Okabe-Ito palette.
#' @param na_color The colour of the missing-cell band. Default `"grey90"`.
#' @param frame Whether to draw a box around each panel, and axis lines.
#'   Default `FALSE`.
#' @param width,height Size of a new graphics device, in inches, opened
#'   only in an interactive session (the other dimension defaults to 7).
#'   Default `NULL` (draw on the current device).
#' @param main The plot title. With several panels on one figure it is
#'   drawn above them; with several panels and `combined = FALSE` it is not
#'   drawn. Default `NULL`.
#' @param show_n Whether to add the number of sequences, `"(n = N)"`, to the
#'   title of each panel. Default `TRUE`.
#' @param time_label The x-axis label. Default `"Time"`; `NULL` or `""`
#'   draws none.
#' @param xlab Another name for `time_label`; when given, it wins.
#' @param y_label The y-axis label. Default `NULL`: `"Proportion"` or
#'   `"Count"`, following `scale`.
#' @param ylab Another name for `y_label`; when given, it wins.
#' @param tick Label every `tick`-th time point. Default `NULL`, which
#'   thins the labels to at most about 15.
#' @param ncol,nrow The panel grid. Default `NULL`: with neither given,
#'   `ncol = ceiling(sqrt(G))` and `nrow = ceiling(G / ncol)` for `G`
#'   panels; with one given, the other follows. With `combined = FALSE`
#'   the grid does not affect the layout, but `ncol * nrow` must still be
#'   at least the number of panels or an error is raised.
#' @param cols For a wide data frame: the columns holding the time points.
#'   Needed when they are not a numbered run (`T1`, `T2`, ...), when one of
#'   them is called `action`, or when any other column sits beside them.
#'   Given for any other input, it raises an error of class
#'   `nestimate_bad_input`. Default `NULL`.
#' @param combined Whether to draw all panels on one figure. `TRUE`
#'   (default) arranges them with [graphics::layout()]; `FALSE` draws each
#'   panel on its own page with its own legend, for knitr (`fig.show =
#'   "asis"`) or for saving each panel to its own file. A single panel
#'   ignores it.
#' @param legend Where the legend goes: `"right"` (default), `"bottom"` or
#'   `"none"`.
#' @param legend_size The legend text size. Default `NULL`, which scales it
#'   with the device width, between 0.65 and 1.2.
#' @param legend_title A title for the legend. Default `NULL` (none).
#' @param legend_ncol The number of legend columns. Default `NULL`: one
#'   column on the right, one row at the bottom.
#' @param legend_border The colour of the swatch borders. Default `NA`
#'   (none).
#' @param legend_bty The legend box: `"n"` (default) draws none, `"o"`
#'   draws one.
#'
#' @return Invisibly, a list describing the drawn figure:
#'   \describe{
#'     \item{counts}{A named list, one entry per panel, each a state by time
#'       point matrix of counts; rows are named by `levels`, columns are the
#'       time points kept.}
#'     \item{proportions}{The same matrices with each column divided by its
#'       total (an empty time point stays 0).}
#'     \item{levels}{The state labels in plotting order, with `"NA"`
#'       appended when `na = TRUE`.}
#'     \item{palette}{The fill colours, one per entry of `levels`.}
#'     \item{groups}{The panel labels (`"all"` when ungrouped), in the order
#'       of `counts` and `proportions`.}
#'   }
#' @seealso \code{\link{sequence_plot}}, \code{\link{build_clusters}}
#' @examples
#' distribution_plot(as.data.frame(trajectories))
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
                              legend_bty     = "n",
                              cols           = NULL) {

  scale  <- match.arg(scale)
  geom   <- match.arg(geom)
  legend <- match.arg(legend)
  if (!is.null(xlab)) time_label <- xlab
  if (!is.null(ylab)) y_label    <- ylab
  .refuse_onehot_network(x, "distribution_plot()")
  .refuse_binary_matrix(x, "distribution_plot()", takes_cols = TRUE)
  if (is.data.frame(x)) {
    .check_declared_columns(x, "distribution_plot()", cols = cols)
    if (is.null(cols)) {
      .refuse_undeclared_roles(x, "distribution_plot()",
                               wide = "`cols`", whole_frame = TRUE,
                               takes_action = FALSE)
    } else {
      x <- x[cols]
    }
  } else if (!is.null(cols)) {
    stop(errorCondition("`cols` applies to a wide data frame.",
                        class = "nestimate_bad_input", call = NULL))
  }

  .assert(
    "`na` must be TRUE or FALSE" =
      is.logical(na) && length(na) == 1L && !is.na(na),
    "`frame` must be TRUE or FALSE" =
      is.logical(frame) && length(frame) == 1L && !is.na(frame),
    "`combined` must be TRUE or FALSE" =
      is.logical(combined) && length(combined) == 1L && !is.na(combined),
    "`trim_clusterwise` must be TRUE or FALSE" =
      is.logical(trim_clusterwise) && length(trim_clusterwise) == 1L &&
      !is.na(trim_clusterwise)
  )

  if (interactive() && (!is.null(width) || !is.null(height))) {
    grDevices::dev.new(width  = if (!is.null(width))  width  else 7,
                       height = if (!is.null(height)) height else 7,
                       noRStudioGD = TRUE)
  }

  # After any dev.new(), pick the device-sized adaptive legend_size.
  legend_size <- .auto_legend_size(legend_size)
  .assert("`legend_size` must be a single positive number" =
            is.numeric(legend_size) && length(legend_size) == 1L &&
            legend_size > 0)

  # A palette attached with set_state_colors() is the object's own default;
  # an explicit `state_colors` still wins for this one figure.
  state_colors <- state_colors %||% .stored_state_colors(x)

  # Extract data and group from various input types. Panels are "clusters"
  # only when the grouping is the input's own clustering, not a `group` given.
  group_given <- !is.null(group)
  extracted <- .extract_seqplot_input(x, group)
  x <- extracted$data
  group <- extracted$group

  if (!is.null(group)) {
    .assert("`group` must have one value per sequence (row of `x`)" =
              length(group) == nrow(x))
    group <- as.factor(group)
  }

  enc        <- .encode_states(x)
  full_codes <- enc$codes
  levels_all <- enc$levels
  K_core     <- length(levels_all)
  K          <- K_core + (if (na) 1L else 0L)
  # One pooled set of columns keeps every panel the same width (aligned axes);
  # trim_clusterwise judges the columns within each group. The row
  # half of `trim` ends each sequence earlier, counted over the kept time
  # points, before the proportions are taken.
  tr         <- .trim_axes(trim)
  global_keep <- .trim_keep(full_codes, tr$cols)
  groups <- if (is.null(group)) factor(rep("all", nrow(full_codes))) else group

  .report_unused_colors(state_colors, levels_all)
  palette       <- .state_palette(state_colors, K_core, levels_all)
  full_palette  <- if (na) c(palette, na_color) else palette
  legend_labels <- if (na) c(levels_all, "NA") else levels_all

  group_levels <- levels(groups)
  G <- length(group_levels)

  count_list <- lapply(group_levels, function(g) {
    sub <- full_codes[groups == g, , drop = FALSE]
    keep <- if (isTRUE(trim_clusterwise)) .trim_keep(sub, tr$cols) else
      global_keep
    sub <- .trim_rows_blank(.trim_columns(sub, keep), tr$rows)
    tab <- vapply(seq_len(ncol(sub)), function(t) {
      col <- sub[, t]
      out <- tabulate(col, nbins = K_core)
      if (na) out <- c(out, sum(is.na(col)))
      out
    }, numeric(K))
    dimnames(tab) <- list(legend_labels, colnames(sub))
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
  y_max <- .distribution_y_max(plot_mats, scale)
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
  .assert("`ncol` * `nrow` must be at least the number of groups" =
            ncol * nrow >= G)

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
    ticks <- .tick_positions(nc, tick, colnames(mat))
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

#' Upper y limit of a distribution plot
#'
#' Areas and bars are stacked, so a count axis must reach the tallest stack,
#' not the largest single cell.
#' @param plot_mats List of state-by-position matrices.
#' @param scale `"proportion"` or `"count"`.
#' @return A single number.
#' @noRd
.distribution_y_max <- function(plot_mats, scale) {
  if (scale == "proportion") return(1)
  max(vapply(plot_mats, \(m) max(colSums(m)), numeric(1L)))
}
