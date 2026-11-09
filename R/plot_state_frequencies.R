# ==============================================================================
# plot_state_frequencies: marimekko + colored-bar state-frequency plotter
#
# Native Nestimate plotter for state (node) frequency distributions across
# the group-bearing classes (netobject, netobject_group, mcml, htna).
#
# Distinct from tna::plot_frequencies (different name, different geometry).
# Defaults to a marimekko (mosaic) layout where column widths reflect
# group totals and segment heights reflect within-group state proportions.
# ==============================================================================


# ------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------

# Drop-in safe state_frequencies() that ignores NA / void markers. Returns
# a data.frame(state, count, proportion) sorted by count descending.
.freq_count_states <- function(values) {
  values <- as.character(values)
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0L) {
    return(data.frame(state = character(0), count = integer(0),
                      proportion = numeric(0), stringsAsFactors = FALSE))
  }
  tbl <- sort(table(values), decreasing = TRUE)
  data.frame(state = names(tbl),
             count = as.integer(tbl),
             proportion = as.numeric(tbl) / sum(tbl),
             stringsAsFactors = FALSE)
}


# Tidy a single (group, data) pair into a freq data.frame. Returns NULL when
# the group has no countable states.
.freq_one_group <- function(group_label, data) {
  if (is.null(data)) return(NULL)
  if (is.data.frame(data) || is.matrix(data)) {
    values <- as.vector(as.matrix(data))
  } else if (is.list(data)) {
    values <- unlist(data, use.names = FALSE)
  } else {
    values <- data
  }
  fr <- .freq_count_states(values)
  if (nrow(fr) == 0L) return(NULL)
  fr$group <- as.character(group_label)
  fr[, c("group", "state", "count", "proportion")]
}


# ------------------------------------------------------------------------------
# Per-class extractors. Each returns data.frame(group, state, count, proportion).
# Public generic state_distribution() at the bottom of the file forwards here.
# ------------------------------------------------------------------------------

.freq_df_netobject <- function(x) {
  .refuse_onehot_network(x, "state_distribution()")
  if (is.null(x$data)) {
    stop(errorCondition("netobject has no `$data` to compute frequencies from.",
                        class = "nestimate_bad_input", call = NULL))
  }

  # If $node_groups present, group states by node_groups$group; otherwise
  # treat the whole data as a single group ("all").
  ng <- x$node_groups
  if (is.null(ng) || nrow(ng) == 0L) {
    return(.freq_one_group("all", x$data))
  }

  # Compute frequency over $data, then attach group via node->group lookup.
  fr <- .freq_count_states(as.vector(as.matrix(x$data)))
  if (nrow(fr) == 0L) return(fr[, integer(0)])
  lookup <- setNames(as.character(ng$group), as.character(ng$node))
  fr$group <- unname(lookup[fr$state])
  fr$group[is.na(fr$group)] <- "(ungrouped)"
  # Re-aggregate per (group, state). State is unique per row already because
  # .freq_count_states returns one row per unique state, but a state could
  # belong to one group only -- keep as-is. Recompute within-group proportion.
  totals <- ave(fr$count, fr$group, FUN = sum)
  fr$proportion <- fr$count / totals
  fr[, c("group", "state", "count", "proportion")]
}


.freq_df_htna <- function(x) {
  # htna inherits from netobject; node_groups is always populated by build_htna.
  .freq_df_netobject(x)
}


.freq_df_mcml <- function(x, include_macro = FALSE) {
  if (is.null(x$clusters)) {
    stop(errorCondition("mcml object has no `$clusters` field.",
                        class = "nestimate_bad_input", call = NULL))
  }
  cluster_names <- names(x$clusters)
  parts <- lapply(cluster_names, function(cn) {
    cl <- x$clusters[[cn]]
    .freq_one_group(cn, cl$data)
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L) {
    stop(errorCondition(
      "mcml clusters have no `$data` to compute frequencies from.",
      class = "nestimate_bad_input", call = NULL))
  }
  out <- do.call(rbind, parts)

  # Apply the labels remap stored by `.apply_node_labels` if present, so
  # raw codes in $data (e.g. "c_evaluation") are shown as their human-
  # readable labels (e.g. "Evaluation").
  map <- attr(x, "labels_map")
  if (!is.null(map)) {
    hits <- out$state %in% names(map)
    out$state[hits] <- unname(map[out$state[hits]])
  }

  if (isTRUE(include_macro)) {
    macro <- aggregate(out$count, by = list(state = out$state), FUN = sum)
    macro$group <- "macro"
    macro$proportion <- macro$x / sum(macro$x)
    macro <- data.frame(group = macro$group, state = macro$state,
                        count = macro$x, proportion = macro$proportion,
                        stringsAsFactors = FALSE)
    out <- rbind(macro, out)
  }
  out
}


.freq_df_netobject_group <- function(x) {
  .refuse_onehot_network(x, "state_distribution()")
  group_names <- names(x)
  if (is.null(group_names)) group_names <- paste0("Group ", seq_along(x))
  parts <- lapply(seq_along(x), function(i) {
    .freq_one_group(group_names[i], x[[i]]$data)
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L) {
    stop(errorCondition("netobject_group has no `$data` in any constituent.",
                        class = "nestimate_bad_input", call = NULL))
  }
  do.call(rbind, parts)
}


# ------------------------------------------------------------------------------
# State ordering and color assignment
# ------------------------------------------------------------------------------

.order_states <- function(states, counts, sort_states) {
  if (sort_states == "none") return(unique(states))
  if (sort_states == "alpha") return(sort(unique(states)))
  agg <- aggregate(counts, by = list(state = states), FUN = sum)
  agg <- agg[order(-agg$x), , drop = FALSE]
  agg$state
}


# ------------------------------------------------------------------------------
# Mosaic of a tidy data frame (internal; reached through mosaic_plot())
# ------------------------------------------------------------------------------

# Marimekko geometry: column widths proportional to the per-column total of
# `weight`; within each column, segments stack to height 1 in proportion to
# each row's share of that column's total. Called by mosaic_plot.data.frame().
.mosaic_df_plot <- function(data,
                            x,
                            y,
                            weight,
                            fill        = "y",
                            colors      = NULL,
                            show_labels = TRUE,
                            label_size  = 3.5,
                            x_label     = NULL,
                            y_label     = NULL) {
  .assert(is.data.frame(data),
            is.character(x), length(x) == 1L,
            is.character(y), length(y) == 1L,
            is.character(weight), length(weight) == 1L)
  needed <- c(x, y, weight)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0L) {
    stop(errorCondition(paste0("mosaic_plot: missing columns: ",
         paste(missing_cols, collapse = ", ")),
         class = "nestimate_bad_input", call = NULL))
  }

  # Build cross-tab. Column = x, row = y. Order respects existing factor
  # levels if present; otherwise alphabetical.
  xv <- if (is.factor(data[[x]])) data[[x]] else factor(data[[x]])
  yv <- if (is.factor(data[[y]])) data[[y]] else factor(data[[y]])
  tab <- tapply(data[[weight]], list(xv, yv), FUN = sum, default = 0)
  tab[is.na(tab)] <- 0
  x_levels <- rownames(tab)
  y_levels <- colnames(tab)

  # Column widths from row sums, segment heights from within-row proportions.
  col_totals <- rowSums(tab)
  total      <- sum(col_totals)
  if (total <= 0) {
    stop(errorCondition("mosaic_plot: total weight is zero -- nothing to draw.",
                        class = "nestimate_bad_input", call = NULL))
  }
  widths_cum <- c(0, cumsum(col_totals)) / total

  # Build rect coordinates: one row per (x_level, y_level) cell.
  rects <- do.call(rbind, lapply(seq_along(x_levels), function(i) {
    row    <- tab[i, , drop = TRUE]
    rs     <- sum(row)
    if (rs <= 0) return(NULL)
    heights_cum <- c(0, cumsum(row / rs))
    data.frame(
      x_level = rep(x_levels[i], length(y_levels)),
      y_level = y_levels,
      xmin    = widths_cum[i],
      xmax    = widths_cum[i + 1L],
      ymin    = heights_cum[seq_len(length(y_levels))],
      ymax    = heights_cum[seq_len(length(y_levels)) + 1L],
      count   = as.numeric(row),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
  rects <- rects[rects$ymax > rects$ymin, , drop = FALSE]

  # A `fill` that names no column would only fail when ggplot builds the plot.
  if (!identical(fill, "y") &&
      !(is.character(fill) && length(fill) == 1L && fill %in% names(data))) {
    stop(errorCondition(paste0(
      "mosaic_plot: `fill` must be \"y\" or the name of one column of the ",
      "data; got ", paste(format(fill), collapse = ", "), "."),
      class = "nestimate_bad_input", call = NULL))
  }

  # Attach extra fill column if requested
  if (!identical(fill, "y") && fill %in% names(data)) {
    extra <- aggregate(data[[fill]],
                       by = list(x_level = as.character(xv),
                                 y_level = as.character(yv)),
                       FUN = function(z) z[1L])
    names(extra)[3L] <- fill
    rects <- merge(rects, extra,
                   by = c("x_level", "y_level"), all.x = TRUE)
  }

  # Map fill aesthetic
  fill_var <- if (identical(fill, "y")) "y_level" else fill
  rects[[fill_var]] <- if (identical(fill, "y")) {
    factor(rects$y_level, levels = y_levels)
  } else {
    rects[[fill_var]]
  }

  # Build palette
  pal <- if (identical(fill, "y")) {
    .state_palette(colors, length(y_levels), y_levels)
  } else {
    NULL
  }

  # Column midpoint x positions for axis ticks
  mid_x <- (widths_cum[-1L] + widths_cum[-length(widths_cum)]) / 2

  p <- ggplot2::ggplot(rects,
                       ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                    ymin = .data$ymin, ymax = .data$ymax,
                                    fill = .data[[fill_var]])) +
    ggplot2::geom_rect(color = "white", linewidth = 0.3) +
    ggplot2::scale_x_continuous(breaks = mid_x, labels = x_levels,
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0),
                                labels = function(v) sprintf("%d%%", round(v * 100))) +
    ggplot2::labs(x = x_label %||% x, y = y_label %||% "Proportion",
                  fill = if (identical(fill, "y")) y else fill) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(angle = 30, hjust = 1)
    )

  if (!is.null(pal)) {
    p <- p + ggplot2::scale_fill_manual(values = pal, drop = FALSE)
  }

  if (isTRUE(show_labels)) {
    rects$pct <- rects$ymax - rects$ymin
    rects$x_mid <- (rects$xmin + rects$xmax) / 2
    rects$y_mid <- (rects$ymin + rects$ymax) / 2
    rects$lab   <- sprintf("%.1f%%", 100 * rects$pct)
    rects$lab[rects$pct < 0.04] <- ""
    fill_cols <- if (is.null(pal)) NULL else {
      pal[match(as.character(rects[[fill_var]]), y_levels)]
    }
    p <- p + .geom_fit_label(rects, label_size, fill_colors = fill_cols)
  }
  p
}


# Tiny null-coalesce helper (avoid rlang dep)
`%||%` <- function(a, b) if (is.null(a)) b else a


# ------------------------------------------------------------------------------
# mosaic_plot: tna-equivalent chi-square mosaic for netobjects
# (.geom_fit_label moved to R/plot-utils.R alongside other shared helpers)
# ------------------------------------------------------------------------------

#' Mosaic Plot of a Network's Transition or Co-occurrence Counts
#'
#' @description
#' `mosaic_plot()` draws a mosaic display of a contingency table and fills
#' each tile with its standardised residual under independence: blue tiles
#' hold more cases than independence predicts, red tiles fewer. It reads the
#' transition or co-occurrence counts of a `netobject` or `htna` network, the
#' state counts of each group of a `netobject_group`, the cluster-level or
#' per-cluster counts of an `mcml`, or a `table` or count `matrix` given
#' directly. A data frame in long form is drawn as a plain marimekko of
#' shares instead, coloured by the stacked variable (or the `fill`
#' column), without residuals.
#'
#' In a network's mosaic each column is a target state, its width
#' proportional to the transitions into it, and its segments are the source
#' states, stacked in proportion to their share of those transitions. For a
#' `netobject_group` the columns are the groups and the segments the states,
#' counted from each group's sequences, so the mosaic shows how the state
#' mix differs between groups. For a `table` or `matrix`, its rows are the
#' columns of the mosaic and its columns the stacked segments. An `htna`
#' network has its states ordered by actor, with each state's label set in
#' bold in its actor's colour.
#'
#' The residuals describe where the counts depart from independence; they
#' are not tests of single cells, and transitions taken from the same
#' sequence are not independent draws. To test edges or group differences,
#' use [permutation()]; for a chi-square or Fisher test with
#' Cram\enc{é}{e}r's V on two variables of a data frame, use
#' [mosaic_analysis()].
#'
#' @details
#' A mosaic needs counts. A network whose weights are whole numbers
#' (`method = "frequency"` or `"co_occurrence"`) is drawn from its weights;
#' otherwise from the transition counts stored with it at estimation, and
#' failing that from first-order transitions recounted from its sequence
#' data. A `netobject_group` is counted from each group's sequences, or
#' from the column sums of its weights when they are counts. A single
#' network with none of these, such as an association network (`cor`,
#' `pcor`, `glasso`, `ising`, `mgm`), raises an error. A `netobject_group`
#' must hold sequence networks: the group method counts whatever its
#' networks' data hold, so a group of association networks is drawn
#' without an error, each distinct numeric value counted as a state. An
#' `mcml` is counted the same way as a network, layer by layer; a cluster
#' with a single state has no table to draw: when it holds self-transitions
#' it is left out of the per-cluster figure with a warning, and when it holds
#' no transition at all the call stops with an error. When a row or
#' column of the table is empty, the asymptotic residuals are undefined
#' and every tile is drawn at 0. The geometry and the asymptotic
#' residuals follow `tna::plot_mosaic()`; the function has its own name so
#' the two packages can be attached together.
#'
#' @param x The object to draw: a `netobject`, `netobject_group`, `mcml` or
#'   `htna` network, a contingency `table`, a numeric count `matrix`, or a
#'   data frame in long form (see `x_var`, `y_var` and `weight`).
#' @param x_var,y_var,weight For a data frame: the column whose levels form
#'   the mosaic's columns (each as wide as its total `weight`), the column
#'   whose levels are stacked within each column (in proportion to their
#'   share), and the column holding each row's weight, such as a count.
#'   Levels keep their factor order, or sort alphabetically. `x_var`,
#'   `y_var`, `weight`, `fill`, `colors` and `show_labels` given with any
#'   other input raise an error of class `nestimate_bad_input`, as does
#'   `label_size` outside `style = "flat"`.
#' @param fill For a data frame: `"y"` (default) colours the segments by
#'   `y_var`; the name of another column maps that column to the fill
#'   instead (the first value per cell).
#' @param colors For a data frame with `fill = "y"`: fill colours, either an
#'   unnamed vector in level order (at least one per level) or a named
#'   lookup such as `c(plan = "#0072B2")` that overrides only the levels
#'   named. Default `NULL`, the Okabe-Ito palette, recycled.
#' @param show_labels,label_size For a data frame: whether to print each
#'   segment's percentage inside it (default `TRUE`; segments under 4% stay
#'   unlabelled), and the label size (default `3.5`). With
#'   `style = "flat"`, `label_size` also sizes the flat mosaic's tile
#'   labels.
#' @param level For an `mcml`: `"macro"` (default) draws one mosaic of the
#'   cluster-by-cluster counts; `"clusters"` draws one mosaic per cluster,
#'   as panels of one plot that share the residual scale.
#' @param xlab,ylab Axis titles. The network and `mcml` methods draw none
#'   by default (`NULL`); the `table` and `matrix` methods default to
#'   `"Row"` and `"Column"`; a data frame defaults to the `x_var` name and
#'   `"Proportion"`. Not drawn with `style = "flat"`.
#' @param range The lower and upper limits of the residual colour scale, two
#'   increasing numbers; residuals beyond them take the end colours. Default
#'   `NULL`, which fits the symmetric range `c(-M, M)`, with `M` the largest
#'   absolute residual (at least 1), so no tile is clipped. `c(-4, 4)`
#'   matches tna's display.
#' @param top_angle,left_angle Rotation, in degrees, of the column labels on
#'   top and the segment labels on the left. Default `NULL`: 90 when an axis
#'   has more than three levels, 0 otherwise. Not used with `style = "flat"`
#'   or for the per-cluster `mcml` panels.
#' @param residuals How the standardised residuals are computed.
#'   `"permutation"` (default) shuffles one variable's labels against the
#'   other `n_perm` times, which keeps both margins, and gives each cell
#'   `(observed - mean) / sd` over the shuffles; a cell that never varies
#'   gets 0. It does not rely on large expected counts, so it suits sparse
#'   tables. `"asymptotic"` gives the adjusted residuals of
#'   [stats::chisq.test()], `(O - E) / sqrt(E (1 - p_row) (1 - p_col))`, as
#'   vcd and tna do.
#' @param n_perm The number of shuffles for `residuals = "permutation"`, at
#'   least 2. Default `500`; use 1000 or more for stable values in the
#'   tails.
#' @param seed A seed for the shuffles, so the residuals repeat exactly. The
#'   session's own random stream is left as it was. Default `NULL`. Ignored with
#'   `residuals = "asymptotic"`.
#' @param ncol For an `mcml` with `level = "clusters"`: the number of panel
#'   columns. Default `2`. Accepted by the `netobject_group` method, whose
#'   single mosaic does not use it.
#' @param values Whether to print each tile's residual (one decimal) inside
#'   it. Default `FALSE`. In the classic style the text is white on
#'   strongly coloured tiles (absolute residual above 1.5) and dark grey
#'   otherwise; with `style = "flat"` it is white on tiles whose absolute
#'   residual exceeds 0.55 of the largest, and tiles too small for a label
#'   stay blank.
#' @param style `"classic"` (default) draws black-bordered tiles with axis
#'   labels. `"flat"` draws the mosaic of [mosaic_analysis()]: white gutters
#'   between tiles, labels beside them, no axes, and a ColorBrewer red-blue
#'   residual ramp; `values = TRUE` then prints the residuals in the tiles,
#'   and its styling arguments (`tile_label`, `col_label_side`,
#'   `row_label_side`, `legend_position`, `legend_size`, `label_size`, ...)
#'   go in `...`. An `mcml` with `level = "clusters"` and more than one
#'   panel is drawn classic, with a message; an `htna` drawn flat loses its
#'   actor-coloured labels, with a warning.
#' @param ... With `style = "flat"`: styling arguments for the flat mosaic
#'   of a network, table or matrix; an unknown name raises an error that
#'   lists the valid ones. The data-frame method takes no `style` and
#'   raises an error of class `nestimate_bad_input` for any argument it
#'   does not take; the `mcml` method does the same with
#'   `style = "classic"`; the other methods ignore it with a warning of
#'   class `nestimate_ignored_argument`.
#'
#' @return A `ggplot` object, one rectangle per cell (empty cells are left
#'   out of the flat and data-frame mosaics), filled by the standardised
#'   residual (by category for a data frame). Print it to draw
#'   it, or add ggplot2 layers and themes. An `mcml` with
#'   `level = "clusters"` gives one faceted `ggplot` with a panel per
#'   cluster and one shared fill scale. An object of an unsupported class
#'   raises an error.
#' @references
#' Friendly, M. (1994). Mosaic displays for multi-way contingency tables.
#' \emph{Journal of the American Statistical Association}, 89(425),
#' 190--200. \doi{10.1080/01621459.1994.10476460}
#' @export
#' @examples
#' data(group_regulation_long, package = "Nestimate")
#' net <- build_network(group_regulation_long, method = "frequency",
#'                      format = "long", actor = "Actor", action = "Action",
#'                      order = "Time")
#' mosaic_plot(net, seed = 1)
#'
#' # A tidy data frame: columns by group, segments by state
#' df <- data.frame(
#'   group = rep(c("A", "B", "C"), each = 3),
#'   state = rep(c("s1", "s2", "s3"), 3),
#'   count = c(10, 5, 2,  7, 8, 3,  4, 6, 12)
#' )
#' mosaic_plot(df, x_var = "group", y_var = "state", weight = "count")
mosaic_plot <- function(x, ...) UseMethod("mosaic_plot")

#' @export
#' @rdname mosaic_plot
mosaic_plot.default <- function(x, ...) {
  stop(errorCondition(paste0("mosaic_plot: no method for class ",
       paste(class(x), collapse = "/"), ".\n",
       "Supported: netobject, netobject_group, mcml, htna, table, matrix, ",
       "data.frame."),
       class = "nestimate_bad_input", call = NULL))
}

#' @export
#' @rdname mosaic_plot
mosaic_plot.netobject <- function(x,
                                  xlab = NULL,
                                  ylab = NULL,
                                  range = NULL,
                                  top_angle = NULL,
                                  left_angle = NULL,
                                  residuals = c("permutation", "asymptotic"),
                                  n_perm = 500L,
                                  seed = NULL,
                                  values = FALSE,
                                  style = c("classic", "flat"),
                                  ...) {
  .mosaic_plot_impl(x, source_class = "netobject",
                    xlab = xlab, ylab = ylab, range = range,
                    top_angle = top_angle, left_angle = left_angle,
                    residuals = match.arg(residuals),
                    n_perm = n_perm, seed = seed, values = values,
                    style = match.arg(style), ...)
}

#' @export
#' @rdname mosaic_plot
mosaic_plot.htna <- function(x,
                             xlab = NULL,
                             ylab = NULL,
                             range = NULL,
                             top_angle = NULL,
                             left_angle = NULL,
                             residuals = c("permutation", "asymptotic"),
                             n_perm = 500L,
                             seed = NULL,
                             values = FALSE,
                             style = c("classic", "flat"),
                             ...) {
  .mosaic_plot_impl(x, source_class = "htna",
                    xlab = xlab, ylab = ylab, range = range,
                    top_angle = top_angle, left_angle = left_angle,
                    residuals = match.arg(residuals),
                    n_perm = n_perm, seed = seed, values = values,
                    style = match.arg(style), ...)
}

#' @export
#' @rdname mosaic_plot
mosaic_plot.mcml <- function(x,
                             level = c("macro", "clusters"),
                             xlab = NULL,
                             ylab = NULL,
                             range = NULL,
                             top_angle = NULL,
                             left_angle = NULL,
                             residuals = c("permutation", "asymptotic"),
                             n_perm = 500L,
                             seed = NULL,
                             ncol = 2L,
                             values = FALSE,
                             style = c("classic", "flat"),
                             ...) {
  style <- match.arg(style)
  .mosaic_reject_df_args(style, ...)
  # The classic path takes no extra dots, so a stray arg is a typo -> reject.
  # The flat path legitimately forwards styling args to the renderer (validated
  # there), so only guard unused dots when classic.
  if (identical(style, "classic")) {
    .mosaic_plot_check_unused_dots("mosaic_plot.mcml", ...)
  }
  .mosaic_plot_impl(x, source_class = "mcml",
                    level = match.arg(level),
                    xlab = xlab, ylab = ylab, range = range,
                    top_angle = top_angle, left_angle = left_angle,
                    residuals = match.arg(residuals),
                    n_perm = n_perm, seed = seed, ncol = ncol,
                    values = values, style = style, ...)
}

.mosaic_plot_check_unused_dots <- function(method, ...) {
  dots <- list(...)
  if (!length(dots)) {
    return(invisible(TRUE))
  }
  dot_names <- names(dots)
  dot_names[!nzchar(dot_names)] <- paste0("..", which(!nzchar(dot_names)))
  stop(errorCondition(
    paste0(
      method, "() got unsupported argument",
      if (length(dots) == 1L) ": " else "s: ",
      paste(dot_names, collapse = ", ")
    ),
    class = "nestimate_bad_input", call = NULL
  ))
}

# An argument that reaches `...` and is read by nothing is named in a
# warning; the call still returns what it returned without the argument.
.warn_ignored_dots <- function(fn, ...) {
  n_dots <- ...length()
  if (n_dots == 0L) {
    return(invisible(FALSE))
  }
  dot_names <- ...names() %||% rep("", n_dots)
  dot_names[is.na(dot_names)] <- ""
  unnamed <- !nzchar(dot_names)
  dot_names[unnamed] <- paste0("..", which(unnamed))
  warning(warningCondition(
    paste0(fn, " ignored ",
           if (n_dots == 1L) "an argument it does not take: "
           else "arguments it does not take: ",
           paste(dot_names, collapse = ", "), "."),
    class = "nestimate_ignored_argument", call = NULL))
  invisible(TRUE)
}

#' @export
#' @rdname mosaic_plot
mosaic_plot.netobject_group <- function(x,
                                        xlab = NULL,
                                        ylab = NULL,
                                        range = NULL,
                                        top_angle = NULL,
                                        left_angle = NULL,
                                        residuals = c("permutation",
                                                      "asymptotic"),
                                        n_perm = 500L,
                                        seed = NULL,
                                        ncol = 2L,
                                        values = FALSE,
                                        style = c("classic", "flat"),
                                        ...) {
  .mosaic_plot_impl(x, source_class = "netobject_group",
                    xlab = xlab, ylab = ylab, range = range,
                    top_angle = top_angle, left_angle = left_angle,
                    residuals = match.arg(residuals),
                    n_perm = n_perm, seed = seed, ncol = ncol,
                    values = values, style = match.arg(style), ...)
}

#' @export
#' @rdname mosaic_plot
mosaic_plot.table <- function(x,
                              xlab = "Row",
                              ylab = "Column",
                              range = NULL,
                              top_angle = NULL,
                              left_angle = NULL,
                              residuals = c("permutation", "asymptotic"),
                              n_perm = 500L,
                              seed = NULL,
                              values = FALSE,
                              style = c("classic", "flat"),
                              ...) {
  .mosaic_plot_tab(x, xlab = xlab, ylab = ylab, range = range,
                   top_angle = top_angle, left_angle = left_angle,
                   residuals = match.arg(residuals),
                   n_perm = n_perm, seed = seed, values = values,
                   style = match.arg(style), ...)
}

#' @export
#' @rdname mosaic_plot
mosaic_plot.matrix <- function(x, ...) mosaic_plot.table(as.table(x), ...)

#' @export
#' @rdname mosaic_plot
mosaic_plot.data.frame <- function(x,
                                   x_var,
                                   y_var,
                                   weight,
                                   fill = "y",
                                   colors = NULL,
                                   show_labels = TRUE,
                                   label_size = 3.5,
                                   xlab = NULL,
                                   ylab = NULL,
                                   ...) {
  .mosaic_plot_check_unused_dots("mosaic_plot.data.frame", ...)
  .mosaic_df_plot(x, x = x_var, y = y_var, weight = weight, fill = fill,
                  colors = colors, show_labels = show_labels,
                  label_size = label_size, x_label = xlab, y_label = ylab)
}


# ---- Shared worker for the 4 data-bearing classes ----------------------------
#
# Pulls one or more weight matrices out of x according to source_class, runs
# them through .mosaic_plot_tab(), and combines panels for multi-matrix cases.
# Single source of truth for the data-bearing flow; the per-class S3 methods
# are one-line forwarders.

.mosaic_plot_impl <- function(x, source_class,
                              level = "macro",
                              xlab = NULL, ylab = NULL,
                              range = NULL,
                              top_angle = NULL, left_angle = NULL,
                              residuals = "permutation",
                              n_perm = 500L, seed = NULL,
                              ncol = 2L, values = FALSE,
                              style = c("classic", "flat"), ...) {
  style  <- match.arg(style)
  .mosaic_reject_df_args(style, ...)
  panels <- .mosaic_panels(x, source_class, level)
  xlab <- xlab %||% panels$xlab_default
  ylab <- ylab %||% panels$ylab_default

  # mcml level="clusters" is the only multi-panel case: one mosaic per
  # cluster, faceted into a single ggplot with one shared fill legend.
  # netobject_group is rendered as a single (group x state) mosaic, so it
  # falls through to the single-panel path below.
  if (length(panels$tabs) > 1L && identical(source_class, "mcml")) {
    if (identical(style, "flat")) {
      message("mosaic_plot: style = 'flat' is not supported for multi-panel ",
              "mcml (level = 'clusters'); using the classic faceted style.")
    }
    return(.mosaic_facet_plot(panels$tabs, panels$titles,
                              panel_states = panels$panel_states,
                              xlab = xlab, ylab = ylab, range = range,
                              residuals = residuals,
                              n_perm = n_perm, seed = seed, ncol = ncol,
                              values = values))
  }

  plots <- lapply(seq_along(panels$tabs), function(i) {
    ag <- if (is.null(panels$axis_groups)) NULL else panels$axis_groups[[i]]
    p <- .mosaic_plot_tab(panels$tabs[[i]],
                          xlab = xlab, ylab = ylab, range = range,
                          top_angle = top_angle, left_angle = left_angle,
                          residuals = residuals,
                          n_perm = n_perm, seed = seed,
                          axis_groups = ag, values = values,
                          style = style, ...)
    if (!is.null(panels$titles)) p + ggplot2::ggtitle(panels$titles[[i]]) else p
  })

  if (length(plots) == 1L) return(plots[[1L]])
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    return(gridExtra::arrangeGrob(grobs = plots, ncol = ncol))
  }
  warning(warningCondition(
    "Install 'gridExtra' to combine grouped mosaics; returning a list.",
    class = "nestimate_missing_package", call = NULL))
  plots
}


# Faceted multi-panel mosaic: one ggplot with one strip per panel and a
# single shared fill legend. Used when the input is an mcml level="within"
# (k clusters, each its own state-by-state mosaic). Each cluster's state set
# is different so we cannot use a shared axis scale -- labels are drawn
# inside the panel via geom_text at the top edge (state names) and via
# strip text (cluster name).
.mosaic_facet_plot <- function(tabs, titles, panel_states = NULL,
                               xlab, ylab, range,
                               residuals, n_perm, seed, ncol,
                               values = FALSE) {
  panel_dfs <- lapply(seq_along(tabs), function(i) {
    d <- .mosaic_rect_data(tabs[[i]], residuals = residuals,
                           n_perm = n_perm, seed = seed)
    d$panel <- titles[[i]]
    d
  })
  d <- do.call(rbind, panel_dfs)
  d$panel <- factor(d$panel, levels = unlist(titles))

  # Per-panel label filter: if panel_states is supplied (netobject_group with
  # harmonised vocabulary), keep only labels for states actually present in
  # the group; otherwise keep labels for any non-zero-extent column / row.
  keep_for_panel <- function(labels, ext, states_i) {
    if (!is.null(states_i)) labels %in% states_i else ext > 1e-6
  }

  # Per-panel x-axis labels at top (above the mosaic).
  x_lab_df <- do.call(rbind, lapply(seq_along(tabs), function(i) {
    d_i <- panel_dfs[[i]]
    w   <- attr(d_i, "widths")
    rls <- attr(d_i, "row_labels")
    keep <- keep_for_panel(rls, diff(w),
                           if (is.null(panel_states)) NULL else panel_states[[i]])
    if (!any(keep)) return(NULL)
    data.frame(x = ((w[-length(w)] + w[-1L]) / 2)[keep],
               y = 1.04, label = rls[keep], panel = titles[[i]],
               stringsAsFactors = FALSE)
  }))
  x_lab_df$panel <- factor(x_lab_df$panel, levels = unlist(titles))

  # Per-panel y-axis labels at left, aligned to the first non-zero-width
  # column's stacking (col 1 may be zero-width when harmonised, in which
  # case its heights are the fallback c(0, 1/m, 2/m, ...) which doesn't
  # describe anything real).
  y_lab_df <- do.call(rbind, lapply(seq_along(tabs), function(i) {
    d_i  <- panel_dfs[[i]]
    h    <- attr(d_i, "heights")
    cls  <- attr(d_i, "col_labels")
    w_widths <- attr(d_i, "widths")
    m    <- length(cls)

    first_real_col <- which(diff(w_widths) > 1e-6)[1L]
    if (is.na(first_real_col)) return(NULL)
    h_col      <- h[, first_real_col]
    col_offset <- (first_real_col - 1L) * m * 0.0025
    y_pos      <- (h_col[-length(h_col)] + h_col[-1L]) / 2 + col_offset

    keep <- keep_for_panel(cls, diff(h_col),
                           if (is.null(panel_states)) NULL else panel_states[[i]])
    if (!any(keep)) return(NULL)
    data.frame(x = -0.04,
               y = y_pos[keep],
               label = cls[keep],
               panel = titles[[i]],
               stringsAsFactors = FALSE)
  }))
  y_lab_df$panel <- factor(y_lab_df$panel, levels = unlist(titles))

  p <- ggplot2::ggplot(d,
    ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                 ymin = .data$ymin, ymax = .data$ymax,
                 fill = .data$stdres)) +
    ggplot2::geom_rect(colour = "black", linewidth = 0.4,
                       show.legend = TRUE) +
    ggplot2::scale_fill_gradient2(
      name   = "Standardized\nresidual",
      oob    = scales::oob_squish,
      low    = "#D33F6A",
      high   = "#4A6FE3",
      limits = .mosaic_residual_limits(d$stdres, range),
      breaks = .mosaic_residual_breaks(d$stdres, range)
    ) +
    ggplot2::geom_text(data = x_lab_df, inherit.aes = FALSE,
                       ggplot2::aes(x = .data$x, y = .data$y,
                                    label = .data$label),
                       size = 3, angle = 90, hjust = 0, vjust = 0.5) +
    ggplot2::geom_text(data = y_lab_df, inherit.aes = FALSE,
                       ggplot2::aes(x = .data$x, y = .data$y,
                                    label = .data$label),
                       size = 3, hjust = 1, vjust = 0.5) +
    ggplot2::facet_wrap(~ panel, ncol = ncol,
                        strip.position = "bottom") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid     = ggplot2::element_blank(),
      axis.text      = ggplot2::element_blank(),
      axis.ticks     = ggplot2::element_blank(),
      axis.line      = ggplot2::element_blank(),
      panel.spacing  = ggplot2::unit(2.5, "lines"),
      strip.text     = ggplot2::element_text(
        face = "bold", size = 11,
        margin = ggplot2::margin(t = 6, b = 4)
      ),
      strip.placement = "outside",
      plot.margin    = ggplot2::margin(t = 36, r = 12, b = 16, l = 36)
    ) +
    ggplot2::labs(x = xlab, y = ylab)

  if (isTRUE(values)) {
    cent <- d
    cent$xcent <- (cent$xmin + cent$xmax) / 2
    cent$ycent <- (cent$ymin + cent$ymax) / 2
    p <- p + ggplot2::geom_text(
      data = cent, inherit.aes = FALSE,
      ggplot2::aes(x = .data$xcent, y = .data$ycent,
                   label = sprintf("%.1f", .data$stdres),
                   colour = abs(.data$stdres) > 1.5),
      size = 3, show.legend = FALSE
    ) + ggplot2::scale_colour_manual(
      values = c(`TRUE` = "white", `FALSE` = "grey20"), guide = "none"
    )
  }
  p
}


# Extract the contingency tables (and panel titles) for each data-bearing
# class. Returns a list with $tabs (list of `table` objects in the orientation
# expected by .mosaic_plot_tab), $titles (NULL for single-panel cases), and
# $xlab_default / $ylab_default (so each class can name its axes naturally).
.mosaic_panels <- function(x, source_class, level = "macro") {
  if (source_class %in% c("netobject", "htna")) {
    # Order-1 transition mosaic: chi-square needs integer transition counts.
    # If $weights is already integer (method = "frequency" / "co_occurrence")
    # use it directly; otherwise recount transitions from the raw $data
    # sequences so the mosaic works on any estimation method (relative,
    # glasso, cor, ...). Falls through to the strict integer guard only
    # when neither path is available.
    if (isTRUE(x$method %in% .association_methods)) {
      stop(errorCondition(sprintf(paste0(
        "mosaic_plot(): this \"%s\" network holds variables, not ",
        "transitions between states. Build it with a transition method or ",
        "pass the sequences."), x$method),
        class = "nestimate_bad_input", call = NULL))
    }
    w <- .mosaic_count_or_stop(x)
    axis_groups <- NULL
    # htna carries actor membership in $node_groups (data.frame node, group).
    # Reorder rows/cols so each actor block is contiguous, and pass the
    # ordered group vector down to the renderer for annotation strips.
    if (identical(source_class, "htna") &&
        !is.null(x$node_groups) && nrow(x$node_groups) > 0L) {
      nodes <- rownames(w)
      if (!is.null(nodes)) {
        lookup <- setNames(as.character(x$node_groups$group),
                           as.character(x$node_groups$node))
        g <- unname(lookup[nodes])
        g[is.na(g)] <- "(ungrouped)"
        ord <- order(g, nodes)
        w <- w[ord, ord]
        axis_groups <- setNames(g[ord], nodes[ord])
      }
    }
    return(list(tabs = list(as.table(t(w))),
                titles = NULL,
                xlab_default = NULL,
                ylab_default = NULL,
                axis_groups = list(axis_groups)))
  }
  if (identical(source_class, "netobject_group")) {
    # Single mosaic of the joint (group x state) contingency table -- an
    # order-0 marginal: how often each state appears within each group.
    # When raw $data is available (typical for sequence-derived netobjects)
    # we count states from the source, so the test works on any estimation
    # method (relative, glasso, cor, ...) -- the chi-square is on state
    # marginals and never needs integer $weights. Falls back to
    # colSums($weights) with the integer-only guard when $data is absent.
    .refuse_onehot_network(x, "mosaic_plot()")
    group_names <- names(x) %||% paste0("Group ", seq_along(x))
    weight_states <- unlist(lapply(x, function(g) rownames(g$weights)))
    data_states   <- unlist(lapply(x, function(g) {
      if (is.null(g$data)) return(NULL)
      as.character(as.vector(as.matrix(g$data)))
    }))
    data_states <- data_states[!is.na(data_states) & nzchar(data_states)]
    all_states  <- sort(unique(c(weight_states, data_states)))
    cnts <- vapply(seq_along(x), function(i) {
      full <- setNames(numeric(length(all_states)), all_states)
      if (!is.null(x[[i]]$data)) {
        fr <- .freq_count_states(as.vector(as.matrix(x[[i]]$data)))
        full[fr$state] <- fr$count
      } else {
        w  <- .mosaic_weights_or_stop(x[[i]]$weights, x[[i]]$method,
                                      ctx = sprintf("group '%s'",
                                                    group_names[i]))
        cs <- colSums(w)
        full[names(cs)] <- cs
      }
      full
    }, numeric(length(all_states)))

    # vapply returned an (n_states x n_groups) matrix; transpose so rows
    # are groups and pass through .mosaic_plot_tab unchanged (no extra t()).
    m <- t(cnts)
    dimnames(m) <- list(group_names, all_states)
    return(list(tabs = list(as.table(m)),
                titles = NULL,
                xlab_default = NULL,
                ylab_default = NULL))
  }
  if (identical(source_class, "mcml")) {
    # build_mcml() / cluster_summary() return $macro$weights (cluster x
    # cluster aggregate) and $clusters[[k]]$weights (per-cluster matrix).
    # mcml weights may be row-normalised (type = "tna"); the mosaic needs
    # integer transition counts. .mosaic_count_or_stop() recovers them from
    # the layer's $data when the stored weights aren't counts, so the mosaic
    # works regardless of how the mcml was estimated.
    if (identical(level, "macro")) {
      if (is.null(x$macro) || is.null(x$macro$weights)) {
        stop(errorCondition("mosaic_plot.mcml: x$macro$weights is missing.",
                          class = "nestimate_bad_input", call = NULL))
      }
      w <- .mosaic_count_or_stop(
        list(weights = x$macro$weights, data = x$macro$data,
             method = x$meta$type), ctx = "macro")
      return(list(tabs = list(as.table(t(w))),
                  titles = NULL,
                  xlab_default = NULL,
                  ylab_default = NULL))
    }
    if (is.null(x$clusters) || length(x$clusters) == 0L) {
      stop(errorCondition(
        "mosaic_plot.mcml: x$clusters is empty (compute_within = FALSE?).",
        class = "nestimate_bad_input", call = NULL))
    }
    titles <- names(x$clusters) %||%
      paste0("Cluster ", seq_along(x$clusters))
    built <- lapply(seq_along(x$clusters), function(i) {
      cl <- x$clusters[[i]]
      # One state gives no table whether or not it has self-transitions, so
      # it is set aside before the count check (which stops on zero counts).
      if (is.matrix(cl$weights) &&
          (nrow(cl$weights) < 2L || ncol(cl$weights) < 2L)) {
        return(NULL)
      }
      w <- .mosaic_count_or_stop(
        list(weights = cl$weights, data = cl$data, method = x$meta$type),
        ctx = sprintf("cluster[%s]", titles[i]))
      # A single-state cluster has a 1x1 matrix -- no 2-D contingency table
      # to draw -- so it is dropped from the faceted figure (with a note)
      # rather than aborting every other cluster's panel.
      if (nrow(w) < 2L || ncol(w) < 2L) NULL else as.table(t(w))
    })
    keep <- !vapply(built, is.null, logical(1L))
    if (!any(keep)) {
      stop(errorCondition(
        "mosaic_plot.mcml: no cluster has >= 2 states to draw a mosaic.",
        class = "nestimate_bad_input", call = NULL))
    }
    if (!all(keep)) {
      warning(warningCondition(
        paste0("mosaic_plot.mcml: skipping single-state cluster(s): ",
               paste(titles[!keep], collapse = ", ")),
        class = "nestimate_skipped_cluster", call = NULL))
    }
    return(list(tabs = built[keep], titles = titles[keep],
                xlab_default = NULL,
                ylab_default = NULL))
  }
  stop(errorCondition(
    paste0("mosaic_plot: unknown source_class '", source_class, "'."),
    class = "nestimate_bad_input", call = NULL))
}


# Resolve a transition count matrix for a single netobject. Mosaic residuals
# are defined for count tables. If $weights is already integer-valued
# (method = "frequency" / "co_occurrence") use it directly. Otherwise, when
# the raw $data sequences are available, recount order-1 transition counts
# from them so the mosaic works on any estimation method (relative, glasso,
# cor, ...) -- mirroring the documented fallback honored by the
# netobject_group path. Falls through to the strict integer guard only when
# $data is also absent.
.mosaic_count_or_stop <- function(x, ctx = NULL) {
  w <- x$weights
  is_count_mat <- function(m) {
    is.matrix(m) && is.numeric(m) && all(is.finite(m)) &&
      all(m >= 0) && all(abs(m - round(m)) <= 1e-8) && sum(m) > 0
  }
  # 1. $weights are already integer counts (method = "frequency" /
  #    "co_occurrence") -- use them directly.
  if (is_count_mat(w)) return(w)
  # 2. The estimator's own stored transition counts. Every build_network
  #    netobject preserves $frequency_matrix, which reflects the actual
  #    count construction (begin_state/end_state/concat/weighted params,
  #    one-hot WTNA, ...). Prefer it over a naive raw recount so the mosaic
  #    agrees with the network and does not break on one-hot data (Codex
  #    review).
  fm <- x$frequency_matrix
  if (is_count_mat(fm)) return(fm)
  # 3. Last resort: recount plain order-1 transitions from raw $data. This
  #    is what lets a row-normalised network (tna / relative / glasso) still
  #    produce a count-based mosaic -- and is the path mcml clusters take,
  #    since their stored $weights are probabilities under type = "tna".
  if (!is.null(x$data) &&
      (is.data.frame(x$data) || is.matrix(x$data))) {
    alphabet <- rownames(w) %||% colnames(w)
    counts <- .count_transitions(as.data.frame(x$data), format = "wide",
                                 id = NULL, alphabet = alphabet)
    return(.mosaic_weights_or_stop(counts, "frequency", ctx))
  }
  .mosaic_weights_or_stop(w, x$method, ctx)
}

# Validate that a weight matrix is mosaic-able (finite, non-negative, integer-
# valued, non-zero total). `ctx` is an optional sub-label used in the error
# message so multi-panel calls report which group / cluster failed.
.mosaic_weights_or_stop <- function(w, method, ctx = NULL) {
  pre <- if (is.null(ctx)) "mosaic_plot: " else paste0("mosaic_plot [", ctx, "]: ")
  if (!is.matrix(w) || !is.numeric(w) || !all(is.finite(w))) {
    stop(errorCondition(
      paste0(pre, "$weights must be a finite numeric matrix."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (any(w < 0) || any(abs(w - round(w)) > 1e-8)) {
    stop(errorCondition(
      paste0(
        pre,
        "only defined for integer-valued weight matrices ",
        "(method = 'frequency' or 'co_occurrence'). Got method = ",
        method %||% "<unknown>",
        "."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (sum(w) <= 0) {
    stop(errorCondition(paste0(pre, "total weight is zero -- nothing to draw."),
                        class = "nestimate_bad_input", call = NULL))
  }
  w
}

# Choose symmetric colour-scale limits: explicit `range` if supplied, else
# auto-fit to the actual stdres range so no signal is squished. Floors at
# +/-1 to keep the legend readable on near-independent tables.
.mosaic_residual_limits <- function(stdres, range) {
  if (!is.null(range)) {
    .assert(is.numeric(range), length(range) == 2L, range[1L] < range[2L])
    return(range)
  }
  m <- max(abs(stdres), na.rm = TRUE)
  m <- max(m, 1)
  c(-m, m)
}

# Pick at most 5 breaks symmetric around 0, rounded to a tidy step (1, 2, 4,
# 5, or 10) so the colour bar is scannable at any data scale.
.mosaic_residual_breaks <- function(stdres, range) {
  lim <- .mosaic_residual_limits(stdres, range)
  m <- lim[2L]
  step <- if      (m <= 4)  1
          else if (m <= 8)  2
          else if (m <= 16) 4
          else if (m <= 25) 5
          else              10
  s <- seq(0, m, by = step)
  unique(c(-rev(s), s))
}

# Build a y-axis scale for the mosaic. All labels shown at the cell
# midpoint of the leftmost column.
.mosaic_y_scale <- function(d, col_labels, ...) {
  breaks_all <- d$ycent[d$xmin == 0]
  o <- order(breaks_all)
  ggplot2::scale_y_continuous(
    breaks = breaks_all[o],
    labels = col_labels[o],
    expand = c(0.01, 0)
  )
}

# Internal: vectorized port of tna::plot_mosaic_(). Builds a marimekko data
# frame from a contingency table and renders it with chi-square stdres fill.
# `tab` rows are the x-axis categories; columns are the within-column stack.
# Permutation-based standardized residuals. Shuffles one variable's labels
# against the other (preserves both marginals under the independence null),
# tabulates n_perm times, and returns a per-cell empirical z-score
# (O - mean_perm) / sd_perm. Vectorized: each iteration is a single
# tabulate() over a column-major linear index, stacked into a (n*m) x B
# matrix; row means and row sds in closed form.
.mosaic_perm_stdres <- function(tab, n_perm = 500L, seed = NULL) {
  .assert(n_perm >= 2L)
  if (!is.null(seed)) .local_seed(seed)
  n <- nrow(tab); m <- ncol(tab)
  cells <- expand.grid(row_i = seq_len(n), col_j = seq_len(m))
  freqs <- as.vector(tab)
  long_row <- rep(cells$row_i, times = freqs)
  long_col <- rep(cells$col_j, times = freqs)

  perm_counts <- vapply(seq_len(n_perm), function(b) {
    pc <- sample.int(length(long_col))
    tabulate(long_row + (long_col[pc] - 1L) * n, nbins = n * m)
  }, numeric(n * m))

  mean_perm <- rowMeans(perm_counts)
  sd_perm   <- sqrt(rowMeans((perm_counts - mean_perm)^2) *
                      n_perm / (n_perm - 1L))
  obs <- as.vector(tab)
  z <- (obs - mean_perm) / sd_perm
  z[!is.finite(z)] <- 0
  matrix(z, n, m, dimnames = dimnames(tab))
}

.mosaic_rect_data <- function(tab, residuals = "permutation",
                              n_perm = 500L, seed = NULL) {
  n <- nrow(tab); m <- ncol(tab)
  if (n < 1L || m < 1L) {
    stop(errorCondition(
      "mosaic_plot: contingency table must have >= 1 row and column.",
      class = "nestimate_bad_input", call = NULL))
  }
  rs <- rowSums(tab)
  total <- sum(rs)
  if (total <= 0) {
    stop(errorCondition("mosaic_plot: total weight is zero -- nothing to draw.",
         class = "nestimate_bad_input", call = NULL))
  }
  widths <- c(0, cumsum(rs)) / total
  heights <- vapply(seq_len(n), function(i) {
    if (rs[i] <= 0) return(c(0, seq_len(m) / m))
    c(0, cumsum(as.numeric(tab[i, ]) / rs[i]))
  }, numeric(m + 1L))

  i_idx <- rep(seq_len(n), each = m)
  j_idx <- rep(seq_len(m), times = n)
  row_offset <- (i_idx - 1L) * n * 0.0025
  col_offset <- (j_idx - 1L) * m * 0.0025

  row_labels <- rownames(tab) %||% as.character(seq_len(n))
  col_labels <- colnames(tab) %||% as.character(seq_len(m))

  stdres <- if (identical(residuals, "permutation")) {
    .mosaic_perm_stdres(tab, n_perm = n_perm, seed = seed)
  } else {
    suppressWarnings(stats::chisq.test(tab))$stdres
  }
  if (is.null(stdres)) {
    stdres <- matrix(0, n, m)
  }
  # A cell of an empty row or column has no expected count, so no residual:
  # draw it as 0 and keep every other cell's residual.
  stdres[!is.finite(stdres)] <- 0

  d <- data.frame(
    xmin   = widths[i_idx] + row_offset,
    xmax   = widths[i_idx + 1L] + row_offset,
    ymin   = heights[cbind(j_idx, i_idx)] + col_offset,
    ymax   = heights[cbind(j_idx + 1L, i_idx)] + col_offset,
    freq   = as.numeric(tab[cbind(i_idx, j_idx)]),
    row    = row_labels[i_idx],
    col    = col_labels[j_idx],
    stdres = as.numeric(stdres[cbind(i_idx, j_idx)]),
    stringsAsFactors = FALSE
  )
  d$xcent <- (d$xmin + d$xmax) / 2
  d$ycent <- (d$ymin + d$ymax) / 2
  attr(d, "widths")     <- widths
  attr(d, "heights")    <- heights
  attr(d, "row_labels") <- row_labels
  attr(d, "col_labels") <- col_labels
  # Keep the cell-aligned residual matrix so the flat renderer can index it
  # positionally instead of rebuilding geometry and rejoining by string key.
  attr(d, "stdres")     <- stdres
  d
}


.mosaic_plot_tab <- function(tab, xlab, ylab, range = NULL,
                             top_angle = NULL, left_angle = NULL,
                             residuals = "permutation", n_perm = 500L,
                             seed = NULL, axis_groups = NULL,
                             values = FALSE, style = "classic", ...) {
  .mosaic_reject_df_args(style, ...)
  d <- .mosaic_rect_data(tab, residuals = residuals,
                         n_perm = n_perm, seed = seed)
  # "flat" routes to the ported Saqrmisc renderer: white gutters, theme_void,
  # in-tile labels. `values = TRUE` maps to in-tile standardized residuals,
  # unless the caller passes an explicit `tile_label` via `...`.
  if (identical(style, "flat")) {
    if (!is.null(axis_groups)) {
      warning(warningCondition(paste0(
        "mosaic_plot: group-coloured axis strips are not drawn in ",
        "style = 'flat'; use style = 'classic' for the htna group ",
        "annotation."),
        class = "nestimate_ignored_argument", call = NULL))
    }
    dots <- list(...)
    if (is.null(dots$tile_label)) {
      dots$tile_label <- if (isTRUE(values)) "residual" else "none"
    }
    .mosaic_flat_check_args(names(dots))
    return(do.call(.mosaic_flat_draw,
                   c(list(res = attr(d, "stdres"), tab = tab, range = range),
                     dots)))
  }
  # Classic style reads nothing from `...`.
  .warn_ignored_dots("mosaic_plot(style = \"classic\")", ...)
  widths     <- attr(d, "widths")
  heights    <- attr(d, "heights")
  row_labels <- attr(d, "row_labels")
  col_labels <- attr(d, "col_labels")
  n <- nrow(tab); m <- ncol(tab)

  top_a  <- top_angle  %||% (if (n > 3) 90 else 0)
  left_a <- left_angle %||% (if (m > 3) 90 else 0)

  p <- ggplot2::ggplot(d,
    ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                 ymin = .data$ymin, ymax = .data$ymax,
                 fill = .data$stdres)) +
    ggplot2::geom_rect(color = "black", linewidth = 0.4,
                       show.legend = TRUE) +
    ggplot2::scale_fill_gradient2(
      name   = "Standardized\nresidual",
      oob    = scales::oob_squish,
      low    = "#D33F6A",
      high   = "#4A6FE3",
      limits = .mosaic_residual_limits(d$stdres, range),
      breaks = .mosaic_residual_breaks(d$stdres, range)
    ) +
    ggplot2::scale_x_continuous(
      breaks   = unique(d$xcent),
      labels   = row_labels,
      position = "top",
      expand   = c(0.01, 0)
    ) +
    .mosaic_y_scale(d, col_labels, left_angle = left_a) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      panel.grid    = ggplot2::element_blank(),
      axis.ticks    = ggplot2::element_blank(),
      axis.line     = ggplot2::element_blank(),
      axis.text.x   = ggplot2::element_text(
        angle = top_a,
        hjust = if (top_a == 0) 0.5 else 0,
        vjust = if (top_a == 0) 0   else 0.5
      ),
      axis.text.y   = ggplot2::element_text(
        angle = left_a,
        hjust = if (left_a == 0) 1   else 0.5,
        vjust = if (left_a == 0) 0.4 else 0.5
      )
    ) +
    ggplot2::labs(x = xlab, y = ylab)

  if (!is.null(axis_groups)) {
    p <- .mosaic_add_group_strips(p, tab, widths, heights,
                                  axis_groups = axis_groups,
                                  top_a = top_a, left_a = left_a)
  }
  if (isTRUE(values)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(x = .data$xcent, y = .data$ycent,
                   label = sprintf("%.1f", .data$stdres),
                   colour = abs(.data$stdres) > 1.5),
      inherit.aes = FALSE, size = 3, show.legend = FALSE
    ) + ggplot2::scale_colour_manual(
      values = c(`TRUE` = "white", `FALSE` = "grey20"), guide = "none"
    )
  }
  p
}


# Validate flat-styling argument names against the renderer's formals, so an
# unknown / mistyped arg (or a mosaic_plot-only arg like `values`) fails with a
# clear, listed message instead of an opaque base-R "unused argument" error.
# The arguments of mosaic_plot.data.frame() mean nothing for a network,
# table or matrix; classic style would otherwise ignore them silently.
# `label_size` also sizes the tile labels of the flat renderer.
.mosaic_reject_df_args <- function(style, ...) {
  df_args <- c("x_var", "y_var", "weight", "fill", "colors", "show_labels",
               if (!identical(style, "flat")) "label_size")
  given <- intersect(names(list(...)), df_args)
  if (length(given)) {
    stop(errorCondition(paste0(
      "`", paste(given, collapse = "`, `"), "` ",
      if (length(given) == 1L) "applies" else "apply",
      " to data-frame input only; a network, table or matrix is drawn ",
      "from its own counts."),
      class = "nestimate_bad_input", call = NULL))
  }
  invisible(TRUE)
}

.mosaic_flat_check_args <- function(nms) {
  nms  <- nms[nzchar(nms)]
  valid <- setdiff(names(formals(.mosaic_flat_draw)), c("res", "tab"))
  bad  <- setdiff(nms, valid)
  if (length(bad)) {
    stop(errorCondition(
      paste0(
        "mosaic flat style: unsupported styling argument",
        if (length(bad) > 1L) "s" else "",
        ": ",
        paste(bad, collapse = ", "),
        ".\n  Valid styling arguments: ",
        paste(valid, collapse = ", "),
        "."),
      class = "nestimate_bad_input", call = NULL))
  }
  invisible(TRUE)
}

# ColorBrewer RdBu (11): red = negative residual, blue = positive. The flat
# mosaic's shared default palette (matches the Saqrmisc style).
.mosaic_rdbu <- c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7",
                  "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac",
                  "#053061")

# Flat ggplot mosaic (variable-width columns, white gutters, theme_void),
# ported from the Saqrmisc style: RdBu residual ramp (`palette`). `res` is the
# cell-aligned standardized-residual matrix (same shape/order as `tab`); `tab`
# supplies the counts and marginals. Reused by mosaic_plot(style = "flat") and
# mosaic_analysis(). Geometry is a clean [0,1] grid (no inter-tile offsets --
# the white gutters supply the separation), so side labels align and the
# residual fill is indexed positionally (no string-key join, no NA on misses).
.mosaic_flat_draw <- function(res, tab, title = "",
                              tile_label = c("none", "count", "percent",
                                             "residual", "category"),
                              pct_base = c("total", "row", "column"),
                              col_label_side = c("top", "bottom", "both", "none"),
                              row_label_side = c("left", "right", "both", "none"),
                              col_label_angle = NULL, row_label_angle = NULL,
                              show_legend = TRUE, legend_position = "right",
                              legend_size = 0.7,
                              legend_title = "Standardized\nresidual",
                              label_size = 3.5, range = NULL, palette = NULL,
                              min_label_h = 0.04, min_label_w = 0.03) {
  if (is.null(palette)) palette <- .mosaic_rdbu
  tile_label     <- match.arg(tile_label)
  pct_base       <- match.arg(pct_base)
  col_label_side <- match.arg(col_label_side)
  row_label_side <- match.arg(row_label_side)
  .assert("residual matrix must match the contingency table shape" =
            identical(dim(res), dim(tab)))

  rn <- rownames(tab); cn <- colnames(tab)
  R  <- nrow(tab); C <- ncol(tab)
  grand   <- sum(tab)
  row_tot <- rowSums(tab)            # first-variable marginal (plot columns)
  col_tot <- colSums(tab)            # second-variable marginal (stacked rows)

  # Column width = first-variable marginal; tile height = within-column
  # conditional. Residuals are read straight from `res` by position.
  xw     <- as.numeric(row_tot) / grand
  xright <- cumsum(xw); xleft <- xright - xw

  tiles <- do.call(rbind, lapply(seq_len(R), function(i) {
    h <- as.numeric(tab[i, ]) / row_tot[i]
    data.frame(
      xmin = xleft[i], xmax = xright[i],
      ymin = 1 - cumsum(h), ymax = 1 - c(0, cumsum(h)[-C]),
      freq = as.numeric(tab[i, ]),
      stdres = as.numeric(res[i, ]),
      row = rn[i], col = cn, stringsAsFactors = FALSE)
  }))
  tiles <- tiles[tiles$freq > 0, , drop = FALSE]
  tiles$xcent <- (tiles$xmin + tiles$xmax) / 2
  tiles$ycent <- (tiles$ymin + tiles$ymax) / 2
  tiles$w <- tiles$xmax - tiles$xmin
  tiles$h <- tiles$ymax - tiles$ymin

  pct <- switch(pct_base,
    total  = tiles$freq / grand * 100,
    row    = tiles$freq / as.numeric(row_tot[tiles$row]) * 100,
    column = tiles$freq / as.numeric(col_tot[tiles$col]) * 100)

  raw_txt <- switch(tile_label,
    none     = rep("", nrow(tiles)),
    count    = formatC(tiles$freq, format = "d", big.mark = ","),
    percent  = paste0(formatC(pct, format = "f", digits = 1), "%"),
    residual = formatC(tiles$stdres, format = "f", digits = 1),
    category = tiles$col)

  max_abs <- max(abs(tiles$stdres), 1e-6, na.rm = TRUE)
  tiles$txt_col <- ifelse(!is.na(tiles$stdres) & abs(tiles$stdres) > 0.55 * max_abs,
                          "white", "grey15")
  tiles$lab <- ifelse(tiles$h >= min_label_h & tiles$w >= min_label_w,
                      raw_txt, "")

  # Auto-rotate top/bottom labels to vertical when there are many (narrow)
  # columns, so wide category names cannot overlap. Caller can override.
  col_angle <- if (is.null(col_label_angle)) (if (R > 6L) 90 else 0) else
    col_label_angle
  row_angle <- if (is.null(row_label_angle)) 0 else row_label_angle
  col_vert  <- abs(col_angle) >= 45

  # --- side category-label frames (on the clean [0,1] grid) ---
  col_centers <- (xleft + xright) / 2
  col_keep    <- xw >= min_label_w
  col_sides <- switch(col_label_side, top = "top", bottom = "bottom",
                      both = c("top", "bottom"), none = character(0))
  col_df <- if (length(col_sides) && any(col_keep))
    do.call(rbind, lapply(col_sides, function(s) {
      top <- s == "top"
      data.frame(
        x   = col_centers[col_keep],
        y   = if (top) 1 else 0,
        lab = rn[col_keep],
        vj  = if (col_vert) 0.5 else (if (top) -0.6 else 1.6),
        hj  = if (col_vert) (if (top) 0 else 1) else 0.5,
        stringsAsFactors = FALSE)
    })) else NULL

  pj          <- as.numeric(col_tot) / grand
  row_centers <- 1 - (cumsum(pj) - pj / 2)
  row_keep    <- pj >= min_label_h
  row_sides <- switch(row_label_side, left = "left", right = "right",
                      both = c("left", "right"), none = character(0))
  row_df <- if (length(row_sides) && any(row_keep))
    do.call(rbind, lapply(row_sides, function(s) data.frame(
      y = row_centers[row_keep], x = if (s == "left") 0 else 1,
      lab = cn[row_keep], hj = if (s == "left") 1.15 else -0.15,
      stringsAsFactors = FALSE))) else NULL

  fill_guide <- if (show_legend && legend_position != "none") {
    ggplot2::guide_colourbar(
      barheight = ggplot2::unit(legend_size * 6, "lines"),
      barwidth  = ggplot2::unit(legend_size * 0.8, "lines"),
      ticks.colour = "grey40", frame.colour = NA)
  } else "none"
  leg_pos <- if (show_legend) legend_position else "none"

  p <- ggplot2::ggplot(tiles,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                   ymin = .data$ymin, ymax = .data$ymax)) +
    ggplot2::geom_rect(ggplot2::aes(fill = .data$stdres),
                       colour = "white", linewidth = 1.1) +
    ggplot2::geom_text(
      ggplot2::aes(x = .data$xcent, y = .data$ycent, label = .data$lab,
                   colour = .data$txt_col),
      hjust = 0.5, vjust = 0.5, size = label_size)

  if (!is.null(col_df))
    p <- p + ggplot2::geom_text(data = col_df, inherit.aes = FALSE,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
      vjust = col_df$vj, hjust = col_df$hj, angle = col_angle,
      size = label_size + 0.1, colour = "grey35")

  if (!is.null(row_df))
    p <- p + ggplot2::geom_text(data = row_df, inherit.aes = FALSE,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
      hjust = row_df$hj, vjust = 0.5, angle = row_angle,
      size = label_size + 0.1, colour = "grey45")

  # Reserve margin only on sides that actually carry labels. Vertical column
  # labels need more room than horizontal ones so they clear the tiles.
  col_e   <- if (col_vert) 0.16 else 0.08
  left_e  <- if (row_label_side %in% c("left", "both"))   0.12 else 0.015
  right_e <- if (row_label_side %in% c("right", "both"))  0.12 else 0.015
  top_e   <- if (col_label_side %in% c("top", "both"))    col_e else 0.015
  bot_e   <- if (col_label_side %in% c("bottom", "both")) col_e else 0.015

  p +
    ggplot2::scale_fill_gradientn(
      name = legend_title, colours = palette, oob = scales::oob_squish,
      limits = .mosaic_residual_limits(as.numeric(res), range),
      breaks = .mosaic_residual_breaks(as.numeric(res), range), guide = fill_guide) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(left_e, right_e))) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(bot_e, top_e))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 13) +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = "#fafafa", colour = NA),
      panel.background = ggplot2::element_rect(fill = "#fafafa", colour = NA),
      plot.title = ggplot2::element_text(face = "bold", size = 15, hjust = 0,
                                         margin = ggplot2::margin(b = 12, l = 4)),
      plot.margin = ggplot2::margin(14, 18, 14, 14),
      legend.position = leg_pos,
      legend.title = ggplot2::element_text(size = legend_size * 11),
      legend.text  = ggplot2::element_text(size = legend_size * 10))
}


# Group-aware axis-label colouring for htna (or any class that supplies an
# axis_groups vector). Each contiguous block of rows / columns gets its own
# Okabe-Ito hue, applied to both x-axis (top) and y-axis (left) tick labels;
# bold weight reinforces the demarcation. Nothing else changes -- the panel
# itself is the unaltered mosaic. Reorder is done upstream in .mosaic_panels()
# so adjacent ticks of the same group sit next to each other.
.mosaic_add_group_strips <- function(p, tab, widths, heights,
                                     axis_groups, top_a, left_a) {
  n <- nrow(tab); m <- ncol(tab)
  if (length(axis_groups) != n || length(axis_groups) != m) {
    return(p)
  }
  group_levels <- unique(axis_groups)
  pal <- .okabe_ito[seq_along(group_levels)]
  names(pal) <- group_levels

  x_label_cols <- unname(pal[axis_groups[seq_len(n)]])
  y_label_cols <- unname(pal[axis_groups[seq_len(m)]])

  p + ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      colour = x_label_cols,
      angle  = top_a,
      hjust  = if (top_a == 0) 0.5 else 0,
      vjust  = if (top_a == 0) 0   else 0.5,
      face   = "bold"
    ),
    axis.text.y = ggplot2::element_text(
      colour = y_label_cols,
      angle  = left_a,
      hjust  = if (left_a == 0) 1   else 0.5,
      vjust  = if (left_a == 0) 0.4 else 0.5,
      face   = "bold"
    )
  )
}


# ------------------------------------------------------------------------------
# Internal renderers
# ------------------------------------------------------------------------------

# Hierarchical 2D marimekko: column widths = group total, segments stacked
# vertically with heights = within-group state proportions. Used for mcml,
# where states partition cleanly into clusters so each rectangle's area
# carries unique information.
.plot_marimekko_hierarchical <- function(freq_df, sort_states, colors,
                                          label, label_size,
                                          legend = "bottom",
                                          legend_dir = "auto",
                                          legend_frame = "none") {
  state_levels <- .order_states(freq_df$state, freq_df$count, sort_states)
  freq_df$state <- factor(freq_df$state, levels = state_levels)
  group_levels <- unique(freq_df$group)
  freq_df$group <- factor(freq_df$group, levels = group_levels)

  pal <- .state_palette(colors, length(state_levels), state_levels)
  names(pal) <- state_levels

  # Build cumulative-width / cumulative-height rectangle coordinates.
  tab <- tapply(freq_df$count,
                list(freq_df$group, freq_df$state),
                FUN = sum, default = 0)
  tab[is.na(tab)] <- 0
  group_levels <- rownames(tab)
  state_levels <- colnames(tab)
  col_totals <- rowSums(tab)
  total <- sum(col_totals)
  widths_cum <- c(0, cumsum(col_totals)) / total

  rects <- do.call(rbind, lapply(seq_along(group_levels), function(i) {
    row <- tab[i, , drop = TRUE]
    rs  <- sum(row)
    if (rs <= 0) return(NULL)
    heights_cum <- c(0, cumsum(row / rs))
    data.frame(
      group = rep(group_levels[i], length(state_levels)),
      state = state_levels,
      xmin  = widths_cum[i],
      xmax  = widths_cum[i + 1L],
      ymin  = heights_cum[seq_len(length(state_levels))],
      ymax  = heights_cum[seq_len(length(state_levels)) + 1L],
      count = as.numeric(row),
      proportion = as.numeric(row) / rs,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
  rects <- rects[rects$ymax > rects$ymin, , drop = FALSE]
  rects$state <- factor(rects$state, levels = state_levels)
  rects$x_mid <- (rects$xmin + rects$xmax) / 2
  rects$y_mid <- (rects$ymin + rects$ymax) / 2
  rects$tile_height <- rects$ymax - rects$ymin

  mid_x <- (widths_cum[-1L] + widths_cum[-length(widths_cum)]) / 2

  legend_ncol <- min(length(state_levels), 5L)

  p <- ggplot2::ggplot(rects,
                       ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                    ymin = .data$ymin, ymax = .data$ymax,
                                    fill = .data$state)) +
    ggplot2::geom_rect(color = "white", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = mid_x, labels = group_levels,
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0),
                                labels = function(v) sprintf("%d%%", round(v * 100))) +
    ggplot2::labs(x = NULL, y = "Within-cluster proportion", fill = "State") +
    .legend_layer(legend, length(state_levels), legend_dir) +
    ggplot2::theme_minimal(base_size = 12) +
    .legend_theme(legend, legend_dir, legend_frame) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(angle = 30, hjust = 1)
    )

  if (label != "none") {
    rects$lab <- .format_value_label(rects$state, rects$count, rects$proportion, label)
    rects$tile_w <- rects$xmax - rects$xmin
    rects$tile_h <- rects$ymax - rects$ymin
    rects$angle <- ifelse(rects$tile_h > rects$tile_w, 90, 0)
    p <- p + .geom_fit_label(rects, label_size,
                             fill_colors = pal[as.character(rects$state)])
  }
  p
}


# Squarified slice-and-dice treemap: split along the longer side at each
# step, placing the largest remaining rectangle as a full slice. Produces
# rectangles whose AREAS are exactly proportional to `values`, with a
# reasonable aspect ratio for small N (3-9 states).
.simple_treemap <- function(values, x = 0, y = 0, w = 1, h = 1) {
  n <- length(values)
  if (n == 0L) return(data.frame())
  if (n == 1L) {
    return(data.frame(xmin = x, xmax = x + w,
                      ymin = y, ymax = y + h, idx = 1L))
  }
  total <- sum(values)
  share <- values[1L] / total

  if (w >= h) {
    first_w <- w * share
    rest <- .simple_treemap(values[-1L], x + first_w, y, w - first_w, h)
    rest$idx <- rest$idx + 1L
    rbind(
      data.frame(xmin = x, xmax = x + first_w,
                 ymin = y, ymax = y + h, idx = 1L),
      rest
    )
  } else {
    first_h <- h * share
    rest <- .simple_treemap(values[-1L], x, y + first_h, w, h - first_h)
    rest$idx <- rest$idx + 1L
    rbind(
      data.frame(xmin = x, xmax = x + w,
                 ymin = y, ymax = y + first_h, idx = 1L),
      rest
    )
  }
}


# Format a tile / segment label according to user choice. All formats
# render inline on a single line (never two-line).
# "none"  -> ""
# "prop"  -> "66%"
# "freq"  -> "1,234"
# "both"  -> "1,234 (66%)"
# "state" -> "Average"
# "all"   -> "Average (66%)"
# .legend_layer / .legend_theme moved to R/plot-utils.R alongside the other
# shared layout helpers.


.format_value_label <- function(state, count, proportion, label) {
  state <- as.character(state)
  freq_str <- format(count, big.mark = ",", trim = TRUE)
  switch(label,
    none  = rep("", length(count)),
    prop  = sprintf("%.1f%%", 100 * proportion),
    freq  = freq_str,
    both  = sprintf("%s (%.1f%%)", freq_str, 100 * proportion),
    state = state,
    all   = sprintf("%s (%.1f%%)", state, 100 * proportion))
}


# Build a single-panel treemap plot with its OWN legend, restricted to
# the states actually present in `sub`. Color is drawn from the global
# named palette so the same state always uses the same color across the
# isolated plots that share a `gridExtra::arrangeGrob` layout.
.single_treemap_plot <- function(sub, pal_named, label, label_size,
                                  legend, legend_dir, legend_frame,
                                  title = NULL, panel_frame = FALSE) {
  states_here <- as.character(sub$state)
  rects <- .simple_treemap(sub$count)
  rects$state      <- states_here[rects$idx]
  rects$count      <- sub$count[rects$idx]
  rects$proportion <- sub$proportion[rects$idx]
  rects$x_mid <- (rects$xmin + rects$xmax) / 2
  rects$y_mid <- (rects$ymin + rects$ymax) / 2
  rects$tile_w <- rects$xmax - rects$xmin
  rects$tile_h <- rects$ymax - rects$ymin
  rects$tile_area <- rects$tile_w * rects$tile_h
  rects$state <- factor(rects$state, levels = states_here)

  pal_local <- pal_named[states_here]

  p <- ggplot2::ggplot(rects,
                       ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                    ymin = .data$ymin, ymax = .data$ymax,
                                    fill = .data$state)) +
    ggplot2::geom_rect(color = "white", linewidth = 0.6) +
    ggplot2::scale_fill_manual(values = pal_local, drop = FALSE,
                               breaks = states_here) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    ggplot2::labs(title = title, x = NULL, y = NULL, fill = "State") +
    .legend_layer(legend, length(states_here), legend_dir) +
    ggplot2::theme_minimal(base_size = 12) +
    .legend_theme(legend, legend_dir, legend_frame) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )

  if (isTRUE(panel_frame)) {
    p <- p + ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = NA, color = "grey60", linewidth = 0.4),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
  }

  if (label != "none") {
    rects$lab <- .format_value_label(rects$state, rects$count,
                                     rects$proportion, label)
    rects$angle <- ifelse(rects$tile_h > rects$tile_w, 90, 0)
    p <- p + .geom_fit_label(rects, label_size,
                             fill_colors = pal_named[as.character(rects$state)])
  }
  p
}


# Per-facet renderer: one separate treemap plot per group, each with its
# OWN legend, arranged via gridExtra. Each panel's legend shows only the
# states present in that panel; colors stay consistent across panels via
# a single shared named palette. Returns a gtable.
.plot_per_facet_grid <- function(freq_df, sort_states, colors,
                                  label, label_size,
                                  legend = "right",
                                  legend_dir = "auto",
                                  legend_frame = "border",
                                  combine = TRUE,
                                  ncol = NULL) {
  groups <- unique(as.character(freq_df$group))
  if (length(groups) == 0L) groups <- "all"

  # Per-panel legend placement: right when there are <= 2 panels (each
  # panel still has plenty of horizontal room, e.g. htna AI + Human);
  # bottom otherwise (3+ panels in a combined gtable squeeze each panel
  # too narrow for a right-side legend, which overflows the column).
  if (legend == "per_facet") {
    legend <- if (length(groups) <= 2L) "right" else "bottom"
  }

  # Resolve combine = "auto": when there are 4+ panels, return a list so
  # each panel renders as its own full-size figure under knitr (the chunk
  # `fig.width`/`fig.height` per panel is much more readable than cramming
  # 6+ panels into a 3x2 gtable). 1-3 panels still combine cleanly.
  if (identical(combine, "auto")) {
    combine <- length(groups) <= 3L
  }
  if (isTRUE(combine) && !requireNamespace("gridExtra", quietly = TRUE)) {
    stop(errorCondition(paste0(
      "legend = 'per_facet' with combine = TRUE requires the ",
      "'gridExtra' package. Install it, set combine = FALSE, or pick ",
      "a different legend value."),
      class = "nestimate_missing_package", call = NULL))
  }

  # Per-panel palettes: each group is independent, so colors restart from
  # the first palette entry inside every panel. This is the right choice
  # when state vocabularies are disjoint (htna actors, mcml clusters) --
  # there is no cross-panel matching to preserve, and starting fresh
  # avoids palette recycling once a panel has many states.
  plots <- lapply(groups, function(g) {
    sub <- freq_df[as.character(freq_df$group) == g, , drop = FALSE]
    local_levels <- .order_states(sub$state, sub$count, sort_states)
    sub <- sub[match(local_levels, as.character(sub$state)), , drop = FALSE]
    sub <- sub[!is.na(sub$state) & sub$count > 0, , drop = FALSE]
    if (nrow(sub) == 0L) return(NULL)

    local_pal <- .state_palette(colors, length(local_levels), local_levels)
    names(local_pal) <- local_levels

    .single_treemap_plot(sub, local_pal, label, label_size,
                          legend, legend_dir, legend_frame,
                          title = if (length(groups) > 1L) g else NULL,
                          panel_frame = FALSE)
  })
  plots <- Filter(Negate(is.null), plots)

  if (!isTRUE(combine)) {
    class(plots) <- c("nestimate_facet_list", "list")
    return(plots)
  }

  ncol_arrange <- if (!is.null(ncol)) {
    as.integer(ncol)
  } else if (length(plots) <= 2L) {
    length(plots)
  } else {
    # Per-panel legends are wide, so default to 2 columns for any
    # 3+ panel layout. User can override with ncol = 3 explicitly.
    2L
  }

  gt <- gridExtra::arrangeGrob(grobs = plots, ncol = ncol_arrange)
  class(gt) <- c("nestimate_facet_plot", class(gt))
  gt
}


#' @rdname plot_state_frequencies
#' @export
print.nestimate_facet_plot <- function(x, ...) {
  grid::grid.newpage()
  grid::grid.draw(x)
  invisible(x)
}


#' @rdname plot_state_frequencies
#' @export
print.nestimate_facet_list <- function(x, ...) {
  lapply(x, print)
  invisible(x)
}


# knitr's auto-print captures only the last `grid.newpage` from a single
# evaluated expression. To get one rendered figure per panel we must hook
# into knitr::knit_print and emit each plot through its own knit_print
# call (which knitr tracks individually). Falls back to print() outside
# knit context.
#' @rawNamespace if (getRversion() >= "3.6.0") S3method(knitr::knit_print, nestimate_facet_list)
#' @keywords internal
knit_print.nestimate_facet_list <- function(x, ...) {
  if (!requireNamespace("knitr", quietly = TRUE)) {
    print(x); return(invisible(NULL))
  }
  parts <- lapply(x, function(p) knitr::knit_print(p, ...))
  knitr::asis_output(paste(unlist(parts), collapse = "\n\n"))
}


# Treemap renderer: one panel per group, each panel a squarified treemap
# whose tile areas are proportional to within-group state counts.
# Single-panel when no groups present.
.plot_treemap_panels <- function(freq_df, sort_states, colors,
                                  label, label_size,
                                  legend = "bottom",
                                  legend_dir = "auto",
                                  legend_frame = "none") {
  state_levels <- .order_states(freq_df$state, freq_df$count, sort_states)
  pal <- .state_palette(colors, length(state_levels), state_levels)
  names(pal) <- state_levels

  groups <- unique(as.character(freq_df$group))
  has_groups <- length(groups) > 1L || !identical(groups, "all")

  # Use the SAME global state order in every facet so the leftmost (and
  # second, third, ...) tile always represents the same state across
  # panels. Tile sizes still reflect within-panel proportions, but the
  # topological position of each color stays put -- making cross-facet
  # comparison by color and tile location possible.
  rects_list <- lapply(groups, function(g) {
    sub <- freq_df[as.character(freq_df$group) == g, , drop = FALSE]
    sub <- sub[match(state_levels, as.character(sub$state)), , drop = FALSE]
    sub <- sub[!is.na(sub$state) & sub$count > 0, , drop = FALSE]
    if (nrow(sub) == 0L) return(NULL)
    rects <- .simple_treemap(sub$count)
    rects$state      <- as.character(sub$state[rects$idx])
    rects$count      <- sub$count[rects$idx]
    rects$proportion <- sub$proportion[rects$idx]
    rects$group      <- g
    rects$x_mid      <- (rects$xmin + rects$xmax) / 2
    rects$y_mid      <- (rects$ymin + rects$ymax) / 2
    rects$tile_area  <- (rects$xmax - rects$xmin) * (rects$ymax - rects$ymin)
    rects
  })
  rects <- do.call(rbind, Filter(Negate(is.null), rects_list))
  rects$state <- factor(rects$state, levels = state_levels)

  p <- ggplot2::ggplot(rects,
                       ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                    ymin = .data$ymin, ymax = .data$ymax,
                                    fill = .data$state)) +
    ggplot2::geom_rect(color = "white", linewidth = 0.6) +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    ggplot2::labs(x = NULL, y = NULL, fill = "State") +
    .legend_layer(legend, length(state_levels), legend_dir) +
    ggplot2::theme_minimal(base_size = 12) +
    .legend_theme(legend, legend_dir, legend_frame) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )

  if (has_groups) {
    p <- p + ggplot2::facet_wrap(~ group)
  }

  if (label != "none") {
    rects$lab <- .format_value_label(rects$state, rects$count, rects$proportion, label)
    rects$tile_w <- rects$xmax - rects$xmin
    rects$tile_h <- rects$ymax - rects$ymin
    rects$angle <- ifelse(rects$tile_h > rects$tile_w, 90, 0)
    p <- p + .geom_fit_label(rects, label_size,
                             fill_colors = pal[as.character(rects$state)])
  }
  p
}




# Horizontal bars: state on Y, count/proportion on X, faceted by group.
# Bar length encodes whichever metric you pass via `metric` ("freq" or
# "prop"). The `label` arg controls the inline numeric annotation.
.plot_state_bars <- function(freq_df, sort_states, colors,
                              label, label_size, metric,
                              legend = "bottom",
                              legend_dir = "auto",
                              legend_frame = "none") {
  state_levels <- .order_states(freq_df$state, freq_df$count, sort_states)
  # Reverse so the largest count appears at the top of the y-axis
  freq_df$state <- factor(freq_df$state, levels = rev(state_levels))
  pal <- .state_palette(colors, length(state_levels), state_levels)
  names(pal) <- state_levels

  x_var <- if (metric == "freq") "count" else "proportion"
  x_lab <- if (metric == "freq") "Count" else "Proportion"

  freq_df$lab <- .format_value_label(freq_df$state, freq_df$count, freq_df$proportion, label)

  groups <- unique(freq_df$group)
  has_groups <- length(groups) > 1L || !identical(groups, "all")

  p <- ggplot2::ggplot(freq_df,
                       ggplot2::aes(y = .data$state, x = .data[[x_var]],
                                    fill = .data$state)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE,
                               breaks = state_levels) +
    ggplot2::labs(x = x_lab, y = NULL, fill = "State") +
    .legend_layer(legend, length(state_levels), legend_dir) +
    ggplot2::theme_minimal(base_size = 12) +
    .legend_theme(legend, legend_dir, legend_frame) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

  if (has_groups) {
    p <- p + ggplot2::facet_wrap(~ group, scales = "free_x")
  }

  if (label != "none") {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$lab), hjust = -0.15, size = label_size,
      color = "grey20"
    ) + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.28)))
  }
  p
}


# ------------------------------------------------------------------------------
# Generic + methods (exported)
# ------------------------------------------------------------------------------

#' Plot State Frequency Distributions
#'
#' @description
#' `plot_state_frequencies()` counts how often each state occurs in the
#' sequences behind a network and draws the shares, per group, as treemaps
#' or bars; it returns the chart together with the table of counts. It
#' takes a `netobject`, a `netobject_group` (one group per network), an
#' `mcml` (one group per cluster) or an `htna` network (one group per
#' actor). A single `netobject` is one group, `"all"`, unless it carries
#' node groups, which then form the groups. Missing and empty cells are not
#' counted.
#'
#' With `style = "marimekko"` (default), each group is a treemap whose
#' tiles have areas proportional to the states' shares within that group.
#' The groups of a `netobject` or `netobject_group` are panels of one plot,
#' with the states in the same order and colours in every panel. The actors
#' of an `htna` network and the clusters of an `mcml` usually have states
#' of their own, so each is drawn as a separate treemap with its own legend
#' (`legend = "per_facet"`). An `mcml` drawn with any other legend is a
#' single marimekko instead: one column per cluster, as wide as the
#' cluster's share of all events, with its states stacked by their share
#' within the cluster. With `style = "bars"`, each group is a panel of
#' horizontal bars, one per state.
#'
#' The shares describe how often states occur, not how they follow each
#' other. Two groups can have the same state mix and different
#' transitions; compare the networks for that, with [permutation()].
#' [mosaic_plot()] shows which cells of a `netobject_group`'s
#' group-by-state table depart from independence.
#'
#' @param x The object to count: a `netobject`, `netobject_group`, `mcml`
#'   or `htna` network built from sequences. A network without sequence
#'   data raises an error; in a `netobject_group` or `mcml`, groups without
#'   it are left out. For the `print()` and `plot()`
#'   methods of the result: the `state_freq` object;
#'   for the `print()` methods of a per-panel figure: an object of class
#'   `nestimate_facet_plot` or `nestimate_facet_list`.
#' @param style The chart. `"marimekko"` (default) draws treemaps, or the
#'   cluster marimekko for an `mcml` (see Description); `"bars"` draws
#'   horizontal bars, the most frequent state on top when `sort_states =
#'   "frequency"`, with one panel per group.
#' @param metric For `style = "bars"`: what the bar length shows,
#'   `"prop"` (default), the share within the group, or `"freq"`, the
#'   count. Treemap and marimekko areas always show the share within the
#'   group.
#' @param label The text printed on each tile or bar, on one line:
#'   `"prop"` (default) the share, as in `"66.0%"`; `"freq"` the count, as
#'   in `"1,234"`; `"both"` the two, `"1,234 (66.0%)"`; `"state"` the state
#'   name; `"all"` the name and share, `"Average (66.0%)"`; `"none"` no
#'   text.
#' @param legend Where the colour legend goes: `"bottom"`, `"top"`,
#'   `"right"`, `"left"`, `"none"` or `"per_facet"`. `"auto"` (default)
#'   resolves to `"none"` for bars, whose axis names the states, to
#'   `"per_facet"` for `htna` and `mcml` treemaps, and to `"bottom"`
#'   otherwise. `"per_facet"` draws each group as its own treemap with its
#'   own legend (on the right with up to two panels, at the bottom with
#'   more), and its colours restart in every panel; it falls back to
#'   `"bottom"` when every group has the same states, and to a legend on
#'   the right for bars. It returns a combined figure or a list of plots,
#'   as `combine` sets.
#' @param legend_dir The direction of the legend's keys: `"auto"`
#'   (default) lays them out horizontally for a legend on top or bottom and
#'   vertically on the left or right; `"horizontal"` and `"vertical"`
#'   force one.
#' @param legend_frame `"none"` (default) draws the legend unframed;
#'   `"border"` draws a thin grey rectangle around it.
#' @param sort_states The order of the states in the legend, the tiles and
#'   the bars: `"frequency"` (default), the most frequent first;
#'   `"alpha"`, alphabetical; or `"none"`, the order of the table.
#' @param colors Fill colours for the states, overriding the palette
#'   attached to `x` with [set_state_colors()] for this chart: an unnamed
#'   vector in state order (at least one per state) or a named lookup such
#'   as `c(plan = "#0072B2")` that overrides only the states named. Default
#'   `NULL`: the attached palette, or else the Okabe-Ito palette, recycled.
#' @param label_size The size of the tile and bar labels; with
#'   \pkg{ggfittext} installed, the largest size, as labels shrink to fit
#'   their tiles. Default `3.5`.
#' @param abbreviate Whether to shorten the state names with
#'   [base::abbreviate()], which drops lower-case vowels and then letters
#'   from the end until a name fits, keeping the names distinct. `FALSE`
#'   (default) keeps the full names; `TRUE` shortens them to at least three
#'   characters; a positive whole number sets that minimum length. It
#'   applies to the labels, the legend and the returned table.
#' @param include_macro For an `mcml`: whether to add a `"macro"` group in
#'   front, holding the state counts over all clusters. Default `FALSE`.
#'   Ignored for the other classes.
#' @param combine With per-panel legends: `"auto"` (default) combines up
#'   to three panels into one figure and returns a list of plots, one per
#'   panel, for four or more, which knitr draws as separate figures at the
#'   chunk's size. `TRUE` always combines them; `FALSE` always returns the
#'   list. A combined figure, whether from `TRUE` or from `"auto"` with up
#'   to three panels, needs \pkg{gridExtra} and raises an error without
#'   it; `FALSE` does not.
#' @param node_groups A named character vector that maps states to broader
#'   groups, such as `c(plan = "Regulation", monitor = "Regulation")`. The
#'   chart is then drawn and coloured by these groups, states not named
#'   keeping their own name; the returned table keeps one row per state.
#'   Default `NULL`.
#' @param ncol For a combined figure with per-panel legends: the number of
#'   columns. Default `NULL`: one column per panel with up to two panels,
#'   two columns with more.
#' @param ... Not used: the `netobject`, `netobject_group` and `htna`
#'   methods ignore an argument they do not take with a warning of class
#'   `nestimate_ignored_argument`; the `mcml` method raises an error for
#'   it. For the `print()` and `plot()` methods of the result:
#'   ignored.
#'
#' @return A `state_freq` object. Printing it shows the event totals per
#'   group and the share of each state within its group, then draws the
#'   chart; `plot()` draws the chart alone. [state_distribution()] of the
#'   result returns the table as drawn, one row per group and state;
#'   [state_distribution()] of the network gives
#'   the same counts, with the full state names and, for
#'   an `mcml`, without the macro group unless `include_macro = TRUE`.
#'   The object holds the chart (a
#'   `ggplot`; with per-panel legends either a combined figure, a `gtable`
#'   of class `nestimate_facet_plot`, or a list of `ggplot`s of class
#'   `nestimate_facet_list`), the table (one row per group and state, with columns `group`,
#'   `state`, `count` and `proportion`, the share within the group), and
#'   the `style`, `metric` and source class of the call. The `print()`
#'   methods return their input invisibly, `plot()` returns `NULL`
#'   invisibly. Any other class raises an error, a tna model included;
#'   draw a tna model with `tna::plot_frequencies()`.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   data(group_regulation_long, package = "Nestimate")
#'   nw <- build_network(group_regulation_long,
#'                       method = "relative", format = "long",
#'                       actor = "Actor", action = "Action",
#'                       order = "Time", group = "Course")
#'   res <- plot_state_frequencies(nw)
#'   print(res)            # tidy frequency table in the console
#'   plot(res)             # ggplot chart
#'   head(state_distribution(res))
#' }
#' @order 1
#' @export
plot_state_frequencies <- function(x, ...) {
  UseMethod("plot_state_frequencies")
}


#' @export
#' @rdname plot_state_frequencies
plot_state_frequencies.netobject <- function(x,
                                              style = "marimekko",
                                              metric = "prop",
                                              label = "prop",
                                              legend = "auto",
                                              legend_dir = "auto",
                                              legend_frame = "none",
                                              sort_states = "frequency",
                                              colors = NULL,
                                              label_size = 3.5,
                                              abbreviate = FALSE,
                                              include_macro = FALSE,
                                              combine = "auto",
                                              ncol = NULL,
                                              node_groups = NULL,
                                              ...) {
  .plot_state_frequencies_impl(
    x, source_class = "netobject", hierarchical = FALSE,
    style = style, metric = metric, label = label, legend = legend,
    legend_dir = legend_dir, legend_frame = legend_frame,
    sort_states = sort_states, colors = colors, label_size = label_size,
    abbreviate = abbreviate, include_macro = include_macro,
    combine = combine, ncol = ncol, node_groups = node_groups, ...)
}

#' @export
#' @rdname plot_state_frequencies
plot_state_frequencies.htna <- function(x,
                                         style = "marimekko",
                                         metric = "prop",
                                         label = "prop",
                                         legend = "auto",
                                         legend_dir = "auto",
                                         legend_frame = "none",
                                         sort_states = "frequency",
                                         colors = NULL,
                                         label_size = 3.5,
                                         abbreviate = FALSE,
                                         include_macro = FALSE,
                                         combine = "auto",
                                         ncol = NULL,
                                         node_groups = NULL,
                                         ...) {
  .plot_state_frequencies_impl(
    x, source_class = "htna", hierarchical = FALSE,
    style = style, metric = metric, label = label, legend = legend,
    legend_dir = legend_dir, legend_frame = legend_frame,
    sort_states = sort_states, colors = colors, label_size = label_size,
    abbreviate = abbreviate, include_macro = include_macro,
    combine = combine, ncol = ncol, node_groups = node_groups, ...)
}

#' @export
#' @rdname plot_state_frequencies
plot_state_frequencies.mcml <- function(x,
                                         style = "marimekko",
                                         metric = "prop",
                                         label = "prop",
                                         legend = "auto",
                                         legend_dir = "auto",
                                         legend_frame = "none",
                                         sort_states = "frequency",
                                         colors = NULL,
                                         label_size = 3.5,
                                         abbreviate = FALSE,
                                         include_macro = FALSE,
                                         combine = "auto",
                                         ncol = NULL,
                                         node_groups = NULL,
                                         ...) {
  .plot_state_freq_check_unused_dots("plot_state_frequencies.mcml", ...)
  .plot_state_frequencies_impl(
    x, source_class = "mcml", hierarchical = TRUE,
    style = style, metric = metric, label = label, legend = legend,
    legend_dir = legend_dir, legend_frame = legend_frame,
    sort_states = sort_states, colors = colors, label_size = label_size,
    abbreviate = abbreviate, include_macro = include_macro,
    combine = combine, ncol = ncol, node_groups = node_groups, ...)
}

#' @export
#' @rdname plot_state_frequencies
plot_state_frequencies.netobject_group <- function(x,
                                                    style = "marimekko",
                                                    metric = "prop",
                                                    label = "prop",
                                                    legend = "auto",
                                                    legend_dir = "auto",
                                                    legend_frame = "none",
                                                    sort_states = "frequency",
                                                    colors = NULL,
                                                    label_size = 3.5,
                                                    abbreviate = FALSE,
                                                    include_macro = FALSE,
                                                    combine = "auto",
                                                    ncol = NULL,
                                                    node_groups = NULL,
                                                          ...) {
  .plot_state_frequencies_impl(
    x, source_class = "netobject_group", hierarchical = FALSE,
    style = style, metric = metric, label = label, legend = legend,
    legend_dir = legend_dir, legend_frame = legend_frame,
    sort_states = sort_states, colors = colors, label_size = label_size,
    abbreviate = abbreviate, include_macro = include_macro,
    combine = combine, ncol = ncol, node_groups = node_groups, ...)
}

#' @export
#' @rdname plot_state_frequencies
plot_state_frequencies.default <- function(x, ...) {
  cls <- paste(class(x), collapse = "/")
  stop(errorCondition(
    paste0(
      "plot_state_frequencies(): no method for class '",
      cls,
      "'.\n",
      "Supported: netobject, netobject_group, mcml, htna.\n",
      "For tna objects, use tna::plot_frequencies(x)."),
    class = "nestimate_bad_input", call = NULL))
}


# Single worker for all four dispatch methods. Validates args once, pulls
# the tidy frame via state_distribution(), runs the existing renderer
# pipeline, wraps both into a state_freq object.
.plot_state_frequencies_impl <- function(x, source_class, hierarchical,
                                          style       = c("marimekko",
                                                          "bars"),
                                          metric      = c("prop", "freq"),
                                          label       = c("prop", "freq",
                                                          "both", "state",
                                                          "all", "none"),
                                          legend      = c("auto", "bottom",
                                                          "right", "top",
                                                          "left", "none",
                                                          "per_facet"),
                                          legend_dir  = c("auto",
                                                          "horizontal",
                                                          "vertical"),
                                          legend_frame = c("none", "border"),
                                          sort_states = c("frequency",
                                                          "alpha", "none"),
                                          colors      = NULL,
                                          label_size  = 3.5,
                                          abbreviate  = FALSE,
                                          include_macro = FALSE,
                                          combine     = "auto",
                                          ncol        = NULL,
                                          node_groups = NULL,
                                          ...) {
  style        <- match.arg(style)
  metric       <- match.arg(metric)
  label        <- match.arg(label)
  legend       <- match.arg(legend)
  legend_dir   <- match.arg(legend_dir)
  legend_frame <- match.arg(legend_frame)
  sort_states  <- match.arg(sort_states)
  # A palette attached with set_state_colors() is the object's own default;
  # an explicit `colors` still wins for this one figure.
  colors       <- colors %||% .stored_state_colors(x)
  # Nothing below reads `...` (the mcml method has already refused it).
  .warn_ignored_dots("plot_state_frequencies()", ...)

  if (identical(legend, "auto")) {
    legend <- if (style == "bars") {
      "none"
    } else if (source_class %in% c("htna", "mcml")) {
      "per_facet"
    } else {
      "bottom"
    }
  }

  freq_df <- if (identical(source_class, "mcml")) {
    state_distribution(x, include_macro = include_macro)
  } else {
    state_distribution(x)
  }

  freq_df <- .abbreviate_states(freq_df, abbreviate)

  # Demote per_facet to a shared bottom legend when every group has the
  # same state vocabulary -- repeating the same legend in every panel is
  # pure redundancy.
  if (identical(legend, "per_facet") &&
      .vocab_is_shared(freq_df)) {
    legend <- "bottom"
  }

  p <- .render_freq(freq_df, style, hierarchical = hierarchical,
                    metric, label, sort_states, colors, label_size,
                    legend, legend_dir, legend_frame,
                    combine = combine, ncol = ncol,
                    node_groups = node_groups)

  .new_state_freq(plot = p, table = freq_df,
                  style = style, metric = metric,
                  source_class = source_class)
}


# Returns TRUE when every group in freq_df contains the same set of
# states (any order). Used to detect per_facet calls that would render
# the same legend k times.
# .vocab_is_shared moved to R/plot-utils.R.


# Abbreviate state names with base R `abbreviate()` so duplicates after
# truncation get extra letters automatically. `abbreviate = FALSE` (default)
# is a no-op; `TRUE` uses minlength = 3; a numeric value sets minlength.
# Aggregates the freq_df after truncation in case different long names
# collapse to the same short name despite the uniqueness guarantee
# (rare; happens with single-letter prefixes or after manual relabel).
.abbreviate_states <- function(freq_df, abbreviate) {
  if (isFALSE(abbreviate) || is.null(abbreviate)) return(freq_df)
  minlen <- if (isTRUE(abbreviate)) 3L else as.integer(abbreviate[[1L]])
  .assert(minlen >= 1L)
  short <- abbreviate(as.character(freq_df$state), minlength = minlen)
  freq_df$state <- unname(short)
  agg <- aggregate(count ~ group + state, data = freq_df, FUN = sum)
  totals <- ave(agg$count, agg$group, FUN = sum)
  agg$proportion <- agg$count / totals
  agg[, c("group", "state", "count", "proportion")]
}


# Collapse the per-state fill onto a smaller categorical legend keyed by
# semantic group, per the `node_groups` @param. `node_groups` is a named
# character vector mapping state (node) labels to group labels. States not
# present in the mapping keep their own name (no states are dropped or
# invented). Counts are re-aggregated per (group, mapped-state) and the
# within-group proportion recomputed -- the same pattern as
# .abbreviate_states. `NULL` (default) is a no-op.
.collapse_states_to_groups <- function(freq_df, node_groups) {
  if (is.null(node_groups)) return(freq_df)
  .assert(is.character(node_groups), length(node_groups) >= 1L,
            !is.null(names(node_groups)))
  lookup <- node_groups
  s <- as.character(freq_df$state)
  mapped <- unname(lookup[s])
  freq_df$state <- ifelse(is.na(mapped), s, mapped)
  agg <- aggregate(count ~ group + state, data = freq_df, FUN = sum)
  totals <- ave(agg$count, agg$group, FUN = sum)
  agg$proportion <- agg$count / totals
  agg[, c("group", "state", "count", "proportion")]
}


# ------------------------------------------------------------------------------
# state_distribution(): public tidy frequency extractor
# ------------------------------------------------------------------------------

#' Per-Class State Distribution as a Tidy Data Frame
#'
#' @description
#' `state_distribution()` counts how often each state occurs in the
#' sequences behind a network and returns one row per group and state,
#' with the count and the share within the group. The groups follow the
#' class: a single `netobject` is one group, `"all"`, unless it carries node
#' groups, which then form the groups (states in none are `"(ungrouped)"`);
#' an `htna` network has one group per actor; a `netobject_group` one per
#' network, named as in the group (or `"Group 1"`, `"Group 2"`, ...); an
#' `mcml` one per cluster, with the node labels it was given in place of
#' the raw codes. Missing and empty cells are not counted.
#'
#' The table is the one behind [plot_state_frequencies()], and the per-group
#' counterpart of [state_frequencies()], which counts the states of raw
#' sequence data.
#'
#' @param x A `netobject`, `netobject_group`, `mcml` or `htna` network
#'   built from sequences, or the `state_freq` result of
#'   [plot_state_frequencies()], whose table is returned as drawn. A network without sequence data raises an error;
#'   in a `netobject_group` or `mcml`, groups without it are left out, and
#'   an error is raised only when none has any. Other classes raise an
#'   error.
#' @param include_macro For an `mcml`: whether to add a `"macro"` group in
#'   front, holding the state counts over all clusters. Default `FALSE`.
#' @param ... Not used. The `mcml` method raises an error for any argument
#'   it does not take; the other methods ignore it.
#'
#' @return A data frame with one row per group and state and the columns
#'   `group` (character), `state` (character), `count` (integer) and
#'   `proportion` (the state's share of its group's events). Within a group
#'   the rows run from the most to the least frequent state; the `"macro"`
#'   group is in alphabetical order.
#' @export
#' @examples
#' data(group_regulation_long, package = "Nestimate")
#' net <- build_network(group_regulation_long, method = "frequency",
#'                      format = "long", actor = "Actor", action = "Action",
#'                      order = "Time", group = "Course")
#' state_distribution(net)
state_distribution <- function(x, ...) UseMethod("state_distribution")

#' @export
#' @rdname state_distribution
state_distribution.netobject <- function(x, ...) .freq_df_netobject(x)

#' @export
#' @rdname state_distribution
state_distribution.htna <- function(x, ...) .freq_df_htna(x)

#' @export
#' @rdname state_distribution
state_distribution.mcml <- function(x, include_macro = FALSE, ...) {
  .plot_state_freq_check_unused_dots("state_distribution.mcml", ...)
  .freq_df_mcml(x, include_macro = include_macro)
}

.plot_state_freq_check_unused_dots <- function(method, ...) {
  dots <- list(...)
  if (!length(dots)) {
    return(invisible(TRUE))
  }
  dot_names <- names(dots)
  dot_names[!nzchar(dot_names)] <- paste0("..", which(!nzchar(dot_names)))
  stop(errorCondition(
    paste0(
      method, "() got unsupported argument",
      if (length(dots) == 1L) ": " else "s: ",
      paste(dot_names, collapse = ", ")
    ),
    class = "nestimate_bad_input", call = NULL
  ))
}

#' @export
#' @rdname state_distribution
state_distribution.netobject_group <- function(x, ...) {
  .freq_df_netobject_group(x)
}

#' @export
#' @rdname state_distribution
state_distribution.state_freq <- function(x, ...) x$table

#' @export
#' @rdname state_distribution
state_distribution.default <- function(x, ...) {
  cls <- paste(class(x), collapse = "/")
  stop(errorCondition(
    paste0(
      "state_distribution(): no method for class '",
      cls,
      "'.\n",
      "Supported: netobject, netobject_group, mcml, htna."),
    class = "nestimate_bad_input", call = NULL))
}


# ------------------------------------------------------------------------------
# state_freq S3 class: holds plot + tidy table together
# ------------------------------------------------------------------------------

.new_state_freq <- function(plot, table, style, metric, source_class) {
  out <- list(
    plot         = plot,
    table        = table,
    style        = style,
    metric       = metric,
    source_class = source_class
  )
  class(out) <- "state_freq"
  out
}

#' @section The state_freq object:
#' In a knitr document, printing the object puts the table in the text and
#' the chart in a figure of the chunk's size; a list of per-panel plots
#' gives one figure per panel.
#'
#' @param digits For `print()` of the result: decimal places of the
#'   percentages. Default `1`.
#' @param max_states For `print()` of the result: the most states shown
#'   per group (default `20`); the rest are folded into one `"(+k more)"`
#'   row. [state_distribution()] gives the full table.
#' @name state_freq
#' @rdname plot_state_frequencies
NULL

# Draw whatever kind of plot object a state_freq carries: a single ggplot,
# a gtable / grob (per_facet, hierarchical marimekko), or a
# nestimate_facet_list (one full-size ggplot per panel, e.g. 4+ panels or
# include_macro). Centralised so print / plot / knit_print stay in sync.
.draw_state_freq_plot <- function(p) {
  if (is.null(p)) return(invisible(NULL))
  if (inherits(p, "ggplot") || inherits(p, "nestimate_facet_list")) {
    print(p)
  } else if (inherits(p, "gtable") || grid::is.grob(p)) {
    grid::grid.newpage()
    grid::grid.draw(p)
  } else {
    print(p)
  }
  invisible(NULL)
}

#' @rdname plot_state_frequencies
#' @export
print.state_freq <- function(x, digits = 1, max_states = 20L, ...) {
  tbl <- x$table
  groups <- unique(as.character(tbl$group))
  total_events <- sum(tbl$count)
  n_states <- length(unique(tbl$state))
  pct <- function(p) sprintf(paste0("%.", digits, "f%%"), 100 * p)

  cat(sprintf(
    "State frequencies (style = %s, source = %s)\n",
    x$style, x$source_class
  ))
  cat(sprintf(
    "  Total events: %s  |  Groups: %d  |  States: %d\n\n",
    format(total_events, big.mark = ","), length(groups), n_states
  ))

  group_totals <- tapply(tbl$count, tbl$group, sum)
  group_totals <- group_totals[groups]
  totals_lines <- .cluster_table_lines(list(
    group  = groups,
    events = format(as.integer(group_totals), big.mark = ","),
    share  = pct(as.numeric(group_totals) / total_events)
  ))
  cat("Per-group totals\n")
  cat(paste(paste0("  ", totals_lines), collapse = "\n"), "\n\n")

  cat("Per-state proportions (within group)\n")
  parts <- lapply(groups, function(g) {
    sub <- tbl[as.character(tbl$group) == g, , drop = FALSE]
    sub <- sub[order(-sub$count), , drop = FALSE]
    if (nrow(sub) > max_states) {
      kept <- sub[seq_len(max_states), , drop = FALSE]
      kept$state <- as.character(kept$state)
      kept <- rbind(kept, data.frame(
        group = g, state = sprintf("(+%d more)", nrow(sub) - max_states),
        count = sum(sub$count[-seq_len(max_states)]),
        proportion = sum(sub$proportion[-seq_len(max_states)]),
        stringsAsFactors = FALSE
      ))
      sub <- kept
    }
    sub
  })
  combined <- do.call(rbind, parts)
  body_lines <- .cluster_table_lines(list(
    group = as.character(combined$group),
    state = as.character(combined$state),
    count = format(as.integer(combined$count), big.mark = ","),
    share = pct(combined$proportion)
  ))
  cat(paste(paste0("  ", body_lines), collapse = "\n"), "\n")

  # Render the chart to the active graphics device so the table + plot
  # appear together (console: opens a window; knitr: embeds inline).
  .draw_state_freq_plot(x$plot)
  invisible(x)
}

#' @rdname plot_state_frequencies
#' @export
plot.state_freq <- function(x, ...) {
  .draw_state_freq_plot(x$plot)
  invisible(NULL)
}

# Kept so that code written against earlier versions still runs; the
# documented verb is state_distribution().
#' @export
#' @noRd
as.data.frame.state_freq <- function(x, ...) x$table


#' @rawNamespace if (getRversion() >= "3.6.0") S3method(knitr::knit_print, state_freq)
#' @keywords internal
knit_print.state_freq <- function(x, options = list(), ...) {
  if (!requireNamespace("knitr", quietly = TRUE)) {
    print(x); return(invisible(NULL))
  }
  # Capture ONLY the textual summary (strip the plot so print() doesn't draw).
  table_text <- paste(utils::capture.output({
    old <- list(plot = x$plot); x$plot <- NULL
    print(x)
    x$plot <- old$plot
  }), collapse = "\n")
  # Draw the chart to the active chunk device so knitr captures it as a
  # properly sized figure. A gtable (per_facet / hierarchical) must NOT be
  # passed through paste()/unlist() -- that explodes the grob tree into tens
  # of thousands of text fragments and never draws the plot.
  .draw_state_freq_plot(x$plot)
  knitr::asis_output(paste0("\n```\n", table_text, "\n```\n"))
}


# Shared dispatcher.
# `hierarchical = TRUE` (mcml) -> 2D cumulative marimekko.
# `hierarchical = FALSE` (everything else) -> per-panel squarified treemap
# (single panel when only one group exists).
.render_freq <- function(freq_df, style, hierarchical,
                          metric, label, sort_states, colors, label_size,
                          legend, legend_dir, legend_frame,
                          combine = TRUE, ncol = NULL,
                          node_groups = NULL) {
  if (nrow(freq_df) == 0L) {
    stop(errorCondition("No state observations available to plot.",
                        class = "nestimate_bad_input", call = NULL))
  }

  # node_groups: collapse the per-state fill onto the semantic groups
  # before any renderer sees the frame, so bars / treemap / hierarchical /
  # per-facet all colour by group (smaller categorical legend) per the
  # @param. NULL is a no-op.
  collapsed <- !is.null(node_groups)
  freq_df <- .collapse_states_to_groups(freq_df, node_groups)
  # When the fill now encodes group rather than state, relabel the fill
  # legend accordingly for ggplot-returning renderers.
  relabel_fill <- function(p) {
    if (collapsed && inherits(p, "ggplot")) {
      p <- p + ggplot2::labs(fill = "Group")
    }
    p
  }

  # Per-facet legend mode applies to the marimekko/treemap path only.
  # For style = "bars", per_facet doesn't make sense (bars already
  # facet-wrap with shared legend), so fall back to a right-side legend.
  if (identical(legend, "per_facet")) {
    if (style == "bars") {
      return(relabel_fill(.plot_state_bars(freq_df, sort_states, colors,
                              label, label_size, metric,
                              "right", legend_dir, legend_frame)))
    }
    return(.plot_per_facet_grid(freq_df, sort_states, colors,
                                 label, label_size,
                                 legend = "per_facet",
                                 legend_dir = legend_dir,
                                 legend_frame = legend_frame,
                                 combine = combine,
                                 ncol = ncol))
  }

  if (style == "bars") {
    return(relabel_fill(.plot_state_bars(freq_df, sort_states, colors,
                            label, label_size, metric,
                            legend, legend_dir, legend_frame)))
  }

  if (isTRUE(hierarchical)) {
    relabel_fill(.plot_marimekko_hierarchical(freq_df, sort_states, colors,
                                  label, label_size,
                                  legend, legend_dir, legend_frame))
  } else {
    relabel_fill(.plot_treemap_panels(freq_df, sort_states, colors,
                          label, label_size,
                          legend, legend_dir, legend_frame))
  }
}
