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
# ------------------------------------------------------------------------------

.freq_df_netobject <- function(x) {
  if (is.null(x$data)) {
    stop("netobject has no `$data` to compute frequencies from.",
         call. = FALSE)
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
    stop("mcml object has no `$clusters` field.", call. = FALSE)
  }
  cluster_names <- names(x$clusters)
  parts <- lapply(cluster_names, function(cn) {
    cl <- x$clusters[[cn]]
    .freq_one_group(cn, cl$data)
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L) {
    stop("mcml clusters have no `$data` to compute frequencies from.",
         call. = FALSE)
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
  group_names <- names(x)
  if (is.null(group_names)) group_names <- paste0("Group ", seq_along(x))
  parts <- lapply(seq_along(x), function(i) {
    .freq_one_group(group_names[i], x[[i]]$data)
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L) {
    stop("netobject_group has no `$data` in any constituent.",
         call. = FALSE)
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
# Mosaic primitive (exported)
# ------------------------------------------------------------------------------

#' Draw a Marimekko / Mosaic Plot from a Tidy Data Frame
#'
#' Low-level rectangle-coordinate builder for marimekko (mosaic) plots.
#' Column widths are proportional to the per-column total of \code{weight};
#' within each column, segments stack to height 1 with sub-heights
#' proportional to each row's share of that column's total.
#'
#' Used internally by \code{\link{plot_state_frequencies}}; exposed so that
#' other plot methods (e.g. permutation-residual visualisations) can reuse
#' the same geometry by supplying a different fill column.
#'
#' @param data A data.frame in long form. Must contain the columns named
#'   in \code{x}, \code{y}, and \code{weight}.
#' @param x Column name for the X (column) variable.
#' @param y Column name for the Y (segment) variable.
#' @param weight Column name for the cell weight (e.g. count).
#' @param fill Either \code{"y"} (color by Y category, e.g. state -- default)
#'   or the name of another column to map to fill (e.g. a residual column
#'   for diverging color).
#' @param colors Optional character vector of fill colors. When
#'   \code{fill = "y"}, length must be at least the number of distinct y
#'   levels. Defaults to recycled Okabe-Ito.
#' @param show_labels If \code{TRUE}, draw within-segment percentage labels.
#' @param label_size Numeric size for segment labels.
#' @param x_label,y_label Optional axis labels.
#'
#' @return A \code{ggplot} object.
#' @import ggplot2
#' @export
#' @examples
#' df <- data.frame(
#'   group = rep(c("A", "B", "C"), each = 3),
#'   state = rep(c("s1", "s2", "s3"), 3),
#'   count = c(10, 5, 2,  7, 8, 3,  4, 6, 12)
#' )
#' plot_mosaic(df, x = "group", y = "state", weight = "count")
plot_mosaic <- function(data,
                        x,
                        y,
                        weight,
                        fill        = "y",
                        colors      = NULL,
                        show_labels = TRUE,
                        label_size  = 3,
                        x_label     = NULL,
                        y_label     = NULL) {
  stopifnot(is.data.frame(data),
            is.character(x), length(x) == 1L,
            is.character(y), length(y) == 1L,
            is.character(weight), length(weight) == 1L)
  needed <- c(x, y, weight)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0L) {
    stop("plot_mosaic: missing columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
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
    stop("plot_mosaic: total weight is zero -- nothing to draw.", call. = FALSE)
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
    .state_palette(colors, length(y_levels))
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
    rects$lab   <- sprintf("%.0f%%", 100 * rects$pct)
    # Drop labels for very thin segments to avoid overlap
    rects$lab[rects$pct < 0.04] <- ""
    p <- p + ggplot2::geom_text(
      data = rects,
      ggplot2::aes(x = .data$x_mid, y = .data$y_mid, label = .data$lab),
      inherit.aes = FALSE, size = label_size, color = "grey15"
    )
  }
  p
}


# Tiny null-coalesce helper (avoid rlang dep)
`%||%` <- function(a, b) if (is.null(a)) b else a


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

  pal <- .state_palette(colors, length(state_levels))
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
    rects$lab[rects$tile_height < 0.06] <- ""
    rects$tile_w <- rects$xmax - rects$xmin
    rects$tile_h <- rects$ymax - rects$ymin
    # Rotate label 90 degrees when the tile is taller than wide so long
    # text (state names, "Average (66%)") follows the tile orientation
    # instead of overflowing its sides.
    rects$angle <- ifelse(rects$tile_h > rects$tile_w, 90, 0)
    p <- p + ggplot2::geom_text(
      data = rects,
      ggplot2::aes(x = .data$x_mid, y = .data$y_mid, label = .data$lab,
                   angle = .data$angle),
      inherit.aes = FALSE, size = label_size, color = "grey15"
    )
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
# Build the guide_legend layer for a given (position, direction) pair.
# legend_dir = "auto" derives from position (bottom/top -> horizontal,
# left/right -> vertical). "horizontal" / "vertical" force the layout.
.legend_layer <- function(legend, n_states, legend_dir = "auto") {
  if (identical(legend, "none")) {
    return(ggplot2::guides(fill = "none"))
  }
  effective_dir <- if (legend_dir == "auto") {
    if (legend %in% c("bottom", "top")) "horizontal" else "vertical"
  } else legend_dir

  if (effective_dir == "horizontal") {
    ncol_legend <- min(n_states, 5L)
    ggplot2::guides(fill = ggplot2::guide_legend(
      direction = "horizontal", ncol = ncol_legend, byrow = TRUE
    ))
  } else {
    ggplot2::guides(fill = ggplot2::guide_legend(
      direction = "vertical", ncol = 1L
    ))
  }
}


# Theme additions controlling legend position, internal layout, and an
# optional border ("frame") around the legend box.
.legend_theme <- function(legend, legend_dir = "auto", legend_frame = "none") {
  if (identical(legend, "none")) {
    return(ggplot2::theme(legend.position = "none"))
  }
  effective_dir <- if (legend_dir == "auto") {
    if (legend %in% c("bottom", "top")) "horizontal" else "vertical"
  } else legend_dir

  bg_rect <- if (identical(legend_frame, "border")) {
    ggplot2::element_rect(color = "grey40", fill = "white", linewidth = 0.4)
  } else {
    ggplot2::element_blank()
  }

  ggplot2::theme(
    legend.position   = legend,
    legend.box        = effective_dir,
    legend.title      = ggplot2::element_text(face = "bold"),
    legend.background = bg_rect,
    legend.box.background = bg_rect,
    legend.margin     = if (identical(legend_frame, "border"))
                          ggplot2::margin(4, 6, 4, 6) else ggplot2::margin()
  )
}


.format_value_label <- function(state, count, proportion, label) {
  state <- as.character(state)
  freq_str <- format(count, big.mark = ",", trim = TRUE)
  switch(label,
    none  = rep("", length(count)),
    prop  = sprintf("%.0f%%", 100 * proportion),
    freq  = freq_str,
    both  = sprintf("%s (%.0f%%)", freq_str, 100 * proportion),
    state = state,
    all   = sprintf("%s (%.0f%%)", state, 100 * proportion))
}


# Build a single-panel treemap plot with its OWN legend, restricted to
# the states actually present in `sub`. Color is drawn from the global
# named palette so the same state always uses the same color across the
# isolated plots that share a `gridExtra::arrangeGrob` layout.
.single_treemap_plot <- function(sub, pal_named, label, label_size,
                                  legend, legend_dir, legend_frame,
                                  title = NULL) {
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
    ggplot2::theme_minimal(base_size = 11) +
    .legend_theme(legend, legend_dir, legend_frame) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )

  if (label != "none") {
    rects$lab <- .format_value_label(rects$state, rects$count,
                                     rects$proportion, label)
    rects$lab[rects$tile_area < 0.05] <- ""
    rects$angle <- ifelse(rects$tile_h > rects$tile_w, 90, 0)
    p <- p + ggplot2::geom_text(
      data = rects,
      ggplot2::aes(x = .data$x_mid, y = .data$y_mid, label = .data$lab,
                   angle = .data$angle),
      inherit.aes = FALSE, size = label_size, color = "grey15"
    )
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
  if (isTRUE(combine) && !requireNamespace("gridExtra", quietly = TRUE)) {
    stop("legend = 'per_facet' with combine = TRUE requires the ",
         "'gridExtra' package. Install it, set combine = FALSE, or pick ",
         "a different legend value.", call. = FALSE)
  }
  if (legend == "per_facet") legend <- "right"

  groups <- unique(as.character(freq_df$group))
  if (length(groups) == 0L) groups <- "all"

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

    local_pal <- .state_palette(colors, length(local_levels))
    names(local_pal) <- local_levels

    .single_treemap_plot(sub, local_pal, label, label_size,
                          legend, legend_dir, legend_frame,
                          title = if (length(groups) > 1L) g else NULL)
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
  } else if (length(plots) <= 4L) {
    2L
  } else {
    3L
  }

  gt <- gridExtra::arrangeGrob(grobs = plots, ncol = ncol_arrange)
  class(gt) <- c("nestimate_facet_plot", class(gt))
  gt
}


#' @export
#' @keywords internal
print.nestimate_facet_plot <- function(x, ...) {
  grid::grid.newpage()
  grid::grid.draw(x)
  invisible(x)
}


#' @export
#' @keywords internal
print.nestimate_facet_list <- function(x, ...) {
  for (p in x) print(p)
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
  pal <- .state_palette(colors, length(state_levels))
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
    rects$lab[rects$tile_area < 0.05] <- ""
    rects$tile_w <- rects$xmax - rects$xmin
    rects$tile_h <- rects$ymax - rects$ymin
    rects$angle <- ifelse(rects$tile_h > rects$tile_w, 90, 0)
    p <- p + ggplot2::geom_text(
      data = rects,
      ggplot2::aes(x = .data$x_mid, y = .data$y_mid, label = .data$lab,
                   angle = .data$angle),
      inherit.aes = FALSE, size = label_size, color = "grey15"
    )
  }
  p
}


# Standardized-residual mosaic: state on Y (rows), group on X (columns),
# row heights proportional to state totals across all groups, column
# widths proportional to group totals. Cells are colored by the
# Pearson/standardized residual of the state x group contingency table:
#   r_{ij} = (O - E) / sqrt(E * (1 - row_p) * (1 - col_p))
# matching the formula used in sequence_compare. Diverging palette
# centered at zero so under-represented (red) and over-represented
# (blue) cells read symmetrically.
.plot_state_residuals <- function(freq_df, sort_states, label, label_size,
                                   legend = "right",
                                   legend_dir = "auto",
                                   legend_frame = "none",
                                   node_groups = NULL,
                                   limit = 4) {
  groups <- unique(as.character(freq_df$group))
  if (length(groups) < 2L) {
    stop("style = 'residual' requires at least 2 groups for the ",
         "contingency table.", call. = FALSE)
  }
  state_levels <- .order_states(freq_df$state, freq_df$count, sort_states)

  counts <- tapply(freq_df$count,
                   list(factor(freq_df$state, levels = state_levels),
                        factor(freq_df$group, levels = groups)),
                   FUN = sum, default = 0L)
  counts[is.na(counts)] <- 0
  N <- sum(counts)
  row_totals <- rowSums(counts)
  col_totals <- colSums(counts)
  row_props <- row_totals / N
  col_props <- col_totals / N
  expected <- outer(row_totals, col_totals) / N
  denom <- sqrt(expected * outer(1 - row_props, 1 - col_props))
  resid <- (counts - expected) / denom
  resid[!is.finite(resid)] <- 0

  # Aligned-column mosaic: row heights from state marginals, column
  # widths from group marginals (shared across all rows). Columns line
  # up vertically across the chart so group dividers are continuous.
  # Cells colored by residual, not by joint area.
  row_heights <- row_totals / N
  col_widths  <- col_totals / N
  y_top  <- c(1, 1 - cumsum(row_heights))
  x_left <- c(0, cumsum(col_widths))

  cells <- expand.grid(state = state_levels, group = groups,
                       stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  cells$count    <- as.numeric(counts[cbind(cells$state, cells$group)])
  cells$expected <- as.numeric(expected[cbind(cells$state, cells$group)])
  cells$residual <- as.numeric(resid[cbind(cells$state, cells$group)])
  ri <- match(cells$state, state_levels)
  ci <- match(cells$group, groups)
  cells$xmin <- x_left[ci]
  cells$xmax <- x_left[ci + 1L]
  # Inset each cell vertically by a small fraction of its band height
  # so adjacent rows show a visible white separator.
  band_h <- y_top[ri] - y_top[ri + 1L]
  inset <- pmin(band_h * 0.08, 0.012)
  cells$ymax <- y_top[ri] - inset
  cells$ymin <- y_top[ri + 1L] + inset
  cells$x_mid <- (cells$xmin + cells$xmax) / 2
  cells$y_mid <- (cells$ymin + cells$ymax) / 2

  abs_max <- if (is.null(limit)) max(abs(cells$residual)) else limit
  if (!is.finite(abs_max) || abs_max == 0) abs_max <- 1
  # Clamp residuals at the limit so colors saturate instead of going
  # blank when one cell dominates (the rare extreme cell).
  cells$residual_clipped <- pmax(pmin(cells$residual, abs_max), -abs_max)

  state_mids <- (y_top[-length(y_top)] + y_top[-1L]) / 2
  group_mids <- (x_left[-length(x_left)] + x_left[-1L]) / 2

  p <- ggplot2::ggplot(cells,
                       ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                    ymin = .data$ymin, ymax = .data$ymax,
                                    fill = .data$residual_clipped)) +
    ggplot2::geom_rect(color = "grey25", linewidth = 0.4) +
    ggplot2::scale_fill_gradient2(
      low = "#B2425A", mid = "#FFFFFF", high = "#5B7DB1",
      midpoint = 0, limits = c(-abs_max, abs_max),
      name = "Standardized\nresidual"
    ) +
    ggplot2::scale_x_continuous(breaks = group_mids, labels = groups,
                                position = "top", expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = state_mids, labels = state_levels,
                                expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 13) +
    .legend_theme(legend, legend_dir, legend_frame) +
    ggplot2::theme(
      panel.grid    = ggplot2::element_blank(),
      axis.text.x   = ggplot2::element_text(face = "bold", size = 13),
      axis.text.y   = ggplot2::element_text(size = 12),
      axis.ticks    = ggplot2::element_blank(),
      plot.background  = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )

  if (!identical(label, "none")) {
    cells$lab <- switch(label,
      prop = sprintf("%.0f%%", 100 * cells$count / N),
      freq = format(cells$count, big.mark = ",", trim = TRUE),
      both = sprintf("%s (%.1f)", format(cells$count, big.mark = ",", trim = TRUE),
                     cells$residual),
      state = as.character(cells$state),
      all  = sprintf("%.1f", cells$residual),
      "")
    cells$txt_color <- ifelse(abs(cells$residual_clipped) > abs_max * 0.85,
                               "white", "grey15")
    p <- p + ggplot2::geom_text(
      data = cells,
      ggplot2::aes(x = .data$x_mid, y = .data$y_mid, label = .data$lab,
                   color = .data$txt_color),
      inherit.aes = FALSE, size = label_size
    ) +
      ggplot2::scale_color_identity()
  }

  # Optional bracket annotations on the LEFT grouping consecutive rows
  # by a node-level grouping (e.g. `htna$node_groups`). Drawn in the
  # negative-x region using line segments + horizontal text labels.
  if (!is.null(node_groups)) {
    ng_map <- if (is.data.frame(node_groups)) {
      setNames(as.character(node_groups[[2]]), as.character(node_groups[[1]]))
    } else if (is.list(node_groups)) {
      unlist(lapply(names(node_groups), function(g)
        setNames(rep(g, length(node_groups[[g]])), node_groups[[g]])))
    } else {
      node_groups  # already named character
    }
    grp_per_state <- ng_map[state_levels]
    if (any(!is.na(grp_per_state))) {
      grp_levels <- unique(grp_per_state[!is.na(grp_per_state)])
      brackets <- do.call(rbind, lapply(grp_levels, function(gl) {
        idx <- which(grp_per_state == gl)
        data.frame(
          group = gl,
          y_top = y_top[min(idx)],
          y_bot = y_top[max(idx) + 1L],
          stringsAsFactors = FALSE
        )
      }))
      brackets$y_mid <- (brackets$y_top + brackets$y_bot) / 2
      bx <- -0.04
      tick <- 0.01
      p <- p +
        ggplot2::geom_segment(
          data = brackets, inherit.aes = FALSE,
          ggplot2::aes(x = bx, xend = bx,
                       y = .data$y_top, yend = .data$y_bot),
          color = "grey20", linewidth = 0.4
        ) +
        ggplot2::geom_segment(
          data = brackets, inherit.aes = FALSE,
          ggplot2::aes(x = bx, xend = bx + tick,
                       y = .data$y_top, yend = .data$y_top),
          color = "grey20", linewidth = 0.4
        ) +
        ggplot2::geom_segment(
          data = brackets, inherit.aes = FALSE,
          ggplot2::aes(x = bx, xend = bx + tick,
                       y = .data$y_bot, yend = .data$y_bot),
          color = "grey20", linewidth = 0.4
        ) +
        ggplot2::geom_text(
          data = brackets, inherit.aes = FALSE,
          ggplot2::aes(x = bx - 0.01, y = .data$y_mid, label = .data$group),
          hjust = 1, size = label_size + 0.2, color = "grey10"
        ) +
        ggplot2::coord_cartesian(xlim = c(-0.32, 1), clip = "off") +
        ggplot2::theme(plot.margin = ggplot2::margin(8, 8, 8, 8))
    }
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
  pal <- .state_palette(colors, length(state_levels))
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
    ) + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18)))
  }
  p
}


# ------------------------------------------------------------------------------
# Generic + methods (exported)
# ------------------------------------------------------------------------------

#' Plot State Frequency Distributions
#'
#' Visualise state (node) frequency distributions across groups for any
#' Nestimate object that carries sequence data: a single \code{netobject},
#' a \code{netobject_group}, an \code{mcml} model, or an \code{htna} network.
#'
#' The marimekko layout is dispatched per class:
#' \itemize{
#'   \item For \code{mcml}, where states partition cleanly into clusters,
#'     the chart is a hierarchical 2D marimekko: cluster columns of width
#'     proportional to cluster total, segments stacked vertically with
#'     heights proportional to within-cluster state proportions.
#'   \item For all other classes (\code{netobject}, \code{netobject_group},
#'     \code{htna}), each group is rendered as its own panel containing a
#'     squarified treemap: each state becomes a rectangular tile whose
#'     AREA is exactly proportional to the state's share within that group.
#'     Single-panel when no groups exist; faceted when groups are present.
#' }
#'
#' The bar style produces horizontal bars (state on the y-axis), faceted by
#' group when groups exist. All variants use the Okabe-Ito palette.
#'
#' @param x A \code{netobject}, \code{netobject_group}, \code{mcml}, or
#'   \code{htna} object.
#' @param style Either \code{"marimekko"} (default) or \code{"bars"}.
#' @param metric For \code{style = "bars"}: which value the bar length
#'   encodes -- \code{"prop"} (default) or \code{"freq"}. Treemap and
#'   hierarchical-marimekko areas always encode proportion within group.
#' @param label Inline tile / bar annotation. All formats render on a
#'   single line.
#'   \itemize{
#'     \item \code{"prop"} (default) -- proportion only, e.g. \code{"66\%"}
#'     \item \code{"freq"} -- count only, e.g. \code{"1,234"}
#'     \item \code{"both"} -- count + proportion, e.g. \code{"1,234 (66\%)"}
#'     \item \code{"state"} -- state name only, e.g. \code{"Average"}
#'     \item \code{"all"} -- state + proportion, e.g. \code{"Average (66\%)"}
#'     \item \code{"none"} -- no inline labels
#'   }
#' @param legend Legend position. \code{"bottom"} (default), \code{"top"},
#'   \code{"right"}, \code{"left"}, \code{"none"} (hide; pair with
#'   \code{label = "state"} or \code{"all"} so state names show on tiles),
#'   or \code{"per_facet"} -- each group renders as its own panel with
#'   an isolated legend showing only the states present in that panel.
#'   For \code{htna} this gives each actor (AI, Human) its own legend;
#'   for \code{mcml} each cluster gets its own. Requires the
#'   \pkg{gridExtra} package and returns a \code{gtable}.
#' @param legend_dir Legend internal layout: \code{"auto"} (default --
#'   horizontal for top/bottom, vertical for left/right), or force
#'   \code{"horizontal"} or \code{"vertical"} regardless of position.
#' @param legend_frame \code{"none"} (default) for an unframed legend, or
#'   \code{"border"} to draw a thin grey rectangle around the legend
#'   ("legend enclosed in a square").
#' @param sort_states One of \code{"frequency"} (default -- most frequent
#'   first), \code{"alpha"}, or \code{"none"}.
#' @param colors Optional character vector overriding the default
#'   Okabe-Ito state palette. Length must be at least the number of unique
#'   states.
#' @param label_size Numeric size of inline labels.
#' @param include_macro For \code{mcml} only: prepend a \code{"macro"}
#'   reference column showing aggregate state frequencies across all
#'   clusters. Default \code{FALSE}.
#' @param combine For \code{legend = "per_facet"} only. When
#'   \code{TRUE} (default), per-panel ggplots are arranged into a single
#'   gtable via \pkg{gridExtra}. When \code{FALSE}, returns a list of
#'   ggplots that are printed one-per-figure by knitr (so each panel
#'   uses the full chunk \code{fig.width} / \code{fig.height}).
#' @param ncol For \code{legend = "per_facet"} with \code{combine = TRUE}:
#'   number of columns in the grid arrangement. \code{NULL} (default)
#'   picks 1, 2, or 3 columns based on the number of panels.
#' @param ... Reserved for future use.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   data(group_regulation_long, package = "Nestimate")
#'   nw <- build_network(group_regulation_long,
#'                       method = "relative", format = "long",
#'                       actor = "Actor", action = "Action",
#'                       order = "Time", group = "Course")
#'   plot_state_frequencies(nw)
#'   plot_state_frequencies(nw, style = "bars")
#' }
#' }
#' @export
plot_state_frequencies <- function(x,
                                    style       = c("marimekko", "bars", "mosaic"),
                                    metric      = c("prop", "freq"),
                                    label       = c("prop", "freq", "both",
                                                    "state", "all", "none"),
                                    legend      = c("bottom", "right",
                                                    "top", "left", "none",
                                                    "per_facet"),
                                    legend_dir  = c("auto", "horizontal", "vertical"),
                                    legend_frame = c("none", "border"),
                                    sort_states = c("frequency", "alpha", "none"),
                                    colors      = NULL,
                                    label_size  = 3,
                                    include_macro = FALSE,
                                    combine     = TRUE,
                                    ncol        = NULL,
                                    node_groups = NULL,
                                    ...) {
  UseMethod("plot_state_frequencies")
}


#' @export
#' @rdname plot_state_frequencies
plot_state_frequencies.netobject <- function(x,
                                              style       = c("marimekko", "bars", "mosaic"),
                                              metric      = c("prop", "freq"),
                                              label       = c("prop", "freq", "both",
                                                              "state", "all", "none"),
                                              legend      = c("bottom", "right",
                                                              "top", "left", "none",
                                                              "per_facet"),
                                              legend_dir  = c("auto", "horizontal", "vertical"),
                                              legend_frame = c("none", "border"),
                                              sort_states = c("frequency", "alpha", "none"),
                                              colors      = NULL,
                                              label_size  = 3,
                                              include_macro = FALSE,
                                              combine     = TRUE,
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
  freq_df <- .freq_df_netobject(x)
  .render_freq(freq_df, style, hierarchical = FALSE,
               metric, label, sort_states, colors, label_size,
               legend, legend_dir, legend_frame,
               combine = combine, ncol = ncol,
               node_groups = node_groups)
}


#' @export
#' @rdname plot_state_frequencies
plot_state_frequencies.htna <- function(x,
                                         style       = c("marimekko", "bars", "mosaic"),
                                         metric      = c("prop", "freq"),
                                         label       = c("prop", "freq", "both",
                                                         "state", "all", "none"),
                                         legend      = c("per_facet", "bottom",
                                                         "right", "top", "left",
                                                         "none"),
                                         legend_dir  = c("auto", "horizontal", "vertical"),
                                         legend_frame = c("border", "none"),
                                         sort_states = c("frequency", "alpha", "none"),
                                         colors      = NULL,
                                         label_size  = 3,
                                         include_macro = FALSE,
                                         combine     = TRUE,
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
  freq_df <- .freq_df_htna(x)
  .render_freq(freq_df, style, hierarchical = FALSE,
               metric, label, sort_states, colors, label_size,
               legend, legend_dir, legend_frame,
               combine = combine, ncol = ncol,
               node_groups = node_groups)
}


#' @export
#' @rdname plot_state_frequencies
plot_state_frequencies.mcml <- function(x,
                                         style       = c("marimekko", "bars", "mosaic"),
                                         metric      = c("prop", "freq"),
                                         label       = c("prop", "freq", "both",
                                                         "state", "all", "none"),
                                         legend      = c("per_facet", "bottom",
                                                         "right", "top", "left",
                                                         "none"),
                                         legend_dir  = c("auto", "horizontal", "vertical"),
                                         legend_frame = c("none", "border"),
                                         sort_states = c("frequency", "alpha", "none"),
                                         colors      = NULL,
                                         label_size  = 3,
                                         include_macro = FALSE,
                                         combine     = TRUE,
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
  freq_df <- .freq_df_mcml(x, include_macro = include_macro)
  .render_freq(freq_df, style, hierarchical = TRUE,
               metric, label, sort_states, colors, label_size,
               legend, legend_dir, legend_frame,
               combine = combine, ncol = ncol,
               node_groups = node_groups)
}


#' @export
#' @rdname plot_state_frequencies
plot_state_frequencies.netobject_group <- function(x,
                                                    style       = c("marimekko", "bars", "mosaic"),
                                                    metric      = c("prop", "freq"),
                                                    label       = c("prop", "freq", "both",
                                                                    "state", "all", "none"),
                                                    legend      = c("bottom", "right",
                                                                    "top", "left", "none",
                                                                    "per_facet"),
                                                    legend_dir  = c("auto", "horizontal", "vertical"),
                                                    legend_frame = c("none", "border"),
                                                    sort_states = c("frequency", "alpha", "none"),
                                                    colors      = NULL,
                                                    label_size  = 3,
                                                    include_macro = FALSE,
                                                    combine     = TRUE,
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
  freq_df <- .freq_df_netobject_group(x)
  .render_freq(freq_df, style, hierarchical = FALSE,
               metric, label, sort_states, colors, label_size,
               legend, legend_dir, legend_frame,
               combine = combine, ncol = ncol,
               node_groups = node_groups)
}


#' @export
#' @rdname plot_state_frequencies
plot_state_frequencies.default <- function(x,
                                            style       = c("marimekko", "bars", "mosaic"),
                                            metric      = c("prop", "freq"),
                                            label       = c("prop", "freq", "both",
                                                            "state", "all", "none"),
                                            legend      = c("bottom", "right",
                                                            "top", "left", "none"),
                                            sort_states = c("frequency", "alpha", "none"),
                                            colors      = NULL,
                                            label_size  = 3,
                                            include_macro = FALSE,
                                            ...) {
  cls <- paste(class(x), collapse = "/")
  stop("plot_state_frequencies(): no method for class '", cls, "'.\n",
       "Supported: netobject, netobject_group, mcml, htna.\n",
       "For tna objects, use tna::plot_frequencies(x).",
       call. = FALSE)
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
    stop("No state observations available to plot.", call. = FALSE)
  }

  # Per-facet legend mode applies to the marimekko/treemap path only.
  # For style = "bars", per_facet doesn't make sense (bars already
  # facet-wrap with shared legend), so fall back to a right-side legend.
  if (identical(legend, "per_facet")) {
    if (style == "bars") {
      return(.plot_state_bars(freq_df, sort_states, colors,
                              label, label_size, metric,
                              "right", legend_dir, legend_frame))
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
    return(.plot_state_bars(freq_df, sort_states, colors,
                            label, label_size, metric,
                            legend, legend_dir, legend_frame))
  }

  if (style == "mosaic") {
    return(.plot_state_residuals(freq_df, sort_states, label, label_size,
                                  legend, legend_dir, legend_frame,
                                  node_groups = node_groups))
  }

  if (isTRUE(hierarchical)) {
    .plot_marimekko_hierarchical(freq_df, sort_states, colors,
                                  label, label_size,
                                  legend, legend_dir, legend_frame)
  } else {
    .plot_treemap_panels(freq_df, sort_states, colors,
                          label, label_size,
                          legend, legend_dir, legend_frame)
  }
}
