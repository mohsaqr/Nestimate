# ==============================================================================
# sequence_plot(periods = ) - one row per actor, the time axis split into
# equal-width period blocks (days, weeks, ...). An actor's steps in a period
# either stretch to fill the block (fill = "stretch", each actor its own cell
# width) or sit in fixed slots as wide as the busiest actor-period
# (fill = "pad", empty slots NA). Both layouts share one coordinate system:
# period p spans [p - 1, p] and step j of n is [p-1 + (j-1)/d, p-1 + j/d],
# with d = n (stretch) or the common slot count (pad).
# ==============================================================================

#' @noRd
.period_bad_input <- function(msg) {
  stop(errorCondition(msg, class = "nestimate_bad_input", call = NULL))
}

# Wide sequences + per-sequence metadata from long data (one row per action),
# via prepare(). With `periods` each actor-period is one sequence (no
# gap-splitting); without, prepare()'s own session rules apply, as in
# build_network(). Only the named columns go in, so prepare() has no other
# metadata to aggregate.
#' @noRd
.seqplot_from_long <- function(x, actor, action, time, order, periods) {
  if (!is.data.frame(x)) {
    .period_bad_input(sprintf(paste0(
      "`action =` reads long data, so `x` must be a data frame; got '%s'."),
      class(x)[1L]))
  }
  cols <- unique(c(actor, action, time, order, periods))
  missing_cols <- setdiff(cols, names(x))
  if (length(missing_cols) > 0L) {
    .period_bad_input(sprintf("Column(s) not in `x`: %s.",
                              paste(missing_cols, collapse = ", ")))
  }
  prep <- prepare(x[cols], actor = actor, action = action, time = time,
                  order = order, session = periods,
                  time_threshold = if (is.null(periods)) 900 else FALSE)
  list(data = prep$sequence_data, metadata = prep$meta_data, actor = actor)
}

# The same three pieces from a netobject built with `actor =` and `session =`.
#' @noRd
.seqplot_from_network <- function(x, periods) {
  if (!inherits(x, "netobject")) {
    .period_bad_input(sprintf(paste0(
      "`periods` needs long data (with `actor =` and `action =`) or a ",
      "netobject from build_network() with `actor =` and `session =`; ",
      "got '%s'."), class(x)[1L]))
  }
  list(data = x$data, metadata = x$metadata, actor = x$build_args$actor)
}

# Check that the source has one sequence per actor-period and return the
# sequences with their actor and period.
#' @noRd
.period_source <- function(src, periods) {
  actor <- src$actor
  meta  <- src$metadata
  if (is.null(src$data) || is.null(meta) || is.null(actor) ||
      length(actor) != 1L || !actor %in% names(meta)) {
    .period_bad_input(paste0(
      "`periods` needs sequences split by actor: give `actor =` (for a ",
      "network, build it with `actor =` and `session = \"", periods, "\"`)."))
  }
  if (!periods %in% names(meta)) {
    .period_bad_input(sprintf(
      "`periods = \"%s\"` is not a column of the network's metadata (%s).",
      periods, paste(setdiff(names(meta), ".session_label"), collapse = ", ")))
  }
  per <- meta[[periods]]
  who <- meta[[actor]]
  if (anyNA(per) || anyNA(who)) {
    .period_bad_input(sprintf("`%s` or `%s` has missing values.",
                              periods, actor))
  }
  key <- paste(who, per, sep = "\r")
  dup <- duplicated(key)
  if (any(dup)) {
    first <- strsplit(key[dup][1L], "\r", fixed = TRUE)[[1L]]
    .period_bad_input(sprintf(paste0(
      "%d actor-period pairs hold more than one sequence (e.g. %s = %s, ",
      "%s = %s). Build with `session = \"%s\"` and `time_threshold = FALSE` ",
      "so each period is one sequence."),
      sum(dup), actor, first[1L], periods, first[2L], periods))
  }
  list(data = src$data, actor = who, period = per)
}

# Cell geometry for both fill modes. Period order: factor levels (so an
# unused level shows as an empty block), else sorted unique values.
#' @noRd
.period_layout <- function(data, actor, period, fill) {
  enc   <- .encode_states(data)
  codes <- enc$codes
  period_levels <- if (is.factor(period)) levels(period) else
    as.character(sort(unique(period)))
  actors <- if (is.factor(actor)) levels(droplevels(actor)) else
    as.character(sort(unique(actor)))

  idx <- which(!is.na(codes), arr.ind = TRUE)
  idx <- idx[order(idx[, "row"], idx[, "col"]), , drop = FALSE]
  r       <- idx[, "row"]
  step    <- stats::ave(r, r, FUN = seq_along)
  n_steps <- rowSums(!is.na(codes))
  width   <- as.integer(max(n_steps))
  denom   <- if (fill == "stretch") n_steps[r] else width
  p <- match(as.character(period), period_levels)[r]
  a <- match(as.character(actor), actors)[r]
  code <- codes[idx]

  cells <- data.frame(actor = a, period = p, step = step,
                      xleft  = p - 1 + (step - 1) / denom,
                      xright = p - 1 + step / denom,
                      code = code)

  # Aligned grid: actor x (period, slot), the pad layout as a matrix.
  grid <- matrix(NA_integer_, length(actors), length(period_levels) * width,
                 dimnames = list(actors, paste(rep(period_levels, each = width),
                                               seq_len(width), sep = ".")))
  grid[cbind(a, (p - 1L) * width + step)] <- code

  # Each actor's whole trajectory (periods in order, steps joined), used
  # only to order the rows.
  o   <- order(a, p, step)
  pos <- stats::ave(a[o], a[o], FUN = seq_along)
  trajectory <- matrix(NA_integer_, length(actors), max(pos))
  trajectory[cbind(a[o], pos)] <- code[o]

  list(cells = cells, grid = grid, trajectory = trajectory,
       actors = actors, period_levels = period_levels, width = width,
       levels = enc$levels)
}

#' @noRd
.sequence_plot_periods <- function(src, periods, fill, type, sort, tree,
                                   dendrogram_width, k, k_color, k_line_width,
                                   row_gap, state_colors, na_color,
                                   cell_border, frame, main, show_n,
                                   time_label, tick, legend, legend_size,
                                   legend_title, legend_ncol, legend_border,
                                   legend_bty) {
  stopifnot(
    "`periods` must be a single column name" =
      is.character(periods) && length(periods) == 1L && !is.na(periods),
    "`row_gap` must be a single number in [0, 1)" =
      is.numeric(row_gap) && length(row_gap) == 1L && row_gap >= 0 &&
      row_gap < 1
  )
  src <- .period_source(src, periods)
  lay <- .period_layout(src$data, src$actor, src$period, fill)
  levels_all <- lay$levels
  n_rows <- length(lay$actors)
  n_per  <- length(lay$period_levels)

  ord  <- .row_order(lay$trajectory, sort, tree, levels_all)
  tree <- if (type == "heatmap") attr(ord, "tree") else NULL
  if (length(ord) != n_rows) {
    .period_bad_input(sprintf("`tree` has %d leaves but there are %d actors.",
                              length(ord), n_rows))
  }
  .report_unused_colors(state_colors, levels_all)
  palette <- .state_palette(state_colors, length(levels_all), levels_all)
  gap <- if (type == "heatmap") 0 else row_gap

  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op), add = TRUE)
  oma <- .legend_oma_size(levels_all, legend, legend_size,
                          legend_ncol, legend_title)
  graphics::par(oma = c(if (legend == "bottom") oma[["oma_b"]] else 0.3,
                        0.3, 0.3,
                        if (legend == "right")  oma[["oma_r"]] else 0.3))
  if (!is.null(tree)) {
    graphics::layout(matrix(c(1L, 2L), nrow = 1L),
                     widths = c(dendrogram_width, 4))
  }
  mar_bottom <- if (!is.null(time_label) && nzchar(time_label)) 3 else 1
  mar_top    <- if (!is.null(main) || isTRUE(show_n)) 2.5 else 1.5
  if (!is.null(tree)) {
    graphics::par(mar = c(mar_bottom, 0.5, mar_top, 0.3))
    graphics::plot(tree, horiz = TRUE, axes = FALSE, yaxs = "i",
                   ylim = c(0.5, n_rows + 0.5),
                   leaflab = "none", xlab = "", ylab = "")
  }
  graphics::par(mar = c(mar_bottom, if (!is.null(tree)) 0 else 2,
                        mar_top, 1))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, n_per), ylim = c(0.5, n_rows + 0.5),
                        xaxs = "i", yaxs = "i")
  graphics::rect(0, 0.5, n_per, n_rows + 0.5, col = na_color, border = NA)

  cells <- lay$cells
  y <- match(cells$actor, ord)
  graphics::rect(xleft   = cells$xleft,
                 ybottom = y - 0.5 + gap / 2,
                 xright  = cells$xright,
                 ytop    = y + 0.5 - gap / 2,
                 col     = palette[cells$code],
                 border  = if (is.na(cell_border)) NA else cell_border,
                 lwd     = 0.6)
  # Period dividers, drawn over the cells so a block edge is always visible.
  if (n_per > 1L) {
    graphics::abline(v = seq_len(n_per - 1L), col = "grey55", lwd = 1)
  }
  if (isTRUE(frame)) graphics::box()

  # Ticks mark period boundaries; labels sit at block centres. With
  # fill = "pad" every slot is a fixed position, so it gets a minor tick,
  # up to 60 slots in all (past that the ticks merge into a solid band).
  graphics::axis(1, at = 0:n_per, labels = FALSE,
                 lwd = if (isTRUE(frame)) 1 else 0, lwd.ticks = 1)
  if (fill == "pad" && lay$width > 1L && n_per * lay$width <= 60L) {
    slots <- outer((seq_len(lay$width) - 0.5) / lay$width,
                   seq_len(n_per) - 1L, `+`)
    graphics::axis(1, at = as.vector(slots), labels = FALSE, lwd = 0,
                   lwd.ticks = 0.6, tcl = -0.2)
  }
  labs <- .tick_positions(n_per, tick, lay$period_levels)
  graphics::axis(1, at = labs$at - 0.5, labels = labs$labels, tick = FALSE,
                 mgp = c(1, 0.5, 0))
  if (!is.null(time_label) && nzchar(time_label)) {
    graphics::mtext(time_label, side = 1, line = 1.7)
  }

  if (!is.null(k) && !is.null(tree)) {
    stopifnot(is.numeric(k), length(k) == 1L, k >= 2L, k < n_rows)
    mem <- stats::cutree(stats::as.hclust(tree), k = as.integer(k))[ord]
    boundaries <- which(diff(mem) != 0L)
    if (length(boundaries) > 0L) {
      graphics::abline(h = boundaries + 0.5, col = k_color, lwd = k_line_width)
    }
  }
  t_line <- .title_line(main, show_n, n_rows)
  if (!is.null(t_line)) graphics::mtext(t_line, side = 3, line = 1, font = 2)
  if (legend != "none") {
    .draw_legend_in_oma(levels_all, palette, legend, legend_size,
                        legend_ncol, legend_title, legend_border, legend_bty)
  }

  invisible(list(ord = ord, codes = lay$grid, palette = palette,
                 levels = levels_all, periods = lay$period_levels,
                 fill = fill))
}
