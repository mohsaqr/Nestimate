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
.seqplot_from_long <- function(x, actor, action, time, order, periods,
                               group = NULL) {
  if (!is.data.frame(x)) {
    .period_bad_input(sprintf(paste0(
      "`action =` reads long data, so `x` must be a data frame; got '%s'."),
      class(x)[1L]))
  }
  cols <- unique(c(actor, action, time, order, periods, group))
  missing_cols <- setdiff(cols, names(x))
  if (length(missing_cols) > 0L) {
    .period_bad_input(sprintf("Column(s) not in `x`: %s.",
                              paste(missing_cols, collapse = ", ")))
  }
  # prepare() would reject these too, but as an "actor/session" error; name
  # the argument the caller actually used.
  id_cols <- c(actor, periods)
  has_na <- id_cols[vapply(id_cols, function(col) anyNA(x[[col]]), logical(1))]
  if (length(has_na) > 0L) {
    .period_bad_input(sprintf(paste0(
      "Column(s) %s have missing values; every row needs an actor%s. ",
      "Drop or relabel these rows first."),
      paste(sprintf("`%s`", has_na), collapse = " and "),
      if (is.null(periods)) "" else " and a period"))
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
.period_source <- function(src, periods, group = NULL) {
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
  list(data = src$data, actor = who, period = per,
       group = .period_group(meta, who, group))
}

# With `periods`, `group` names a metadata column (a covariate of the long
# data) that is constant within actor; the plot gets one panel per value.
# Returns one value per sequence, aligned with `who`, or NULL.
#' @noRd
.period_group <- function(meta, who, group) {
  if (is.null(group)) return(NULL)
  if (!is.character(group) || length(group) != 1L || is.na(group)) {
    .period_bad_input(paste0(
      "With `periods`, `group` names one column of the long data (or of the ",
      "network's metadata) that is constant within actor, not a vector."))
  }
  if (!group %in% names(meta)) {
    .period_bad_input(sprintf(
      "`group = \"%s\"` is not a column of the data (%s).", group,
      paste(setdiff(names(meta), c(".session_label", ".session_id")),
            collapse = ", ")))
  }
  g <- meta[[group]]
  g <- if (is.factor(g)) droplevels(g) else as.vector(g)
  if (anyNA(g)) {
    .period_bad_input(sprintf("`%s` has missing values.", group))
  }
  n_per_actor <- tapply(as.character(g), as.character(who),
                        function(v) length(unique(v)))
  if (any(n_per_actor > 1L)) {
    .period_bad_input(sprintf(paste0(
      "`group = \"%s\"` must be constant within actor; %d actor(s) have ",
      "more than one value (e.g. %s)."), group, sum(n_per_actor > 1L),
      names(n_per_actor)[which(n_per_actor > 1L)[1L]]))
  }
  g
}

# The period axis. Factor levels (an unused level is an empty block); every
# day from the first to the last for a Date; every whole number from the
# smallest to the largest for a day counter; otherwise the sorted observed
# values. So a day on which nobody acted stays on the axis as white space.
#' @noRd
# `n_max` (a count cut from `trim`) bounds the expansion so a wide range is
# never materialised beyond what is drawn; past 100,000 blocks the axis is
# unreadable and the request is refused.
.period_levels <- function(period, n_max = NULL) {
  if (is.factor(period)) return(levels(period))
  gridded <- inherits(period, "Date") ||
    (is.numeric(period) && all(period == round(period)))
  if (!gridded) return(as.character(sort(unique(period))))
  first <- min(period)
  last  <- max(period)
  n <- as.numeric(last - first) + 1
  if (!is.null(n_max)) {
    n    <- min(n, n_max)
    last <- first + (n - 1)
  }
  if (n > 1e5) {
    .period_bad_input(sprintf(paste0(
      "`periods` spans %s blocks from %s to %s; give `trim` a count of ",
      "periods to draw, or a label or factor column."),
      format(n, big.mark = ","), format(first), format(last)))
  }
  as.character(seq(first, last, by = 1))
}

# One row per actor, that actor's states in period-then-step order, NA
# padded to the longest; an actor with no cell is an all-NA row.
#' @noRd
.period_trajectory <- function(a, p, step, code, n_actors) {
  o   <- order(a, p, step)
  pos <- stats::ave(a[o], a[o], FUN = seq_along)
  trajectory <- matrix(NA_integer_, n_actors, max(pos, 1L))
  trajectory[cbind(a[o], pos)] <- code[o]
  trajectory
}

# Cell geometry for both fill modes.
#' @noRd
.period_layout <- function(data, actor, period, fill, trim = NULL) {
  enc   <- .encode_states(data)
  codes <- enc$codes
  tr <- .trim_axes(trim)
  period_levels <- .period_levels(
    period, n_max = if (!is.null(tr$cols) && tr$cols >= 1) tr$cols)
  actors <- if (is.factor(actor)) levels(droplevels(actor)) else
    as.character(sort(unique(actor)))

  idx <- which(!is.na(codes), arr.ind = TRUE)
  idx <- idx[order(idx[, "row"], idx[, "col"]), , drop = FALSE]
  r    <- idx[, "row"]
  step <- stats::ave(r, r, FUN = seq_along)
  p    <- match(as.character(period), period_levels)[r]
  a    <- match(as.character(actor), actors)[r]
  code <- codes[idx]
  # A period past a count-bounded axis is outside the drawn range.
  if (anyNA(p)) {
    inside <- !is.na(p)
    r <- r[inside]; step <- step[inside]; p <- p[inside]; a <- a[inside]
    code <- code[inside]
  }

  # The time half of `trim` cuts the period axis as it cuts time points
  # elsewhere: a count keeps the first n periods, a fraction keeps the
  # periods up to that quantile of the actors' last active period.
  if (!is.null(tr$cols)) {
    present <- matrix(NA_integer_, length(actors), length(period_levels))
    present[cbind(a, p)] <- 1L
    cut <- .trim_cut(present, tr$cols)
    if (cut < length(period_levels)) {
      keep <- p <= cut
      if (!any(keep)) {
        .period_bad_input(sprintf(
          "`trim = %s` keeps %d period(s) with no action in them.",
          format(tr$cols), cut))
      }
      r <- r[keep]; step <- step[keep]; p <- p[keep]; a <- a[keep]
      code <- code[keep]
      period_levels <- period_levels[seq_len(cut)]
    }
  }
  # The row half drops the least active actors, counted over the kept
  # periods; the remaining actors are renumbered so `a` stays an index
  # into `actors`.
  if (!is.null(tr$rows)) {
    keep_actor <- .trim_rows_keep(
      .period_trajectory(a, p, step, code, length(actors)), tr$rows)
    if (!all(keep_actor)) {
      kept    <- which(keep_actor)
      in_keep <- a %in% kept
      r <- r[in_keep]; step <- step[in_keep]; p <- p[in_keep]
      code <- code[in_keep]
      a <- match(a[in_keep], kept)
      actors <- actors[kept]
    }
  }

  n_steps <- rowSums(!is.na(codes))
  width   <- as.integer(max(n_steps[unique(r)]))
  denom   <- if (fill == "stretch") n_steps[r] else width

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
  trajectory <- .period_trajectory(a, p, step, code, length(actors))

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
                                   legend_bty, trim = NULL, group = NULL,
                                   ncol = NULL, nrow = NULL, combined = TRUE) {
  stopifnot(
    "`periods` must be a single column name" =
      is.character(periods) && length(periods) == 1L && !is.na(periods),
    "`row_gap` must be a single number in [0, 1)" =
      is.numeric(row_gap) && length(row_gap) == 1L && row_gap >= 0 &&
      row_gap < 1
  )
  src <- .period_source(src, periods, group)
  lay <- .period_layout(src$data, src$actor, src$period, fill, trim)
  levels_all <- lay$levels
  n_rows <- length(lay$actors)
  n_per  <- length(lay$period_levels)
  .report_unused_colors(state_colors, levels_all)
  palette <- .state_palette(state_colors, length(levels_all), levels_all)
  gap <- if (type == "heatmap") 0 else row_gap

  # One panel per group value (index type), else one panel of all actors.
  # `actor_group` is one value per actor, in `lay$actors` order.
  actor_group <- if (is.null(src$group)) {
    factor(rep("all", n_rows))
  } else {
    g <- src$group[match(lay$actors, as.character(src$actor))]
    # droplevels: a row trim can empty a level, which must not be a panel.
    if (is.factor(src$group)) droplevels(factor(g, levels = levels(src$group)))
    else factor(g)
  }
  group_levels <- levels(actor_group)
  G <- length(group_levels)

  # Rows of each panel in plotted order (indices into lay$actors). The
  # dendrogram, when drawn, orders all actors at once.
  panels <- lapply(group_levels, function(g) {
    rows <- which(actor_group == g)
    ord  <- .row_order(lay$trajectory[rows, , drop = FALSE], sort,
                       if (G == 1L) tree else NULL, levels_all)
    structure(rows[ord], tree = attr(ord, "tree"))
  })
  tree <- if (type == "heatmap" && G == 1L) attr(panels[[1L]], "tree") else
    NULL
  if (G == 1L && length(panels[[1L]]) != n_rows) {
    .period_bad_input(sprintf(paste0(
      "`tree` has %d leaves but there are %d actors",
      if (!is.null(.trim_axes(trim)$rows)) {
        " (the row half of `trim` drops actors; it cannot be combined with `tree`)."
      } else "."), length(panels[[1L]]), n_rows))
  }

  if (is.null(ncol) && is.null(nrow)) {
    ncol <- as.integer(ceiling(sqrt(G)))
    nrow <- as.integer(ceiling(G / ncol))
  } else if (is.null(ncol)) {
    ncol <- as.integer(ceiling(G / nrow))
  } else if (is.null(nrow)) {
    nrow <- as.integer(ceiling(G / ncol))
  }
  stopifnot("`ncol * nrow` must hold every group" = ncol * nrow >= G)

  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op), add = TRUE)
  oma <- .legend_oma_size(levels_all, legend, legend_size,
                          legend_ncol, legend_title)
  set_oma <- function() {
    graphics::par(oma = c(if (legend == "bottom") oma[["oma_b"]] else 0.3,
                          0.3, 0.3,
                          if (legend == "right")  oma[["oma_r"]] else 0.3))
  }
  set_oma()
  use_layout <- combined && G > 1L
  if (!is.null(tree)) {
    graphics::layout(matrix(c(1L, 2L), nrow = 1L),
                     widths = c(dendrogram_width, 4))
  } else if (use_layout) {
    graphics::layout(matrix(c(seq_len(G), integer(nrow * ncol - G)),
                            nrow = nrow, ncol = ncol, byrow = TRUE))
  }
  mar_bottom <- if (!is.null(time_label) && nzchar(time_label)) 3 else 1
  mar_top    <- if (!is.null(main) || isTRUE(show_n) || G > 1L) 2.5 else 1.5
  if (!is.null(tree)) {
    graphics::par(mar = c(mar_bottom, 0.5, mar_top, 0.3))
    graphics::plot(tree, horiz = TRUE, axes = FALSE, yaxs = "i",
                   ylim = c(0.5, n_rows + 0.5),
                   leaflab = "none", xlab = "", ylab = "")
  }

  lapply(seq_len(G), function(g_idx) {
    ord <- as.integer(panels[[g_idx]])
    if (!combined && g_idx > 1L) set_oma()
    graphics::par(mar = c(mar_bottom, if (!is.null(tree)) 0 else 2,
                          mar_top, 1))
    .period_panel(lay, ord, palette, gap, na_color, cell_border, frame,
                  fill, tick, time_label)
    if (!is.null(k) && !is.null(tree)) {
      stopifnot(is.numeric(k), length(k) == 1L, k >= 2L, k < n_rows)
      mem <- stats::cutree(stats::as.hclust(tree), k = as.integer(k))[ord]
      boundaries <- which(diff(mem) != 0L)
      if (length(boundaries) > 0L) {
        graphics::abline(h = boundaries + 0.5, col = k_color,
                         lwd = k_line_width)
      }
    }
    t_line <- if (G > 1L) {
      .title_line(group_levels[g_idx], show_n, length(ord))
    } else .title_line(main, show_n, n_rows)
    if (!is.null(t_line)) {
      graphics::mtext(t_line, side = 3, line = if (G > 1L) 0.5 else 1,
                      font = 2, cex = if (G > 1L) 0.9 else 1)
    }
    if (!combined && legend != "none") {
      .draw_legend_in_oma(levels_all, palette, legend, legend_size,
                          legend_ncol, legend_title, legend_border,
                          legend_bty)
    }
  })
  if (combined && G > 1L && !is.null(main)) {
    graphics::mtext(main, side = 3, line = 1.5, font = 2, outer = TRUE)
  }
  if (combined && legend != "none") {
    .draw_legend_in_oma(levels_all, palette, legend, legend_size,
                        legend_ncol, legend_title, legend_border, legend_bty)
  }

  out <- list(ord = unlist(lapply(panels, as.integer)), codes = lay$grid,
              palette = palette, levels = levels_all,
              periods = lay$period_levels, fill = fill)
  if (!is.null(src$group)) {
    out$groups <- group_levels
    out$orders <- lapply(panels, as.integer)
  }
  invisible(out)
}

# Draw one panel: the actors `ord` (indices into lay$actors, top to bottom
# in plotted order) over the full period axis, with the period axis.
#' @noRd
.period_panel <- function(lay, ord, palette, gap, na_color, cell_border,
                          frame, fill, tick, time_label) {
  n_per <- length(lay$period_levels)
  n_row <- length(ord)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, n_per), ylim = c(0.5, n_row + 0.5),
                        xaxs = "i", yaxs = "i")
  graphics::rect(0, 0.5, n_per, n_row + 0.5, col = na_color, border = NA)

  cells <- lay$cells[lay$cells$actor %in% ord, , drop = FALSE]
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
  invisible(NULL)
}
