# magnitude_difference() — compare the probability view (TNA) of a
# transition network against the frequency view (FTNA) on a single dataset.
#
# The two matrices encode the same data in different units: raw counts
# (FTNA) versus row-conditional probabilities (TNA). Row-normalization
# promotes rare-source transitions and demotes dense-source ones. This
# module quantifies that per-edge discrepancy ("magnitude difference")
# on a common scale, with print and plot methods.

# Bare column names referenced inside ggplot2::aes() in the plot helpers;
# declared so R CMD check does not flag "no visible binding for global
# variable" (the package-wide convention, mirroring R/Nestimate-package.R).
utils::globalVariables(c(
  "x", "y", "hjust", "from_x", "from_y", "to_x", "to_y", "grp", "lx", "ly",
  "x_lab", "y_lab", "text_angle", "lab_col", "color"
))

# Five metrics, all computed after the two weight matrices are placed on a
# common scale (see `scale`). Each passed |Spearman rho vs. bootstrap
# stability| > 0.5 in the 1000-simulation benchmark.
.magdiff_metrics <- c("abs_diff", "chord_dist", "atanh_diff",
                      "geom_norm_diff", "cv_inflation")

# Per-edge discrepancy between reference (FTNA) and estimate (TNA) vectors.
.magdiff_compute <- function(ref, est, metric) {
  guard <- 1e-3
  switch(metric,
    abs_diff       = abs(est - ref),
    chord_dist     = abs(est / sqrt(sum(est^2)) - ref / sqrt(sum(ref^2))),
    atanh_diff     = abs(atanh(pmin(est, 1 - 1e-6)) -
                          atanh(pmin(ref, 1 - 1e-6))),
    geom_norm_diff = {
      gm_e <- exp(mean(log(est + guard))); gm_r <- exp(mean(log(ref + guard)))
      abs(est / gm_e - ref / gm_r)
    },
    cv_inflation   = abs(est / max(stats::sd(est), guard) -
                          ref / max(stats::sd(ref), guard))
  )
}

# Min-max scale a matrix's values to [0, 1].
.magdiff_minmax <- function(M) {
  v <- as.vector(M); rng <- range(v, finite = TRUE)
  if (diff(rng) == 0) return(M)
  out <- M; out[] <- (v - rng[1]) / diff(rng); out
}

# Rank-min-max: rank each matrix's values, scale ranks to [0, 1]. Puts FTNA
# (counts) and TNA (rates) on the same comparable grid despite different
# native units.
.magdiff_rank_minmax <- function(M) {
  v <- as.vector(M)
  rks <- rank(v, ties.method = "average")
  out <- M
  out[] <- (rks - 1) / (length(v) - 1)
  out
}

# Rescale a matrix linearly into the [min, max] range of a reference matrix.
# Maps FTNA (counts) into TNA's probability range so |TNA - FTNA_rescaled|
# reads in probability units. Default scaling for magnitude_difference().
.magdiff_rescale_to <- function(M, ref) {
  v <- as.vector(M)
  ref_v <- as.vector(ref)
  rng_M <- range(v, finite = TRUE)
  rng_R <- range(ref_v, finite = TRUE)
  if (diff(rng_M) == 0) return(M)
  out <- M
  out[] <- rng_R[1] + (v - rng_M[1]) / diff(rng_M) * diff(rng_R)
  out
}

#' Magnitude difference between the frequency and probability views
#'
#' Quantifies, per edge, how much row-normalization moves a transition
#' network between its two natural summaries: raw transition counts
#' (frequency / FTNA, `build_network(method = "frequency")`) and
#' row-conditional probabilities (TNA, `build_network(method = "relative")`).
#' The two matrices rank edges differently — an edge that is large in counts
#' can be modest in probability, and a rare-source edge can dominate its row
#' in probability. The per-edge discrepancy on a common scale is the
#' *magnitude difference*.
#'
#' @param data Long- or wide-format event log (`data.frame`).
#' @param actor,action,time Column names in `data` (long format). `time` may
#'   be `NULL`.
#' @param metric Discrepancy metric. One of `"abs_diff"` (default; absolute
#'   difference, simplest and most robust), `"chord_dist"` (chord distance on
#'   the unit sphere), `"atanh_diff"` (Fisher z-style on a bounded scale),
#'   `"geom_norm_diff"` (geometric-mean normalized; amplifies small edges),
#'   `"cv_inflation"` (per-vector SD-standardized then absolute difference).
#' @param scale How the two weight matrices are placed on a common scale
#'   before differencing. `"tna_range"` (default) rescales FTNA linearly into
#'   TNA's `[min, max]` range and leaves TNA untouched, so the difference is
#'   in TNA probability units. `"rank_minmax"` converts each matrix's values
#'   to ranks scaled to `[0, 1]` (ordinal). `"minmax"` scales each matrix's
#'   raw values to `[0, 1]` separately (asymmetric — TNA's max and FTNA's max
#'   map to the same value despite differing native ranges). `"none"` uses
#'   raw weights.
#' @param format Input format passed through to [build_network()];
#'   `"auto"` (default) treats the data as wide when `action` is not a column.
#'
#' @return An object of class `"magnitude_difference"`: a list with
#'   `$edges` (per-edge `data.frame` with columns `from`, `to`, `ftna`,
#'   `tna`, `signed` = `tna - ftna`, and `value` = the chosen metric),
#'   `$metric`, `$scale`, `$weights_ftna`, `$weights_tna`, and `$states`.
#'
#' @examples
#' data(group_regulation_long, package = "Nestimate")
#' fit <- magnitude_difference(group_regulation_long,
#'                             actor = "Actor", action = "Action",
#'                             time = "Time")
#' print(fit)
#' # Edges most promoted by row-normalization (rare-source transitions):
#' head(fit$edges[order(-fit$edges$signed), c("from", "to", "ftna", "tna")])
#' \donttest{
#' plot(fit)                       # stacked polar portrait
#' plot(fit, type = "circular")    # chord-style diagram
#' }
#' @seealso [build_network()], [compare_model()]
#' @export
magnitude_difference <- function(data, actor = "Actor", action = "Action",
                                 time = NULL,
                                 metric = c("abs_diff", "chord_dist",
                                            "atanh_diff", "geom_norm_diff",
                                            "cv_inflation"),
                                 scale = c("tna_range", "rank_minmax",
                                           "minmax", "none"),
                                 format = c("auto", "long", "wide")) {
  metric <- match.arg(metric)
  scale  <- match.arg(scale)
  format <- match.arg(format)
  stopifnot(is.data.frame(data))

  # Build both views. Wide-format input (one sequence per row, columns =
  # positions) is detected by the absence of an `action` column.
  is_wide <- format == "wide" ||
             (format == "auto" && !(action %in% names(data)))
  if (is_wide) {
    args_f <- list(data = data, method = "frequency", format = "wide")
    args_t <- list(data = data, method = "relative",  format = "wide")
  } else {
    args_f <- list(data = data, actor = actor, action = action,
                   method = "frequency")
    args_t <- list(data = data, actor = actor, action = action,
                   method = "relative")
    if (!is.null(time)) { args_f$time <- time; args_t$time <- time }
  }
  net_f <- do.call(build_network, args_f)
  net_t <- do.call(build_network, args_t)

  # tna_range RESCALES FTNA into TNA's [min, max] range; TNA stays raw.
  # Every other mode scales both matrices the same way.
  if (scale == "tna_range") {
    W_t <- net_t$weights
    W_f <- .magdiff_rescale_to(net_f$weights, ref = W_t)
  } else {
    W_f <- switch(scale,
                  rank_minmax = .magdiff_rank_minmax(net_f$weights),
                  minmax      = .magdiff_minmax(net_f$weights),
                  none        = net_f$weights)
    W_t <- switch(scale,
                  rank_minmax = .magdiff_rank_minmax(net_t$weights),
                  minmax      = .magdiff_minmax(net_t$weights),
                  none        = net_t$weights)
  }
  states <- rownames(W_f)

  ed <- expand.grid(from = states, to = states, KEEP.OUT.ATTRS = FALSE,
                    stringsAsFactors = FALSE)
  ed$ftna   <- as.vector(W_f)
  ed$tna    <- as.vector(W_t)
  ed$signed <- ed$tna - ed$ftna                 # direction (TNA - FTNA)
  ed$value  <- .magdiff_compute(ed$ftna, ed$tna, metric)

  out <- list(
    edges        = ed,
    metric       = metric,
    scale        = scale,
    weights_ftna = W_f,
    weights_tna  = W_t,
    states       = states
  )
  class(out) <- "magnitude_difference"
  out
}

#' @describeIn magnitude_difference Print a compact summary of the per-edge
#'   magnitude-difference distribution.
#' @param x A `magnitude_difference` object.
#' @param ... Passed to plotting helpers (ignored by `print`).
#' @return `print` invisibly returns `x`.
#' @export
print.magnitude_difference <- function(x, ...) {
  cat("magnitude_difference object\n")
  cat("  metric:", x$metric, "  scale:", x$scale, "\n")
  cat("  states (", length(x$states), "):  ",
      paste(x$states, collapse = ", "), "\n", sep = "")
  cat("  edges:", nrow(x$edges), "\n")
  cat("  magnitude-difference summary (", x$metric, "):\n", sep = "")
  print(summary(x$edges$value))
  invisible(x)
}

# --- diverging palette for the circular (chord) plot ----------------------
# Okabe-Ito blue <-> vermillion, with light/dark anchors.
.magdiff_blue_dark   <- "#003D63"; .magdiff_blue        <- "#0072B2"
.magdiff_blue_pale   <- "#B4D5EC"; .magdiff_orange      <- "#D55E00"
.magdiff_orange_dark <- "#9C3A00"; .magdiff_orange_pale <- "#FAD7B4"

.magdiff_circular <- function(ed, value_col = "signed", abs_col = "value",
                              title = NULL,
                              self_loops = TRUE,
                              label_offset = 1.15,
                              edge_curve = 0.18,
                              min_show = 0.01) {
  ed$value <- ed[[value_col]]
  ed$abs   <- ed[[abs_col]]
  ed <- ed[ed$abs >= min_show * max(ed$abs, na.rm = TRUE), , drop = FALSE]

  states <- sort(unique(c(ed$from, ed$to)))
  n <- length(states)
  ang <- seq(pi / 2, pi / 2 - 2 * pi, length.out = n + 1)[-(n + 1)]
  nodes <- data.frame(node = states,
                      x = cos(ang), y = sin(ang),
                      lx = cos(ang) * label_offset,
                      ly = sin(ang) * label_offset,
                      hjust = ifelse(cos(ang) > 0.01, 0,
                                ifelse(cos(ang) < -0.01, 1, 0.5)),
                      stringsAsFactors = FALSE)
  ed$from_x <- nodes$x[match(ed$from, nodes$node)]
  ed$from_y <- nodes$y[match(ed$from, nodes$node)]
  ed$to_x   <- nodes$x[match(ed$to,   nodes$node)]
  ed$to_y   <- nodes$y[match(ed$to,   nodes$node)]
  ed$is_self <- ed$from == ed$to
  ed <- ed[order(ed$abs), , drop = FALSE]
  non_self <- ed[!ed$is_self, , drop = FALSE]
  self_e   <- ed[ ed$is_self, , drop = FALSE]
  max_abs <- max(abs(ed$value), na.rm = TRUE)
  if (max_abs <= 0) max_abs <- 1

  p <- ggplot2::ggplot() +
    ggplot2::annotate("path",
                       x = cos(seq(0, 2 * pi, length.out = 200)),
                       y = sin(seq(0, 2 * pi, length.out = 200)),
                       color = "grey92", linewidth = 0.4) +
    ggplot2::geom_curve(data = non_self,
      ggplot2::aes(x = from_x, y = from_y, xend = to_x, yend = to_y,
                    color = value, linewidth = abs, alpha = abs),
      curvature = edge_curve, lineend = "round",
      arrow = ggplot2::arrow(length = ggplot2::unit(0.018, "npc"),
                              type = "closed"))

  if (self_loops && nrow(self_e) > 0) {
    self_e$cx <- self_e$from_x * 1.10
    self_e$cy <- self_e$from_y * 1.10
    self_e$r  <- 0.06
    theta <- seq(0, 2 * pi, length.out = 60)
    loop_paths <- do.call(rbind, lapply(seq_len(nrow(self_e)), function(i) {
      data.frame(grp = i,
                 x = self_e$cx[i] + self_e$r[i] * cos(theta),
                 y = self_e$cy[i] + self_e$r[i] * sin(theta),
                 value = self_e$value[i], abs = self_e$abs[i])
    }))
    p <- p + ggplot2::geom_path(data = loop_paths,
      ggplot2::aes(x = x, y = y, group = grp,
                    color = value, linewidth = abs, alpha = abs))
  }

  p +
    ggplot2::geom_point(data = nodes, ggplot2::aes(x = x, y = y),
                        size = 4, color = "grey30", fill = "white",
                        shape = 21, stroke = 0.6) +
    ggplot2::geom_text(data = nodes,
      ggplot2::aes(x = lx, y = ly, label = node, hjust = hjust),
      size = 3.6, family = "sans") +
    ggplot2::scale_color_gradientn(
      name = "signed\nmagnitude\ndifference",
      colors = c(.magdiff_blue_dark, .magdiff_blue, .magdiff_blue_pale,
                 "grey92", .magdiff_orange_pale, .magdiff_orange,
                 .magdiff_orange_dark),
      values = scales::rescale(c(-max_abs, -max_abs * 0.6,
                                  -max_abs * 0.2, 0,
                                   max_abs * 0.2,  max_abs * 0.6, max_abs)),
      limits = c(-max_abs, max_abs),
      guide  = ggplot2::guide_colorbar(barheight = 9, barwidth = 0.7,
                                        ticks.colour = "grey50")) +
    ggplot2::scale_linewidth(range = c(0.25, 2.4), guide = "none") +
    ggplot2::scale_alpha(range = c(0.15, 1), guide = "none") +
    ggplot2::coord_equal(xlim = c(-1.35, 1.35), ylim = c(-1.35, 1.35),
                         clip = "off") +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(legend.position = "right",
                    plot.title = ggplot2::element_text(hjust = 0.5, size = 13),
                    plot.margin = ggplot2::margin(8, 8, 8, 8))
}

# --- stacked-bar polar portrait -------------------------------------------
# Each from-state owns a sector; within it, each to-state is a wedge.
#   r_inner          - base of the bar
#   to_r(min(F, T))  - grey "shared" portion top (the two views agree here)
#   to_r(max(F, T))  - colored tip = |F - T|, the magnitude difference
# So one wedge encodes the rescaled FTNA, the TNA probability, and the
# discrepancy at once.

# Sequential green -> yellow -> orange -> red -> dark-red ramp (Okabe-Ito
# based) keyed to magnitude. `threshold` defaults to 10% of the max.
.magdiff_ramp <- function(mag, threshold = NULL) {
  m <- max(mag, na.rm = TRUE)
  if (is.null(threshold)) threshold <- 0.10 * m
  list(
    cols      = c("#009E73", "#F0E442", "#E69F00", "#D55E00", "#7A1F00"),
    breaks    = c(0, threshold, 0.25 * m, 0.50 * m, m),
    threshold = threshold,
    max_mag   = m
  )
}

.magdiff_stacked <- function(ed, title = NULL, threshold = NULL,
                             r_inner = 0.45, r_outer = 1.00,
                             gap_rad = 0.06, label_offset = 0.04,
                             label_size = 3.4, src_label_size = 4.1) {
  # ed expects: from, to, ftna, tna, value (the metric)
  ed <- ed[ed$ftna > 0 | ed$tna > 0, , drop = FALSE]
  ed <- ed[order(ed$from, ed$to), , drop = FALSE]
  ed$mag <- abs(ed$value)

  from_nodes <- sort(unique(ed$from))
  n_from <- length(from_nodes)
  node_col <- stats::setNames(.okabe_ito[((seq_len(n_from) - 1L)
                                    %% length(.okabe_ito)) + 1L],
                       from_nodes)

  # sector start / size — vectorized via cumsum (no for loop). unname() so
  # the from-state names do not propagate into per-sector data.frames (which
  # would warn "row names were found from a short variable").
  edge_counts <- unname(vapply(from_nodes,
                        function(n) sum(ed$from == n), integer(1)))
  available <- 2 * pi - gap_rad * n_from
  sector_sz <- (edge_counts / sum(edge_counts)) * available
  starts_offset <- c(0, cumsum(sector_sz + gap_rad)[-n_from])
  sector_start  <- pi / 2 - starts_offset
  sector_mid    <- sector_start - sector_sz / 2

  # per-edge angle + slot width — lapply over sectors
  geom_per_sector <- lapply(seq_along(from_nodes), function(i) {
    idx <- which(ed$from == from_nodes[i])
    n_e <- length(idx); s <- sector_start[i]; sz <- sector_sz[i]
    pad <- sz * 0.06
    if (n_e == 1L) {
      data.frame(idx = idx, angle = s - sz / 2, slot_w = sz - 2 * pad)
    } else {
      step <- (sz - 2 * pad) / n_e
      data.frame(idx = idx,
                 angle  = seq(s - pad - step / 2, s - sz + pad + step / 2,
                              length.out = n_e),
                 slot_w = step * 0.78)
    }
  })
  geom_df <- do.call(rbind, geom_per_sector)
  geom_df <- geom_df[order(geom_df$idx), ]
  ed$angle  <- geom_df$angle
  ed$slot_w <- geom_df$slot_w

  # Bar anatomy:
  #   r_inner ─ to_r(min(F,T))         = grey shared portion
  #   to_r(min(F,T)) ─ to_r(max(F,T))  = colored magnitude-difference tip
  scale_max <- max(c(ed$ftna, ed$tna))
  to_r <- function(v) r_inner + (v / scale_max) * (r_outer - r_inner)
  ed$r_grey_top <- to_r(pmin(ed$ftna, ed$tna))
  ed$r_color_hi <- to_r(pmax(ed$ftna, ed$tna))

  ramp <- .magdiff_ramp(ed$mag, threshold)

  ed$id <- seq_len(nrow(ed))
  make_wedge <- function(theta_c, w, r1, r2, id, value) {
    if (r2 <= r1) return(NULL)
    arc <- seq(theta_c - w / 2, theta_c + w / 2, length.out = 16)
    xs <- c(r1 * cos(arc), r2 * cos(rev(arc)))
    ys <- c(r1 * sin(arc), r2 * sin(rev(arc)))
    data.frame(id = id, value = value, x = xs, y = ys)
  }
  grey_polys <- do.call(rbind, lapply(ed$id, function(i)
    make_wedge(ed$angle[i], ed$slot_w[i],
                r_inner, ed$r_grey_top[i], i, NA_real_)))
  col_polys  <- do.call(rbind, lapply(ed$id, function(i)
    make_wedge(ed$angle[i], ed$slot_w[i],
                ed$r_grey_top[i], ed$r_color_hi[i], i, ed$mag[i])))

  theta_seq  <- seq(0, 2 * pi, length.out = 300)
  ring_inner <- data.frame(x = r_inner * cos(theta_seq),
                           y = r_inner * sin(theta_seq))
  ring_outer <- data.frame(x = r_outer * cos(theta_seq),
                           y = r_outer * sin(theta_seq))

  label_r <- r_outer + label_offset
  deg <- ed$angle * 180 / pi
  flip <- cos(ed$angle) < 0
  ed$x_lab <- label_r * cos(ed$angle)
  ed$y_lab <- label_r * sin(ed$angle)
  ed$text_angle <- ifelse(flip, deg + 180, deg)
  ed$hjust      <- ifelse(flip, 1, 0)
  ed$lab_col    <- node_col[ed$from]

  src_r <- r_inner * 0.82
  src_df <- data.frame(node = from_nodes,
                       angle = sector_mid,
                       x_lab = src_r * cos(sector_mid),
                       y_lab = src_r * sin(sector_mid),
                       color = node_col[from_nodes])
  src_df$text_angle <- ifelse(cos(sector_mid) < 0,
                               sector_mid * 180 / pi + 90,
                               sector_mid * 180 / pi - 90)

  ggplot2::ggplot() +
    ggplot2::geom_path(data = ring_inner,
                        ggplot2::aes(x = x, y = y),
                        colour = "grey88", linewidth = 0.3) +
    ggplot2::geom_path(data = ring_outer,
                        ggplot2::aes(x = x, y = y),
                        colour = "grey88", linewidth = 0.3) +
    ggplot2::geom_polygon(data = grey_polys,
                          ggplot2::aes(x = x, y = y, group = id),
                          fill = "grey75", colour = "white", linewidth = 0.15) +
    ggplot2::geom_polygon(data = col_polys,
                          ggplot2::aes(x = x, y = y, group = id, fill = value),
                          colour = "white", linewidth = 0.15) +
    ggplot2::geom_text(data = ed,
                       ggplot2::aes(x = x_lab, y = y_lab, label = to,
                                     angle = text_angle, hjust = hjust,
                                     colour = lab_col),
                       size = label_size, show.legend = FALSE) +
    ggplot2::geom_text(data = src_df,
                       ggplot2::aes(x = x_lab, y = y_lab, label = node,
                                     angle = text_angle, colour = color),
                       size = src_label_size, fontface = "bold",
                       show.legend = FALSE) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_gradientn(
      name = "|magnitude\ndifference|",
      colors = ramp$cols,
      values = scales::rescale(ramp$breaks),
      limits = c(0, ramp$max_mag),
      breaks = ramp$breaks,
      labels = c("0  stable",
                 sprintf("%.2f  small",    ramp$breaks[2]),
                 sprintf("%.2f  moderate", ramp$breaks[3]),
                 sprintf("%.2f  large",    ramp$breaks[4]),
                 sprintf("%.2f  extreme",  ramp$breaks[5])),
      guide  = ggplot2::guide_colorbar(barwidth = 0.7, barheight = 9,
                                        ticks.colour = "grey50")) +
    ggplot2::coord_equal(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4),
                          clip = "off") +
    ggplot2::ggtitle(title) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(legend.position = "right",
                    plot.title = ggplot2::element_text(hjust = 0.5,
                                                         face = "bold"),
                    plot.margin = ggplot2::margin(8, 8, 8, 8))
}

#' @describeIn magnitude_difference Plot the per-edge magnitude difference as
#'   a polar portrait. `type = "stacked"` (default) draws one sector per
#'   from-state with stacked wedges (grey base = shared value, colored tip =
#'   magnitude difference); `type = "circular"` draws a chord-style diagram
#'   with signed differences on a diverging blue-orange scale.
#' @param type Plot style, `"stacked"` (default) or `"circular"`.
#' @param min_show For `type = "circular"`, drop edges whose magnitude is
#'   below this fraction of the maximum.
#' @param title Plot title. `NULL` generates one from the metric and scale.
#' @return `plot` returns a `ggplot` object.
#' @export
plot.magnitude_difference <- function(x, type = c("stacked", "circular"),
                                      min_show = 0.01,
                                      title = NULL, ...) {
  type <- match.arg(type)
  if (is.null(title))
    title <- sprintf("Magnitude difference - %s (%s scaling)",
                      x$metric, x$scale)
  if (type == "stacked") {
    .magdiff_stacked(x$edges, title = title, ...)
  } else {
    .magdiff_circular(x$edges, title = title, min_show = min_show, ...)
  }
}
