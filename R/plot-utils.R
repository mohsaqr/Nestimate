# ==============================================================================
# Internal plotting helpers shared by sequence_plot() and distribution_plot().
# Not exported. Keep this file thin: only logic used in 2+ places belongs here.
# ==============================================================================

.okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")


# Encode a wide sequence matrix (states as character/factor/integer) into
# integer codes and return the sorted level set. Returns a list with:
#   codes      - integer matrix, same shape as input, NA preserved
#   levels     - sort.int(unique(non-NA))
.encode_states <- function(x) {
  stopifnot(is.data.frame(x) || is.matrix(x))
  if (is.data.frame(x)) x <- as.matrix(x)
  levels <- sort.int(unique(as.vector(x)), na.last = NA)
  if (length(levels) == 0L) {
    stop("`x` contains no non-NA values.", call. = FALSE)
  }
  codes <- matrix(match(x, levels),
                  nrow = nrow(x), ncol = ncol(x),
                  dimnames = dimnames(x))
  list(codes = codes, levels = levels)
}


# Compute only the dissimilarity matrix for a sequence of states using the
# same dispatch as cluster_data(), but without running PAM / hclust /
# silhouette. Used by sequence_plot() when sort = any distance metric.
.sequence_distance <- function(x, dissimilarity = "lcs") {
  enc <- .encode_sequences(x)
  use_stringdist <- requireNamespace("stringdist", quietly = TRUE) &&
    enc$n_states < 62L &&
    !(dissimilarity == "hamming")
  if (use_stringdist) {
    .dissimilarity_matrix_stringdist(enc, dissimilarity, lambda = 0)
  } else {
    .dissimilarity_matrix_r(enc, dissimilarity, lambda = 0)
  }
}


# Translate the legend knobs into outer-margin sizes (in text-lines). Uses
# the device's actual char-width / line-height ratio so long labels or
# titles don't clip. Returns a 2-vector c(oma_bottom, oma_right); all
# other edges get 0.3 elsewhere.
.legend_oma_size <- function(labels, position, size, ncol, title) {
  label_chars <- max(nchar(as.character(labels)))
  title_chars <- if (!is.null(title)) nchar(title) else 0
  wide_chars  <- max(label_chars, title_chars)
  cw_per_line <- graphics::par("cin")[1] / graphics::par("csi")
  rows <- if (!is.null(ncol) && ncol > 0) {
    ceiling(length(labels) / ncol)
  } else if (position == "bottom") 1L else length(labels)
  cols <- if (!is.null(ncol) && ncol > 0) {
    ncol
  } else if (position == "bottom") length(labels) else 1L
  # text line factor 1.6 (was 1.25): reserves enough height that cex = 1.2
  # descenders + padding fit inside the oma strip without clipping at the
  # device bottom edge.
  h <- rows * size * 1.6 + (if (!is.null(title)) 1.2 else 0.6)
  w <- cols * ((wide_chars + 3) * cw_per_line * size * 1.15) + 2.5
  c(oma_b = if (position == "bottom") h else 0.3,
    oma_r = if (position == "right")  w else 0.3)
}


# Draw a legend anchored to the inner-figure edge, growing outward into
# the reserved oma strip. Called after all plot panels are drawn.
.draw_legend_in_oma <- function(labels, colors, position, size, ncol,
                                title, border, bty) {
  args <- list(legend = as.character(labels),
               fill   = colors,
               border = border,
               bty    = bty,
               cex    = size,
               pt.cex = size * 1.3,
               title  = title,
               xpd    = NA)
  if (is.null(ncol)) {
    args$horiz <- position == "bottom"
  } else {
    args$ncol <- ncol
  }
  omd <- graphics::par("omd")
  gap <- 0.005
  anchor <- switch(position,
                   bottom = list(x = 0.5,            y = omd[3] - gap,
                                 xj = 0.5,           yj = 1),
                   right  = list(x = omd[2] + gap,   y = 0.5,
                                 xj = 0,             yj = 0.5))
  x <- graphics::grconvertX(anchor$x, from = "ndc", to = "user")
  y <- graphics::grconvertY(anchor$y, from = "ndc", to = "user")
  do.call(graphics::legend,
          c(list(x = x, y = y, xjust = anchor$xj, yjust = anchor$yj), args))
}


# Every-Nth tick positions. When `tick` is NULL, thin automatically so at
# most ~15 labels are drawn.
.tick_positions <- function(n_cols, tick, colnames_ = NULL) {
  by <- if (!is.null(tick)) {
    stopifnot(is.numeric(tick), length(tick) == 1L, tick >= 1)
    as.integer(tick)
  } else {
    max(1L, as.integer(ceiling(n_cols / 15)))
  }
  at  <- seq(1L, n_cols, by = by)
  lab <- if (!is.null(colnames_)) colnames_[at] else at
  list(at = at, labels = lab)
}


# Adaptive legend text size. NULL = auto-scale from device width so the
# legend visually matches the plot: ~0.7 on a 5-inch figure, ~1.0 on a
# 10-inch figure, clamped to [0.65, 1.2] to stay readable.
.auto_legend_size <- function(user_size) {
  if (!is.null(user_size)) return(user_size)
  din_w <- graphics::par("din")[1]
  max(0.65, min(1.2, din_w / 10))
}


# Build the state palette from an optional user-supplied vector, falling
# back to recycled Okabe-Ito. Returns a character vector length = n_states.
.state_palette <- function(user_colors, n_states) {
  if (is.null(user_colors)) return(rep_len(.okabe_ito, n_states))
  stopifnot(is.character(user_colors), length(user_colors) >= n_states)
  user_colors[seq_len(n_states)]
}
