# ---- Sequence Pattern Comparison ----

#' Compare Subsequence Patterns Between Groups
#'
#' @description
#' Extracts all k-gram patterns (subsequences of length k) from sequences
#' in each group, compares their frequencies, and optionally runs a
#' permutation test to identify statistically significant differences.
#'
#' @param x A \code{netobject_group} (from grouped \code{build_network}),
#'   a \code{netobject} (requires \code{group}), or a wide-format
#'   \code{data.frame} (requires \code{group}).
#' @param group Character or vector. Column name or vector of group labels.
#'   Not needed for \code{netobject_group}.
#' @param sub Integer vector. Pattern lengths to analyze. Default: \code{2:4}.
#' @param min_freq Integer. Minimum frequency in each group for a pattern
#'   to be included. Default: 5.
#' @param test Logical. Run permutation test? Default: \code{FALSE}.
#' @param iter Integer. Permutation iterations. Default: 1000.
#' @param adjust Character. P-value correction method (see
#'   \code{\link[stats]{p.adjust}}). Default: \code{"BH"}.
#'
#' @return An object of class \code{"net_sequence_comparison"} containing:
#' \describe{
#'   \item{patterns}{Tidy data.frame: pattern, length, group frequencies,
#'     proportions, and optionally effect_size and p_value.}
#'   \item{groups}{Character vector of group names.}
#'   \item{n_patterns}{Integer. Number of patterns passing min_freq.}
#'   \item{params}{List of sub, min_freq, test, iter, adjust.}
#' }
#'
#' @examples
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:4], 60, TRUE),
#'   V2 = sample(LETTERS[1:4], 60, TRUE),
#'   V3 = sample(LETTERS[1:4], 60, TRUE),
#'   V4 = sample(LETTERS[1:4], 60, TRUE)
#' )
#' grp <- rep(c("X", "Y"), 30)
#' net <- build_network(seqs, method = "relative")
#' res <- compare_sequences(net, group = grp, sub = 2:3)
#'
#' @importFrom data.table data.table rbindlist setDT setorderv
#' @importFrom stats p.adjust sd
#' @export
compare_sequences <- function(x, group = NULL, sub = 2:4,
                              min_freq = 5L, test = FALSE,
                              iter = 1000L, adjust = "BH") {

  min_freq <- as.integer(min_freq)
  iter <- as.integer(iter)
  sub <- as.integer(sub)
  stopifnot(
    length(sub) >= 1L, all(sub >= 1L),
    min_freq >= 1L,
    iter >= 1L
  )
  adjust <- match.arg(adjust, stats::p.adjust.methods)

  # ---- Extract sequences and group labels ----
  parsed <- .cs_parse_input(x, group)
  seq_mat <- parsed$seq_mat   # character matrix: rows=sequences, cols=timesteps
  group_vec <- parsed$group   # character vector, length = nrow(seq_mat)
  groups <- sort(unique(group_vec))
  n_groups <- length(groups)
  stopifnot("Need at least 2 groups" = n_groups >= 2L)

  n <- nrow(seq_mat)
  p <- ncol(seq_mat)
  sub <- sub[sub <= p]
  stopifnot("No valid pattern lengths (sub > sequence length)" = length(sub) >= 1L)

  # ---- Extract all k-gram patterns ----
  pat_list <- .cs_extract_patterns(seq_mat, sub)

  # ---- Count per group ----
  group_idx <- lapply(groups, function(g) which(group_vec == g))
  names(group_idx) <- groups

  results <- rbindlist(lapply(seq_along(sub), function(si) {
    k <- sub[si]
    pat_mat <- pat_list[[si]]  # rows = positions, cols = sequences
    unique_pats <- sort(unique(pat_mat[nzchar(pat_mat)]))
    if (length(unique_pats) == 0L) return(NULL)

    # Count per group via tabulate
    counts <- vapply(groups, function(g) {
      vals <- pat_mat[, group_idx[[g]], drop = FALSE]
      idx <- match(vals, unique_pats, nomatch = 0L)
      tabulate(idx, nbins = length(unique_pats))
    }, integer(length(unique_pats)))

    if (!is.matrix(counts)) counts <- matrix(counts, ncol = n_groups)
    colnames(counts) <- groups

    # Proportions (computed on ALL patterns, before min_freq filter)
    col_totals <- colSums(counts)
    props <- sweep(counts, 2, col_totals, `/`)
    props[!is.finite(props)] <- NA_real_

    # Filter by min_freq in ALL groups
    row_min <- apply(counts, 1, min)
    keep <- row_min >= min_freq
    if (!any(keep)) return(NULL)
    counts <- counts[keep, , drop = FALSE]
    props <- props[keep, , drop = FALSE]
    unique_pats <- unique_pats[keep]

    # Build tidy rows
    freq_df <- as.data.frame(counts, stringsAsFactors = FALSE)
    names(freq_df) <- paste0("freq_", groups)
    prop_df <- as.data.frame(props, stringsAsFactors = FALSE)
    names(prop_df) <- paste0("prop_", groups)

    out <- data.table(pattern = unique_pats, length = k)
    out <- cbind(out, freq_df, prop_df)

    # ---- Permutation test ----
    if (test) {
      stat_obs <- .cs_pattern_statistic(counts)
      stat_perm <- matrix(0, nrow = nrow(counts), ncol = iter)
      stat_ge <- numeric(nrow(counts))

      for (i in seq_len(iter)) {
        perm_group <- group_vec[sample(n)]
        perm_idx <- lapply(groups, function(g) which(perm_group == g))
        perm_counts <- vapply(seq_along(groups), function(gi) {
          vals <- pat_mat[, perm_idx[[gi]], drop = FALSE]
          idx <- match(vals, unique_pats, nomatch = 0L)
          tabulate(idx, nbins = length(unique_pats))
        }, integer(length(unique_pats)))
        if (!is.matrix(perm_counts)) {
          perm_counts <- matrix(perm_counts, ncol = n_groups)
        }
        stat_perm[, i] <- .cs_pattern_statistic(perm_counts)
        stat_ge <- stat_ge + (stat_perm[, i] >= stat_obs)
      }

      perm_mean <- rowMeans(stat_perm)
      perm_sd <- apply(stat_perm, 1, sd)
      perm_sd[perm_sd == 0] <- 1  # avoid division by zero
      out$effect_size <- (stat_obs - perm_mean) / perm_sd
      out$p_value <- p.adjust((stat_ge + 1) / (iter + 1), method = adjust)
    }
    out
  }))

  if (is.null(results) || nrow(results) == 0L) {
    stop(sprintf("No patterns with min_freq >= %d in all groups.", min_freq),
         call. = FALSE)
  }

  if (test) {
    results <- results[order(results$p_value, -abs(results$effect_size)), ]
  } else {
    results <- results[order(results$length, results$pattern), ]
  }
  rownames(results) <- NULL

  structure(list(
    patterns   = as.data.frame(results, stringsAsFactors = FALSE),
    groups     = groups,
    n_patterns = nrow(results),
    params     = list(sub = sub, min_freq = min_freq,
                      test = test, iter = iter, adjust = adjust)
  ), class = "net_sequence_comparison")
}


# ---- Input parsing ----

#' @noRd
.cs_parse_input <- function(x, group) {
  # netobject_group: each element has $data
  if (inherits(x, "netobject_group")) {
    group_names <- names(x)
    mats <- lapply(seq_along(x), function(i) {
      d <- x[[i]]$data
      if (is.null(d)) stop("netobject_group element has no $data.", call. = FALSE)
      m <- as.matrix(as.data.frame(d, stringsAsFactors = FALSE))
      # Decode integer-encoded data
      lbl <- rownames(x[[i]]$weights)
      if (!is.null(lbl) && length(lbl) > 0L &&
          (is.integer(m[, 1]) || is.numeric(m[, 1]))) {
        storage.mode(m) <- "character"
        m[] <- lbl[as.integer(m)]
      }
      m
    })
    # Align column count
    max_cols <- max(vapply(mats, ncol, integer(1)))
    mats <- lapply(mats, function(m) {
      if (ncol(m) < max_cols) {
        pad <- matrix(NA_character_, nrow(m), max_cols - ncol(m))
        cbind(m, pad)
      } else m
    })
    seq_mat <- do.call(rbind, mats)
    group_vec <- rep(group_names, vapply(mats, nrow, integer(1)))
    return(list(seq_mat = seq_mat, group = group_vec))
  }

  # netobject or tna: extract $data
  if (inherits(x, "tna") || inherits(x, "netobject") ||
      inherits(x, "cograph_network")) {
    if (inherits(x, "cograph_network") && !inherits(x, "netobject")) {
      x <- .as_netobject(x)
    }
    if (inherits(x, "tna")) {
      d <- x$data
      lbl <- x$labels
    } else {
      d <- x$data
      lbl <- rownames(x$weights)
    }
    if (is.null(d)) stop("Object has no $data.", call. = FALSE)
    seq_mat <- as.matrix(as.data.frame(d, stringsAsFactors = FALSE))
    if (!is.null(lbl) && length(lbl) > 0L &&
        (is.integer(seq_mat[, 1]) || is.numeric(seq_mat[, 1]))) {
      storage.mode(seq_mat) <- "character"
      seq_mat[] <- lbl[as.integer(seq_mat)]
    }
  } else if (is.data.frame(x) || is.matrix(x)) {
    seq_mat <- as.matrix(as.data.frame(x, stringsAsFactors = FALSE))
  } else {
    stop("x must be a netobject_group, netobject, tna, data.frame, or matrix.",
         call. = FALSE)
  }

  # Resolve group
  if (is.null(group)) {
    stop("'group' is required when x is not a netobject_group.", call. = FALSE)
  }
  if (length(group) == 1L && is.character(group) && group %in% colnames(seq_mat)) {
    group_vec <- seq_mat[, group]
    seq_mat <- seq_mat[, colnames(seq_mat) != group, drop = FALSE]
  } else {
    stopifnot("'group' length must match nrow(x)" = length(group) == nrow(seq_mat))
    group_vec <- as.character(group)
  }

  list(seq_mat = seq_mat, group = group_vec)
}


# ---- Pattern extraction (vectorized) ----

#' @noRd
.cs_extract_patterns <- function(m, sub) {
  n <- nrow(m)
  p <- ncol(m)
  mis <- is.na(m)

  lapply(sub, function(k) {
    n_pos <- p - k + 1L
    if (n_pos < 1L) {
      return(matrix("", nrow = 0, ncol = n))
    }
    # Build pattern strings for each position × sequence
    out <- matrix("", nrow = n_pos, ncol = n)
    for (i in seq_len(n_pos)) {
      idx <- i:(i + k - 1L)
      has_na <- .rowSums(mis[, idx, drop = FALSE], m = n, n = k) > 0L
      cols <- m[, idx, drop = FALSE]
      pat <- do.call(paste, c(as.data.frame(cols, stringsAsFactors = FALSE),
                              list(sep = "->")))
      pat[has_na] <- ""
      out[i, ] <- pat
    }
    out
  })
}


# ---- Test statistic ----

#' @noRd
.cs_pattern_statistic <- function(x) {
  n <- sum(x)
  nr <- nrow(x)
  nc <- ncol(x)
  rs <- .rowSums(x, m = nr, n = nc)
  cs <- .colSums(x, m = nr, n = nc)
  expected <- outer(rs, cs) / n
  sqrt(.rowSums((x - expected)^2, m = nr, n = nc))
}


# ---- S3 Methods ----

#' Print Method for net_sequence_comparison
#'
#' @param x A \code{net_sequence_comparison} object.
#' @param ... Additional arguments (ignored).
#' @return The input object, invisibly.
#' @export
print.net_sequence_comparison <- function(x, ...) {
  cat(sprintf("Sequence Comparison  [%d patterns | %d groups: %s]\n",
              x$n_patterns, length(x$groups), paste(x$groups, collapse = ", ")))
  cat(sprintf("  Lengths: %s  |  min_freq: %d",
              paste(x$params$sub, collapse = ", "), x$params$min_freq))
  if (x$params$test) cat(sprintf("  |  permutation: %d iter (%s)", x$params$iter, x$params$adjust))
  cat("\n")

  if (x$n_patterns > 0) {
    top <- utils::head(x$patterns, 10)
    cat("\n")
    print(top, row.names = FALSE)
    if (x$n_patterns > 10) {
      cat(sprintf("  ... and %d more patterns\n", x$n_patterns - 10))
    }
  }
  invisible(x)
}


#' Summary Method for net_sequence_comparison
#'
#' @param object A \code{net_sequence_comparison} object.
#' @param ... Additional arguments (ignored).
#' @return The patterns data.frame, invisibly.
#' @export
summary.net_sequence_comparison <- function(object, ...) {
  print(object$patterns, row.names = FALSE)
  invisible(object$patterns)
}


#' Plot Method for net_sequence_comparison
#'
#' @description
#' Dumbbell plot comparing pattern proportions across groups. Each pattern
#' is a row with dots for each group connected by a segment. When a
#' permutation test was run, significant patterns are highlighted.
#'
#' @param x A \code{net_sequence_comparison} object.
#' @param top Integer. Show top N patterns. Default: 20.
#' @param alpha Numeric. Significance threshold for highlighting.
#'   Default: 0.05.
#' @param ... Additional arguments (ignored).
#' @return A \code{ggplot} object, invisibly.
#'
#' @import ggplot2
#' @export
plot.net_sequence_comparison <- function(x, top = 20L, alpha = 0.05, ...) {
  pat <- x$patterns
  groups <- x$groups
  if (nrow(pat) == 0) {
    message("No patterns to plot.")
    return(invisible(NULL))
  }

  prop_cols <- paste0("prop_", groups)

  # Rank by effect size (if tested) or max proportion difference
  if ("effect_size" %in% names(pat)) {
    pat <- pat[order(-abs(pat$effect_size)), ]
  } else {
    pat$.diff <- apply(pat[, prop_cols, drop = FALSE], 1, function(r) {
      diff(range(r, na.rm = TRUE))
    })
    pat <- pat[order(-pat$.diff), ]
  }
  pat <- utils::head(pat, top)

  # Significance flag
  if ("p_value" %in% names(pat)) {
    pat$significant <- pat$p_value < alpha
  } else {
    pat$significant <- FALSE
  }

  # For the dumbbell segment: need min and max proportion per pattern
  pat$prop_min <- apply(pat[, prop_cols, drop = FALSE], 1, min, na.rm = TRUE)
  pat$prop_max <- apply(pat[, prop_cols, drop = FALSE], 1, max, na.rm = TRUE)

  # Order patterns by the difference (bottom = smallest gap, top = biggest)
  pat$pattern <- factor(pat$pattern,
                         levels = rev(pat$pattern))

  # Segment data
  seg_df <- data.frame(
    pattern  = pat$pattern,
    xmin     = pat$prop_min,
    xmax     = pat$prop_max,
    sig      = pat$significant,
    stringsAsFactors = FALSE
  )

  # Point data (long form)
  long_rows <- lapply(seq_along(groups), function(gi) {
    data.frame(
      pattern    = pat$pattern,
      group      = groups[gi],
      proportion = pat[[prop_cols[gi]]],
      stringsAsFactors = FALSE
    )
  })
  pts <- do.call(rbind, long_rows)

  # Colors
  n_groups <- length(groups)
  if (n_groups == 2L) {
    grp_colors <- c("#2171B5", "#CB181D")
  } else {
    grp_colors <- scales::hue_pal()(n_groups)
  }

  # Back-to-back pyramid with Pearson residual shading
  g1 <- groups[1]; g2 <- groups[2]
  freq1 <- paste0("freq_", g1); freq2 <- paste0("freq_", g2)
  prop1 <- paste0("prop_", g1); prop2 <- paste0("prop_", g2)

  # Pearson residuals: (observed - expected) / sqrt(expected)
  # Expected under independence: row_total * col_total / grand_total
  f1 <- pat[[freq1]]; f2 <- pat[[freq2]]
  row_total <- f1 + f2
  grand_total <- sum(row_total)
  col1_total <- sum(f1); col2_total <- sum(f2)
  exp1 <- row_total * col1_total / grand_total
  exp2 <- row_total * col2_total / grand_total
  resid1 <- (f1 - exp1) / sqrt(pmax(exp1, 1e-10))
  resid2 <- (f2 - exp2) / sqrt(pmax(exp2, 1e-10))

  # Order by total proportion (largest at top)
  total_prop <- pat[[prop1]] + pat[[prop2]]
  ord <- order(total_prop)
  pat <- pat[ord, ]
  resid1 <- resid1[ord]; resid2 <- resid2[ord]

  pyramid <- data.frame(
    pattern = factor(pat$pattern, levels = pat$pattern),
    left    = -pat[[prop1]],
    right   = pat[[prop2]],
    freq_left   = pat[[freq1]],
    freq_right  = pat[[freq2]],
    resid_left  = resid1,
    resid_right = resid2,
    stringsAsFactors = FALSE
  )

  r_max <- max(abs(c(resid1, resid2)), na.rm = TRUE)

  x_max <- max(abs(c(pyramid$left, pyramid$right)), na.rm = TRUE)
  gutter <- x_max * 0.06
  has_pval <- "p_value" %in% names(pat)
  bar_hw <- 0.42

  y_num <- as.numeric(pyramid$pattern)

  # Left bars: xmin..(-gutter), Right bars: (+gutter)..xmax
  rect_df <- rbind(
    data.frame(yidx = y_num, ymin = y_num - bar_hw, ymax = y_num + bar_hw,
               xmin = pyramid$left - gutter, xmax = -gutter,
               freq = pyramid$freq_left, resid = pyramid$resid_left,
               side = "left", stringsAsFactors = FALSE),
    data.frame(yidx = y_num, ymin = y_num - bar_hw, ymax = y_num + bar_hw,
               xmin = gutter, xmax = pyramid$right + gutter,
               freq = pyramid$freq_right, resid = pyramid$resid_right,
               side = "right", stringsAsFactors = FALSE)
  )

  # P-value labels in the gutter
  if (has_pval) {
    p_df <- data.frame(
      yidx = y_num,
      p_label = ifelse(pat$p_value < 0.001, "<.001",
                ifelse(pat$p_value < 0.01,
                       sprintf("%.3f", pat$p_value),
                       sprintf("%.2f", pat$p_value))),
      stringsAsFactors = FALSE
    )
  }

  p <- ggplot() +
    geom_rect(data = rect_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                  fill = resid)) +
    # Freq labels outside each bar
    geom_text(data = rect_df,
              aes(x = ifelse(side == "left", xmin, xmax),
                  y = yidx, label = freq),
              hjust = ifelse(rect_df$side == "left", 1.15, -0.15),
              size = 3, fontface = "bold",
              color = ifelse(abs(rect_df$resid) > r_max * 0.5, "grey25", "grey45")) +
    # P-values in center gutter
    {if (has_pval)
      geom_text(data = p_df,
                aes(x = 0, y = yidx, label = p_label),
                size = 2.5, color = "grey40", fontface = "italic")
    } +
    scale_fill_gradient2(low = "#2171B5", mid = "grey92", high = "#CB181D",
                         midpoint = 0, limits = c(-r_max, r_max),
                         name = "Pearson\nresidual") +
    scale_y_continuous(breaks = y_num,
                       labels = as.character(pat$pattern)) +
    scale_x_continuous(
      labels = function(x) format(abs(x), scientific = FALSE)
    ) +
    labs(x = "Proportion", y = NULL,
         title = "Sequence Pattern Comparison",
         subtitle = paste0(g1, "  \u2190\u2190    \u2192\u2192  ", g2)) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.y          = element_text(family = "mono", size = 8),
      axis.ticks.y         = element_blank(),
      panel.grid.major.y   = element_blank(),
      panel.grid.minor     = element_blank(),
      plot.title           = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle        = element_text(hjust = 0.5, color = "grey40",
                                           size = 11),
      legend.position      = "right",
      plot.margin          = margin(5, 10, 5, 5)
    )

  print(p)
  invisible(p)
}
