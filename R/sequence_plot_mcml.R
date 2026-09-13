# ==============================================================================
# Multichannel sequence plots for MCML objects.
#
# An mcml built from sequences stores, per cluster, the full subject x time
# matrix with non-cluster states masked to NA (`fit$clusters[[k]]$data`). Each
# cluster is therefore a "channel": stacking them gives a TraMineR-style
# multichannel view. Reached via sequence_plot(mcml, type = ...); never called
# directly by users.
# ==============================================================================

utils::globalVariables(c("time", "y", "key", "prop"))

# ---- internal: pull per-channel matrices out of an mcml ---------------------
.mcml_seq_channels <- function(x, trim = NULL, expand = NULL) {
  layers <- x$clusters
  if (is.null(layers) || !length(layers)) {
    stop("mcml has no clusters to plot.", call. = FALSE)
  }
  datas <- lapply(layers, function(z) z$data)
  if (any(vapply(datas, is.null, logical(1L)))) {
    stop("This mcml carries no sequence data (built from a matrix). ",
         "Multichannel sequence plots need an mcml built from sequences.",
         call. = FALSE)
  }
  cluster_names <- names(layers)
  clusters <- lapply(layers, function(z) z$labels)
  tcols    <- colnames(datas[[1L]])
  cmats    <- lapply(datas, function(d) as.matrix(d[, tcols, drop = FALSE]))
  names(cmats) <- cluster_names

  all_states    <- sort(unique(unlist(clusters, use.names = FALSE)))
  state2cluster <- stats::setNames(rep(cluster_names, lengths(clusters)),
                                   unlist(clusters, use.names = FALSE))
  # The summary band IS the macro over time, so it follows the macro's
  # resolution: a state whose cluster was expanded (build_mcml(expand = ))
  # keys as itself, every other state still keys as its cluster. Without this
  # the panel collapses exactly the distinction `expand` was asked to show.
  expanded <- .mcml_resolve_expand(expand, cluster_names)
  if (!is.null(expanded)) {
    keep <- unlist(clusters[names(clusters) %in% expanded], use.names = FALSE)
    state2cluster[keep] <- keep
  }
  # Clusters partition the states, so coalescing recovers the full sequence.
  summary_mat <- Reduce(function(a, b) { a[is.na(a)] <- b[is.na(a)]; a }, cmats)
  summary_cluster_mat <- summary_mat
  summary_cluster_mat[] <- state2cluster[summary_mat]

  # Trim the time axis on the coalesced full sequence (a single masked
  # channel mostly contains NA, so its per-row length understates the real
  # sequence length). Apply the same cut to every channel to keep rows
  # aligned across panels.
  cut <- .trim_cut(summary_mat, trim)
  if (cut < length(tcols)) {
    keep                <- seq_len(cut)
    tcols               <- tcols[keep]
    cmats               <- lapply(cmats, function(m) m[, keep, drop = FALSE])
    summary_mat         <- summary_mat[, keep, drop = FALSE]
    summary_cluster_mat <- summary_cluster_mat[, keep, drop = FALSE]
  }

  # Categories of the summary band: a cluster name, or the member states of a
  # cluster that `expand` opened. Distinct from `cluster_names`, which stays
  # the partition and names the lower channels.
  macro_keys <- unique(unlist(lapply(cluster_names, function(k) {
    if (!is.null(expanded) && k %in% expanded) clusters[[k]] else k
  }), use.names = FALSE))

  list(cluster_names = cluster_names, clusters = clusters, cmats = cmats,
       times = seq_along(tcols), all_states = all_states,
       macro_keys = macro_keys,
       summary_mat = summary_mat, summary_cluster_mat = summary_cluster_mat)
}

# Resolve `expand` against the cluster names: NULL stays NULL, "all"/TRUE
# means every cluster, anything else must name clusters that exist.
.mcml_resolve_expand <- function(expand, cluster_names) {
  if (is.null(expand)) {
    return(NULL)
  }
  if (isTRUE(expand) || identical(expand, "all")) {
    return(cluster_names)
  }
  stopifnot("`expand` must be a character vector, TRUE, \"all\" or NULL" =
              is.character(expand))
  unknown <- setdiff(expand, cluster_names)
  if (length(unknown) > 0L) {
    stop("Unknown cluster(s) in `expand`: ",
         paste(utils::head(unknown, 5L), collapse = ", "),
         ". Available: ", paste(cluster_names, collapse = ", "),
         call. = FALSE)
  }
  expand
}

# ---- internal: state / cluster / faded palettes -----------------------------
.mcml_seq_palettes <- function(ch, state_colors) {
  states <- ch$all_states
  cn     <- ch$macro_keys
  state_pal <- stats::setNames(.state_palette(state_colors, length(states)), states)
  base_cl   <- c("#264653", "#2A9D8F", "#E76F51", "#E9C46A", "#8AB17D",
                 "#5B5F97", "#B5838D", "#6D6875")
  cluster_pal <- stats::setNames(rep_len(base_cl, length(cn)), cn)
  chan_pal <- cluster_pal[ch$cluster_names]
  chan_pal[is.na(names(chan_pal))] <- base_cl[1L]
  faded <- stats::setNames(
    vapply(chan_pal, function(col) {
      m <- 0.22 * (grDevices::col2rgb(col) / 255) + 0.78   # 22% colour, 78% white
      grDevices::rgb(m[1L], m[2L], m[3L])
    }, character(1L)),
    paste0(ch$cluster_names, " (elsewhere)"))
  list(state = state_pal, cluster = cluster_pal, faded = faded)
}

# ---- internal: legend keys grouped by cluster -------------------------------
# Order the fill legend so each cluster name is followed by its own states:
#   Cognitive, discuss, synthesis, ... , Regulation, plan, ... .
# The cluster-name key carries the cluster colour (the same colour used in the
# Summary channel), so it reads as a group header above its member states.
# Used as `breaks` for the fill scale; rendered as a single vertical column.
.mcml_grouped_breaks <- function(ch) {
  unique(unlist(lapply(ch$macro_keys, function(k) c(k, ch$clusters[[k]])),
                use.names = FALSE))
}

# A cluster may legitimately be named after a state it contains (a singleton
# cluster is the natural case). Cluster keys and state keys then collide in
# the single `key` fill scale: the level vector gains a duplicate, which
# `factor()` rejects outright, and the values vector gains a duplicate name.
# Both helpers below drop the later duplicate, so the shared name keeps one
# level and one colour.
.mcml_uniq_levels <- function(...) unique(c(...))

.mcml_uniq_values <- function(...) {
  v <- c(...)
  v[!duplicated(names(v))]
}

# ---- internal: melt a subject x time matrix to long, dropping NA cells ------
.mcml_mat_long <- function(m, channel, times, value) {
  d <- data.frame(
    y       = rep(seq_len(nrow(m)), times = ncol(m)),
    time    = rep(times, each = nrow(m)),
    v       = as.vector(m),
    channel = channel,
    stringsAsFactors = FALSE
  )
  names(d)[names(d) == "v"] <- value
  d[!is.na(d[[value]]), , drop = FALSE]
}

# ---- internal: multichannel index (carpet) ----------------------------------
.mcml_index_plot <- function(ch, pals, main, time_label) {
  cn          <- ch$cluster_names
  chan_levels <- c("Summary", cn)
  times       <- ch$times

  # Order rows by the macro (cluster) sequence so the carpet reads cleanly, and
  # apply the SAME order to every channel so rows stay aligned across panels.
  ord <- do.call(order, c(lapply(times, function(j) ch$summary_cluster_mat[, j]),
                          list(na.last = TRUE)))
  scm <- ch$summary_cluster_mat[ord, , drop = FALSE]
  sm  <- ch$summary_mat[ord, , drop = FALSE]
  cm  <- lapply(ch$cmats, function(m) m[ord, , drop = FALSE])

  summary_tiles <- .mcml_mat_long(scm, "Summary", times, "key")
  ghost_tiles <- do.call(rbind, lapply(cn, function(k) {
    gm <- scm
    gm[!(is.na(cm[[k]]) & !is.na(sm))] <- NA           # keep "active elsewhere"
    d <- .mcml_mat_long(gm, k, times, "key")
    d$key <- paste0(d$key, " (elsewhere)")
    d
  }))
  state_tiles <- do.call(rbind, lapply(cn, function(k)
    .mcml_mat_long(cm[[k]], k, times, "key")))

  tiles <- rbind(summary_tiles, ghost_tiles, state_tiles)
  tiles$channel <- factor(tiles$channel, levels = chan_levels)
  tiles$key     <- factor(tiles$key,
                          levels = .mcml_uniq_levels(ch$all_states, cn,
                                                     paste0(ch$cluster_names, " (elsewhere)")))

  ggplot2::ggplot(tiles, ggplot2::aes(x = time, y = y, fill = key)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~ channel, ncol = 1L, strip.position = "left") +
    ggplot2::scale_fill_manual(
      values = .mcml_uniq_values(pals$state, pals$cluster, pals$faded),
      breaks = .mcml_grouped_breaks(ch), na.value = "white",
      name = "Cluster / State") +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1L, byrow = TRUE)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_reverse(expand = c(0, 0)) +
    ggplot2::labs(x = time_label, y = NULL, title = main) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid        = ggplot2::element_blank(),
      panel.background  = ggplot2::element_rect(fill = "white", colour = "grey80"),
      axis.text.y       = ggplot2::element_blank(),
      axis.ticks.y      = ggplot2::element_blank(),
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold"),
      panel.spacing     = ggplot2::unit(8, "pt"),
      legend.position   = "right")
}

# ---- internal: multichannel distribution (seqdplot) -------------------------
.mcml_dist_plot <- function(ch, pals, na_color, normalize, main, time_label,
                            keep = "all") {
  cn          <- ch$cluster_names
  chan_levels <- c("Summary", cn)
  times       <- ch$times
  N           <- nrow(ch$summary_mat)

  share_at <- function(m, cats) {
    r <- vapply(seq_len(ncol(m)),
                function(j) tabulate(factor(m[, j], levels = cats), nbins = length(cats)),
                integer(length(cats))) / N
    matrix(r, nrow = length(cats), dimnames = list(cats, NULL))   # keep 2-D for k=1
  }
  cl_share <- share_at(ch$summary_cluster_mat, ch$macro_keys)
  band <- function(channel, key, prop)
    data.frame(time = times, channel = channel, key = key, prop = prop,
               stringsAsFactors = FALSE)

  if (isTRUE(normalize)) {
    # seqdplot: each time point sums to 1 within its channel (composition).
    norm_cols <- function(m) { cs <- colSums(m); sweep(m, 2L, ifelse(cs > 0, cs, 1), "/") }
    summ <- norm_cols(cl_share)
    bands <- rbind(
      do.call(rbind, lapply(ch$macro_keys,
                            function(cl) band("Summary", cl, summ[cl, ]))),
      do.call(rbind, lapply(cn, function(k) {
        own <- norm_cols(share_at(ch$cmats[[k]], ch$clusters[[k]]))
        do.call(rbind, lapply(seq_along(ch$clusters[[k]]),
                              function(i) band(k, ch$clusters[[k]][i], own[i, ])))
      })))
    bands$key <- factor(bands$key, levels = .mcml_uniq_levels(ch$all_states, cn))
    fillvals  <- .mcml_uniq_values(pals$state, pals$cluster)
    brks      <- c(cn, ch$all_states)
    y_lab     <- "Composition (sums to 1)"
  } else {
    # prevalence: own states (solid) + other clusters (faded) + NA, to 100%.
    inactive   <- pmax(0, 1 - colSums(cl_share))
    faded_keys <- paste0(cn, " (elsewhere)")
    bands <- rbind(
      do.call(rbind, lapply(cn, function(cl) band("Summary", cl, cl_share[cl, ]))),
      band("Summary", "NA", inactive),
      do.call(rbind, lapply(cn, function(k) {
        own <- share_at(ch$cmats[[k]], ch$clusters[[k]])
        own_b <- do.call(rbind, lapply(seq_along(ch$clusters[[k]]),
                                       function(i) band(k, ch$clusters[[k]][i], own[i, ])))
        oth_b <- do.call(rbind, lapply(setdiff(cn, k),
                                       function(j) band(k, paste0(j, " (elsewhere)"), cl_share[j, ])))
        rbind(own_b, oth_b, band(k, "NA", inactive))
      })))
    bands$key <- factor(bands$key,
                        levels = .mcml_uniq_levels(ch$all_states, cn,
                                                   faded_keys, "NA"))
    fillvals  <- .mcml_uniq_values(pals$state, pals$cluster, pals$faded,
                                   "NA" = na_color)
    brks      <- c(cn, ch$all_states, "NA")
    y_lab     <- "Share of subjects"
  }
  bands$channel <- factor(bands$channel, levels = chan_levels)

  # Restrict to one half when the caller wants the macro drawn on its own.
  if (!identical(keep, "all")) {
    want <- if (identical(keep, "summary")) "Summary" else setdiff(chan_levels, "Summary")
    bands <- bands[bands$channel %in% want, , drop = FALSE]
    bands$channel <- droplevels(bands$channel)
    bands$key     <- droplevels(bands$key)
    fillvals <- fillvals[names(fillvals) %in% levels(bands$key)]
    brks     <- intersect(brks, levels(bands$key))
    legend_name <- if (identical(keep, "summary")) "Cluster" else "State"
    if (is.null(main)) {
      main <- if (identical(keep, "summary")) {
        "Macro: cluster composition"
      } else {
        "Within-cluster: state composition"
      }
    }
  } else {
    legend_name <- "Cluster / State"
  }
  grouped <- if (identical(keep, "all")) {
    c(.mcml_grouped_breaks(ch), if ("NA" %in% brks) "NA")
  } else {
    intersect(c(.mcml_grouped_breaks(ch), "NA"), levels(bands$key))
  }

  ggplot2::ggplot(bands, ggplot2::aes(x = time, y = prop, fill = key)) +
    ggplot2::geom_area(position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::facet_wrap(~ channel, ncol = 1L, strip.position = "left") +
    ggplot2::scale_fill_manual(values = fillvals, breaks = grouped,
                               name = legend_name) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1L, byrow = TRUE)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0),
                                labels = scales::percent_format(accuracy = 1)) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = time_label, y = y_lab, title = main) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor  = ggplot2::element_blank(),
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold"),
      panel.spacing     = ggplot2::unit(8, "pt"),
      legend.position   = "right")
}

# ---- internal: mcml dispatcher (called from sequence_plot) ------------------
.sequence_plot_mcml <- function(x, type, normalize, state_colors, na_color,
                                main, time_label, trim = NULL,
                                panel = "both", expand = NULL) {
  panel <- match.arg(panel, c("both", "summary", "channels"))
  ch   <- .mcml_seq_channels(x, trim, expand)
  pals <- .mcml_seq_palettes(ch, state_colors)
  if (type %in% c("heatmap", "index")) {
    if (!identical(panel, "both")) {
      stop("`panel` applies to type = \"distribution\".", call. = FALSE)
    }
    return(.mcml_index_plot(ch, pals, main, time_label))
  }
  # The macro channel is keyed by cluster and the rest by state, so a shared
  # fill scale can put a cluster and a state on one colour. Drawing a single
  # panel gives that panel its own palette, legend and default title.
  .mcml_dist_plot(ch, pals, na_color, normalize, main, time_label,
                  keep = switch(panel, both = "all", panel))
}
