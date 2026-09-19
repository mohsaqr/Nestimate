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
.mcml_seq_channels <- function(x, trim = NULL, expand = NULL, combine = NULL) {
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

  # `combine` merges clusters into one channel before anything else, so a
  # merged group behaves as a single cluster everywhere below (its own panel,
  # one Summary key, one faded other-states band) and `expand` can name it.
  groups <- .mcml_resolve_combine(combine, cluster_names)
  if (length(groups) > 0L) {
    owner <- stats::setNames(cluster_names, cluster_names)
    owner[unlist(groups, use.names = FALSE)] <- rep(names(groups), lengths(groups))
    new_names <- unique(unname(owner))
    clusters <- stats::setNames(lapply(new_names, function(g) {
      unlist(clusters[names(owner)[owner == g]], use.names = FALSE)
    }), new_names)
    cmats <- stats::setNames(lapply(new_names, function(g) {
      Reduce(function(a, b) { a[is.na(a)] <- b[is.na(a)]; a },
             cmats[names(owner)[owner == g]])
    }), new_names)
    cluster_names <- new_names
  }

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
       macro_keys = macro_keys, expanded = expanded,
       summary_mat = summary_mat, summary_cluster_mat = summary_cluster_mat)
}

# Legend keys for time spent in other clusters: one per cluster
# ("Social (Other states)") or a single pooled key ("Other states"). The text
# is the user's `rest_label`, carried on the channel list.
.mcml_other_key <- function(ch, cluster) {
  paste0(cluster, " (", ch$rest_label, ")", recycle0 = TRUE)
}

# Resolve `combine` against the cluster names: NULL stays NULL; a character
# vector is one group, a list is several. Every group names >= 2 existing
# clusters, no cluster sits in two groups, and a group is labelled by its list
# name or, failing that, "A + B".
.mcml_resolve_combine <- function(combine, cluster_names) {
  if (is.null(combine)) {
    return(NULL)
  }
  groups <- if (is.list(combine)) combine else list(combine)
  stopifnot(
    "`combine` must be a character vector or a list of character vectors" =
      all(vapply(groups, is.character, logical(1L))),
    "each `combine` group must name at least two clusters" =
      all(lengths(groups) >= 2L))
  members <- unlist(groups, use.names = FALSE)
  unknown <- setdiff(members, cluster_names)
  if (length(unknown) > 0L) {
    stop("Unknown cluster(s) in `combine`: ",
         paste(utils::head(unknown, 5L), collapse = ", "),
         ". Available: ", paste(cluster_names, collapse = ", "),
         call. = FALSE)
  }
  if (anyDuplicated(members) > 0L) {
    stop("A cluster appears in more than one `combine` group: ",
         paste(unique(members[duplicated(members)]), collapse = ", "),
         call. = FALSE)
  }
  labels <- names(groups) %||% character(length(groups))
  missing <- !nzchar(labels)
  labels[missing] <- vapply(groups[missing], paste, character(1L),
                            collapse = " + ")
  clash <- intersect(labels, setdiff(cluster_names, members))
  if (length(clash) > 0L) {
    stop("`combine` label clashes with an existing cluster: ",
         paste(clash, collapse = ", "), call. = FALSE)
  }
  stats::setNames(groups, labels)
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
  # A cluster opened by `expand` has no Summary colour; give each such
  # cluster its own unused base colour so its faded band stays distinct.
  chan_pal <- stats::setNames(cluster_pal[ch$cluster_names], ch$cluster_names)
  open_cl  <- is.na(chan_pal)
  spare    <- setdiff(base_cl, chan_pal)
  chan_pal[open_cl] <- rep_len(if (length(spare)) spare else base_cl, sum(open_cl))
  faded <- stats::setNames(
    vapply(chan_pal, function(col) {
      m <- 0.22 * (grDevices::col2rgb(col) / 255) + 0.78   # 22% colour, 78% white
      grDevices::rgb(m[1L], m[2L], m[3L])
    }, character(1L)),
    .mcml_other_key(ch, ch$cluster_names))
  list(state = state_pal, cluster = cluster_pal, faded = faded,
       rest = stats::setNames("grey78", ch$rest_label))
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
.mcml_index_plot <- function(ch, pals, main, time_label, rest = "clusters") {
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
  # Cells where the subject is active in another cluster. Keyed by the
  # partition (not the Summary keys, which are states once `expand` is set):
  # per cluster, pooled into one `rest_label` key, or left blank.
  state2cl <- stats::setNames(rep(cn, lengths(ch$clusters[cn])),
                              unlist(ch$clusters[cn], use.names = FALSE))
  ghost_tiles <- if (identical(rest, "none")) NULL else do.call(rbind, lapply(cn, function(k) {
    gm <- sm
    gm[] <- state2cl[sm]
    gm[!(is.na(cm[[k]]) & !is.na(sm))] <- NA           # keep "active in another cluster"
    d <- .mcml_mat_long(gm, k, times, "key")
    d$key <- if (identical(rest, "pooled")) rep(ch$rest_label, nrow(d)) else .mcml_other_key(ch, d$key)
    d
  }))
  state_tiles <- do.call(rbind, lapply(cn, function(k)
    .mcml_mat_long(cm[[k]], k, times, "key")))

  tiles <- rbind(summary_tiles, ghost_tiles, state_tiles)
  tiles$channel <- factor(tiles$channel, levels = chan_levels)
  tiles$key     <- factor(tiles$key,
                          levels = .mcml_uniq_levels(ch$all_states, ch$macro_keys, cn,
                                                     .mcml_other_key(ch, ch$cluster_names),
                                                     ch$rest_label))

  .mcml_stack_channels(
    tiles, ch,
    values = .mcml_uniq_values(pals$state, pals$cluster, pals$faded, pals$rest),
    layer  = function(d) ggplot2::ggplot(d, ggplot2::aes(x = time, y = y, fill = key)) +
      ggplot2::geom_tile() +
      ggplot2::scale_y_reverse(expand = c(0, 0)),
    y_lab = NULL, main = main, time_label = time_label,
    theme_extra = ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = "grey80"),
      axis.text.y      = ggplot2::element_blank(),
      axis.ticks.y     = ggplot2::element_blank()))
}

# ---- internal: one ggplot per channel, each with its own legend -------------
# A faceted ggplot has a single fill scale, so every cluster's states pile into
# one long legend. Each channel is drawn as its own ggplot instead: the Summary
# legend lists the clusters, a cluster channel lists its own states followed by
# the faded other-cluster bands. Legend order: macro keys, then each cluster's
# states, then the other-cluster keys, then NA.
.mcml_channel_plots <- function(d, ch, values, layer, y_lab, main, time_label,
                                theme_extra = NULL) {
  d$channel <- droplevels(d$channel)
  chans     <- levels(d$channel)
  n         <- length(chans)
  key_order <- c(ch$macro_keys, ch$all_states,
                 .mcml_other_key(ch, ch$cluster_names), ch$rest_label, "NA")
  plots <- lapply(seq_len(n), function(i) {
    di      <- d[d$channel == chans[i], , drop = FALSE]
    di$key  <- droplevels(di$key)
    own     <- if (identical(chans[i], "Summary")) character(0) else ch$clusters[[chans[i]]]
    present <- levels(di$key)
    brks    <- unique(c(intersect(own, present), intersect(key_order, present),
                        setdiff(present, key_order)))
    first <- i == 1L
    last  <- i == n
    layer(di) +
      ggplot2::facet_wrap(~ channel, ncol = 1L, strip.position = "left") +
      ggplot2::scale_fill_manual(values = values[names(values) %in% brks],
                                 breaks = brks, na.value = "white",
                                 name = if (identical(chans[i], "Summary")) "Cluster" else "State") +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1L, byrow = TRUE)) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::labs(x = if (last) time_label else NULL, y = y_lab,
                    title = if (first) main else NULL) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        strip.text.y.left    = ggplot2::element_text(angle = 0, face = "bold"),
        legend.position      = "right",
        legend.justification = "left") +
      theme_extra +
      if (last) NULL else ggplot2::theme(axis.text.x  = ggplot2::element_blank(),
                                         axis.ticks.x = ggplot2::element_blank())
  })
  names(plots) <- chans
  plots
}

# Stack the per-channel ggplots into one figure. A single channel stays a plain
# ggplot. Several are bound as gtables with `size = "max"`, which aligns every
# panel column (strip and legend widths differ between channels) and splits the
# flexible panel height evenly, so a title or axis never shrinks one panel.
.mcml_stack_channels <- function(d, ch, values, layer, y_lab, main, time_label,
                                 theme_extra = NULL) {
  plots <- .mcml_channel_plots(d, ch, values, layer, y_lab, main, time_label,
                               theme_extra)
  if (length(plots) == 1L) {
    return(plots[[1L]])
  }
  stacked <- do.call(rbind, c(lapply(plots, ggplot2::ggplotGrob),
                              list(size = "max")))
  class(stacked) <- c("mcml_sequence_plot", class(stacked))
  attr(stacked, "panels") <- plots   # the per-channel ggplots, for inspection
  stacked
}

#' Draw a stacked multichannel mcml sequence plot
#'
#' Print method for the figure \code{\link{sequence_plot}} returns for an
#' \code{mcml} with more than one channel: one panel per channel (the macro
#' \code{Summary} and one per cluster), each with its own legend.
#'
#' @param x An \code{mcml_sequence_plot} (a \code{gtable}).
#' @param ... Ignored.
#' @return \code{x}, invisibly. Called for the side effect of drawing it on a
#'   new page of the current graphics device.
#' @export
print.mcml_sequence_plot <- function(x, ...) {
  grid::grid.newpage()
  grid::grid.draw(x)
  invisible(x)
}

# ---- internal: multichannel distribution (seqdplot) -------------------------
.mcml_dist_plot <- function(ch, pals, na_color, normalize, main, time_label,
                            keep = "all", rest = "clusters", na = TRUE) {
  cn          <- ch$cluster_names
  chan_levels <- c("Summary", cn)
  times       <- ch$times
  # Denominator per time point: every subject (na = TRUE, ended subjects form
  # the NA band), or only the subjects still running (na = FALSE, no NA band;
  # each time point stacks to 100% of the running sequences).
  N <- if (isTRUE(na)) {
    rep(nrow(ch$summary_mat), ncol(ch$summary_mat))
  } else {
    pmax(colSums(!is.na(ch$summary_mat)), 1)
  }

  share_at <- function(m, cats) {
    r <- vapply(seq_len(ncol(m)),
                function(j) tabulate(factor(m[, j], levels = cats), nbins = length(cats)),
                integer(length(cats)))
    r <- sweep(matrix(r, nrow = length(cats)), 2L, N, "/")
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
    y_lab     <- "Composition (sums to 1)"
  } else {
    # prevalence: own states (solid) + other clusters (faded) + NA, to 100%.
    inactive   <- pmax(0, 1 - colSums(cl_share))
    na_band    <- function(channel) if (isTRUE(na)) band(channel, "NA", inactive)
    faded_keys <- .mcml_other_key(ch, cn)
    # An expanded cluster has no row of its own in `cl_share` (its states are
    # the macro keys), so its faded share in other panels sums its members.
    cluster_share <- function(j) {
      rows <- if (j %in% ch$expanded) ch$clusters[[j]] else j
      colSums(cl_share[rows, , drop = FALSE])
    }
    bands <- rbind(
      do.call(rbind, lapply(ch$macro_keys,
                            function(cl) band("Summary", cl, cl_share[cl, ]))),
      na_band("Summary"),
      do.call(rbind, lapply(cn, function(k) {
        own <- share_at(ch$cmats[[k]], ch$clusters[[k]])
        own_b <- do.call(rbind, lapply(seq_along(ch$clusters[[k]]),
                                       function(i) band(k, ch$clusters[[k]][i], own[i, ])))
        # The rest of the panel: one faded band per other cluster, one pooled
        # `rest_label` band, or nothing (own states only, no NA cap either).
        others <- setdiff(cn, k)
        switch(rest,
          clusters = rbind(own_b,
                           do.call(rbind, lapply(others, function(j)
                             band(k, .mcml_other_key(ch, j), cluster_share(j)))),
                           na_band(k)),
          pooled   = rbind(own_b,
                           band(k, ch$rest_label,
                                Reduce(`+`, lapply(others, cluster_share), 0)),
                           na_band(k)),
          none     = own_b)
      })))
    bands$key <- factor(bands$key,
                        levels = .mcml_uniq_levels(ch$all_states, ch$macro_keys, cn,
                                                   faded_keys, ch$rest_label, "NA"))
    fillvals  <- .mcml_uniq_values(pals$state, pals$cluster, pals$faded,
                                   pals$rest, "NA" = na_color)
    y_lab     <- if (isTRUE(na)) "Share of subjects" else "Share of active subjects"
  }
  bands$channel <- factor(bands$channel, levels = chan_levels)

  # Restrict to one half when the caller wants the macro drawn on its own.
  if (!identical(keep, "all")) {
    want <- if (identical(keep, "summary")) "Summary" else setdiff(chan_levels, "Summary")
    bands <- bands[bands$channel %in% want, , drop = FALSE]
    bands$channel <- droplevels(bands$channel)
    bands$key     <- droplevels(bands$key)
    if (is.null(main)) {
      main <- if (identical(keep, "summary")) {
        "Macro: cluster composition"
      } else {
        "Within-cluster: state composition"
      }
    }
  }

  .mcml_stack_channels(
    bands, ch, values = fillvals,
    layer = function(d) ggplot2::ggplot(d, ggplot2::aes(x = time, y = prop, fill = key)) +
      ggplot2::geom_area(position = ggplot2::position_stack(reverse = TRUE)) +
      ggplot2::scale_y_continuous(expand = c(0, 0),
                                  labels = scales::percent_format(accuracy = 1)) +
      ggplot2::coord_cartesian(ylim = c(0, 1)),
    y_lab = y_lab, main = main, time_label = time_label,
    theme_extra = ggplot2::theme(panel.grid.minor = ggplot2::element_blank()))
}

# ---- internal: mcml dispatcher (called from sequence_plot) ------------------
.sequence_plot_mcml <- function(x, type, normalize, state_colors, na_color,
                                main, time_label, trim = NULL,
                                panel = "both", expand = NULL,
                                combine = NULL, rest = "clusters", na = TRUE,
                                rest_label = "Other states") {
  panel <- match.arg(panel, c("both", "summary", "channels"))
  stopifnot("`rest_label` must be a single non-empty string" =
              is.character(rest_label) && length(rest_label) == 1L &&
              !is.na(rest_label) && nzchar(rest_label))
  ch   <- .mcml_seq_channels(x, trim, expand, combine)
  if (rest_label %in% c(ch$all_states, ch$cluster_names, ch$macro_keys, "NA")) {
    stop("`rest_label` \"", rest_label, "\" is already a state or cluster ",
         "name, or the reserved \"NA\" key.", call. = FALSE)
  }
  ch$rest_label <- rest_label
  pals <- .mcml_seq_palettes(ch, state_colors)
  if (type %in% c("heatmap", "index")) {
    if (!identical(panel, "both")) {
      stop("`panel` applies to type = \"distribution\".", call. = FALSE)
    }
    return(.mcml_index_plot(ch, pals, main, time_label, rest = rest))
  }
  # The macro channel is keyed by cluster and the rest by state, so a shared
  # fill scale can put a cluster and a state on one colour. Drawing a single
  # panel gives that panel its own palette, legend and default title.
  .mcml_dist_plot(ch, pals, na_color, normalize, main, time_label,
                  keep = switch(panel, both = "all", panel), rest = rest,
                  na = na)
}
