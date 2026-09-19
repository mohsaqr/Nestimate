# ==============================================================================
# sequence_plot() - unified entry point for categorical sequence plots.
# type = "heatmap"      - dense carpet with optional dendrogram (single panel).
# type = "index"        - row-gapped carpet without dendrogram; supports group
#                         + ncol/nrow facet grid (like distribution_plot).
# type = "distribution" - dispatches to distribution_plot().
# ==============================================================================

# ---- Input extraction helper ------------------------------------------------
#' Extract sequence data and clustering info from various input types
#'
#' Supports: data.frame, matrix, net_clustering, netobject, netobject_group,
#' net_mmm. Returns a list with data, group (assignments), and clustering info.
#' @noRd
.extract_seqplot_input <- function(x, group = NULL) {
  clustering <- NULL
  data <- NULL


  # --- net_clustering (from build_clusters) ---
  if (inherits(x, "net_clustering")) {
    data <- x$data
    if (is.null(group)) group <- x$assignments
    clustering <- x
  }
  # --- net_mmm (from build_mmm) ---
  else if (inherits(x, "net_mmm")) {
    data <- x$data %||% x$models[[1L]]$data
    if (is.null(group)) group <- x$assignments
    # No distance matrix in MMM; clustering info for reference only
    clustering <- list(assignments = x$assignments, k = x$k)
  }
  # --- netobject_group (from cluster_network / build_network on clustering) ---
  else if (inherits(x, "netobject_group")) {
    cl <- attr(x, "clustering")
    # Prefer the full N-row data carried by the clustering attribute --
    # cluster_network() splits per-cluster $data subsets into each member,
    # so x[[1L]]$data has only n1 rows whereas $assignments has all N.
    # build_network(net_mmm) stashes the same full data to make this invariant.
    if (!is.null(cl) && !is.null(cl$data)) {
      data <- cl$data
      if (is.null(group)) group <- cl$assignments
      clustering <- cl
    } else {
      # Fallback: no clustering attribute (plain group_col split). Rebuild
      # the full sequence frame from per-group subsets and label rows by
      # their group membership.
      parts <- lapply(seq_along(x), function(i) x[[i]]$data)
      if (any(vapply(parts, is.null, logical(1L)))) {
        stop("netobject_group members have no $data; rebuild with ",
             "build_network() carrying data, or pass the original ",
             "clustering object instead.", call. = FALSE)
      }
      part_nrow <- vapply(parts, NROW, integer(1L))
      data <- do.call(rbind, parts)
      if (is.null(group)) {
        labs  <- if (!is.null(names(x))) names(x) else as.character(seq_along(x))
        group <- factor(rep(labs, part_nrow), levels = labs)
      }
    }
  }
  # --- netobject (single network) ---
  else if (inherits(x, "netobject") || inherits(x, "cograph_network")) {
    if (is.null(x$data)) {
      stop("netobject has no $data. Pass the original sequence data.",
           call. = FALSE)
    }
    data <- x$data
  }
  # --- tna model ---
  else if (inherits(x, "tna")) {
    if (is.null(x$data) || is.null(x$labels)) {
      stop("tna object missing $data or $labels.", call. = FALSE)
    }
    # Decode integer matrix back to state names
    decoded <- matrix(x$labels[x$data], nrow = nrow(x$data),
                      ncol = ncol(x$data))
    colnames(decoded) <- colnames(x$data)
    data <- as.data.frame(decoded, stringsAsFactors = FALSE)
  }
  # --- data.frame or matrix (pass through) ---
  else if (is.data.frame(x) || is.matrix(x)) {
    data <- x
  }
  else {
    stop("Unsupported input type: ", paste(class(x), collapse = ", "),
         ". Expected data.frame, matrix, netobject, netobject_group, ",
         "net_clustering, net_mmm, or tna.", call. = FALSE)
  }

  list(data = data, group = group, clustering = clustering)
}

#' Sequence Plot (heatmap, index, or distribution)
#'
#' Single entry point for three categorical-sequence visualisations.
#' \itemize{
#'   \item \code{type = "heatmap"} (default): dense carpet, rows reordered
#'     by \code{sort} / dendrogram (single panel).
#'   \item \code{type = "index"}: same data layout, but rows separated by
#'     thin gaps (no dendrogram). Supports grouping via \code{group} or a
#'     \code{net_clustering}, plus a \code{ncol} x \code{nrow} facet grid.
#'   \item \code{type = "distribution"}: dispatches to
#'     \code{\link{distribution_plot}}.
#' }
#'
#' @param x Wide-format sequence data. Accepts:
#'   \describe{
#'     \item{data.frame / matrix}{Rows = sequences, columns = time points.}
#'     \item{netobject}{Extracts \code{$data}.}
#'     \item{net_clustering}{From \code{\link{build_clusters}}. Uses
#'       \code{$data}, \code{$assignments} for grouping, and \code{$distance}
#'       for dendrogram.}
#'     \item{netobject_group}{From \code{\link{cluster_network}} or
#'       \code{\link{build_network}} on a clustering. Extracts data and
#'       assignments from \code{attr(, "clustering")}.
#'     }
#'     \item{net_mmm}{From \code{\link{build_mmm}}. Uses \code{$data} (falling
#'       back to \code{$models[[1]]$data}) and \code{$assignments}.}
#'     \item{tna}{From the tna package. Decodes integer-encoded sequences.}
#'     \item{mcml}{From \code{\link{build_mcml}} (built from sequences).
#'       Produces a \strong{multichannel} plot: one panel per cluster plus a
#'       macro \code{Summary} panel. \code{type = "heatmap"}/\code{"index"}
#'       draw the carpet (each channel's own states solid, other clusters a
#'       faded wash); \code{type = "distribution"} draws the stacked
#'       distribution (add \code{normalize = TRUE} for a TraMineR-style
#'       \code{seqdplot} where each time point sums to 1). See the section
#'       \emph{Multichannel view of an mcml} for the options that shape it,
#'       and \emph{Value} for what it returns.}
#'   }
#' @param type One of \code{"heatmap"} (default), \code{"index"}, or
#'   \code{"distribution"}.
#' @param sort Row-ordering strategy for heatmap / within-panel for index.
#'   One of \code{"lcs"} (default), \code{"frequency"}, \code{"start"},
#'   \code{"end"}, or any \code{\link{build_clusters}} distance
#'   (\code{"hamming"}, \code{"osa"}, \code{"lv"}, \code{"dl"},
#'   \code{"qgram"}, \code{"cosine"}, \code{"jaccard"}, \code{"jw"}).
#' @param tree Optional \code{hclust}/\code{dendrogram}/\code{agnes}
#'   object to supply row ordering (heatmap only; overrides \code{sort}).
#' @param group Optional grouping vector (length \code{nrow(x)}) producing
#'   one facet per group. Index/distribution only. Ignored for heatmap.
#' @param scale,geom,na Passed to \code{\link{distribution_plot}} when
#'   \code{type = "distribution"}.
#'   For an \code{mcml} (\code{type = "distribution"}), \code{na = FALSE}
#'   drops the \code{NA} (ended) band and shows every time point as shares
#'   of the sequences still running there, so each panel stacks to 100
#'   percent (cluster panels only with \code{rest = "clusters"} or
#'   \code{"pooled"}; a time point where no sequence is running stays empty).
#' @param normalize \code{mcml} + \code{type = "distribution"} only. When
#'   \code{TRUE}, each time point is normalised to sum to 1 within its channel
#'   (TraMineR-style \code{seqdplot} composition); when \code{FALSE} (default)
#'   the stack shows prevalence and is capped with an \code{NA} band.
#' @param expand For an \code{mcml}, names of clusters whose member states
#'   are shown individually in the Summary band; \code{"all"} or \code{TRUE}
#'   expands every cluster. The per-cluster channels are unaffected. Default
#'   \code{NULL} keys the Summary band by cluster.
#' @param combine For an \code{mcml}, clusters to merge into one channel. A
#'   character vector merges one group (e.g.
#'   \code{combine = c("Cognitive", "Affective")}); a list merges several,
#'   and its names label the merged channels (default label
#'   \code{"Cognitive + Affective"}). A merged group acts as one cluster
#'   throughout the figure: one per-cluster panel holding all its states, one
#'   key in the Summary band, and one faded band in the other panels.
#'   \code{expand} is resolved after merging, so it can name the merged
#'   label. Errors on unknown clusters, a group of fewer than two, or a
#'   cluster in two groups. Default \code{NULL} draws the partition as built.
#' @param rest For an \code{mcml}, how a cluster's panel shows the time its
#'   subjects spend in \emph{other} clusters. \code{"clusters"} (default): one
#'   faded band (or wash, in the carpet) per other cluster.
#'   \code{"pooled"}: all other clusters as one grey band labelled
#'   \code{rest_label}.
#'   \code{"none"}: left blank, so the panel shows only its own states; in the
#'   distribution view the \code{NA} (ended) band is dropped too, and the
#'   panel's height at each time point is the share of subjects in that
#'   cluster. Ignored with \code{normalize = TRUE}, which rescales each panel
#'   to its own states. The Summary panel is unaffected.
#' @param rest_label For an \code{mcml}, the legend text for time spent in
#'   other clusters. Default \code{"Other states"}; e.g. \code{"Others"} or
#'   \code{"Rest of states"}. The pooled band (\code{rest = "pooled"}) takes
#'   it as is; the per-cluster bands read \code{"Social (Other states)"}.
#'   Must not equal a state or cluster name.
#' @param panel \code{mcml} + \code{type = "distribution"} only. Which panel
#'   to draw. \code{"both"} (default) stacks the macro \code{Summary}
#'   channel and the per-cluster channels on one figure, each with its own
#'   legend; \code{"summary"} draws the macro channel alone,
#'   keyed and coloured by cluster; \code{"channels"} draws the per-cluster
#'   channels alone, keyed and coloured by state. The macro channel is keyed
#'   by cluster and the rest by state, so a cluster and a state can land on
#'   the same colour -- drawing one panel avoids that and gives it a default
#'   title.
#' @param trim Optional time-axis truncation, to stop a few long
#'   sequences from stretching the plot. Applies to all three types
#'   (including the \code{mcml} multichannel view). \code{NULL} (default)
#'   plots the full width. A fraction in \code{(0, 1)} drops everything
#'   past that quantile of sequence lengths (e.g. \code{trim = 0.95} keeps
#'   the columns covering the shortest 95\% of sequences); a value
#'   \code{>= 1} is an absolute cut (\code{trim = 50} keeps the first 50
#'   time points).
#' @param trim_clusterwise Grouped \code{type = "index"} /
#'   \code{"distribution"} only, and only when \code{trim} is a fraction.
#'   \code{FALSE} (default) computes one cutoff on the pooled data and
#'   applies it to every panel, so all facets share the same width and the
#'   time axes stay aligned. \code{TRUE} crops each group to its own
#'   length quantile, so panels can end up at different widths (ragged
#'   axes). Absolute \code{trim} (\code{>= 1}) ignores this - the column is
#'   the same everywhere either way.
#' @param row_gap Fraction of row height used as vertical gap between
#'   sequences in index plots. \code{0} (default) = dense like heatmap.
#'   Try \code{0.15} for visible separators at low row counts.
#' @param dendrogram_width Width ratio of the dendrogram panel (heatmap).
#' @param k Optional integer. When supplied in \code{type = "heatmap"},
#'   cuts the dendrogram into \code{k} clusters and draws thin horizontal
#'   separators between them in the carpet. Ignored when there is no
#'   dendrogram (e.g. \code{sort = "start"}) or for other types.
#' @param k_color Colour for the cluster separator lines. Default
#'   \code{"white"}.
#' @param k_line_width Line width for the cluster separators. Default
#'   \code{2.5}.
#' @param state_colors Colours for the fill keys. Two forms:
#'   \emph{unnamed} - one colour per state, in level order (states are
#'   ordered as \code{sort(unique(...))}); \emph{named} - a lookup, where
#'   only the keys you name are overridden and every other key keeps its
#'   default. Names this figure does not draw are dropped with a message
#'   naming them, so one project-wide palette can be handed to every plot
#'   and each takes the keys that apply to it.
#'
#'   For an \code{mcml} the named form reaches the whole figure, not just
#'   the states: a cluster name colours its \code{Summary} band, its
#'   channel strip and its faded band in the other panels, and a group
#'   merged by \code{combine} is named by its label (the list name you
#'   gave it, or \code{"A + B"}). \code{rest_label} is a key too. So
#'   \code{state_colors = c(plan = "#0072B2", "Planning + Monitoring" =
#'   "#D55E00")} recolours one state and one combined cluster and leaves
#'   the rest of the palette alone.
#' @param na_color Colour for \code{NA} cells.
#' @param cell_border Cell border colour. \code{NA} (default) = off.
#' @param frame \code{FALSE} (default) draws no box - axis ticks and labels
#'   still appear. \code{TRUE} draws a box around each panel.
#' @param width,height Optional device dimensions in inches. When supplied,
#'   opens a new graphics device via \code{grDevices::dev.new()}. In knitr
#'   chunks use the \code{fig.width} / \code{fig.height} chunk options
#'   instead.
#' @param main Plot title.
#' @param show_n Append \code{"(n = N)"} to the title.
#' @param time_label,xlab X-axis label. \code{xlab} is an alias.
#' @param y_label,ylab Y-axis label (distribution only). \code{ylab} alias.
#' @param tick Show every Nth x-axis label. \code{NULL} = auto.
#' @param ncol,nrow Facet grid dimensions (index + distribution).
#'   Ignored when \code{combined = FALSE}.
#' @param combined Index and distribution types only. When \code{TRUE}
#'   (default), groups are arranged on one figure via
#'   \code{graphics::layout()}. When \code{FALSE}, each group is drawn
#'   on its own page (one full-size figure per group, with its own
#'   legend). Single-group calls (\code{G == 1}) ignore this argument.
#'   Heatmap is always single-figure.
#' @param legend Legend position: \code{"bottom"}, \code{"right"}, or
#'   \code{"none"}. \code{NULL} (default) resolves to \code{"right"} for
#'   every type.
#' @param legend_size Legend text size. \code{NULL} (default) auto-scales
#'   from the device width so the legend looks proportional at 5 in vs
#'   12 in figures (clamped to \code{[0.65, 1.2]}).
#' @param legend_title Optional legend title.
#' @param legend_ncol Number of legend columns.
#' @param legend_border Swatch border colour.
#' @param legend_bty \code{"n"} or \code{"o"}.
#'
#' @section Multichannel view of an mcml:
#' An \code{mcml} built from sequences stores, for every cluster, the full
#' sequence matrix with the other clusters' states blanked out. Each cluster
#' is therefore a \emph{channel}, and \code{sequence_plot()} stacks them:
#' \describe{
#'   \item{\code{Summary}}{The macro sequence: at every time point, the
#'     cluster each subject is in.}
#'   \item{One panel per cluster}{That cluster's own states, plus the time
#'     its subjects spend in the other clusters.}
#' }
#' The options apply in this order, so each one sees the result of the one
#' before:
#' \enumerate{
#'   \item \code{combine} merges clusters into one channel. The merged group
#'     is then one cluster throughout the figure: one panel with all its
#'     states, one Summary key, one band in the other panels. The object
#'     itself is not changed.
#'   \item \code{expand} opens clusters (including a merged group, by its
#'     label) into their member states in the Summary panel only.
#'   \item \code{rest} and \code{rest_label} set how a cluster's panel shows
#'     the other clusters: one faded band per cluster labelled
#'     \code{"<cluster> (<rest_label>)"} (\code{rest = "clusters"}), one grey
#'     band labelled \code{rest_label} (\code{"pooled"}), or nothing
#'     (\code{"none"}).
#'   \item \code{na} (distribution only) keeps the \code{NA} band of
#'     sequences that have ended (\code{TRUE}, shares of all subjects) or
#'     drops it (\code{FALSE}, shares of the subjects still running, so every
#'     panel stacks to 100 percent unless \code{rest = "none"}). \code{normalize = TRUE} instead rescales each
#'     cluster panel to its own states, which ignores \code{rest} and
#'     \code{na}.
#' }
#' \code{panel} draws the Summary or the cluster panels alone, and
#' \code{trim} cuts the time axis for every panel at once.
#'
#' @return An \code{mcml} input returns the multichannel figure: one panel
#'   per channel (the macro \code{Summary} and one per cluster), stacked,
#'   each with its own legend of its own clusters or states. When more than
#'   one channel is drawn this is an \code{mcml_sequence_plot} (a
#'   \code{gtable} whose print method draws it); when \code{panel =
#'   "summary"} leaves a single channel it is a plain \code{ggplot}. Every
#'   other input draws with base graphics and
#'   returns, invisibly, a list whose shape depends on \code{type}:
#'   \describe{
#'     \item{\code{"heatmap"}}{\code{ord} (integer row order actually
#'       plotted), \code{codes} (the integer-encoded, trimmed sequence
#'       matrix), \code{palette}, \code{levels} (state labels, parallel to
#'       \code{palette}), and \code{sort_used} (the ordering strategy
#'       applied, \code{"net_clustering"} when a clustering dendrogram was
#'       used).}
#'     \item{\code{"index"}}{\code{codes}, \code{palette}, \code{levels},
#'       \code{orders} (list of integer row orders, one per panel, indexing
#'       the original rows) and \code{groups} (panel labels).}
#'     \item{\code{"distribution"}}{Whatever
#'       \code{\link{distribution_plot}} returns: \code{counts},
#'       \code{proportions}, \code{levels}, \code{palette}, \code{groups}.}
#'   }
#' @seealso \code{\link{distribution_plot}}, \code{\link{build_clusters}},
#'   \code{\link{build_mcml}}
#' @examples
#' \donttest{
#' sequence_plot(trajectories)
#' sequence_plot(trajectories, type = "index")
#' sequence_plot(trajectories, type = "distribution")
#'
#' # Multichannel MCML view: one channel per cluster + a macro Summary.
#' fit <- build_mcml(
#'   group_regulation_long,
#'   clusters = list(Cognitive  = c("discuss", "synthesis", "consensus", "cohesion"),
#'                   Regulation = c("plan", "monitor", "adapt", "coregulate"),
#'                   Affective  = "emotion"),
#'   actor = "Actor", action = "Action", time = "Time")
#' sequence_plot(fit)                                          # multichannel carpet
#' sequence_plot(fit, type = "distribution")                  # prevalence + NA band
#' sequence_plot(fit, type = "distribution", normalize = TRUE) # seqdplot (sums to 1)
#'
#' # Shape the multichannel view (see the section above).
#' sequence_plot(fit, type = "distribution",
#'               combine = c("Cognitive", "Affective"))        # two channels as one
#' sequence_plot(fit, type = "distribution",
#'               combine = list(Task = c("Cognitive", "Regulation")),
#'               expand = "Task")                              # merge, then open
#' sequence_plot(fit, type = "distribution", panel = "channels",
#'               rest = "pooled", na = FALSE)                  # one grey band, no NA
#' sequence_plot(fit, type = "distribution", panel = "channels",
#'               rest = "pooled", rest_label = "Rest of states")
#' sequence_plot(fit, type = "distribution", panel = "channels",
#'               rest = "none")                                # own states only
#' sequence_plot(fit, rest = "pooled")                         # carpet, pooled wash
#'
#' # Colour by name: one state, one cluster, one combined group. Everything
#' # not named keeps its default colour.
#' sequence_plot(fit, type = "distribution",
#'               combine = list(Task = c("Cognitive", "Regulation")),
#'               state_colors = c(Task = "#0072B2", Affective = "#D55E00",
#'                                emotion = "#CC79A7"))
#' }
#' @export
sequence_plot <- function(x,
                          type             = c("heatmap", "index",
                                               "distribution"),
                          sort             = c("lcs", "frequency",
                                               "start", "end",
                                               "hamming", "osa", "lv", "dl",
                                               "qgram", "cosine",
                                               "jaccard", "jw"),
                          tree             = NULL,
                          group            = NULL,
                          scale            = c("proportion", "count"),
                          geom             = c("area", "bar"),
                          na               = TRUE,
                          normalize        = FALSE,
                          trim             = NULL,
                          panel            = c("both", "summary",
                                               "channels"),
                          expand           = NULL,
                          combine          = NULL,
                          rest             = c("clusters", "pooled", "none"),
                          rest_label       = "Other states",
                          trim_clusterwise = FALSE,
                          row_gap          = 0,
                          dendrogram_width = 1.2,
                          k                = NULL,
                          k_color          = "white",
                          k_line_width     = 2.5,
                          state_colors     = NULL,
                          na_color         = "grey90",
                          cell_border      = NA,
                          frame            = FALSE,
                          width            = NULL,
                          height           = NULL,
                          main             = NULL,
                          show_n           = TRUE,
                          time_label       = "Time",
                          xlab             = NULL,
                          y_label          = NULL,
                          ylab             = NULL,
                          tick             = NULL,
                          ncol             = NULL,
                          nrow             = NULL,
                          combined         = TRUE,
                          legend           = NULL,
                          legend_size      = NULL,
                          legend_title     = NULL,
                          legend_ncol      = NULL,
                          legend_border    = NA,
                          legend_bty       = "n") {

  type <- match.arg(type)
  sort <- match.arg(sort)
  # A palette attached with set_state_colors() is the object's own default;
  # an explicit `state_colors` still wins for this one figure.
  state_colors <- state_colors %||% .stored_state_colors(x)
  stopifnot(is.logical(combined), length(combined) == 1L,
            is.logical(trim_clusterwise), length(trim_clusterwise) == 1L)
  if (is.null(legend)) legend <- "right"
  legend <- match.arg(legend, c("bottom", "right", "none"))
  if (!is.null(xlab)) time_label <- xlab

  # mcml: multichannel sequence / distribution plot (ggplot, faceted by
  # cluster-channel). Self-contained path; returns a ggplot object.
  if (!inherits(x, "mcml") &&
      (!is.null(combine) || !is.null(expand) || !missing(rest) ||
       !missing(rest_label))) {
    stop("`combine`, `expand`, `rest` and `rest_label` apply to an mcml ",
         "(from build_mcml()); got '", class(x)[1L], "'.", call. = FALSE)
  }
  if (inherits(x, "mcml")) {
    return(.sequence_plot_mcml(
      x, type = type, normalize = isTRUE(normalize),
      state_colors = state_colors, na_color = na_color,
      main = main, time_label = time_label, trim = trim,
      panel = match.arg(panel), expand = expand, combine = combine,
      rest = match.arg(rest), na = na, rest_label = rest_label))
  }

  # Open a new device when width/height supplied (interactive use). In
  # knitr, set fig.width / fig.height in the chunk header instead - this
  # call would open a detached device the chunk won't capture.
  if (interactive() && (!is.null(width) || !is.null(height))) {
    grDevices::dev.new(width  = if (!is.null(width))  width  else 7,
                       height = if (!is.null(height)) height else 7,
                       noRStudioGD = TRUE)
  }

  # Adaptive legend_size: measured against the *actual* device dimensions,
  # so the same call produces a proportional legend on a 7-inch knitr figure
  # and on a 12-inch interactive window.
  legend_size <- .auto_legend_size(legend_size)

  if (type == "distribution") {
    return(distribution_plot(
      x, group = group, scale = match.arg(scale), geom = match.arg(geom),
      na = na, trim = trim, trim_clusterwise = trim_clusterwise,
      state_colors = state_colors, na_color = na_color,
      frame = frame,
      main = main, show_n = show_n,
      time_label = time_label, y_label = y_label, ylab = ylab,
      tick = tick, ncol = ncol, nrow = nrow, combined = combined,
      legend = legend, legend_size = legend_size,
      legend_title = legend_title, legend_ncol = legend_ncol,
      legend_border = legend_border, legend_bty = legend_bty))
  }

  if (type == "heatmap") {
    return(.sequence_plot_heatmap(
      x, sort, tree, dendrogram_width, k, k_color, k_line_width,
      state_colors, na_color, cell_border, frame,
      main, show_n, time_label, tick,
      legend, legend_size, legend_title, legend_ncol,
      legend_border, legend_bty, trim))
  }

  .sequence_plot_index(
    x, sort, group, row_gap,
    state_colors, na_color, cell_border, frame,
    main, show_n, time_label, tick, ncol, nrow, combined,
    legend, legend_size, legend_title, legend_ncol,
    legend_border, legend_bty, trim, trim_clusterwise)
}


# ---- Heatmap (single-panel dense carpet with optional dendrogram) ---------
.sequence_plot_heatmap <- function(x, sort, tree, dendrogram_width,
                                   k, k_color, k_line_width,
                                   state_colors, na_color, cell_border, frame,
                                   main, show_n, time_label, tick,
                                   legend, legend_size, legend_title,
                                   legend_ncol, legend_border, legend_bty,
                                   trim = NULL) {

  # Extract data from various input types
  extracted <- .extract_seqplot_input(x)
  clustering <- extracted$clustering
  sort_used <- sort

  # Build dendrogram from clustering distance matrix if available
  if (!is.null(clustering) && is.null(tree) && !is.null(clustering$distance)) {
    hc_method <- if (!is.null(clustering$method) &&
                     clustering$method %in% c("ward.D", "ward.D2", "single",
                                              "complete", "average", "mcquitty",
                                              "median", "centroid")) {
      clustering$method
    } else "ward.D2"
    tree <- stats::hclust(stats::as.dist(clustering$distance), method = hc_method)
    sort_used <- "net_clustering"
  }

  x <- extracted$data

  enc        <- .encode_states(x)
  codes      <- .trim_codes(enc$codes, trim)
  levels_all <- enc$levels
  n_rows     <- nrow(codes); n_cols <- ncol(codes)

  ord  <- .row_order(codes, sort, tree, levels_all)
  tree <- attr(ord, "tree")

  if (length(ord) != n_rows) {
    stop(sprintf("`tree` has %d leaves but `x` has %d rows.",
                 length(ord), n_rows), call. = FALSE)
  }

  .report_unused_colors(state_colors, levels_all)
  palette <- .state_palette(state_colors, length(levels_all), levels_all)
  z       <- t(codes[ord, , drop = FALSE])

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
    # mar[4] = 0.3: leave a sliver of room for the root branch, otherwise
    # the joining line where all clades merge gets clipped at the panel
    # boundary. `ylim = c(0.5, n_rows + 0.5)` aligns leaves with image
    # cells (which span y = 0.5 .. n_rows + 0.5).
    graphics::par(mar = c(mar_bottom, 0.5, mar_top, 0.3))
    graphics::plot(tree, horiz = TRUE, axes = FALSE, yaxs = "i",
                   ylim = c(0.5, n_rows + 0.5),
                   leaflab = "none", xlab = "", ylab = "")
  }
  graphics::par(mar = c(mar_bottom, if (!is.null(tree)) 0 else 2,
                        mar_top, 1))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0.5, n_cols + 0.5),
                        ylim = c(0.5, n_rows + 0.5), yaxs = "i", xaxs = "i")
  graphics::rect(0.5, 0.5, n_cols + 0.5, n_rows + 0.5,
                 col = na_color, border = NA)
  # useRaster = TRUE: render as bitmap instead of one polygon per cell.
  # Massively smaller output, crisper on hi-DPI displays, fine for our
  # uniform grid. Safe: image() itself decides not to raster if the
  # device can't handle it.
  graphics::image(x = seq_len(n_cols), y = seq_len(n_rows),
                  z = z, col = palette, add = TRUE, useRaster = TRUE)
  if (isTRUE(frame)) graphics::box()
  ticks <- .tick_positions(n_cols, tick, colnames(codes))
  graphics::axis(1, at = ticks$at, labels = ticks$labels,
                 mgp = c(1, 0.5, 0), lwd = if (isTRUE(frame)) 1 else 0,
                 lwd.ticks = 1)
  if (!is.null(time_label) && nzchar(time_label)) {
    graphics::mtext(time_label, side = 1, line = 1.7)
  }
  if (!is.na(cell_border)) {
    graphics::abline(v = seq(0.5, n_cols + 0.5, 1),
                     col = cell_border, lwd = 0.3)
    graphics::abline(h = seq(0.5, n_rows + 0.5, 1),
                     col = cell_border, lwd = 0.3)
    if (isTRUE(frame)) graphics::box()
  }

  # Cluster separators: cut the dendrogram (if any) into k groups, find
  # where group membership changes in the plotted row order, draw thin
  # lines at those row boundaries.
  if (!is.null(k) && !is.null(tree)) {
    stopifnot(is.numeric(k), length(k) == 1L, k >= 2L, k < n_rows)
    hc  <- stats::as.hclust(tree)
    mem <- stats::cutree(hc, k = as.integer(k))[ord]
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

  invisible(list(ord = ord, codes = codes, palette = palette,
                 levels = levels_all, sort_used = sort_used))
}


# ---- Index (gapped per-row bars, optional facet grid) ---------------------
.sequence_plot_index <- function(x, sort, group, row_gap,
                                 state_colors, na_color, cell_border, frame,
                                 main, show_n, time_label, tick, ncol, nrow,
                                 combined,
                                 legend, legend_size, legend_title,
                                 legend_ncol, legend_border, legend_bty,
                                 trim = NULL, trim_clusterwise = FALSE) {

  stopifnot(is.numeric(row_gap), length(row_gap) == 1L,
            row_gap >= 0, row_gap < 1)

  # Extract data and group from various input types
  extracted <- .extract_seqplot_input(x, group)
  x <- extracted$data
  group <- extracted$group

  if (!is.null(group)) {
    stopifnot(length(group) == nrow(x))
    group <- as.factor(group)
  }

  enc        <- .encode_states(x)
  full_codes <- enc$codes
  levels_all <- enc$levels
  col_names  <- colnames(full_codes)
  .report_unused_colors(state_colors, levels_all)
  palette    <- .state_palette(state_colors, length(levels_all), levels_all)
  # One global cut keeps every panel the same width (aligned axes). With
  # trim_clusterwise = TRUE each group is cropped to its own length
  # quantile instead, so panels can end up different widths.
  global_cut <- .trim_cut(full_codes, trim)

  groups       <- if (is.null(group)) factor(rep("all", nrow(full_codes))) else group
  group_levels <- levels(groups)
  G            <- length(group_levels)
  group_sizes  <- vapply(group_levels,
                         function(g) sum(groups == g), integer(1))

  # Per-group: crop columns (own cut if clusterwise, else the global cut),
  # then order rows within the cropped view. `order` keeps the global row
  # indices so the returned `orders` still maps back to the input rows.
  panels <- lapply(group_levels, function(g) {
    rows <- which(groups == g)
    sub  <- full_codes[rows, , drop = FALSE]
    cut  <- if (isTRUE(trim_clusterwise)) .trim_cut(sub, trim) else global_cut
    sub  <- sub[, seq_len(cut), drop = FALSE]
    ord  <- .row_order(sub, sort, tree = NULL, levels_all)
    list(sub = sub[ord, , drop = FALSE], order = rows[ord])
  })

  if (is.null(ncol) && is.null(nrow)) {
    ncol <- as.integer(ceiling(sqrt(G)))
    nrow <- as.integer(ceiling(G / ncol))
  } else if (is.null(ncol)) {
    ncol <- as.integer(ceiling(G / nrow))
  } else if (is.null(nrow)) {
    nrow <- as.integer(ceiling(G / ncol))
  }
  stopifnot(ncol * nrow >= G)

  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op), add = TRUE)
  oma <- .legend_oma_size(levels_all, legend, legend_size,
                          legend_ncol, legend_title)
  graphics::par(oma = c(if (legend == "bottom") oma[["oma_b"]] else 0.3,
                        0.3, 0.3,
                        if (legend == "right")  oma[["oma_r"]] else 0.3))
  use_layout <- combined && nrow * ncol > 1L
  if (use_layout) {
    lm <- matrix(c(seq_len(G), integer(nrow * ncol - G)),
                 nrow = nrow, ncol = ncol, byrow = TRUE)
    graphics::layout(lm)
  }

  mar_bottom <- if (!is.null(time_label) && nzchar(time_label)) 3 else 1
  mar_top    <- if (!is.null(main) || isTRUE(show_n) || G > 1L) 2.5 else 1.5

  invisible(lapply(seq_len(G), function(g_idx) {
    g    <- group_levels[g_idx]
    sub  <- panels[[g_idx]]$sub
    nN   <- nrow(sub)
    nc   <- ncol(sub)
    ticks <- .tick_positions(nc, tick, col_names[seq_len(nc)])

    if (!combined && g_idx > 1L) {
      graphics::par(oma = c(if (legend == "bottom") oma[["oma_b"]] else 0.3,
                            0.3, 0.3,
                            if (legend == "right")  oma[["oma_r"]] else 0.3))
    }
    graphics::par(mar = c(mar_bottom, 2, mar_top, 1))
    graphics::plot.new()
    graphics::plot.window(xlim = c(0.5, nc + 0.5),
                          ylim = c(0.5, nN + 0.5), yaxs = "i", xaxs = "i")

    # Fill whole panel with na_color so NA cells show through.
    graphics::rect(0.5, 0.5, nc + 0.5, nN + 0.5,
                   col = na_color, border = NA)

    # Vectorised per-cell rects with a vertical gap.
    ts <- rep(seq_len(nc), each  = nN)
    ss <- rep(seq_len(nN), times = nc)
    cell_col <- palette[as.vector(sub)]
    keep <- !is.na(cell_col)
    graphics::rect(xleft   = ts[keep] - 0.5,
                   ybottom = ss[keep] - 0.5 + row_gap / 2,
                   xright  = ts[keep] + 0.5,
                   ytop    = ss[keep] + 0.5 - row_gap / 2,
                   col     = cell_col[keep],
                   border  = if (is.na(cell_border)) NA else cell_border,
                   lwd     = 0.3)
    if (isTRUE(frame)) graphics::box()
    graphics::axis(1, at = ticks$at, labels = ticks$labels,
                   mgp = c(1, 0.5, 0),
                   lwd = if (isTRUE(frame)) 1 else 0, lwd.ticks = 1)
    if (!is.null(time_label) && nzchar(time_label)) {
      graphics::mtext(time_label, side = 1, line = 1.7)
    }

    panel_title <- if (G > 1L) {
      base <- sprintf("Cluster %s", g)
      if (isTRUE(show_n)) sprintf("%s  (n = %d)", base,
                                  group_sizes[g_idx]) else base
    } else .title_line(main, show_n, nN)
    if (!is.null(panel_title)) {
      graphics::mtext(panel_title, side = 3, line = 0.5, font = 2,
                      cex = if (G > 1L) 0.9 else 1)
    }
    if (!combined && legend != "none") {
      .draw_legend_in_oma(levels_all, palette, legend, legend_size,
                          legend_ncol, legend_title, legend_border, legend_bty)
    }
  }))
  if (combined && G > 1L && !is.null(main)) {
    graphics::mtext(main, side = 3, line = 1.5, font = 2,
                    outer = TRUE, cex = 1)
  }

  if (combined && legend != "none") {
    .draw_legend_in_oma(levels_all, palette, legend, legend_size,
                        legend_ncol, legend_title, legend_border, legend_bty)
  }

  invisible(list(codes = full_codes[, seq_len(global_cut), drop = FALSE],
                 palette = palette, levels = levels_all,
                 orders = lapply(panels, `[[`, "order"),
                 groups = group_levels))
}


# ---- internal: row ordering for a codes matrix ----------------------------
.row_order <- function(codes, sort, tree, levels_all) {
  n_rows <- nrow(codes); n_cols <- ncol(codes)
  if (!is.null(tree)) {
    if (!inherits(tree, "dendrogram")) tree <- stats::as.dendrogram(tree)
    return(structure(stats::order.dendrogram(tree), tree = tree))
  }
  if (sort == "start") {
    ord <- do.call(order, c(lapply(seq_len(n_cols), function(t) codes[, t]),
                            list(na.last = TRUE)))
    return(structure(ord, tree = NULL))
  }
  if (sort == "end") {
    ord <- do.call(order, c(lapply(rev(seq_len(n_cols)),
                                   function(t) codes[, t]),
                            list(na.last = TRUE)))
    return(structure(ord, tree = NULL))
  }
  if (sort == "frequency") {
    counts <- vapply(seq_along(levels_all),
                     function(k) rowSums(codes == k, na.rm = TRUE),
                     numeric(n_rows))
    tr <- stats::as.dendrogram(stats::hclust(stats::dist(counts),
                                             method = "ward.D2"))
    return(structure(stats::order.dendrogram(tr), tree = tr))
  }
  # distance-based
  d  <- .sequence_distance(codes_to_char(codes, levels_all),
                           dissimilarity = sort)
  tr <- stats::as.dendrogram(stats::hclust(stats::as.dist(d),
                                           method = "ward.D2"))
  structure(stats::order.dendrogram(tr), tree = tr)
}


# ---- internal: decode integer codes back to character (for .sequence_distance)
codes_to_char <- function(codes, levels_all) {
  out <- matrix(NA_character_, nrow = nrow(codes), ncol = ncol(codes))
  nz  <- !is.na(codes)
  out[nz] <- levels_all[codes[nz]]
  out
}


# ---- internal: standard title with "(n = N)" suffix ----------------------
.title_line <- function(main, show_n, n) {
  if (!is.null(main)) {
    if (isTRUE(show_n)) sprintf("%s  (n = %d)", main, n) else main
  } else if (isTRUE(show_n)) sprintf("n = %d", n) else NULL
}
