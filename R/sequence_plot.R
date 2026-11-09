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
        stop(errorCondition(
          paste0(
            "netobject_group members have no $data; rebuild with ",
            "build_network() carrying data, or pass the original ",
            "clustering object instead."),
          class = "nestimate_bad_input", call = NULL))
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
      stop(errorCondition(
        "netobject has no $data. Pass the original sequence data.",
        class = "nestimate_bad_input", call = NULL))
    }
    data <- x$data
  }
  # --- tna model ---
  else if (inherits(x, "tna")) {
    if (is.null(x$data) || is.null(x$labels)) {
      stop(errorCondition("tna object missing $data or $labels.",
                          class = "nestimate_bad_input", call = NULL))
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
    stop(errorCondition(
      paste0(
        "Unsupported input type: ",
        paste(class(x), collapse = ", "),
        ". Expected data.frame, matrix, netobject, netobject_group, ",
        "net_clustering, net_mmm, or tna."),
      class = "nestimate_bad_input", call = NULL))
  }

  list(data = data, group = group, clustering = clustering)
}

#' Sequence Plot (heatmap, index, or distribution)
#'
#' @description
#' `sequence_plot()` draws categorical sequences: one row per sequence and
#' one coloured cell per time point, or the share of sequences in each state
#' at every time point. Three views share one call:
#' \itemize{
#'   \item `type = "heatmap"` (default): a dense carpet in a single panel,
#'     rows ordered by `sort` or a dendrogram.
#'   \item `type = "index"`: the same layout with optional thin gaps between
#'     rows and no dendrogram, split into one panel per group (`group`, or
#'     the clusters of a clustering) on an `ncol` x `nrow` grid.
#'   \item `type = "distribution"`: the stacked state distribution over
#'     time, drawn by [distribution_plot()].
#' }
#' It reads wide sequences, event logs (with `action`), networks built from
#' sequences, clusterings, mixture fits and `tna` models; `periods` gives one
#' row per actor with the time axis split into days, weeks or another
#' period. An `mcml` from [build_mcml()] gets a multichannel figure instead:
#' a macro `Summary` panel and one panel per cluster (see the section
#' *Multichannel view of an mcml*). That figure is drawn with ggplot2 and
#' takes only `type`, `trim`, `state_colors`, `na_color`, `main`,
#' `time_label` (or `xlab`) and the mcml options (`normalize`, `na`,
#' `panel`, `expand`, `combine`, `rest`, `rest_label`); every other figure
#' is drawn with base graphics.
#'
#' The plot describes the sequences; it does not model them. The row order
#' of a heatmap comes from `sort` and is a display choice, so blocks of
#' similar rows are a reason to cluster the sequences with
#' [build_clusters()], not a clustering in themselves. `trim` changes only
#' what is drawn, never the data.
#'
#' @param x The sequences to draw:
#'   \describe{
#'     \item{event log}{A data frame with one row per action and `action`
#'       naming the action column (with `actor`, `time` and `order` as
#'       needed). It is converted with [prepare()]: one row per actor, or
#'       per actor-period with `periods`.}
#'     \item{data frame or matrix}{One row per sequence, one column per
#'       time point; for a data frame, `cols` names the position columns
#'       when they are not a numbered run.}
#'     \item{netobject}{A network built from sequences; its sequence data is
#'       drawn. A network without its sequence data (built from a matrix)
#'       raises an error, and a one-hot or association network raises an
#'       error of class `nestimate_bad_input`.}
#'     \item{net_clustering}{From [build_clusters()]: its sequences, with
#'       the cluster assignments as groups and its distance matrix as the
#'       heatmap's dendrogram.}
#'     \item{netobject_group}{From [cluster_network()] or [build_network()]
#'       on a clustering: the full sequences and assignments carried by the
#'       clustering, or, for a plain group split, the groups' sequences
#'       stacked and labelled by group.}
#'     \item{net_mmm}{From [build_mmm()]: its sequences (or, when the fit
#'       stores none, those of its first component model), with the hard
#'       assignments as groups.}
#'     \item{tna}{A `tna` model; its integer-coded sequences are decoded to
#'       their state labels.}
#'     \item{mcml}{From [build_mcml()], built from sequences (an `mcml`
#'       built from a matrix carries none and raises an error). It gives
#'       the multichannel figure: one panel per cluster plus a macro
#'       `Summary` panel. `type = "heatmap"` and `"index"` draw the carpet
#'       (each channel's own states solid, the other clusters a faded
#'       wash); `type = "distribution"` draws the stacked distribution (add
#'       `normalize = TRUE` for a TraMineR-style `seqdplot`, where each time
#'       point sums to 1).}
#'   }
#'   For the `print()` method: an `mcml_sequence_plot`.
#' @param type The view: `"heatmap"` (default), `"index"` or
#'   `"distribution"`.
#' @param sort How rows are ordered, in the heatmap and within each panel of
#'   the index view. `"lcs"` (default) orders them by a Ward (`ward.D2`)
#'   hierarchical clustering of their longest-common-subsequence distances,
#'   and any other [build_clusters()] distance (`"hamming"`, `"osa"`,
#'   `"lv"`, `"dl"`, `"qgram"`, `"cosine"`, `"jaccard"`, `"jw"`) does the
#'   same with that distance. `"frequency"` clusters the rows by how often
#'   each state occurs in them. `"start"` sorts the rows alphabetically by
#'   their states from the first time point on, `"end"` from the last time
#'   point back; these two draw no dendrogram. A clustering input (a
#'   `net_clustering`, or a `netobject_group` from [cluster_network()])
#'   orders the heatmap by a hierarchical clustering of its own distance
#'   matrix, using its linkage when it was clustered hierarchically and
#'   `ward.D2` otherwise (for example with the default `pam`), unless `tree`
#'   is given.
#' @param tree A `hclust`, `dendrogram` or `agnes` object that gives the row
#'   order of the heatmap and is drawn beside it; it overrides `sort`. It
#'   must have one leaf per row. Default `NULL`.
#' @param group The groups to draw in separate panels, one value per
#'   sequence (`nrow(x)` values), for the index and distribution views; the
#'   heatmap ignores it. `NULL` (default) uses the cluster assignments of a
#'   `net_clustering`, `net_mmm` or `netobject_group` input, and draws one
#'   panel otherwise. With `periods`, `group` instead names a column of the
#'   event log (or of the network's metadata) that is constant within
#'   actor, such as a class or a condition: one panel per value, one row
#'   per actor, all panels on the same period axis. With an event log it
#'   works only together with `periods`; otherwise it raises an error of
#'   class `nestimate_bad_input`.
#' @param periods The column that splits each actor's sequence into periods
#'   (days, weeks, ...): a column of the event log, or of the network's
#'   metadata. `x` must be an event log with `actor` and `action`, or a
#'   netobject built with `actor` and `session = periods`, so that each
#'   actor-period pair holds one sequence; anything else raises an error of
#'   class `nestimate_bad_input`. The plot then has one row per actor and
#'   the time axis is split into equal-width blocks, one per period. A
#'   `Date` column gives one block per calendar day and a whole-number
#'   column (a day or week counter) one block per number, from the first to
#'   the last, so a period in which nobody acted is an empty block; factor
#'   levels give their order (an unused level is an empty block); other
#'   columns are sorted. Ticks mark the period boundaries and each label
#'   sits under its block; `tick` thins the labels, and the x-axis label
#'   defaults to the column name. Heatmap and index views only (the
#'   distribution view raises an error of class `nestimate_bad_input`).
#'   `trim` then counts periods, not time points. Default `NULL` plots the
#'   sequences as they are.
#' @param fill With `periods`, how an actor's steps fill a block.
#'   `"stretch"` (default) spreads them over the whole block, so one step
#'   fills it and three steps take a third each; the cell width differs
#'   between actors. `"pad"` gives every block as many fixed slots as the
#'   busiest actor-period and leaves unused slots as `NA`, so a slot is the
#'   same step position for every actor (minor ticks mark the slots when
#'   there are at most 60 in all). With `periods`, empty space is white
#'   unless `na_color` is given, so it is never confused with a state drawn
#'   in grey. Repeats of a state merge into one spell; `cell_border =
#'   "white"` (TraMineR's `border`) outlines each action so they show as
#'   separate, equal shares of the period. Setting `fill` without `periods`
#'   raises an error of class `nestimate_bad_input`; an `mcml` ignores it.
#' @param actor,action,time,order For an event log: column names, as in
#'   [build_network()]. Giving `action` declares `x` an event log (one row
#'   per action); `actor` identifies who acted (one row per actor in the
#'   plot), and `time` and `order` put each actor's actions in order
#'   (otherwise row order is used). Without `periods`, a `time` column also
#'   splits an actor's actions into sessions at gaps over 15 minutes, as
#'   [prepare()] does; with `periods` it does not, so each actor-period is
#'   one sequence. `actor`, `time` or `order` without `action`, a column
#'   that is not in `x`, and a missing actor or period raise an error of
#'   class `nestimate_bad_input`. Default `NULL` each.
#' @param scale,geom,na For the distribution view, passed to
#'   [distribution_plot()]: `scale` is `"proportion"` (default) or
#'   `"count"`, `geom` is `"area"` (default) or `"bar"`, and `na = TRUE`
#'   (default) draws missing cells as a band of their own. For an `mcml`
#'   distribution, `na = FALSE` drops the `NA` (ended) band and shows every
#'   time point as shares of the sequences still running there, so each
#'   panel stacks to 100 percent (cluster panels only with `rest =
#'   "clusters"` or `"pooled"`; a time point where no sequence is running
#'   stays empty).
#' @param normalize For an `mcml` distribution only: `TRUE` normalises each
#'   time point to sum to 1 within its channel (the composition of
#'   TraMineR's `seqdplot`); `FALSE` (default) shows prevalence, capped with
#'   an `NA` band.
#' @param expand For an `mcml`: the clusters whose member states are shown
#'   individually in the `Summary` panel; `"all"` or `TRUE` expands every
#'   cluster. The per-cluster panels are unaffected. Default `NULL` keys the
#'   `Summary` panel by cluster.
#' @param combine For an `mcml`: clusters to merge into one channel. A
#'   character vector merges one group (e.g. `combine = c("Cognitive",
#'   "Affective")`); a list merges several, and its names label the merged
#'   channels (default label `"Cognitive + Affective"`). A merged group acts
#'   as one cluster throughout the figure: one per-cluster panel holding all
#'   its states, one key in the `Summary` panel, and one faded band in the
#'   other panels. `expand` is resolved after merging, so it can name the
#'   merged label. An unknown cluster, a group of fewer than two, or a
#'   cluster in two groups raises an error. Default `NULL` draws the
#'   partition as built.
#' @param rest For an `mcml`: how a cluster's panel shows the time its
#'   subjects spend in *other* clusters. `"clusters"` (default) draws one
#'   faded band (or wash, in the carpet) per other cluster; `"pooled"` draws
#'   all other clusters as one grey band labelled `rest_label`; `"none"`
#'   leaves them blank, so the panel shows only its own states. In the
#'   distribution view `"none"` drops the `NA` (ended) band too, and the
#'   panel's height at each time point is the share of subjects in that
#'   cluster. Ignored with `normalize = TRUE`, which rescales each panel to
#'   its own states. The `Summary` panel is unaffected.
#' @param rest_label For an `mcml`: the legend text for time spent in other
#'   clusters. Default `"Other states"`; `"Others"` or `"Rest of states"`
#'   are alternatives. The pooled band (`rest = "pooled"`) takes it as is;
#'   the per-cluster bands read `"Social (Other states)"`. A label equal to
#'   a state or cluster name, or to `"NA"`, raises an error. `combine`,
#'   `expand`, `rest` and `rest_label` on any input other than an `mcml`
#'   raise an error.
#' @param panel For an `mcml` distribution only: which panels to draw.
#'   `"both"` (default) stacks the macro `Summary` channel and the
#'   per-cluster channels on one figure, each with its own legend;
#'   `"summary"` draws the macro channel alone, keyed and coloured by
#'   cluster; `"channels"` draws the per-cluster channels alone, keyed and
#'   coloured by state. The macro channel is keyed by cluster and the rest
#'   by state, so a cluster and a state can land on the same colour; drawing
#'   one panel avoids that and gives it a default title. With an `mcml` heatmap
#'   or index view, a value other than `"both"` raises an error.
#' @param trim A cut of the time axis (one number) or of the time axis and
#'   the rows (two numbers), to stop a few long sequences from stretching
#'   the plot and a few near-empty ones from padding it. It applies to all
#'   three views, the `mcml` multichannel figure included. Default `NULL`
#'   plots everything.
#'
#'   The first number cuts the time axis. A fraction in (0, 1) judges each
#'   time point on its own activity and keeps the busiest ones that together
#'   hold that share of the observed states, ties going to the earlier one;
#'   the rest are removed wherever they sit and the kept ones close up under
#'   their own labels. So `trim = 0.9` removes the quiet time points holding
#'   the last tenth of the states, which is what a long tail consists of,
#'   and those with no data. A value `>= 1` is an absolute cut (`trim = 50`
#'   keeps the first 50 time points). With `periods` the unit is the period
#'   block: `trim = 60` keeps the first 60 periods and a fraction judges
#'   each period by its actions.
#'
#'   The second number ends each row earlier and never removes a row,
#'   counted over the time points kept by the first: a value `>= 1` keeps
#'   each row's first that many observed states (`trim = c(60, 5)`), a
#'   fraction its first share of them (`trim = c(0.9, 0.9)`). What is left
#'   stays in place; the row is white after that point. `NA` leaves an axis
#'   alone (`trim = c(NA, 5)`). Either way only the drawing changes, not the
#'   data.
#' @param trim_clusterwise For a grouped index or distribution view with a
#'   fractional `trim`: `FALSE` (default) computes one cut on the pooled
#'   data and applies it to every panel, so all panels share the same width
#'   and the time axes stay aligned; `TRUE` judges the time points within
#'   each group, so panels can end up at different widths. An absolute
#'   `trim` (`>= 1`) ignores it, since the cut is the same everywhere.
#' @param row_gap The vertical gap between rows of the index view, as a
#'   fraction of the row height, in \[0, 1). Default `0` (dense, like the
#'   heatmap); `0.15` gives visible separators when there are few rows.
#' @param dendrogram_width The width of the heatmap's dendrogram panel,
#'   against a carpet width of 4: the default `1.2` makes the dendrogram 0.3
#'   times as wide as the carpet. Used only when a dendrogram is drawn.
#' @param k For the heatmap: the number of clusters to cut the dendrogram
#'   into, a single number of at least 2 and less than the number of rows;
#'   thin horizontal lines then separate them in the carpet. Ignored
#'   without a dendrogram (e.g. `sort = "start"`) and in the other views.
#'   Default `NULL` (no separators).
#' @param k_color The colour of the cluster separator lines. Default
#'   `"white"`.
#' @param k_line_width The line width of the cluster separators. Default
#'   `2.5`.
#' @param state_colors Colours for the fill keys, in one of two forms:
#'   *unnamed*, one colour per state in level order (states are ordered as
#'   `sort(unique(...))`); or *named*, a lookup where only the keys you name
#'   are overridden and every other key keeps its default (Okabe-Ito).
#'   Names this figure does not draw are dropped with a message naming
#'   them, so one project-wide palette can be handed to every plot and each
#'   takes the keys that apply to it. Default `NULL` uses the palette
#'   attached to `x` with [set_state_colors()], if any; a palette given
#'   here wins for this figure.
#'
#'   For an `mcml` the named form reaches the whole figure, not only the
#'   states: a cluster name colours its `Summary` band, its channel strip
#'   and its faded band in the other panels, and a group merged by
#'   `combine` is named by its label (the list name you gave it, or `"A +
#'   B"`). `rest_label` is a key too. So `state_colors = c(plan =
#'   "#0072B2", "Planning + Monitoring" = "#D55E00")` recolours one state
#'   and one combined cluster and leaves the rest of the palette alone.
#' @param na_color The colour of missing cells. Default `"grey90"`; with
#'   `periods`, `"white"` unless given.
#' @param cell_border The colour of each cell's border. Default `NA` (no
#'   border).
#' @param frame Whether to draw a box around each panel; axis ticks and
#'   labels appear either way. Default `FALSE`.
#' @param width,height Device dimensions in inches. In an interactive
#'   session, giving either opens a new graphics device with
#'   `grDevices::dev.new()` (the other defaults to 7); otherwise they are
#'   ignored. In knitr chunks, use the `fig.width` and `fig.height` chunk
#'   options instead. Default `NULL`.
#' @param main The plot title. Default `NULL`: no title, apart from
#'   `"n = N"` when `show_n = TRUE`. With several groups, `main` is drawn
#'   once above the panels, and not at all when `combined = FALSE`.
#' @param show_n Whether to append `"(n = N)"` to the title (to each panel
#'   title when grouped). Default `TRUE`.
#' @param time_label,xlab The x-axis label; `xlab` is an alias that wins
#'   when given. Default `"Time"` (with `periods`, the period column's
#'   name).
#' @param y_label,ylab The y-axis label of the distribution view; `ylab` is
#'   an alias. Default `NULL`, which [distribution_plot()] reads as
#'   `"Proportion"` or `"Count"` according to `scale`.
#' @param tick Label every `tick`-th time point on the x-axis. Default
#'   `NULL` thins the labels to about 15.
#' @param ncol,nrow The panel grid of the grouped index and distribution
#'   views. Default `NULL`: with neither given, `ncol = ceiling(sqrt(G))`
#'   for `G` groups and `nrow = ceiling(G / ncol)`; with one given, the
#'   other follows. With `combined = FALSE` the grid does not affect the
#'   layout, but `ncol * nrow` must still be at least the number of groups
#'   or an error is raised.
#' @param combined For the grouped index and distribution views: `TRUE`
#'   (default) arranges the groups on one figure with `graphics::layout()`;
#'   `FALSE` draws each group on its own page, at full size with its own
#'   legend. A single group ignores it, and the heatmap is always one
#'   figure.
#' @param legend The legend position: `"right"`, `"bottom"` or `"none"`.
#'   Default `NULL`, read as `"right"` for every view.
#' @param legend_size The legend text size. Default `NULL` scales it with
#'   the device width so the legend stays in proportion on a 5-inch and a
#'   12-inch figure (between 0.65 and 1.2).
#' @param legend_title The legend title. Default `NULL` (none).
#' @param legend_ncol The number of legend columns. Default `NULL`: one row
#'   for a bottom legend, one column for a right legend.
#' @param legend_border The colour of the legend swatches' borders. Default
#'   `NA` (none).
#' @param legend_bty The legend box: `"n"` (default, none) or `"o"` (boxed).
#' @param cols For a wide data frame: the columns holding the positions.
#'   Needed when they are not a numbered run (`T1`, `T2`, ...), when one of
#'   them is called `action`, or when any other column sits beside them;
#'   given for any other input it raises an error of class
#'   `nestimate_bad_input`, except for an `mcml`, which ignores it.
#'   Default `NULL`.
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
#'     panel stacks to 100 percent unless \code{rest = "none"}).
#'     \code{normalize = TRUE} instead rescales each cluster panel to its
#'     own states, which ignores \code{rest} and \code{na}.
#' }
#' \code{panel} draws the Summary or the cluster panels alone, and
#' \code{trim} cuts the time axis for every panel at once.
#'
#' @return For an `mcml`: the multichannel figure, one panel per channel
#'   (the macro `Summary` and one per cluster), stacked, each with its own
#'   legend of its clusters or states. When more than one channel is drawn
#'   it is an `mcml_sequence_plot`, a `gtable` that its `print()` method
#'   draws on a new page; when `panel = "summary"` leaves a single channel it
#'   is a `ggplot`. The `print()` method returns its input invisibly.
#'
#'   Every other input is drawn with base graphics, and the call returns,
#'   invisibly, a list describing the figure, whose fields depend on the
#'   view:
#'   \describe{
#'     \item{heatmap}{`ord` (the row order plotted, as integer indices),
#'       `codes` (the integer-coded, trimmed sequence matrix), `palette`,
#'       `levels` (the state labels, parallel to `palette`) and `sort_used`
#'       (the ordering applied; `"net_clustering"` when a clustering's
#'       dendrogram was used).}
#'     \item{index}{`codes`, `palette`, `levels`, `orders` (a list of row
#'       orders, one per panel, indexing the original rows) and `groups`
#'       (the panel labels).}
#'     \item{with `periods`}{`ord` (the row order of the actors), `codes`
#'       (an integer matrix, one row per actor and one column per period
#'       slot in the `fill = "pad"` layout, named `"<period>.<slot>"`),
#'       `palette`, `levels`, `periods` (the period labels in plotted
#'       order, after `trim`) and `fill`; with `group`, also `groups` (the
#'       panel labels) and `orders` (one row order per panel, indexing the
#'       actors), `ord` being their concatenation.}
#'     \item{distribution}{What [distribution_plot()] returns: `counts`,
#'       `proportions`, `levels`, `palette` and `groups`.}
#'   }
#' @seealso \code{\link{distribution_plot}}, \code{\link{build_clusters}},
#'   \code{\link{build_mcml}}
#' @examples
#' sequence_plot(trajectories)
#'
#' \donttest{
#' sequence_plot(trajectories, type = "index")
#' sequence_plot(trajectories, type = "distribution")
#'
#' # One row per student, one block per day; each day's actions fill it.
#' daily <- data.frame(
#'   student = rep(c("s1", "s2", "s3"), times = c(5, 9, 6)),
#'   day     = c("d1", "d2", "d2", "d3", "d3",
#'               "d1", "d1", "d2", "d2", "d2", "d3", "d3", "d3", "d3",
#'               "d1", "d1", "d1", "d2", "d3", "d3"),
#'   action  = c("plan", "read", "write", "read", "plan",
#'               "read", "read", "plan", "write", "write", "read", "plan",
#'               "write", "read", "plan", "plan", "write", "read", "write",
#'               "read"))
#' sequence_plot(daily, actor = "student", action = "action", periods = "day")
#' sequence_plot(daily, actor = "student", action = "action", periods = "day",
#'               fill = "pad", type = "index")
#'
#' # The same from a network built with the period as its session.
#' daily_net <- build_network(daily, method = "relative", actor = "student",
#'                            action = "action", session = "day")
#' sequence_plot(daily_net, periods = "day")
#'
#' # Multichannel MCML view: one channel per cluster + a macro Summary.
#' fit <- build_mcml(
#'   group_regulation_long,
#'   clusters = list(Cognitive  = c("discuss", "synthesis", "consensus", "cohesion"),
#'                   Regulation = c("plan", "monitor", "adapt", "coregulate"),
#'                   Affective  = "emotion"),
#'   actor = "Actor", action = "Action", time = "Time")
#' sequence_plot(fit)                                          # multichannel carpet
#'
#' # Shape the multichannel view (see the section above).
#' sequence_plot(fit, type = "distribution",
#'               combine = list(Task = c("Cognitive", "Regulation")),
#'               expand = "Task")                              # merge, then open
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
                          periods          = NULL,
                          fill             = c("stretch", "pad"),
                          actor            = NULL,
                          action           = NULL,
                          time             = NULL,
                          order            = NULL,
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
                          legend_bty       = "n",
                          cols             = NULL) {

  type <- match.arg(type)
  sort <- match.arg(sort)
  # A palette attached with set_state_colors() is the object's own default;
  # an explicit `state_colors` still wins for this one figure.
  state_colors <- state_colors %||% .stored_state_colors(x)
  .assert(is.logical(combined), length(combined) == 1L, !is.na(combined),
            is.logical(trim_clusterwise), length(trim_clusterwise) == 1L)
  if (is.null(legend)) legend <- "right"
  legend <- match.arg(legend, c("bottom", "right", "none"))
  if (!is.null(xlab)) time_label <- xlab

  # mcml: multichannel sequence / distribution plot (ggplot, faceted by
  # cluster-channel). Self-contained path; returns a ggplot object.
  if (!inherits(x, "mcml") &&
      (!is.null(combine) || !is.null(expand) || !missing(rest) ||
       !missing(rest_label))) {
    stop(errorCondition(
      paste0(
        "`combine`, `expand`, `rest` and `rest_label` apply to an mcml ",
        "(from build_mcml()); got '",
        class(x)[1L],
        "'."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (inherits(x, "mcml") && !is.null(periods)) {
    .period_bad_input(paste0("`periods` needs a netobject built with ",
                             "`actor =` and `session =`, not an mcml."))
  }
  if (inherits(x, "mcml")) {
    return(.sequence_plot_mcml(
      x, type = type, normalize = isTRUE(normalize),
      state_colors = state_colors, na_color = na_color,
      main = main, time_label = time_label, trim = trim,
      panel = match.arg(panel), expand = expand, combine = combine,
      rest = match.arg(rest), na = na, rest_label = rest_label))
  }

  if (is.null(periods) && !missing(fill)) {
    .period_bad_input("`fill` applies only with `periods`.")
  }
  if (!is.null(periods) && type == "distribution") {
    .period_bad_input(paste0("`periods` works with type = \"heatmap\" or ",
                             "\"index\"."))
  }

  .refuse_onehot_network(x, "sequence_plot()")
  .refuse_binary_matrix(x, "sequence_plot()", takes_cols = TRUE)
  if (is.data.frame(x) && !inherits(x, "nestimate_data")) {
    .check_declared_columns(x, "sequence_plot()", action = action, actor = actor,
                            time = time, order = order, cols = cols)
    if (!is.null(cols)) x <- x[cols]
  } else if (!is.null(cols)) {
    stop(errorCondition("`cols` applies to a wide data frame.",
                        class = "nestimate_bad_input", call = NULL))
  }

  # Long data (one row per action): `action =` names the action column and
  # prepare() turns it into one sequence per actor (or actor-period).
  long <- NULL
  if (!is.null(action)) {
    if (!is.null(group) && is.null(periods)) {
      .period_bad_input(paste0("`group` is one value per sequence, so it ",
                               "cannot be combined with long data ",
                               "(`action =`) unless `periods` is given, ",
                               "when it names a column of the data."))
    }
    long <- .seqplot_from_long(x, actor = actor, action = action,
                               time = time, order = order, periods = periods,
                               group = if (type == "index") group)
    if (is.null(periods)) x <- long$data
  } else if (!is.null(actor) || !is.null(time) || !is.null(order)) {
    .period_bad_input(paste0("`actor`, `time` and `order` describe long ",
                             "data; name the action column with `action =` ",
                             "too."))
  }

  if (is.data.frame(x) && is.null(action) && is.null(cols) &&
      !inherits(x, "nestimate_data")) {
    .refuse_undeclared_roles(x, "sequence_plot()", taken = group, whole_frame = TRUE)
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
      legend_border = legend_border, legend_bty = legend_bty,
      cols = if (is.data.frame(x)) names(x)))
  }

  if (!is.null(periods)) {
    if (missing(time_label) && is.null(xlab)) time_label <- periods
    # Empty period space is white: the grey NA fill reads as a state when
    # Okabe-Ito deals grey (#999999) to the eighth state.
    if (missing(na_color)) na_color <- "white"
    return(.sequence_plot_periods(
      long %||% .seqplot_from_network(x, periods), periods = periods, fill = match.arg(fill), type = type, sort = sort,
      tree = tree, dendrogram_width = dendrogram_width, k = k,
      k_color = k_color, k_line_width = k_line_width, row_gap = row_gap,
      state_colors = state_colors, na_color = na_color,
      cell_border = cell_border, frame = frame, main = main,
      show_n = show_n, time_label = time_label, tick = tick,
      legend = legend, legend_size = legend_size,
      legend_title = legend_title, legend_ncol = legend_ncol,
      legend_border = legend_border, legend_bty = legend_bty, trim = trim,
      group = if (type == "index") group, ncol = ncol, nrow = nrow,
      combined = combined))
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

  tr         <- .trim_axes(trim)
  enc        <- .encode_states(x)
  codes      <- .trim_codes(enc$codes, tr$cols)
  levels_all <- enc$levels
  n_rows     <- nrow(codes); n_cols <- ncol(codes)

  # Ordered before the row half of `trim`, so a row ending earlier stays put.
  ord  <- .row_order(codes, sort, tree, levels_all)
  tree <- attr(ord, "tree")
  codes <- .trim_rows_blank(codes, tr$rows)

  if (length(ord) != n_rows) {
    stop(errorCondition(sprintf("`tree` has %d leaves but `x` has %d rows.",
                 length(ord), n_rows),
                 class = "nestimate_bad_input", call = NULL))
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
  # breaks pin code k to palette[k]; without them image() rescales to the
  # codes present, so a trimmed carpet missing some states is mis-coloured.
  graphics::image(x = seq_len(n_cols), y = seq_len(n_rows),
                  z = z, col = palette, add = TRUE, useRaster = TRUE,
                  breaks = seq_len(length(palette) + 1L) - 0.5)
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
    .assert(is.numeric(k), length(k) == 1L, k >= 2L, k < n_rows)
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

  .assert(is.numeric(row_gap), length(row_gap) == 1L,
            row_gap >= 0, row_gap < 1)

  # Extract data and group from various input types
  extracted <- .extract_seqplot_input(x, group)
  x <- extracted$data
  group <- extracted$group

  if (!is.null(group)) {
    .assert(length(group) == nrow(x))
    group <- as.factor(group)
  }

  tr         <- .trim_axes(trim)
  enc        <- .encode_states(x)
  full_codes <- enc$codes
  levels_all <- enc$levels
  .report_unused_colors(state_colors, levels_all)
  palette    <- .state_palette(state_colors, length(levels_all), levels_all)
  # One global set of columns keeps every panel the same width (aligned axes).
  # With trim_clusterwise = TRUE each group judges the columns on its own
  # activity instead, so panels can end up different widths.
  global_keep <- .trim_keep(full_codes, tr$cols)
  groups <- if (is.null(group)) factor(rep("all", nrow(full_codes))) else group
  group_levels <- levels(groups)
  G            <- length(group_levels)
  group_sizes  <- vapply(group_levels,
                         function(g) sum(groups == g), integer(1))

  # Per-group: crop columns (own set if clusterwise, else the global set), then
  # order rows within the cropped view, before the row half of `trim` so a row
  # ending earlier stays put. `order` keeps the global row indices so the
  # returned `orders` still maps back to the input rows.
  panels <- lapply(group_levels, function(g) {
    rows <- which(groups == g)
    sub  <- full_codes[rows, , drop = FALSE]
    keep <- if (isTRUE(trim_clusterwise)) .trim_keep(sub, tr$cols) else
      global_keep
    sub  <- .trim_columns(sub, keep)
    ord  <- .row_order(sub, sort, tree = NULL, levels_all)
    sub  <- .trim_rows_blank(sub, tr$rows)
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
  .assert(ncol * nrow >= G)

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
    ticks <- .tick_positions(nc, tick, colnames(sub))

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

  invisible(list(codes = .trim_rows_blank(
                   .trim_columns(full_codes, global_keep), tr$rows),
                 palette = palette, levels = levels_all,
                 orders = lapply(panels, `[[`, "order"),
                 groups = group_levels))
}


# ---- internal: row ordering for a codes matrix ----------------------------
.row_order <- function(codes, sort, tree, levels_all) {
  n_rows <- nrow(codes); n_cols <- ncol(codes)
  # One sequence (a singleton cluster, a single actor) has nothing to order
  # and nothing to cluster.
  if (n_rows < 2L && is.null(tree)) {
    return(structure(seq_len(n_rows), tree = NULL))
  }
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
