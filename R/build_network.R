#' Build a Network
#'
#' @description
#' `build_network()` estimates a network from sequences, an event log,
#' one-hot codes or numeric variables, and returns a `netobject`: the object
#' every other Nestimate verb reads (bootstrap, permutation, centrality,
#' clustering, plotting). One call covers two kinds of network.
#' Transition networks (`"relative"`, `"frequency"`, `"co_occurrence"` and
#' their variants) are directed or undirected graphs of states, where an edge
#' says how often or how likely one state follows another. Association
#' networks (`"cor"`, `"pcor"`, `"glasso"`, `"ising"`, `"mgm"`) are undirected
#' graphs of variables, where an edge is a correlation or a (regularised)
#' partial correlation. The *Estimators* section lists every method and its
#' aliases; estimators added with [register_estimator()] are called the same
#' way.
#'
#' Column roles are declared, never read from names or values. An event log
#' needs `action`, and usually `actor`, `time` or `session`; a wide table
#' needs position columns named as a numbered run (`T1`, `T2`, ...) or listed
#' in `state_cols`; one-hot data needs `codes` or `format = "onehot"`. Input
#' that fits none of these raises an error of class `nestimate_bad_input`
#' instead of being read one way or another.
#'
#' The weights are point estimates. Before reading single edges, see how far
#' they move under resampling with [bootstrap_network()], and compare groups
#' with [permutation()] rather than by eye. An event log that yields a single
#' sequence is reported with a message of class `nestimate_single_sequence`:
#' that network describes one sequence and cannot be resampled.
#'
#' @param data The data to estimate from: a data frame of sequences (one row
#'   per sequence, one column per position), an event log (one row per
#'   event, with `action` set), one-hot code columns, numeric variables (one
#'   row per observation) for the association methods, or a square symmetric
#'   correlation or covariance matrix. A fitted `net_clustering` or `net_mmm`
#'   object also works: its per-cluster networks are built, with the
#'   `scaling`, `threshold`, `start`, `end` and `params` given here, and
#'   returned as a `netobject_group`. The relative networks of a `net_mmm`
#'   are the fitted chains: `scaling` and `threshold` are applied to them,
#'   and any other argument raises an error of class `nestimate_bad_input`.
#' @param method The estimator: a built-in method or alias (see the
#'   *Estimators* section), or a name added with [register_estimator()];
#'   [list_estimators()] lists them. Required, except when `data` is a
#'   `net_clustering` or `net_mmm` object, whose own method is used.
#' @param actor The column naming who each row belongs to (a person, a team).
#'   In an event log, each actor's events form a sequence, split further by
#'   `session` or `time_threshold`. With `level`, it is the id the multilevel
#'   decomposition groups by. Default `NULL`.
#' @param action The column holding each event's state; setting it declares
#'   the data an event log. Column names give no role: `action`, `actor`,
#'   `time` and `session` are read only when set. A column called `Action`
#'   with nothing declared raises an error of class `nestimate_bad_input`,
#'   and an undeclared time or session column beside `action` is named in a
#'   message of class `nestimate_column_role`. Default `NULL`.
#' @param time The timestamp column of an event log. Events are ordered by
#'   it, and a gap longer than `time_threshold` starts a new sequence.
#'   Default `NULL`.
#' @param session The column identifying sessions in an event log; each
#'   session is its own sequence, also without `actor`. Default `NULL`.
#' @param order The column that orders events sharing a timestamp; without
#'   it, those events keep their row order. Default `NULL`.
#' @param codes The one-hot code columns (one per state, 1 when the state is
#'   active); setting them declares one-hot data. Default `NULL`.
#' @param group A column to split the data by: one network is estimated per
#'   value and the result is a `netobject_group`. Default `NULL`.
#' @param format How to read `data`. `"auto"` (default) reads an event log
#'   when `action` is set, one-hot data when `codes` is set, and wide
#'   positions otherwise; it never infers one-hot data from the values.
#'   `"wide"`, `"long"` and `"onehot"` force a reading. The association
#'   methods and the `wtna` methods read 0/1 columns as they always do.
#' @param window_size For one-hot data: how many consecutive rows form one
#'   window (default `3`). Windows on sequence data are set through `params`
#'   instead.
#' @param mode For one-hot data: `"non-overlapping"` (default) or
#'   `"overlapping"` windows; [wtna()] describes how overlapping windows are
#'   counted. Wide and long sequence data ignore it.
#' @param scaling Rescaling applied to the weights after estimation, in the
#'   order given: `"minmax"` maps the non-zero weights into (0, 1],
#'   `"max"` divides by the largest absolute weight, `"rank"` replaces the
#'   non-zero weights by their ranks, and `"normalize"` makes each row's
#'   absolute weights sum to 1. Zero weights stay zero under every option.
#'   Combine them as in `c("rank", "minmax")`. Default `NULL` (none).
#' @param threshold Weights whose absolute value is below this are set to
#'   zero. Default `0` (none). For `method = "mgm"`, `"LW"` or `"none"`
#'   instead selects mgm's coefficient-thresholding rule (a `threshold` in
#'   `params` takes precedence); any other non-numeric value raises an
#'   error.
#' @param level Multilevel decomposition for the undirected association
#'   methods (`cor`, `pcor`, `glasso`, `mgm`): `"between"` estimates the
#'   network of person means, `"within"` the network of deviations from them,
#'   and `"both"` returns the two as a `netobject_ml`. It needs an id column,
#'   given as `actor` or as `params$id`; a directed method raises an error.
#'   Default `NULL` (no decomposition). For `method = "mgm"`, a numeric
#'   vector instead gives one category count per modelled variable (`1` for
#'   continuous ones), after group, actor and metadata columns are set
#'   aside; it is checked and kept for compatibility, and the fit reads the
#'   categories from the data.
#' @param time_threshold The longest gap, in seconds, between two events of
#'   the same sequence in an event log; a longer gap starts a new sequence.
#'   `FALSE` turns the split off, so each actor (or actor-session) is one
#'   sequence. Default `900`.
#' @param timezone The Olson time zone for timestamps without an offset
#'   (passed to [prepare()]); timestamps with an offset, such as `...Z` or
#'   `+02:00`, are converted from it. Default `"UTC"`.
#' @param predictability Whether to compute each node's predictability (the
#'   R-squared of the node on its neighbours) for the undirected association
#'   methods `cor`, `pcor` and `glasso`. It is stored on the network,
#'   and `cograph::splot()` draws it as a ring around each node. Default
#'   `TRUE`.
#' @param state_cols The columns of a wide table that are sequence positions;
#'   every other column of the cleaned input goes to the network's metadata. Declared
#'   columns are positions even when they hold only 0 and 1; one-hot data is
#'   declared with `codes` or `format = "onehot"`. With nothing declared, the
#'   frame must be one column or hold position columns named as a numbered
#'   run (a shared prefix and a number, `T1`, `T2`, ..., or `1`, `2`, ...);
#'   numeric columns beside the run and character or factor columns outside
#'   it are then covariates. Any other frame on a sequence method, and
#'   columns of only 0 and 1, raise an error of class `nestimate_bad_input`:
#'   nothing says whether they are positions or one-hot codes. Default
#'   `NULL`.
#' @param metadata_cols Columns to keep as the network's metadata (covariates), never as
#'   positions; they cannot overlap `state_cols`. A column left out of a
#'   numbered run of positions is named in a message of class
#'   `nestimate_column_role`, and listing it here or in `state_cols` settles
#'   it. Default `NULL`.
#' @param start A start state put before every sequence, so the network shows
#'   where sequences begin (an edge from it to each sequence's first state).
#'   `FALSE` (default) adds none, `TRUE` adds one labelled `"Start"`, and a
#'   single string adds one with that label. Only for the transition methods
#'   `relative`, `frequency`, `co_occurrence`, `attention`, `ngram`, `gap`
#'   and `reverse`; any other method, `wtna` included, raises an error.
#' @param end An end state put in the cell after each sequence's last
#'   observed state, so the network shows where sequences stop. It has no
#'   outgoing edges and no self-loop, unlike [mark_terminal_state()], which
#'   fills all trailing missing cells with an absorbing state. `FALSE`
#'   (default) adds none, `TRUE` adds one labelled `"End"`, and a single
#'   string adds one with that label. Same methods as `start`.
#' @param params Settings passed to the estimator, as a named list, such as
#'   `list(gamma = 0.5)` for `glasso` or `list(n_gram = 3)` for `ngram`. The
#'   list is stored on the network, so [bootstrap_network()] and the other
#'   resampling verbs re-estimate with exactly the same settings. The
#'   transition methods take tna's sequence options `weighted` and `concat`,
#'   and the low-level `begin_state` / `end_state` that `start` and `end`
#'   set. `window_size` (with `window_span`, or `windowed = TRUE` for a size
#'   of 2) counts `relative`, `frequency` and `co_occurrence` between
#'   consecutive non-overlapping blocks of that many positions, as
#'   `tna::tna()` does; other transition methods raise an error of class
#'   `nestimate_bad_input`, and `mode = "overlapping"` a warning of class
#'   `nestimate_ignored_mode`. Column roles given here (`action`, `id`,
#'   `id_col`, `actor`, `time`, `session`, `order`, `cols`, `codes`,
#'   `group`) must name existing columns, and must agree with the same role
#'   given directly.
#' @param labels New display labels for the nodes, applied after estimation:
#'   a two-column data frame (`name`, `label`), a named character vector
#'   `c(name = "label")` or a named list. It rewrites the node labels and the
#'   weight matrix's dimnames; nodes not listed keep their names.
#' @param ... Further arguments for the estimator. For event logs,
#'   [prepare()]'s time options `unix_time_unit`, `is_unix_time` and
#'   `custom_format` go to [prepare()] instead (e.g.
#'   `unix_time_unit = "milliseconds"`). For the `print()` and `summary()`
#'   methods: ignored.
#' @return A `netobject` (class `c("netobject", "cograph_network")`). Print it
#'   for an overview; [extract_edges()] gives one row per edge, [extract_transition_matrix()]
#'   and [extract_initial_probs()] the weights and initial-state
#'   probabilities, and `cograph::splot()` draws it. [summary()] returns the
#'   descriptive metrics, one row per metric (a `netobject_group` gives one
#'   column per group). The `print()` methods return their input invisibly.
#'   The `netobject`'s fields:
#' \describe{
#'   \item{data}{The state columns of the cleaned input data, as a data frame.}
#'   \item{metadata}{Data frame of the non-state columns of the cleaned input
#'     (and, for long input, the per-sequence metadata), or NULL.}
#'   \item{weights}{The estimated network weight matrix.}
#'   \item{nodes}{Data frame with columns \code{id}, \code{label}, \code{name},
#'     \code{x}, \code{y}.}
#'   \item{edges}{Data frame of non-zero edges with integer \code{from}/\code{to}
#'     (node IDs) and numeric \code{weight}.}
#'   \item{directed}{Logical. Whether the network is directed.}
#'   \item{method}{The resolved method name.}
#'   \item{params}{The params list used (for reproducibility).}
#'   \item{scaling}{The scaling applied (or NULL).}
#'   \item{threshold}{The threshold applied.}
#'   \item{n_nodes}{Number of nodes.}
#'   \item{n_edges}{Number of non-zero edges.}
#'   \item{level}{Decomposition level used, or NULL; for `mgm`, the category
#'     count of each variable.}
#'   \item{build_args}{The resolved column/format arguments (\code{actor},
#'     \code{action}, \code{time}, \code{session}, \code{order}, \code{codes},
#'     \code{format}, \code{window_size}, \code{mode}) used for this build.}
#'   \item{meta}{List with \code{source}, \code{layout}, and \code{tna} metadata
#'     (cograph-compatible).}
#'   \item{node_groups}{Node groupings data frame, or NULL.}
#'   \item{predictability}{Named numeric vector of R-squared predictability
#'     values per node, for `cor`, `pcor` and `glasso` when
#'     `predictability = TRUE`; NULL otherwise.}
#' }
#' The estimator's own results are kept beside these fields, depending on
#' the method: for example `frequency_matrix` and `initial` for `relative`,
#' `frequency`, `ngram`, `gap` and `reverse` (only `initial` for
#' `attention`, neither for `co_occurrence`), `cor_matrix` and the sample
#' size `n` for `cor`, `pcor` and `glasso`, `precision_matrix` for `pcor`
#' and `glasso`, and `lambda_selected` for `glasso`.
#'
#' With `level = "both"`, the result is a `netobject_ml` holding the
#' between-person and within-person networks (`between` and `within`, each a
#' `netobject`) and the `method`. `level = "between"` or `"within"` returns
#' a single `netobject` estimated on the decomposed data.
#'
#' With `group`, or when `data` is a `net_clustering` or `net_mmm` object,
#' the result is a `netobject_group`: a named list of `netobject`s, one per
#' group, that records the grouping column in its `group_col` attribute.
#'
#' @section Estimators:
#' \describe{
#'   \item{\code{"relative"}}{Transition probabilities (row-normalised counts); directed. Aliases: \code{"tna"}, \code{"transition"}.}
#'   \item{\code{"frequency"}}{Transition counts; directed. Aliases: \code{"ftna"}, \code{"counts"}.}
#'   \item{\code{"co_occurrence"}}{Co-occurrence within sequences; undirected. Aliases: \code{"cna"}, \code{"wcna"}, \code{"co-occurrence"}.}
#'   \item{\code{"attention"}}{Decay-weighted transitions; directed. Alias: \code{"atna"}.}
#'   \item{\code{"ngram"}}{Adjacent pairs counted once per n-gram window; directed. Aliases: \code{"n-gram"}, \code{"n_gram"}.}
#'   \item{\code{"gap"}}{Pairs up to \code{max_gap + 1} positions apart, weighted \code{1 / distance}; directed.}
#'   \item{\code{"reverse"}}{Reply network, the transpose of \code{"frequency"}; directed.}
#'   \item{\code{"wtna"}}{Window-based transitions on one-hot data; directed. Alias: \code{"wtna_transition"}.}
#'   \item{\code{"wtna_cooccurrence"}}{Window-based co-occurrence on one-hot data; undirected.}
#'   \item{\code{"cor"}}{Pearson correlations; undirected. Aliases: \code{"corr"}, \code{"correlation"}.}
#'   \item{\code{"pcor"}}{Partial correlations, unregularised; undirected. Alias: \code{"partial"}.}
#'   \item{\code{"glasso"}}{EBIC graphical lasso partial correlations; undirected. Aliases: \code{"ebicglasso"}, \code{"regularized"}.}
#'   \item{\code{"ising"}}{Ising model (L1-penalised logistic regression); undirected. Alias: \code{"isingfit"}.}
#'   \item{\code{"mgm"}}{Mixed graphical model; undirected. Aliases: \code{"mixed"}, \code{"mixed_graphical"}.}
#' }
#'
#' \code{"ngram"}, \code{"gap"} and \code{"reverse"} are \code{tna::build_model()}'s
#' types \code{"n-gram"}, \code{"gap"} and \code{"reverse"}. Their settings go in
#' \code{params}: \code{list(n_gram = 2)}, \code{list(max_gap = 1)} and
#' \code{list(weighted = FALSE)} are the defaults. All three return raw weights;
#' add \code{scaling = "normalize"} for row probabilities.
#'
#' @details
#' `build_network()` resolves the method alias, checks the declared column
#' arguments before the format is resolved, and takes the estimator from the
#' registry. For an association method with `level`, it first decomposes
#' the data into person means (`"between"`) or deviations from them
#' (`"within"`). It then calls the estimator with `data` and `params`,
#' applies `scaling` and then `threshold` to the weight matrix, and extracts
#' the edges.
#'
#' In an event log, `action` alone makes all rows one sequence in row order;
#' `session`, or a `time` gap above `time_threshold`, splits them without an
#' `actor`. When the log yields a single sequence, a message of class
#' `nestimate_single_sequence` says so, because such a network cannot be
#' validated by bootstrap or the other confirmatory tests.
#'
#' For `method = "mgm"`, missing continuous values are filled with the
#' observed column mean and missing categorical values with the modal level.
#' This keeps the rows and is \pkg{psychnets}' documented
#' `na_method = "pairwise"` policy; it is not pairwise deletion.
#'
#' For a network with non-negative weights, the `summary()` values equal
#' those of `tna::summary.tna()` from tna 1.5.4 on.
#'
#' @examples
#' seqs <- data.frame(V1 = c("A","B","C","A"), V2 = c("B","C","A","B"))
#' net <- build_network(seqs, method = "relative")
#' net
#' \donttest{
#' # Transition network (relative probabilities)
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
#'   V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' print(net)
#'
#' # Association network (glasso)
#' freq_data <- convert_sequence_format(seqs, format = "frequency")
#' net_glasso <- build_network(freq_data, method = "glasso",
#'                              params = list(gamma = 0.5, nlambda = 50))
#'
#' # With scaling
#' net_scaled <- build_network(seqs, method = "relative",
#'                              scaling = c("rank", "minmax"))
#' }
#'
#' @references
#' Saqr, M., L\enc{ó}{o}pez-Pernas, S., T\enc{ö}{o}rm\enc{ä}{a}nen, T.,
#' Kaliisa, R., Misiejuk, K., & Tikka, S. (2025). Transition network
#' analysis: A novel framework for modeling, visualizing, and identifying the
#' temporal patterns of learners and learning processes. In
#' \emph{Proceedings of the 15th International Learning Analytics and
#' Knowledge Conference (LAK25)}. ACM. \doi{10.1145/3706468.3706513}
#'
#' Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse
#' covariance estimation with the graphical lasso. \emph{Biostatistics},
#' 9(3), 432--441. \doi{10.1093/biostatistics/kxm045}
#'
#' Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria
#' for Gaussian graphical models. \emph{Advances in Neural Information
#' Processing Systems}, 23, 604--612.
#'
#' van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F., Boschloo,
#' L., Schoevers, R. A., & Waldorp, L. J. (2014). A new method for
#' constructing networks from binary data. \emph{Scientific Reports}, 4,
#' 5918. \doi{10.1038/srep05918}
#'
#' Haslbeck, J. M. B., & Waldorp, L. J. (2020). mgm: Estimating time-varying
#' mixed graphical models in high-dimensional data. \emph{Journal of
#' Statistical Software}, 93(8), 1--46. \doi{10.18637/jss.v093.i08}
#'
#' @seealso \code{\link{register_estimator}}, \code{\link{list_estimators}},
#'   \code{\link{bootstrap_network}}
#'
#' @importFrom stats aggregate ave cor complete.cases var
#' @importFrom utils capture.output
#' @export
build_network <- function(data,
                          method,
                          actor = NULL,
                          action = NULL,
                          time = NULL,
                          session = NULL,
                          order = NULL,
                          codes = NULL,
                          group = NULL,
                          format = "auto",
                          window_size = 3L,
                          mode = c("non-overlapping", "overlapping"),
                          scaling = NULL,
                          threshold = 0,
                          level = NULL,
                          time_threshold = 900,
                          timezone = "UTC",
                          predictability = TRUE,
                          state_cols = NULL,
                          metadata_cols = NULL,
                          start = FALSE,
                          end = FALSE,
                          params = list(),
                          labels = NULL,
                          ...) {
  # --- Coerce sequence matrices (character / logical / integer) to data.frame ---
  # Numeric square matrices are left alone - they may already be transition
  # matrices for downstream estimators that accept them.
  if (is.matrix(data) && !is.numeric(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  # A compare_mmm() table is a data frame of fit statistics, not sequences;
  # read as wide data its log-likelihoods and BICs became states.
  if (inherits(data, "mmm_compare")) {
    k_best <- data[["k"]][which.min(data[["BIC"]])]
    stop(errorCondition(sprintf(
      "Build from one fit: build_network(build_mmm(data, k = %d)).",
      k_best), class = "nestimate_bad_input", call = NULL))
  }

  # The network-shaping formals a caller set explicitly; a fitted clustering
  # rebuilds its per-cluster networks with them (they are formals, so they
  # never reach `...`, and were silently dropped).
  shaping <- intersect(c("scaling", "threshold", "start", "end", "params"),
                       names(match.call()))
  shaping_args <- mget(shaping, envir = environment())

  # --- Early dispatch for net_clustering objects ---
  if (inherits(data, "net_clustering")) {
    if (missing(method)) method <- data$network_method %||% "relative"
    return(do.call(.build_network_clustering,
                   c(list(data, method = method), shaping_args, list(...))))
  }

  # --- Early dispatch for net_mmm objects ---
  if (inherits(data, "net_mmm")) {
    if (missing(method)) method <- data$network_method %||% "relative"
    resolved <- .resolve_method_alias(method)
    dots <- c(shaping_args, list(...))
    if (resolved != "relative") {
      # Re-build per-component networks from hard assignments using requested method
      raw_data    <- data$data %||% data$models[[1L]]$data
      assignments <- data$assignments
      k_comp      <- data$k
      # Merge stored build_args with caller's ...; caller takes precedence
      call_args <- if (!is.null(data$build_args)) modifyList(data$build_args, dots) else dots
      empty <- setdiff(seq_len(k_comp), assignments)
      if (length(empty) > 0L) {
        stop(errorCondition(sprintf(paste(
          "build_network(): no sequence is assigned to cluster %s, so its",
          "\"%s\" network cannot be built. Refit the mixture with fewer",
          "clusters or more starts."), paste(empty, collapse = ", "), method),
          class = "nestimate_empty_cluster", call = NULL))
      }
      nets <- lapply(seq_len(k_comp), function(m) {
        sub <- .keep_window_attrs(raw_data[assignments == m, , drop = FALSE],
                                  raw_data)
        net <- do.call(build_network, c(list(data = sub, method = method), call_args))
        net <- .attach_set_aside_metadata(net, data$metadata, assignments == m)
        # Inject EM-fitted initials only for directed sequence methods
        if (resolved %in% c("relative", "frequency", "attention")) {
          net$initial <- data$models[[m]]$initial
        }
        net
      })
      names(nets) <- paste0("Cluster ", seq_len(k_comp))
      class(nets) <- "netobject_group"
      return(.attach_mmm_clustering(nets, data, full_data = raw_data))
    }
    # params identical to the fitted components' own are not a change;
    # scaling and threshold only reshape the fitted weights, so they apply
    # to the EM networks as build_network() applies them (scaling, then
    # threshold).
    if (identical(dots$params, data$models[[1L]]$params)) dots$params <- NULL
    reshape <- dots[intersect(names(dots), c("scaling", "threshold"))]
    dots <- dots[setdiff(names(dots), c("scaling", "threshold"))]
    if (length(dots)) {
      stop(errorCondition(paste0(
        "Unsupported argument(s) for build_network(net_mmm) with the ",
        "pre-built relative MMM networks: ",
        paste(names(dots), collapse = ", "), ". ",
        "Use a non-relative `method` to rebuild per-cluster networks, or ",
        "set these options before fitting the MMM."),
        class = "nestimate_bad_input", call = NULL))
    }
    # Default: wrap pre-built "relative" models (retain $initial from EM)
    nets <- Map(function(net, m) {
      net <- .reshape_weights(net, reshape$scaling, reshape$threshold %||% 0)
      .attach_set_aside_metadata(net, data$metadata, data$assignments == m)
    }, data$models, seq_along(data$models))
    if (is.null(names(nets))) names(nets) <- paste0("Cluster ", seq_along(nets))
    class(nets) <- "netobject_group"
    return(.attach_mmm_clustering(nets, data,
                                  full_data = data$data %||% data$models[[1L]]$data))
  }

  .assert(is.character(method), length(method) == 1)
  .assert(is.list(params))

  # ---- mgm's `level`/`threshold` collide with build_network's own
  # same-named formals. build_network's `level` is the "between"/"within"/
  # "both" multilevel-decomposition enum; its `threshold` is a numeric
  # edge-weight cutoff. mgm instead documents an integer `level` vector and
  # a "LW"/"none" `threshold` enum. Route mgm's variants into `params`
  # (the estimator's own argument channel) so the documented mgm values are
  # reachable via build_network(method = "mgm", ...) directly, and error
  # clearly for invalid values. Other methods keep build_network's
  # numeric-`threshold` / enum-`level` semantics unchanged.
  resolved_method <- .resolve_method_alias(method)
  if (identical(resolved_method, "mgm")) {
    if (!is.numeric(threshold)) {
      if (length(threshold) == 1L && threshold %in% c("LW", "none")) {
        if (is.null(.param_get(params, "threshold"))) {
          params$threshold <- threshold
        }
        threshold <- 0
      } else {
        stop(errorCondition(paste0(
          "For method = \"mgm\", `threshold` must be \"LW\" or \"none\" ",
          "(coefficient thresholding rule); got: ",
          paste(threshold, collapse = ", "), "."),
          class = "nestimate_bad_input", call = NULL
        ))
      }
    }
    if (!is.null(level) && is.numeric(level)) {
      # Route now; validate only once the final modeled columns are known.
      if (is.null(.param_get(params, "level"))) params$level <- level
      level <- NULL
    }
  }

  .assert(is.numeric(threshold), length(threshold) == 1, threshold >= 0)
  mode <- match.arg(mode)

  if (!is.null(state_cols)) {
    .assert(is.character(state_cols))
  }
  if (!is.null(metadata_cols)) {
    .assert(is.character(metadata_cols))
  }
  if (!is.null(state_cols) && !is.null(metadata_cols)) {
    overlap <- intersect(state_cols, metadata_cols)
    if (length(overlap)) {
      stop(errorCondition(paste0(
        "state_cols and metadata_cols overlap: ",
        paste(overlap, collapse = ", "), ". A column must be either state ",
        "or metadata, not both."),
        class = "nestimate_bad_input", call = NULL
      ))
    }
  }
  if (is.data.frame(data)) {
    missing_state <- setdiff(state_cols,    names(data))
    missing_meta  <- setdiff(metadata_cols, names(data))
    if (length(missing_state)) {
      stop(errorCondition(paste0("state_cols not found in data: ",
           paste(missing_state, collapse = ", ")),
           class = "nestimate_bad_input", call = NULL))
    }
    if (length(missing_meta)) {
      stop(errorCondition(paste0("metadata_cols not found in data: ",
           paste(missing_meta, collapse = ", ")),
           class = "nestimate_bad_input", call = NULL))
    }
  }

  # Merge ... into params (... takes precedence)
  dots <- list(...)
  if (length(dots) > 0) {
    params <- modifyList(params, dots)
  }

  # Resolve method aliases early (needed for format detection)
  method <- .resolve_method_alias(method)

  # ---- start/end boundary markers ----
  # Public surface for the internal begin_state/end_state transition options.
  # TRUE inserts a default-labelled marker; a single string uses that label;
  # FALSE/NULL inserts nothing. `start` prepends one token to every sequence
  # (a pure source); `end` places one token right after each sequence's last
  # observed state (a pure sink). Only the transition estimators honour these,
  # so anything else errors loudly rather than silently ignoring the request.
  .resolve_boundary <- function(x, default_label, arg) {
    if (is.null(x) || isFALSE(x)) return(NULL)
    if (isTRUE(x)) return(default_label)
    if (is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)) return(x)
    stop(errorCondition(paste0(
      "`", arg, "` must be TRUE/FALSE or a single non-empty string label."),
      class = "nestimate_bad_input", call = NULL))
  }
  begin_label <- .resolve_boundary(start, "Start", "start")
  end_label   <- .resolve_boundary(end,   "End",   "end")
  if (!is.null(begin_label) || !is.null(end_label)) {
    boundary_methods <- c("relative", "frequency", "co_occurrence", "attention",
                          "ngram", "gap", "reverse")
    if (!method %in% boundary_methods) {
      stop(errorCondition(paste0(
        "`start`/`end` boundary markers are only supported for the ",
        "transition methods (", paste(boundary_methods, collapse = ", "),
        "); got method = \"", method, "\"."),
        class = "nestimate_bad_input", call = NULL))
    }
    if (!is.null(begin_label)) params$begin_state <- begin_label
    if (!is.null(end_label))   params$end_state   <- end_label
  }

  if (is.data.frame(data)) {
    canonical <- .canonicalize_build_network_params(
      params,
      actor = actor, action = action, time = time, session = session,
      order = order, codes = codes, group = group, format = format
    )
    actor <- canonical$actor
    action <- canonical$action
    time <- canonical$time
    session <- canonical$session
    order <- canonical$order
    codes <- canonical$codes
    group <- canonical$group
    format <- canonical$format

    .validate_build_network_columns(
      data,
      actor = actor, action = action, time = time, session = session,
      order = order, codes = codes, group = group,
      state_cols = state_cols, metadata_cols = metadata_cols
    )
    .validate_build_network_params(data, params, actor = actor,
                                   action = action, time = time,
                                   codes = codes)
  }

  # ---- Column roles are declared, never read from names ----
  # The wtna methods read one-hot codes; declared positions would be ignored
  if (method %in% c("wtna", "wtna_cooccurrence") &&
      (!is.null(state_cols) || !is.null(.param_get(params, "cols")))) {
    stop(errorCondition(
      sprintf(paste0("build_network(): method \"%s\" reads one-hot codes, and ",
                     "`state_cols` declares positions. Use \"relative\" or ",
                     "\"frequency\" for positions, or `codes` to name the one-hot ",
                     "columns."), method),
      class = "nestimate_bad_input", call = NULL))
  }
  if (is.data.frame(data)) {
    .refuse_duplicate_names(data, "build_network()")
    .stop_undeclared_identifier(data, method, action, actor, time, session,
                                order, codes, group, state_cols, metadata_cols,
                                params$cols, format)
    .stop_undeclared_roles(data, method, action, actor, time, session, order,
                           codes, group, state_cols, metadata_cols, params$cols,
                           format)
    .unused_role_notice(data, action, time, session,
                        c(action, actor, time, session, group, order,
                          metadata_cols, state_cols), order = order)
  }

  # ---- Group dispatch: per-group networks ----
  # Resolve wide sequence roles before collecting a shared group alphabet.
  # Child estimation and metadata extraction must use the same columns.
  transition_methods <- c("relative", "frequency", "co_occurrence", "attention",
                          "ngram", "gap", "reverse")
  if (method %in% transition_methods && is.data.frame(data) &&
      is.null(action) && is.null(codes) && !identical(format, "onehot")) {
    roles_declared <- !is.null(params$cols) || !is.null(state_cols) ||
      !is.null(metadata_cols)
    params$cols <- .build_sequence_cols(data, params$cols, state_cols,
                                        metadata_cols,
                                        exclude = c(group, actor, session,
                                          .param_get(params, "id"),
                                          .param_get(params, "id_col")))
    if (!roles_declared && !is.null(params$cols)) {
      .set_aside_notice(data, params$cols,
                        c(group, actor, session, .param_get(params, "id"),
                          .param_get(params, "id_col")))
    }
  }
  if (!is.null(group)) {
    .assert(is.character(group))
    if (method %in% transition_methods && is.data.frame(data) &&
        is.null(.param_get(params, "alphabet"))) {
      if (!is.null(action) && action %in% names(data)) {
        vals <- .clean_states(as.character(data[[action]]))
      } else {
        exclude <- c(group, actor, session, metadata_cols)
        state_names <- if (!is.null(params$cols)) {
          intersect(params$cols, names(data))
        } else {
          setdiff(names(data), exclude)
        }
        vals <- .clean_states(unlist(lapply(data[state_names], as.character),
                                     use.names = FALSE))
      }
      params$alphabet <- sort(unique(c(
        vals[!is.na(vals)],
        .param_get(params, "begin_state"),
        .param_get(params, "end_state")
      )))
    }
    if (length(group) == 1L) {
      grp_key <- data[[group]]
    } else {
      grp_key <- interaction(data[, group, drop = FALSE], sep = "-",
                             drop = TRUE)
    }
    if (anyNA(grp_key)) {
      stop(errorCondition(sprintf(paste0(
        "`group` has %d row(s) with a missing value. Remove them or give ",
        "them a group label."), sum(is.na(grp_key))),
        class = "nestimate_bad_input", call = NULL))
    }
    grp_levels <- unique(grp_key)
    # Drop group column(s) from sub-data so they don't become state cols
    drop_cols <- intersect(group, names(data))
    # Strip group from params so the recursive per-group call does not
    # try to re-canonicalize a `group` column that has been dropped from
    # `sub`. The parent call already consumed it.
    child_params <- if (is.list(params)) params[setdiff(names(params), "group")] else params
    nets <- lapply(grp_levels, function(g) {
      sub <- data[grp_key == g, , drop = FALSE]
      if (length(drop_cols) > 0L) {
        sub <- sub[, setdiff(names(sub), drop_cols), drop = FALSE]
      }
      build_network(
        sub, method = method, actor = actor, action = action,
        time = time, session = session, order = order, codes = codes,
        group = NULL, format = format, window_size = window_size,
        mode = mode, scaling = scaling, threshold = threshold,
        level = level, time_threshold = time_threshold,
        timezone = timezone, predictability = predictability,
        state_cols = state_cols, metadata_cols = metadata_cols,
        labels = labels, params = child_params, ...
      )
    })
    names(nets) <- as.character(grp_levels)
    attr(nets, "group_col") <- group
    class(nets) <- "netobject_group"
    return(nets)
  }

  # ---- Auto-detect input format ----
  is_onehot <- FALSE
  onehot_codes <- NULL
  if (format == "auto" && is.data.frame(data)) {
    if (!is.null(codes)) {
      # Explicit codes = one-hot
      is_onehot <- TRUE
      format <- "onehot"
    } else if (!is.null(action) && action %in% names(data)) {
      # Has action column = long event log
      format <- "long"
    } else {
      # One-hot is never read from the values (owner, 2026-11-07): only
      # `method = "wtna"` declares it, by its name. Every other method reads
      # 0/1 columns as positions (sequence methods, which refuse them unless
      # declared) or as variables (association methods).
      # Declared roles are not codes: the same exclusion `.resolve_codes()` uses.
      exclude <- c(actor, session, metadata_cols, group)
      check_cols <- setdiff(names(data), exclude)
      declared <- !is.null(state_cols) || !is.null(.param_get(params, "cols")) ||
        !is.null(.marked_positions(data))
      if (!declared && method %in% c("wtna", "wtna_cooccurrence") &&
          length(check_cols) >= 2L && all(vapply(
        data[, check_cols, drop = FALSE],
        function(x) (is.numeric(x) || is.logical(x)) &&
          all(x[!is.na(x)] %in% c(0, 1)),
        logical(1)
      ))) {
        is_onehot <- TRUE
        format <- "onehot"
      } else {
        format <- "wide"
      }
    }
  }

  # ---- Long format: prepare event log data ----
  if (format == "long" && !is.null(action) && is.data.frame(data) &&
      action %in% names(data)) {
    prep_args <- list(data = data, action = action)
    if (!is.null(actor)) prep_args$actor <- actor
    if (!is.null(time)) prep_args$time <- time
    if (!is.null(session)) prep_args$session <- session # nocov start
    if (!is.null(order)) prep_args$order <- order # nocov end
    prep_args$time_threshold <- time_threshold
    prep_args$timezone <- timezone
    # prepare()'s time-parsing options arrive through `...`; they belong to
    # prepare(), not to the estimator (which would ignore them silently).
    time_opts <- intersect(c("custom_format", "is_unix_time", "unix_time_unit"),
                           names(params))
    prep_args[time_opts] <- params[time_opts]
    params[time_opts] <- NULL

    prepared <- do.call(prepare, prep_args)
    data <- prepared$sequence_data
    if (nrow(data) == 1L) .single_sequence_note()
    format <- "wide"
    params$format <- "wide"
  } else {
    prepared <- NULL
  }

  # ---- One-hot format: dispatch to wtna/cna path ----
  if (format == "onehot" || is_onehot) {
    # These count positions of states, and would read 0, 1 and the ids as states
    if (method %in% c("attention", "ngram", "gap", "reverse")) {
      stop(errorCondition(
        sprintf(paste0("build_network(): method \"%s\" reads positions of states and ",
                       "has no reading of one-hot data. Use \"relative\", ",
                       "\"frequency\", \"co_occurrence\" or \"wtna\"."), method),
        class = "nestimate_bad_input", call = NULL))
    }
    resolved_codes <- .resolve_codes(
      as.data.frame(data), codes,
      exclude = c(actor, session, metadata_cols, group)
    )

    # Warning for one-hot without windowing or session
    if (window_size <= 1L && is.null(session)) {
      warning(
        "One-hot data without windowing or sessions produces sparse networks. ",
        "Consider setting window_size > 1 or providing session.",
        call. = FALSE
      )
    }

    # Build grouping from actor + session
    grp_col <- NULL
    if (!is.null(actor) || !is.null(session)) {
      grp_parts <- c(actor, session)
      grp_col <- grp_parts
    }

    params$codes <- resolved_codes
    onehot_codes <- resolved_codes
    params$window_size <- window_size
    params$mode <- mode
    if (!is.null(grp_col)) params$actor <- grp_col
    # Let the estimator handle it
  } else {
    # Wide sequence data: pass format through
    if (!"format" %in% names(params)) params$format <- format
    # Explicit column roles decide what the estimator reads. Without this a
    # declared metadata column was counted as a sequence position (its values
    # became states) and was only moved to $metadata after estimation.
    if (is.data.frame(data) &&
        (method %in% transition_methods || !is.null(state_cols) ||
         !is.null(metadata_cols))) {
      params$cols <- .build_sequence_cols(data, params$cols, state_cols,
                                          metadata_cols,
                                          exclude = c(actor, session,
                                            .param_get(params, "id"),
                                            .param_get(params, "id_col")),
                                          infer_numeric_metadata =
                                            method %in% transition_methods)
    }
    if (is.data.frame(data) && method %in% transition_methods) {
      windowing <- .sequence_window(data, params, method, mode)
      data <- windowing$data
      params <- windowing$params
    }
  }

  # Declared metadata columns are not variables of an association network
  # (cor, pcor, glasso, ising, mgm); without this they entered as nodes.
  # Dropped before the sequence check below, so a declared character
  # covariate cannot turn numeric data into "sequences".
  if (method %in% c("cor", "pcor", "glasso", "ising", "mgm") &&
      is.data.frame(data) && !is.null(metadata_cols)) {
    data <- data[, setdiff(names(data), metadata_cols), drop = FALSE]
  }
  # ---- Auto-convert sequences to frequencies for association methods ----
  # Association methods (cor, pcor, glasso) need numeric data. Categorical
  # columns are sequence positions when there are at least two of them and
  # fewer than two numeric columns: they are converted to per-row state
  # frequencies and a numeric covariate is left out. Numeric variables with
  # a stray character column (gender) are not sequences: the estimator uses
  # the numeric columns.
  assoc_methods <- c("cor", "pcor", "glasso")
  if (method %in% assoc_methods && is.data.frame(data)) {
    exclude <- intersect(c(actor, session), names(data))
    check_cols <- setdiff(names(data), exclude)
    is_cat <- vapply(data[check_cols],
                     function(x) is.character(x) || is.factor(x), logical(1))
    is_num <- vapply(data[check_cols], is.numeric, logical(1))
    if (sum(is_cat) >= 2L && sum(is_num) < 2L) {
      seq_data <- data[c(exclude, check_cols[is_cat])]
      freq_data <- convert_sequence_format(
        seq_data, seq_cols = check_cols[is_cat],
        id_col = if (length(exclude) == 0L) character(0) else exclude,
        format = "frequency"
      )
      # Drop rid and ID columns, keep only numeric frequency columns
      drop <- c("rid", exclude)
      data <- freq_data[, setdiff(names(freq_data), drop), drop = FALSE]
      params$format <- "wide"
    }
  }

  # `actor` names the ID column in wide and association data as well as in
  # long data: hand it to the estimator under the name it takes (`id` for
  # sequence estimators, `id_col` for association estimators), otherwise the
  # ID values become states or the ID column a variable.
  estimator <- get_estimator(method)
  # (All-binary association data is detected as one-hot above; that only
  # matters to sequence estimators, so `id_col` is passed regardless.)
  if (!is.null(actor) && is.null(prepared) && is.data.frame(data) &&
      all(actor %in% names(data))) {
    est_args <- names(formals(estimator$fn))
    one_hot <- is_onehot || identical(format, "onehot")
    if ("id_col" %in% est_args && is.null(.param_get(params, "id_col"))) {
      params$id_col <- actor
    } else if (!one_hot && "id" %in% est_args &&
               is.null(.param_get(params, "id"))) {
      params$id <- actor
    }
  }

  # ---- Multilevel decomposition ----
  id_col <- .param_get(params, "id") %||% .param_get(params, "id_col") %||%
    actor

  # Validate level parameter
  if (!is.null(level)) {
    level <- match.arg(level, c("between", "within", "both"))
    if (is.null(id_col)) {
      stop(errorCondition("'id' is required when 'level' is specified.",
                          class = "nestimate_bad_input", call = NULL))
    }
    if (!is.data.frame(data)) {
      stop(errorCondition(
        "'data' must be a data frame when 'level' is specified.",
        class = "nestimate_bad_input", call = NULL))
    }
  }

  .check_scaling(scaling)

  # Multilevel decomposition is only defined for undirected association
  # methods. For directed transition estimators the decomposition is a
  # numeric no-op, so accepting `level` silently produced a pooled network
  # falsely labelled "[between-person]" (and a fake netobject_ml with
  # byte-identical between/within). Error clearly instead of mislabelling.
  if (!is.null(level) && isTRUE(estimator$directed)) {
    stop(errorCondition(paste0(
      "'level' (multilevel decomposition) is only supported for undirected ",
      "association methods (cor, pcor, glasso, mgm). Method '", method,
      "' is a directed estimator and has no between-/within-person ",
      "decomposition."),
      class = "nestimate_bad_input", call = NULL
    ))
  }

  # level = "both": recursive dispatch
  if (identical(level, "both")) {
    between <- build_network(
      data, method = method, actor = actor, params = params,
      scaling = scaling, threshold = threshold, level = "between",
      predictability = predictability, labels = labels
    )
    within_net <- build_network(
      data, method = method, actor = actor, params = params,
      scaling = scaling, threshold = threshold, level = "within",
      predictability = predictability, labels = labels
    )
    result <- list(between = between, within = within_net, method = method)
    class(result) <- "netobject_ml"
    return(result)
  }

  # Multilevel decomposition for association methods
  if (!is.null(level) && !estimator$directed) {
    data <- .decompose_multilevel(data, id_col = id_col, level = level)
  }

  if (identical(method, "mgm") && is.data.frame(data)) {
    # Keep identifiers until any requested multilevel decomposition has
    # used them. MGM itself has no id_col formal and must not model them.
    identifiers <- c(actor, .param_get(params, "id"), .param_get(params, "id_col"))
    data <- data[, setdiff(names(data), identifiers), drop = FALSE]
    mgm_level <- .param_get(params, "level")
    if (!is.null(mgm_level) &&
        (!is.numeric(mgm_level) || length(mgm_level) != ncol(data))) {
      stop(errorCondition(paste0(
        "For method = \"mgm\", numeric `level` must be an integer vector ",
        "of length ", ncol(data), " (one level per modeled variable, 1 ",
        "for continuous; group, actor and metadata columns are excluded)."),
        class = "nestimate_bad_input", call = NULL))
    }
  }

  # Call estimator
  est_result <- do.call(estimator$fn, c(list(data = data), params))

  # Validate estimator output
  if (!is.list(est_result) ||
      is.null(est_result$matrix) ||
      is.null(est_result$nodes) ||
      is.null(est_result$directed)) {
    stop(errorCondition(paste0(
      "Estimator '", method,
      "' must return a list with 'matrix', 'nodes', and 'directed'."),
      class = "nestimate_bad_input", call = NULL))
  }

  net_matrix <- est_result$matrix
  nodes <- est_result$nodes
  directed <- est_result$directed

  # Apply scaling
  if (!is.null(scaling)) {
    net_matrix <- .apply_scaling(
      net_matrix, scaling, include_zeros = .scaling_includes_zeros(method)
    )
  }

  # Apply threshold
  if (threshold > 0) {
    net_matrix[abs(net_matrix) < threshold] <- 0
  }

  # Extract edges
  edges <- .extract_edges_from_matrix(net_matrix, directed = directed)

  # Split cleaned_data into state-only $data and $metadata
  # A column is a state column if all its non-void/non-NA values are in nodes
  raw_data <- est_result$cleaned_data
  metadata <- NULL
  # Use prepared$meta_data when available (long format path)
  if (!is.null(prepared) && !is.null(prepared$meta_data)) {
    md <- prepared$meta_data
    # Drop internal .session_id column
    md <- md[, setdiff(names(md), ".session_id"), drop = FALSE]
    if (ncol(md) > 0L) metadata <- md
  }
  if (is.data.frame(raw_data)) {
    # If the user passed explicit overrides that name columns present in
    # raw_data, those take priority over the values-in-nodes heuristic.
    auto_onehot_state <- if (is.null(state_cols)) {
      intersect(onehot_codes, names(raw_data))
    } else {
      character(0)
    }
    auto_wide_state <- if (is.null(state_cols) &&
                           method %in% transition_methods &&
                           !identical(format, "onehot")) {
      intersect(params$cols, names(raw_data))
    } else {
      character(0)
    }
    auto_variable_state <- if (is.null(state_cols) &&
                               !method %in% transition_methods &&
                               all(nodes %in% names(raw_data)) &&
                               any(vapply(raw_data[intersect(nodes, names(raw_data))],
                                 function(col) {
                                   vals <- .clean_states(as.character(col))
                                   any(!is.na(vals) & !vals %in% nodes)
                                 }, logical(1)))) {
      intersect(nodes, names(raw_data))
    } else {
      character(0)
    }
    # A window with nothing in it still holds its place in the window blocks.
    auto_window_state <- if (is.null(state_cols) &&
                             isTRUE(attr(raw_data, "windowed"))) {
      grep("^W[0-9]+_T[0-9]+$", names(raw_data), value = TRUE)
    } else {
      character(0)
    }
    user_state <- union(intersect(state_cols, names(raw_data)),
                        union(auto_variable_state,
                              union(auto_window_state,
                                    union(auto_onehot_state, auto_wide_state))))
    user_meta  <- intersect(metadata_cols, names(raw_data))

    if (length(user_state) > 0L) {
      # Explicit state_cols: those ARE state; rest of raw_data goes to metadata.
      resolved_state <- user_state
      resolved_extra <- setdiff(names(raw_data), user_state)
    } else {
      is_state_col <- vapply(raw_data, function(col) {
        vals <- .clean_states(as.character(col))
        vals <- vals[!is.na(vals)]
        length(vals) > 0L && all(vals %in% nodes)
      }, logical(1))
      # Honor metadata_cols by forcing those FALSE in the auto-detection.
      if (length(user_meta) > 0L) is_state_col[user_meta] <- FALSE
      resolved_state <- names(raw_data)[is_state_col]
      resolved_extra <- names(raw_data)[!is_state_col]
    }

    if (length(resolved_extra) > 0L) {
      extra_df <- raw_data[, resolved_extra, drop = FALSE]
      metadata <- if (is.null(metadata)) extra_df else cbind(metadata, extra_df)
      raw_data <- .keep_window_attrs(raw_data[, resolved_state, drop = FALSE],
                                     raw_data)
    }
    # Clean void/missing markers in character/factor state columns
    if (is.data.frame(raw_data)) {
      raw_data[] <- lapply(raw_data, function(col) {
        if (is.character(col) || is.factor(col)) {
          .clean_states(as.character(col))
        } else {
          col
        }
      })
    }
  }

  # Build unified netobject / cograph_network
  nodes_df <- data.frame(
    id = seq_along(nodes),
    label = nodes,
    name = nodes,
    stringsAsFactors = FALSE
  )

  result <- list(
    # the network holds its roles, so a verb that rebuilds from it needs no guess
    data = if (is.data.frame(raw_data) && method %in% transition_methods &&
               !identical(format, "onehot")) .mark_positions(raw_data) else raw_data,
    metadata = metadata,
    weights = net_matrix,
    nodes = nodes_df,
    edges = edges,
    directed = directed,
    method = method,
    params = params,
    scaling = scaling,
    threshold = threshold,
    n_nodes = length(nodes),
    n_edges = nrow(edges),
    level = level,
    build_args = list(
      actor = actor, action = action, time = time, session = session,
      order = order, codes = codes, format = format,
      window_size = window_size, mode = mode
    ),
    meta = list(
      source = "nestimate",
      layout = NULL,
      tna = list(method = method)
    ),
    node_groups = NULL
  )

  # Carry over method-specific extras (exclude cleaned_data - identical to $data)
  extras <- setdiff(names(est_result),
                    c("matrix", "nodes", "directed", "cleaned_data"))
  for (key in extras) {
    result[[key]] <- est_result[[key]]
  }

  # Auto-compute predictability (R^2) for undirected association methods.
  # Stored as a named numeric vector for backward compatibility with downstream
  # consumers (e.g. cograph::splot pie ring). Users calling predictability()
  # directly get the full tidy data.frame with R2 + RMSE.
  if (isTRUE(predictability) && !directed) {
    class(result) <- c("netobject", "cograph_network")
    pred_df <- tryCatch(predictability(result), error = function(e) NULL)
    if (!is.null(pred_df)) {
      vec <- pred_df$R2
      names(vec) <- pred_df$node
      result$predictability <- vec
    }
  }

  result <- structure(result, class = c("netobject", "cograph_network"))
  if (!is.null(labels)) result <- .apply_node_labels(result, labels)
  result
}


#' @noRd
.canonicalize_build_network_params <- function(params,
                                               actor = NULL,
                                               action = NULL,
                                               time = NULL,
                                               session = NULL,
                                               order = NULL,
                                               codes = NULL,
                                               group = NULL,
                                               format = "auto") {
  .take_param <- function(current, param_name, formal_name = param_name) {
    value <- .param_get(params, param_name)
    if (is.null(value)) {
      return(current)
    }
    if (is.null(current)) {
      return(value)
    }
    if (!identical(as.character(current), as.character(value))) {
      stop(errorCondition(paste0(
        "'", formal_name, "' and 'params$", param_name,
        "' specify different columns: ",
        paste(as.character(current), collapse = ", "), " vs ",
        paste(as.character(value), collapse = ", ")),
        class = "nestimate_bad_input", call = NULL
      ))
    }
    current
  }

  params_format <- .param_get(params, "format")
  if (!is.null(params_format)) {
    if (!is.character(params_format) || length(params_format) != 1L) {
      stop(errorCondition("'params$format' must be a single character value.",
                          class = "nestimate_bad_input", call = NULL))
    }
    params_format <- match.arg(params_format, c("auto", "wide", "long", "onehot"))
    if (!identical(format, "auto") && !identical(format, params_format)) {
      stop(errorCondition(paste0(
        "'format' and 'params$format' specify different values: ",
        format, " vs ", params_format),
        class = "nestimate_bad_input", call = NULL
      ))
    }
    if (identical(format, "auto")) {
      format <- params_format
    }
  }

  list(
    actor = .take_param(.take_param(actor, "id", formal_name = "actor"),
                        "actor"),
    action = .take_param(action, "action"),
    time = .take_param(time, "time"),
    session = .take_param(session, "session"),
    order = .take_param(order, "order"),
    codes = .take_param(codes, "codes"),
    group = .take_param(group, "group"),
    format = format
  )
}


#' Window a sequence network from `params`
#'
#' `params = list(window_size =, window_span =, windowed =)` count transitions
#' between consecutive blocks of positions, as tna() does. The window travels as
#' the data attributes .prepare_onehot() sets, so the estimators and the
#' resampling verbs read one definition; windows the data already records win.
#' @noRd
.sequence_window <- function(data, params, method, mode = NULL) {
  size <- .param_get(params, "window_size")
  span <- .param_get(params, "window_span")
  if (is.null(size) && is.null(span) &&
      !isTRUE(.param_get(params, "windowed"))) {
    return(list(data = data, params = params))
  }
  if (!method %in% c("relative", "frequency", "co_occurrence")) {
    stop(errorCondition(
      sprintf(paste0("`params$window_size` applies to the methods relative, ",
                     "frequency and co_occurrence, not to \"%s\"."), method),
      class = "nestimate_bad_input", call = NULL))
  }
  is_positive_whole <- function(x) {
    is.null(x) || (is.numeric(x) && length(x) == 1L && is.finite(x) &&
                     x >= 1 && x == round(x))
  }
  if (!is_positive_whole(size) || !is_positive_whole(span)) {
    stop(errorCondition(
      "`params$window_size` and `params$window_span` must be positive whole numbers.",
      class = "nestimate_bad_input", call = NULL))
  }
  # Blocks are the only window sequence positions have; `mode` is one-hot only.
  if (identical(.param_get(params, "mode"), "overlapping") ||
      identical(mode, "overlapping")) {
    warning(warningCondition(
      paste0("`mode = \"overlapping\"` applies to one-hot data; sequence ",
             "positions are counted in non-overlapping blocks."),
      class = "nestimate_ignored_mode", call = NULL))
  }
  params$mode <- NULL
  if (!isTRUE(attr(data, "windowed"))) {
    attr(data, "windowed") <- TRUE
    attr(data, "window_size") <- as.integer(size %||% 2L)
    attr(data, "window_span") <- as.integer(span %||% 1L)
  }
  list(data = data, params = params)
}


#' @noRd
.param_get <- function(params, name, default = NULL) {
  if (!is.list(params) || !name %in% names(params)) {
    return(default)
  }
  params[[name]]
}


#' @noRd
.validate_build_network_columns <- function(data, ...) {
  args <- list(...)
  cols <- names(data)

  missing_by_arg <- lapply(names(args), function(arg) {
    value <- args[[arg]]
    if (is.null(value)) {
      return(character(0))
    }
    if (!is.character(value)) {
      stop(errorCondition(paste0(
        "'", arg, "' must be a character vector of column name(s)."),
        class = "nestimate_bad_input", call = NULL))
    }
    setdiff(value, cols)
  })
  names(missing_by_arg) <- names(args)

  missing_by_arg <- missing_by_arg[lengths(missing_by_arg) > 0L]
  if (length(missing_by_arg) == 0L) {
    return(invisible(TRUE))
  }

  first_arg <- names(missing_by_arg)[1L]
  missing_cols <- missing_by_arg[[1L]]
  available <- paste(cols, collapse = ", ")
  stop(errorCondition(paste0(
    "'", first_arg, "' column",
    if (length(missing_cols) == 1L) " " else "s ",
    "not found in data: ", paste(missing_cols, collapse = ", "),
    ". Available columns: ", available),
    class = "nestimate_bad_input", call = NULL
  ))
}


#' @noRd
.validate_build_network_params <- function(data, params,
                                           actor = NULL,
                                           action = NULL,
                                           time = NULL,
                                           codes = NULL) {
  column_params <- intersect(
    names(params),
    c("action", "id", "time", "cols", "codes", "actor", "id_col",
      "session", "order", "group")
  )
  missing_by_arg <- lapply(column_params, function(arg) {
    value <- params[[arg]]
    if (is.null(value)) {
      return(character(0))
    }
    if (!is.character(value)) {
      stop(errorCondition(paste0(
        "'params$", arg, "' must be a character vector of column name(s)."),
        class = "nestimate_bad_input", call = NULL))
    }
    setdiff(value, names(data))
  })
  names(missing_by_arg) <- column_params
  missing_by_arg <- missing_by_arg[lengths(missing_by_arg) > 0L]
  if (length(missing_by_arg) > 0L) {
    first_arg <- names(missing_by_arg)[1L]
    missing_cols <- missing_by_arg[[1L]]
    stop(errorCondition(paste0(
      "'params$", first_arg, "' column",
      if (length(missing_cols) == 1L) " " else "s ",
      "not found in data: ", paste(missing_cols, collapse = ", "),
      ". Available columns: ", paste(names(data), collapse = ", ")),
      class = "nestimate_bad_input", call = NULL
    ))
  }

  .check_param_conflict <- function(formal_value, param_name,
                                    formal_name = param_name) {
    if (is.null(formal_value) || is.null(params[[param_name]])) {
      return(invisible(TRUE))
    }
    if (!identical(as.character(formal_value), as.character(params[[param_name]]))) {
      stop(errorCondition(paste0(
        "'", formal_name, "' and 'params$", param_name,
        "' specify different columns: ",
        paste(as.character(formal_value), collapse = ", "), " vs ",
        paste(as.character(params[[param_name]]), collapse = ", ")),
        class = "nestimate_bad_input", call = NULL
      ))
    }
    invisible(TRUE)
  }

  .check_param_conflict(action, "action")
  .check_param_conflict(time, "time")
  .check_param_conflict(codes, "codes")
  .check_param_conflict(actor, "actor")
  .check_param_conflict(actor, "id", formal_name = "actor")
  .check_param_conflict(actor, "id_col", formal_name = "actor")
}


# ---- S3 methods ----

#' @rdname build_network
#' @param x For the `print()` methods: a `netobject`, `netobject_group` or
#'   `netobject_ml`, or the result of `summary()` on one.
#' @param digits For `print()` of a `netobject_group`: decimal places of the
#'   weight range. Default `3`.
#' @param object For `summary()`: a `netobject` or `netobject_group`.
#' @param combined For `summary()` of a `netobject_group`: `TRUE` (default)
#'   joins the groups' metrics into one data frame with a column per group;
#'   `FALSE` returns a named list of the per-group tables.
#' @export
print.netobject <- function(x, ...) {
  method_labels <- c(
    relative      = "Transition Network (relative probabilities)",
    frequency     = "Transition Network (frequency counts)",
    co_occurrence = "Co-occurrence Network",
    glasso        = "Partial Correlation Network (EBICglasso)",
    pcor          = "Partial Correlation Network (unregularised)",
    cor           = "Correlation Network",
    ising         = "Ising Model Network",
    attention     = "Attention Network (decay-weighted transitions)",
    wtna          = "Window TNA (transitions)"
  )

  label <- if (x$method %in% names(method_labels)) {
    method_labels[[x$method]]
  } else {
    sprintf("Network (method: %s)", x$method)
  }

  dir_label   <- if (x$directed) " [directed]" else " [undirected]"
  level_label <- if (!is.null(x$level)) sprintf(" [%s-person]", x$level) else ""
  cat(label, dir_label, level_label, "\n", sep = "")

  if (!is.null(x$n)) cat(sprintf("  Sample size: %d\n", x$n))

  # ---- Weight summary (one line) ----
  mat <- x$weights
  if (x$directed) {
    nz <- mat[mat != 0]
  } else {
    nz <- mat[mat != 0 & row(mat) <= col(mat)]
  }
  if (length(nz) > 0) {
    is_assoc <- x$method %in% c("cor", "pcor", "glasso", "ising")
    if (is_assoc) {
      cat(sprintf("  Weights: [%.3f, %.3f]  |  +%d / -%d edges\n",
                  min(nz), max(nz), sum(nz > 0), sum(nz < 0)))
    } else {
      cat(sprintf("  Weights: [%.3f, %.3f]  |  mean: %.3f\n",
                  min(nz), max(nz), mean(nz)))
    }
  }

  # ---- Weight matrix ----
  cat("\n  Weight matrix:\n")
  digits <- if (all(nz == floor(nz))) 0L else 3L
  mat_r  <- round(mat, digits)
  labels <- x$nodes$label
  dimnames(mat_r) <- list(labels, labels)
  formatted <- capture.output(print(mat_r))
  cat(paste0("  ", formatted, collapse = "\n"), "\n")

  # ---- Initial probabilities ----
  if (!is.null(x$initial) && length(x$initial) > 0) {
    cat("\n  Initial probabilities:\n")
    init  <- x$initial
    ord   <- order(init, decreasing = TRUE)
    bar_w <- 40L
    max_v <- max(init)
    vapply(ord, function(i) {
      bars <- if (max_v > 0) strrep("\u2588", round(init[i] / max_v * bar_w)) else ""
      cat(sprintf("  %-12s  %.3f  %s\n", names(init)[i], init[i], bars))
      invisible("")
    }, character(1L))
  }

  # ---- Predictability (R\u00b2) ----
  if (!is.null(x$predictability) && length(x$predictability) > 0) {
    cat("\n  Predictability (R\u00b2):\n")
    pred <- x$predictability
    ord  <- order(pred, decreasing = TRUE)
    bar_w <- 40L
    max_v <- max(pred)
    vapply(ord, function(i) {
      bars <- if (max_v > 0) strrep("\u2588", round(pred[i] / max_v * bar_w)) else ""
      cat(sprintf("  %-12s  %.3f  %s\n", names(pred)[i], pred[i], bars))
      invisible("")
    }, character(1L))
  }

  # ---- Method-specific params ----
  if (x$method == "glasso" && !is.null(x$gamma)) {
    cat(sprintf("\n  Gamma: %.2f  |  Lambda: %.4f\n", x$gamma, x$lambda_selected))
  }
  if (x$method == "ising") {
    cat(sprintf("\n  Gamma: %.2f  |  Rule: %s\n", x$gamma, x$rule))
    if (!is.null(x$thresholds)) {
      thr_rng <- range(x$thresholds)
      cat(sprintf("  Thresholds: [%.3f, %.3f]\n", thr_rng[1], thr_rng[2]))
    }
  }
  if (!is.null(x$scaling))  cat(sprintf("\n  Scaling: %s\n",   paste(x$scaling, collapse = " -> ")))
  if (x$threshold > 0)      cat(sprintf("  Threshold: %g\n",  x$threshold))

  invisible(x)
}


#' @rdname build_network
#' @section Methods:
#' Printing a `netobject` shows the kind of network, the range of its
#' weights (with the numbers of positive and negative edges for `cor`,
#' `pcor`, `glasso` and `ising`, the mean weight otherwise), the weight
#' matrix under the node labels, the initial probabilities and the
#' predictability values when present, the sample size for `cor`, `pcor`,
#' `glasso` and `ising`, the tuning values of `glasso` and `ising`, and any
#' scaling or threshold. Printing a `netobject_group` gives a header naming
#' its source (a clustering from `build_network()` on a `net_clustering` or
#' `net_mmm` fit, or [cluster_network()], or a plain split by the grouping
#' column) and one row per group with its node and edge counts and weight
#' range; when the groups are clusters, it also shows the number and share
#' of sequences in each, in the layout of the [build_clusters()] and
#' [build_mmm()] prints. A `netobject_ml` prints the node, edge and sample
#' counts of its between- and within-person networks.
#'
#' `summary()` gives the descriptive metrics `tna::summary.tna()` reports:
#' node and edge counts, density, mean shortest-path distance, the mean and
#' SD of out- and in-strength, the mean and SD of out-degree, out- and
#' in-degree centralisation (Freeman) and reciprocity. A self-loop counts as
#' an edge, an undirected edge is counted once, and density is capped at 1.
#' Strengths sum absolute weights. Distances take each edge's absolute
#' weight as its length, not its inverse, and pairs with no path between
#' them are left out of the mean. Degree centralisation ignores self-loops,
#' and reciprocity is the share of off-diagonal edges whose reverse is also
#' present (1 for an undirected network). Mean distance and reciprocity are
#' `NaN` when no pair qualifies. Printing a summary shows the counts as
#' integers and the other metrics to four significant digits.
#' @export
print.netobject_group <- function(x, digits = 3L, ...) {
  digits <- as.integer(digits)
  grps   <- names(x)
  if (is.null(grps) || any(grps == "")) {
    grps <- paste0("Group ", seq_along(x))
  }
  k <- length(grps)
  cl <- attr(x, "clustering")

  # ---- Header --------------------------------------------------------
  if (inherits(cl, "net_mmm_clustering")) {
    cat(sprintf("Group Networks (%d clusters from MMM)\n", k))
  } else if (inherits(cl, "net_clustering")) {
    cat(sprintf("Group Networks (%d clusters via %s / %s)\n",
                k, cl$method, cl$dissimilarity))
  } else {
    group_col <- attr(x, "group_col")
    if (!is.null(group_col)) {
      cat(sprintf("Group Networks (%d groups, group_col: %s)\n",
                  k, group_col))
    } else {
      cat(sprintf("Group Networks (%d groups)\n", k))
    }
  }

  # ---- Per-group table ----------------------------------------------
  nodes_v <- integer(k)
  edges_v <- integer(k)
  wmin_v  <- numeric(k)
  wmax_v  <- numeric(k)
  for (i in seq_len(k)) {
    net <- x[[i]]
    nodes_v[i] <- as.integer(net$n_nodes %||% 0L)
    edges_v[i] <- as.integer(net$n_edges %||% 0L)
    mat <- net$weights
    if (!is.null(mat)) {
      nz <- if (isTRUE(net$directed)) mat[mat != 0] else
        mat[mat != 0 & row(mat) <= col(mat)]
      if (length(nz) > 0L) {
        wmin_v[i] <- min(nz); wmax_v[i] <- max(nz)
      } else {
        wmin_v[i] <- NA_real_; wmax_v[i] <- NA_real_
      }
    } else {
      wmin_v[i] <- NA_real_; wmax_v[i] <- NA_real_
    }
  }

  weights_str <- ifelse(
    is.na(wmin_v),
    "--",
    sprintf(paste0("[%.", digits, "f, %.", digits, "f]"), wmin_v, wmax_v)
  )

  cat("\n")
  cols <- list(
    Group   = grps,
    Nodes   = sprintf("%d", nodes_v),
    Edges   = sprintf("%d", edges_v),
    Weights = weights_str
  )

  # If a clustering attribute provides assignments, append N + %.
  if (!is.null(cl) && !is.null(cl$assignments)) {
    sizes <- as.integer(tabulate(cl$assignments, nbins = k))
    n_total <- sum(sizes)
    if (n_total > 0L) {
      cols[["N"]] <- .fmt_size_pct(sizes, n_total)
    }
  }

  cat(paste(.cluster_table_lines(cols), collapse = "\n"), "\n", sep = "")

  invisible(x)
}


#' @rdname build_network
#' @export
print.netobject_ml <- function(x, ...) {
  cat(sprintf("Multilevel Network (method: %s)\n", x$method))
  cat("-- Between-person --\n")
  b <- x$between
  cat(sprintf("  Nodes: %d  |  Edges: %d\n", b$n_nodes, b$n_edges))
  if (!is.null(b$n)) {
    cat(sprintf("  Sample size: %d (unique persons)\n", b$n))
  }
  cat("-- Within-person --\n")
  w <- x$within
  cat(sprintf("  Nodes: %d  |  Edges: %d\n", w$n_nodes, w$n_edges))
  if (!is.null(w$n)) {
    cat(sprintf("  Sample size: %d (observations)\n", w$n))
  }
  invisible(x)
}


# ---- Predictability ----

#' Compute Node Predictability
#'
#' @description
#' `predictability()` computes how well each node of an association network
#' is predicted by the others: the proportion of its variance explained
#' (R\eqn{^2}), with the root mean squared error (RMSE) of that prediction,
#' following Haslbeck & Waldorp (2018). It reads networks estimated with
#' `method = "glasso"`, `"pcor"` or `"cor"`, single, multilevel
#' (`netobject_ml`) or grouped (`netobject_group`). [build_network()]
#' stores the same R\eqn{^2} values on these networks when
#' `predictability = TRUE` (the default), and `cograph::splot()` draws them
#' as a ring around each node.
#'
#' R\eqn{^2} says how much of a node's variance its connections in the
#' network account for; a low value leaves most of it to variables outside
#' the network or to noise. It is computed from the matrices estimated on
#' the sample, so it describes that sample. For `"glasso"` it reflects the
#' regularised estimate, and for `"cor"` it depends on which edges survive
#' `threshold`.
#'
#' @details
#' For `"glasso"` and `"pcor"`, R\eqn{^2} comes from the estimated precision
#' matrix \eqn{\Omega} (the inverse of the correlation matrix, regularised
#' for `"glasso"`):
#' \deqn{R^2_j = 1 - 1 / \Omega_{jj}}
#' For `"cor"`, it is the multiple R\eqn{^2} of each node regressed on its
#' neighbours (the nodes it shares a non-zero edge with), computed from the
#' stored correlation matrix: 0 for a node without neighbours, the squared
#' correlation for a node with one. Values are clipped to \[0, 1\].
#'
#' RMSE predicts each standardised variable from the others (for `"cor"`,
#' from its neighbours) with the regression weights implied by the same
#' matrix, on the complete rows of the data, and reports the error on the
#' variable's original scale. For `"cor"`, a node without neighbours has an
#' RMSE equal to its standard deviation (prediction by the mean), and a
#' singular correlation block of its neighbours gives an R\eqn{^2} of 0 and
#' an RMSE of `NA`.
#'
#' @param object The network: a `netobject` estimated with `method =
#'   "glasso"`, `"pcor"` or `"cor"`, or a `netobject_ml` or `netobject_group`
#'   of them. A network from any other method raises an error.
#' @param data The variables the network was estimated from, used only for
#'   the `RMSE` column; R\eqn{^2} comes from the matrices stored on the
#'   network. Its columns are matched by node label or, failing that, by
#'   node name. Default `NULL`, which uses the data stored on the network;
#'   without either, or when the columns do not match, `RMSE` is `NA`. The
#'   `netobject_ml` and `netobject_group` methods take no `data` (one passed
#'   to them is ignored) and use each network's own.
#' @param ... Not used.
#'
#' @return A data frame with one row per node and the columns `node` (the
#'   node label), `R2` (between 0 and 1) and `RMSE` (`NA` when no data is
#'   available). A `netobject_ml` gives a list of two such data frames,
#'   `between` and `within`, and a `netobject_group` a named list of them,
#'   one per group.
#'
#' @references
#' Haslbeck, J. M. B., & Waldorp, L. J. (2018). How well do network models
#' predict observations? On the importance of predictability in network models.
#' \emph{Behavior Research Methods}, 50(2), 853--861.
#' \doi{10.3758/s13428-017-0910-x}
#'
#' @examples
#' set.seed(42)
#' mat <- matrix(rnorm(60), ncol = 4)
#' colnames(mat) <- LETTERS[1:4]
#' net <- build_network(as.data.frame(mat), method = "glasso")
#' predictability(net)
#'
#' @export
predictability <- function(object, ...) {
  UseMethod("predictability")
}


#' @rdname predictability
#' @export
predictability.netobject <- function(object, data = NULL, ...) {
  labels <- object$nodes$label
  p <- length(labels)

  # ---- R^2 ----
  if (!is.null(object$precision_matrix)) {
    omega_diag <- diag(object$precision_matrix)
    r2 <- 1 - 1 / omega_diag
  } else if (!is.null(object$cor_matrix)) {
    S   <- object$cor_matrix
    net <- object$weights
    r2 <- vapply(seq_len(p), function(j) {
      neighbors <- which(net[j, ] != 0)
      if (length(neighbors) == 0L) return(0)
      if (length(neighbors) == 1L) return(S[neighbors, j]^2)
      r_vec <- S[neighbors, j]
      R_nn  <- S[neighbors, neighbors]
      tryCatch(as.numeric(crossprod(r_vec, solve(R_nn, r_vec))),
               error = function(e) 0)
    }, numeric(1))
  } else {
    stop(errorCondition(paste0(
      "predictability() requires a precision or correlation matrix ",
      "(methods: glasso, pcor, cor). Method '", object$method,
      "' does not support predictability."),
      class = "nestimate_bad_input", call = NULL))
  }
  r2 <- pmin(pmax(r2, 0), 1)

  # ---- RMSE ----
  if (is.null(data)) data <- object$data
  rmse <- rep(NA_real_, p)

  if (is.data.frame(data) || is.matrix(data)) {
    # Columns carry the node labels (relabelled $data) or, for data the
    # caller passes, the original node names.
    cols <- if (all(labels %in% colnames(data))) labels else object$nodes$name
    X <- if (all(cols %in% colnames(data))) {
      as.matrix(data[, cols, drop = FALSE])
    } else NULL
    if (!is.null(X) && is.numeric(X)) {
      X  <- X[complete.cases(X), , drop = FALSE]
      mu <- colMeans(X)
      sd <- apply(X, 2, stats::sd)
      Z  <- scale(X, center = mu, scale = sd)

      if (!is.null(object$precision_matrix)) {
        Omega <- object$precision_matrix
        rmse <- vapply(seq_len(p), function(j) {
          beta_z <- -Omega[-j, j] / Omega[j, j]
          z_hat  <- Z[, -j, drop = FALSE] %*% beta_z
          y_hat  <- mu[j] + sd[j] * z_hat
          sqrt(mean((X[, j] - y_hat)^2))
        }, numeric(1))
      } else {
        S   <- object$cor_matrix
        net <- object$weights
        rmse <- vapply(seq_len(p), function(j) {
          nbrs <- which(net[j, ] != 0)
          if (length(nbrs) == 0L) return(sd[j])
          beta_z <- tryCatch(solve(S[nbrs, nbrs, drop = FALSE], S[nbrs, j]),
                             error = function(e) NULL)
          if (is.null(beta_z)) return(NA_real_)
          z_hat <- Z[, nbrs, drop = FALSE] %*% beta_z
          y_hat <- mu[j] + sd[j] * z_hat
          sqrt(mean((X[, j] - y_hat)^2))
        }, numeric(1))
      }
    }
  }

  data.frame(node = labels, R2 = r2, RMSE = rmse,
             stringsAsFactors = FALSE, row.names = NULL)
}


#' @rdname predictability
#' @export
predictability.netobject_ml <- function(object, ...) {
  list(
    between = predictability(object$between),
    within  = predictability(object$within)
  )
}


#' @rdname predictability
#' @export
predictability.netobject_group <- function(object, ...) {
  lapply(object, predictability)
}

#' Resolve the wide sequence columns before estimation and group dispatch
#'
#' @param data A wide data frame.
#' @param cols,state_cols Explicit estimator or public sequence selectors.
#' @param metadata_cols Columns declared to be metadata.
#' @param exclude Identifier and grouping columns.
#' @param infer_numeric_metadata Whether to infer covariates: numeric columns
#'   of mixed categorical/numeric input, and categorical columns whose values
#'   lie outside the shared state alphabet (see
#'   \code{.categorical_covariates()}).
#' @return A character selector, or NULL for the estimator's default.
#' @noRd
.build_sequence_cols <- function(data, cols = NULL, state_cols = NULL,
                                  metadata_cols = NULL, exclude = NULL,
                                  infer_numeric_metadata = TRUE) {
  stopifnot(is.data.frame(data))
  if (!is.null(cols)) return(cols)
  if (!is.null(state_cols)) return(state_cols)
  candidates <- setdiff(names(data), c(metadata_cols, exclude))
  if (!is.null(metadata_cols)) return(candidates)
  # positions a verb already resolved keep their roles, whatever they are called;
  # a role declared in this call comes first
  if (!is.null(.marked_positions(data))) return(candidates)
  if (isTRUE(infer_numeric_metadata)) {
    is_num <- vapply(data[candidates], is.numeric, logical(1))
    is_cat <- vapply(data[candidates],
                     function(x) is.character(x) || is.factor(x), logical(1))
    all_binary <- all(vapply(data[candidates], .is_binary_column, logical(1)))
    if (any(is_num) && !all_binary) {
      # numeric-coded states: when the largest run of numbered position
      # names is entirely numeric, it is the sequence and the other columns
      # (an id, a grade, a session label) are covariates, not the other way round
      numeric_run <- .numeric_position_run(data[candidates])
      if (length(numeric_run) > 0L) return(numeric_run)
    }
    inferred <- if (any(is_num) && any(is_cat)) candidates[!is_num] else
      candidates
    covariates <- .categorical_covariates(data[inferred])
    if (length(covariates) > 0L || !identical(inferred, candidates)) {
      return(setdiff(inferred, covariates))
    }
  }
  # Returning NULL would let the estimator include the excluded columns
  # again when its own default selector does not know the session/group role.
  if (any(exclude %in% names(data))) return(candidates)
  NULL
}

# Numeric or logical, and only 0, 1 or missing.
#' @noRd
.is_binary_column <- function(x) {
  (is.numeric(x) || is.logical(x)) && all(x[!is.na(x)] %in% c(0, 1))
}

# Names give no role, so a frame a sequence method would read with nothing
# declared has to be one it can read safely. Association methods read numeric
# columns and pass; a registered estimator reads what it likes.
#' @noRd
.stop_undeclared_roles <- function(data, method, action, actor, time, session,
                                   order, codes, group, state_cols,
                                   metadata_cols, cols, format) {
  # the wtna methods read one-hot data by definition, so naming one declares it
  sequence_methods <- c("relative", "frequency", "co_occurrence", "attention",
                        "ngram", "gap", "reverse")
  if (!method %in% sequence_methods || !is.null(action) ||
      !is.null(state_cols) || !is.null(cols) || !is.null(codes) ||
      !identical(format, "auto")) {
    return(invisible(NULL))
  }
  .refuse_undeclared_roles(
    data, "build_network()",
    taken = c(actor, time, session, order, group, metadata_cols),
    wide = "`state_cols` or `metadata_cols`",
    rest_declared = !is.null(metadata_cols),
    binary_advice = if (method %in% c("relative", "frequency", "co_occurrence")) {
      paste0("Set `format = \"onehot\"` or `codes` for one-hot data, ",
             "or `state_cols` for positions.")
    } else {
      sprintf(paste0("Set `state_cols` if they are positions. Method \"%s\" has no ",
                     "reading of one-hot data; \"relative\", \"frequency\", ",
                     "\"co_occurrence\" and \"wtna\" have."), method)
    }
  )
}

# Reading an event log as positions turns its ids and times into states, without
# a word, and no value tells a log from a table. So a frame is read as positions
# with nothing declared only when it holds one column or a numbered run of positions
# (T1, T2, ...). Anything else, any column called `action`, and columns of only 0 and
# 1, stop and ask. This refuses; it never assigns a role from a name.
#' @noRd
.refuse_undeclared_roles <- function(data, verb, taken = NULL, wide = "`cols`",
                                     via_network = FALSE, whole_frame = FALSE,
                                     rest_declared = FALSE, takes_action = TRUE,
                                     advice = NULL, binary_advice = NULL,
                                     reads_onehot = FALSE) {
  .refuse_duplicate_names(data, verb)
  if (!is.null(.marked_positions(data))) return(invisible(NULL))
  rest <- setdiff(names(data), taken)
  if (length(rest) == 0L) return(invisible(NULL))
  # Columns of 0 and 1 could be one-hot codes or positions, and nothing in them says
  # which. No reading is guessed.
  if (all(vapply(data[rest], .is_binary_column, logical(1)))) {
    stop(errorCondition(
      sprintf("%s: these columns hold only 0 and 1, so they could be one-hot codes or positions: %s. %s",
              verb, paste(sprintf("`%s`", utils::head(rest, 6L)), collapse = ", "),
              binary_advice %||% if (via_network && reads_onehot) {
                paste0("Declare them with build_network(format = \"onehot\") or ",
                       "build_network(state_cols = ), and pass the network.")
              } else if (via_network) {
                paste0("If they are positions, declare them with ",
                       "build_network(state_cols = ) and pass the network. One-hot ",
                       "codes are not sequences of states, which this verb reads.")
              } else {
                sprintf("Set %s if they are positions.", wide)
              }),
      class = "nestimate_bad_input", call = NULL
    ))
  }
  hit <- rest[tolower(rest) == "action"]
  frame <- data[rest]
  # One column has nothing beside it to mistake for a state. Otherwise a verb that
  # sets covariates aside needs a run among the columns, and one that reads every
  # column needs the frame to be nothing else.
  readable <- rest_declared || length(rest) == 1L || if (whole_frame) {
    .all_one_position_run(rest)
  } else {
    length(.categorical_positions(frame)) > 0L ||
      length(.numeric_position_run(frame)) > 0L
  }
  if (length(hit) == 0L && readable) return(invisible(NULL))
  problem <- if (length(hit) > 0L) {
    sprintf("column names are not read, and `%s` is not declared", hit[1L])
  } else {
    sprintf(paste0("no column role is declared, and these columns are not a ",
                   "numbered run of positions (T1, T2, ...): %s"),
            paste(sprintf("`%s`", utils::head(rest, 6L)), collapse = ", "))
  }
  advice <- advice %||% if (via_network) {
    paste0("Declare them with build_network(action = ) for an event log or ",
           "build_network(state_cols = ) for a wide table, and pass the network.")
  } else if (takes_action) {
    sprintf("Set `action` for an event log, or %s for a wide table.", wide)
  } else {
    sprintf("Set %s to the columns holding the positions.", wide)
  }
  stop(errorCondition(
    sprintf("%s: %s. %s", verb, problem, advice),
    class = "nestimate_bad_input", call = NULL
  ))
}

# A verb that returns positions has had them declared or checked already, so its
# result reads back in without declaring them twice, whatever they are called. The
# mark holds the columns, so a column added afterwards is not trusted. A numbered
# run needs no mark, and stays a plain data frame.
#' @noRd
.mark_positions <- function(data) {
  if (!.all_one_position_run(names(data))) {
    attr(data, "nestimate_positions") <- names(data)
  }
  data
}

# The position columns in their current order, or NULL when the frame no longer
# holds exactly the marked ones.
#' @noRd
.marked_positions <- function(data) {
  marked <- attr(data, "nestimate_positions", exact = TRUE)
  nms <- names(data)
  # setequal() ignores repeats, so a duplicated name would let a column in
  if (is.character(marked) && !anyDuplicated(nms) &&
      length(marked) == length(nms) && setequal(marked, nms)) nms
}

# A numeric or logical matrix of only 0 and 1 could hold one-hot codes or
# positions, and it has no columns to declare; a character matrix holds states.
#' @noRd
.refuse_binary_matrix <- function(x, verb, takes_cols = FALSE) {
  if (is.matrix(x) && (is.numeric(x) || is.logical(x)) &&
      all(x[!is.na(x)] %in% c(0, 1))) {
    stop(errorCondition(
      sprintf(paste0("%s: this matrix holds only 0 and 1, so it could hold one-hot ",
                     "codes or positions. If they are positions, give them as a ",
                     "character matrix%s."), verb,
              if (takes_cols) " or a data frame with `cols`" else ""),
      class = "nestimate_bad_input", call = NULL))
  }
  invisible(NULL)
}

# Two columns with one name cannot be told apart by any role, and different
# readers would keep different ones of them.
#' @noRd
.refuse_duplicate_names <- function(data, verb) {
  if (!is.data.frame(data)) return(invisible(NULL))
  dup <- unique(names(data)[duplicated(names(data))])
  if (length(dup) == 0L) return(invisible(NULL))
  stop(errorCondition(
    sprintf("%s: column names must be unique; repeated: %s.", verb,
            paste(sprintf("`%s`", dup), collapse = ", ")),
    class = "nestimate_bad_input", call = NULL))
}

# A one-hot network stores 0 and 1 per code, and an association network its
# variables, so reading either as positions gives meaningless states. Verbs with
# no reading of them stop and say so, for each network of a group too.
#' @noRd
.refuse_onehot_network <- function(x, verb) {
  if (inherits(x, "netobject_group")) {
    lapply(x, .refuse_onehot_network, verb = verb)
    return(invisible(NULL))
  }
  if (!inherits(x, c("netobject", "cograph_network"))) return(invisible(NULL))
  if (.is_onehot_params(x$params)) {
    stop(errorCondition(
      sprintf(paste0("%s: this is a one-hot network, whose data hold 0 and 1 per code, ",
                     "not positions of states. Build it from sequences."), verb),
      class = "nestimate_bad_input", call = NULL))
  }
  if (isTRUE(x$method %in% .association_methods)) {
    stop(errorCondition(
      sprintf(paste0("%s: this \"%s\" network holds variables, not sequences of ",
                     "states. Build it with a transition method or pass the ",
                     "sequences."), verb, x$method),
      class = "nestimate_bad_input", call = NULL))
  }
  invisible(NULL)
}

# Methods whose data are variables (columns measured on each row), not states.
#' @noRd
.association_methods <- c("cor", "pcor", "glasso", "ising", "mgm")

# Every name numbered, with one prefix: T1, T2, ...; 1, 2, ...; W0_T1, W1_T1, ...
#' @noRd
.all_one_position_run <- function(nms) {
  prefix <- sub("[0-9]+$", "", nms)
  prefix[grepl("^W[0-9]+_T[0-9]+$", nms)] <- "W_T"
  all(grepl("[0-9]+$", nms)) && length(unique(prefix)) == 1L
}

# A role that names a column the data lacks is a mistake, not a request to read
# the frame some other way.
#' @noRd
.check_declared_columns <- function(data, verb, ...) {
  .refuse_duplicate_names(data, verb)
  declared <- list(...)
  declared <- declared[lengths(declared) > 0L]
  missing_cols <- vapply(names(declared), function(role) {
    absent <- setdiff(declared[[role]], names(data))
    if (length(absent) == 0L) "" else
      sprintf("`%s` (%s)", role, paste(absent, collapse = ", "))
  }, character(1))
  missing_cols <- missing_cols[nzchar(missing_cols)]
  if (length(missing_cols) == 0L) return(invisible(NULL))
  stop(errorCondition(
    sprintf("%s: column(s) not found in `data`: %s.", verb,
            paste(missing_cols, collapse = "; ")),
    class = "nestimate_bad_input", call = NULL
  ))
}

# 0/1 columns beside a numeric column that is not a run of positions: an id left
# undeclared. Read as wide, its values would become the network's states.
#' @noRd
.stop_undeclared_identifier <- function(data, method, action, actor, time,
                                        session, order, codes, group,
                                        state_cols, metadata_cols, cols,
                                        format) {
  sequence_methods <- c("relative", "frequency", "co_occurrence", "attention",
                        "ngram", "gap", "reverse", "wtna", "wtna_cooccurrence")
  if (!method %in% sequence_methods || !is.null(action) ||
      !is.null(state_cols) || !is.null(cols) || !is.null(codes) ||
      !identical(format, "auto")) {
    return(invisible(NULL))
  }
  rest <- data[setdiff(names(data),
                       c(actor, time, session, order, group, metadata_cols))]
  if (ncol(rest) == 0L || !all(vapply(rest, is.numeric, logical(1)))) {
    return(invisible(NULL))
  }
  binary <- vapply(rest, .is_binary_column, logical(1))
  if (sum(binary) < 2L || all(binary) ||
      length(.numeric_position_run(rest)) > 0L) {
    return(invisible(NULL))
  }
  stop(errorCondition(
    sprintf(paste0("build_network(): %s beside 0/1 columns would become states. ",
                   "Set `actor` or `session` if it identifies the units, or ",
                   "`codes`, `state_cols` or `metadata_cols`."),
            paste(sprintf("`%s`", names(rest)[!binary]), collapse = ", ")),
    class = "nestimate_bad_input", call = NULL
  ))
}

#' @noRd
.column_role_message <- function(text) {
  cond <- simpleMessage(paste0("build_network(): ", text, "\n"))
  class(cond) <- c("nestimate_column_role", class(cond))
  message(cond)
  invisible(NULL)
}

# In an event log an undeclared time or session column changes the result, so
# the call says it went unused.
#' @noRd
.unused_role_notice <- function(data, action, time, session, declared,
                                order = NULL) {
  if (is.null(action)) return(invisible(NULL))
  lower <- tolower(names(data))
  unused <- c(
    if (is.null(time)) names(data)[lower == "time"],
    if (is.null(session)) names(data)[lower %in% c("session", "session_id")]
  )
  # a column declared under any other role is in use, whatever it is called
  unused <- setdiff(unused, declared)
  if (length(unused) == 0L) return(invisible(NULL))
  ordered <- !is.null(time) || !is.null(order)
  effects <- c(if (!ordered) "events are read in row order",
               if (is.null(session)) "no session column splits them",
               if (ordered && !is.null(session) && is.null(time)) "`order` sets the events' order")
  roles <- c(if (is.null(time)) "`time`", if (is.null(session)) "`session`")
  .column_role_message(sprintf(
    "%s not used: %s. Set %s to use them, or list them in `metadata_cols`.",
    paste(sprintf("`%s`", unused), collapse = ", "),
    paste(effects, collapse = " and "), paste(roles, collapse = " or ")))
}

# A character column could hold states, so leaving one out of the sequence is
# said aloud. A numeric column beside character states could not.
#' @noRd
.set_aside_notice <- function(data, cols, exclude) {
  set_aside <- setdiff(names(data), c(cols, exclude))
  if (length(set_aside) == 0L) return(invisible(NULL))
  .column_role_message(sprintf(
    paste0("setting %s aside as metadata, outside the numbered position ",
           "columns. Name the column in `state_cols` to count it as a state."),
    paste(sprintf("`%s`", set_aside), collapse = ", ")))
}

#' The unique largest run of numbered position names, when it is numeric
#'
#' Same naming rule as .categorical_covariates(): a shared prefix followed
#' by a number. Returns the run's column names when the run has at least two
#' columns, is larger than any other run, and every column in it is
#' numeric or entirely missing, with at least one numeric column; otherwise
#' character(0).
#' @noRd
.numeric_position_run <- function(data) {
  nms <- names(data)
  numbered <- grepl("[0-9]+$", nms)
  if (!any(numbered)) return(character(0))
  prefix <- sub("[0-9]+$", "", nms)
  runs <- table(prefix[numbered])
  if (max(runs) < 2L || sum(runs == max(runs)) > 1L) return(character(0))
  run <- nms[numbered & prefix == names(runs)[which.max(runs)]]
  numeric <- vapply(data[run], is.numeric, logical(1))
  missing <- vapply(data[run], function(x) all(is.na(x)), logical(1))
  if (any(numeric) && all(numeric | missing)) run else character(0)
}

#' Categorical columns that are covariates, not sequence positions
#'
#' Wide sequence positions are named as a run: a shared prefix followed by a
#' number (T1, T2, ...; V1, V2, ...; 1, 2, ...). When the categorical
#' columns contain one such run that is larger than any other, its columns
#' are the positions and every other categorical column (an outcome, a
#' course, a group label) is a covariate. Nothing is reassigned when no run
#' of at least two columns exists or when two runs tie. Names, not values,
#' decide: a position whose states occur nowhere else (a constant first
#' state, a terminal state) looks exactly like a covariate by its values.
#'
#' @param data Candidate state columns of a wide data frame.
#' @return Character vector of covariate column names (possibly empty).
#' @noRd
.categorical_covariates <- function(data) {
  positions <- .categorical_positions(data)
  if (length(positions) == 0L) return(character(0))
  is_cat <- vapply(data, \(x) is.character(x) || is.factor(x) ||
                     all(is.na(x)), logical(1))
  setdiff(names(data)[is_cat], positions)
}

# The unique largest run of numbered categorical columns; character(0) if none.
#' @noRd
.categorical_positions <- function(data) {
  # an all-missing column (logical NA) is an unobserved position, not a
  # covariate, so it must not break the run of position names
  is_cat <- vapply(data, \(x) is.character(x) || is.factor(x) ||
                     all(is.na(x)), logical(1))
  if (sum(is_cat) < 2L) return(character(0))
  nms <- names(data)
  numbered <- grepl("[0-9]+$", nms)
  prefix <- sub("[0-9]+$", "", nms)
  # W0_T1, W1_T1, ... from .prepare_onehot() are one run of positions: windows
  # of different lengths must not read as competing runs.
  prefix[grepl("^W[0-9]+_T[0-9]+$", nms)] <- "W_T"
  runs <- table(prefix[is_cat & numbered])
  if (length(runs) == 0L || max(runs) < 2L ||
      sum(runs == max(runs)) > 1L) {
    return(character(0))
  }
  nms[is_cat & numbered & prefix == names(runs)[which.max(runs)]]
}

# Validate build_network()'s `scaling`.
#' @noRd
.check_scaling <- function(scaling) {
  if (!is.null(scaling)) {
    valid_scaling <- c("minmax", "max", "rank", "normalize")
    bad <- setdiff(scaling, valid_scaling)
    if (length(bad) > 0) {
      stop(errorCondition(paste0(
        "Unknown scaling method(s): ", paste(bad, collapse = ", "),
        ". Options: ", paste(valid_scaling, collapse = ", ")),
        class = "nestimate_bad_input", call = NULL))
    }
  }
  invisible(NULL)
}

# Scale and threshold a fitted network's weights as build_network() does
# (scaling first, then |w| < threshold set to 0) and rebuild its edges.
#' @noRd
.reshape_weights <- function(net, scaling = NULL, threshold = 0) {
  .check_scaling(scaling)
  .assert("`threshold` must be a single number of at least 0" =
              is.numeric(threshold) && length(threshold) == 1L &&
              threshold >= 0)
  if (is.null(scaling) && threshold == 0) return(net)
  W <- net$weights
  if (!is.null(scaling)) {
    W <- .apply_scaling(W, scaling,
                        include_zeros = .scaling_includes_zeros(net$method))
  }
  if (threshold > 0) W[abs(W) < threshold] <- 0
  net$weights <- W
  net$edges <- .extract_edges_from_matrix(W, directed = net$directed)
  net$n_edges <- nrow(net$edges)
  net$scaling <- scaling
  net$threshold <- threshold
  net
}
