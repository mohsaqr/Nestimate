#' Window-based Transition Network Analysis
#'
#' @description
#' `wtna()` builds a network from one-hot data, where each row is a time
#' point and each code column holds 1 when that code is active, by counting
#' codes across windows of consecutive rows. It returns a directed
#' transition network (codes active in one window, followed by codes active
#' in the next), an undirected co-occurrence network (codes active in the
#' same window), or both. With `actor`, each actor's rows are counted on
#' their own and the counts summed, so no window spans two actors. The
#' networks are the ones `build_network(method = "wtna")` and
#' `build_network(method = "wtna_cooccurrence")` estimate.
#'
#' The weights count pairs of active codes, not events: a row with several
#' active codes adds to several edges at once, and a larger window pairs
#' more rows. Compare networks built with the same `window_size` and
#' `mode`. How far the weights move under resampling of actors is checked
#' with [bootstrap_network()], and two networks are compared with
#' [permutation()]; both need `actor`, since rows without it form a single
#' unit.
#'
#' @param data A data frame or matrix whose rows are time points, in order,
#'   holding the one-hot code columns and, optionally, the `actor` column.
#'   A column named in `codes` may hold only 0, 1 or missing; a missing
#'   cell counts as inactive, and any other value raises an error of class
#'   `nestimate_bad_input`. Under the default `codes = NULL` such a column
#'   is not read as a code.
#' @param method The network to build: `"transition"` (default), directed,
#'   from the codes of one window to those of the next; `"cooccurrence"`,
#'   undirected, between the codes of the same window; or `"both"`, which
#'   returns the two networks together as a `wtna_mixed`.
#' @param type The weights: `"frequency"` (default) keeps the counts, and
#'   `"relative"` divides each row by its sum (a row without counts stays
#'   zero). On a co-occurrence network `"relative"` gives an asymmetric
#'   matrix, the co-occurrence conditional on the row's code, so keep
#'   `"frequency"` when symmetric counts are needed.
#' @param codes The one-hot columns, as names, as column positions, or as a
#'   single `"first:last"` string naming a range of adjacent columns.
#'   Default `NULL`, every numeric or logical column holding only 0, 1 and
#'   missing values, except `actor`. At least two codes are needed, and a
#'   name that is not a column stops `wtna()`.
#' @param window_size The number of consecutive rows in a window, a whole
#'   number of at least 1. Default `3`. With `1` there is no windowing:
#'   transitions are counted from each row to the next, and co-occurrence
#'   within each row, and `mode` is ignored.
#' @param mode How windows are formed when `window_size` is above 1.
#'   `"non-overlapping"` (default) cuts the rows into consecutive, separate
#'   blocks of `window_size` rows. `"overlapping"` first replaces the rows
#'   by one row per sliding window of `window_size` rows (step 1), marking
#'   a code active when it is active in any of them, so n rows give
#'   n - `window_size` + 1 windows, and then counts those in blocks of
#'   `window_size`, as `tna::import_onehot(window_type = "sliding")` does.
#' @param actor The column, or columns, naming who each row belongs to.
#'   The counts are computed per actor (per combination of values, for
#'   several columns) and summed. Default `NULL`, all rows as one sequence.
#'   A missing identifier raises an error of class `nestimate_bad_input`,
#'   and a name that is not a column stops `wtna()`. With several columns,
#'   the initial-state distribution groups the rows by the first column
#'   only.
#'
#' @return For `method = "transition"` or `"cooccurrence"`: a `netobject`
#'   (class `c("netobject", "cograph_network")`, see [build_network()]) with
#'   `method` set to `"wtna_transition"` or `"wtna_cooccurrence"`. `print()`
#'   shows the weights and, for transitions, the initial probabilities;
#'   [extract_edges()] gives one row per edge, [extract_initial_probs()]
#'   the initial-state distribution of a transition network,
#'   `cograph::splot()` draws it,
#'   [bootstrap_network()] resamples its actors, and [permutation()]
#'   compares it with a second network by permuting actors. It is directed
#'   only for transitions, holds the input as given in `data`,
#'   and records the windowing settings (`type`, `window_size`, `mode`,
#'   `codes`, `actor`) in `params`. A transition network also stores the
#'   initial-state distribution in `initial`: for each actor (the first
#'   `actor` column when there are several), the codes active in its first
#'   row with an active code share its weight equally, and the shares are
#'   averaged over actors (without `actor`, the first row with an active
#'   code in the whole data).
#'
#'   For `method = "both"`: a `wtna_mixed` object, a list of the
#'   `transition` and `cooccurrence` networks (each a `netobject` as above)
#'   and `method = "wtna_both"`. Its `print()` gives the node and edge
#'   counts of each network and returns its input invisibly;
#'   [bootstrap_network()] and [permutation()] test both networks.
#'
#' @details
#' \strong{Transitions}: with `window_size = 1`, `crossprod(X[-n, ], X[-1, ])`
#' counts how often code i is active in row t and code j in row t + 1. With
#' larger windows, every row of a window is paired with every row of the
#' next window, which is the outer product of the two windows' column sums;
#' a final, shorter window is counted as the successor of the last full
#' one.
#'
#' \strong{Co-occurrence}: with `window_size = 1`, `crossprod(X)` counts the
#' codes active in the same row, so the diagonal is how often each code is
#' active. With larger windows, every row of a window is paired with every
#' row of the same window, itself included, which is the outer product of
#' the window's column sums.
#'
#' \strong{Windows}: in non-overlapping windows a code counts once for each
#' row it is active in. These counts follow tna's windowed algorithm.
#'
#' \strong{Per actor}: with `actor`, the counts are computed per group and
#' summed.
#'
#' @examples
#' oh <- matrix(c(1,0,0, 0,1,0, 0,0,1, 1,0,0), nrow = 4, byrow = TRUE,
#'              dimnames = list(NULL, c("A","B","C")))
#' w <- wtna(oh)
#'
#' \donttest{
#' # Simple one-hot data
#' df <- data.frame(
#'   A = c(1, 0, 1, 0, 1),
#'   B = c(0, 1, 0, 1, 0),
#'   C = c(0, 0, 1, 0, 0)
#' )
#'
#' # Transition network
#' net <- wtna(df)
#' print(net)
#'
#' # Both networks
#' nets <- wtna(df, method = "both")
#' print(nets)   # the transition and the co-occurrence network
#'
#' # With windowing
#' net <- wtna(df, window_size = 2, mode = "non-overlapping")
#' }
#'
#' @seealso \code{\link{build_network}}
#'
#' @export
wtna <- function(data,
                 method = c("transition", "cooccurrence", "both"),
                 type = c("frequency", "relative"),
                 codes = NULL,
                 window_size = 3L,
                 mode = c("non-overlapping", "overlapping"),
                 actor = NULL) {
  method <- match.arg(method)
  type <- match.arg(type)
  mode <- match.arg(mode)
  # Mirror .prepare_onehot()'s window_size guard so an invalid value (non-
  # integer, zero, negative) errors cleanly instead of producing a distinct
  # garbage matrix (2.7) or silently falling to the non-windowed branch (0/-5).
  .assert("`window_size` must be a single whole number of at least 1" =
            length(window_size) == 1L && is.numeric(window_size) &&
            !is.na(window_size) && window_size == as.integer(window_size))
  window_size <- as.integer(window_size)
  .assert("`window_size` must be a single whole number of at least 1" =
            window_size >= 1L)

  df <- as.data.frame(data)
  codes <- .resolve_codes(df, codes, exclude = actor)

  .assert("At least two one-hot code columns are needed" = length(codes) >= 2L)

  if (is.null(actor)) {
    X_raw <- .wtna_code_matrix(df, codes)
    weights <- .wtna_compute_weights(X_raw, method, window_size, mode)
  } else {
    .assert("`actor` must name columns of `data`" = all(actor %in% names(df)))
    weights <- .wtna_compute_by_actor(df, codes, window_size, mode, actor, method)
  }

  # Initial state probs for transition networks (directed only)
  initial <- if (method %in% c("transition", "both"))
    .wtna_initial_probs(df, codes, actor) else NULL

  if (method == "both") {
    result <- list(
      transition   = .wtna_finalize(weights$transition,   type, codes, data, "transition",
                                    initial = initial, window_size = window_size,
                                    mode = mode, actor = actor),
      cooccurrence = .wtna_finalize(weights$cooccurrence, type, codes, data, "cooccurrence",
                                    window_size = window_size, mode = mode, actor = actor),
      method = "wtna_both"
    )
    class(result) <- "wtna_mixed"
    return(result)
  }

  .wtna_finalize(weights, type, codes, data, method, initial = initial,
                 window_size = window_size, mode = mode, actor = actor)
}


# ---- Private helpers ----

#' Prepare one-hot code matrix
#'
#' Missing one-hot indicators are treated as inactive cells. This keeps sparse
#' indicator exports from propagating NA through cross-product counts.
#' @noRd
.wtna_code_matrix <- function(df, codes) {
  bad <- codes[!vapply(df[codes], function(x) all(x[!is.na(x)] %in% c(0, 1)),
                       logical(1))]
  if (length(bad) > 0L) {
    stop(errorCondition(
      sprintf(paste0("one-hot column(s) must hold only 0, 1 or missing: `%s` ",
                     "(found %s). Recode them first, for example `x > 0` for ",
                     "counts or `x == \"Yes\"` for labels."),
              paste(bad, collapse = "`, `"),
              paste(utils::head(setdiff(unique(as.vector(unlist(df[bad]))),
                                        c(0, 1, NA)), 3L), collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }
  X <- as.matrix(df[, codes, drop = FALSE])
  X[is.na(X)] <- 0L
  storage.mode(X) <- "integer"
  X
}


#' Compute transition counts
#'
#' For window_size <= 1: consecutive crossprod (t->t+1).
#' For window_size > 1: pairwise between-window transitions matching tna's
#' windowed algorithm - every position in window_i paired with every position
#' in window_{i+1}.
#' @noRd
.wtna_transitions <- function(X, window_size = 3L, mode = "non-overlapping") {
  n <- nrow(X)
  k <- ncol(X)
  if (n < 2L) return(matrix(0, k, k))

  if (window_size <= 1L) {
    return(crossprod(X[-n, , drop = FALSE], X[-1L, , drop = FALSE]))
  }

  weights <- matrix(0, k, k)
  X <- .wtna_window_rows(X, window_size, mode)
  n <- nrow(X)
  if (n < 2L) return(weights)

  # Match tna's compute_transitions_windowed: pair every position in
  # window_i with every position in window_{i+1}
  divides <- n %% window_size == 0L
  q <- n %/% window_size - 1L * divides
  # Loop kept: a running sum over windows; one k x k accumulator instead of
  # a stack of n_windows matrices.
  for (i in seq_len(q)) {
    j_idx <- seq((i - 1L) * window_size + 1L, i * window_size)
    k_idx <- seq(i * window_size + 1L, min(n, (i + 1L) * window_size))
    # sum_j sum_ki tcrossprod(X[j,], X[ki,]) = tcrossprod(colSums, colSums)
    s_curr <- colSums(X[j_idx, , drop = FALSE])
    s_next <- colSums(X[k_idx, , drop = FALSE])
    weights <- weights + tcrossprod(s_curr, s_next)
  }
  weights
}


#' The rows the windowed counts are taken over
#'
#' Overlapping: one row per sliding window, marking the codes active anywhere in
#' it, as tna::import_onehot(window_type = "sliding") does. Both modes then
#' count those rows in blocks of window_size, so "overlapping" means the same
#' here as in tna.
#' @noRd
.wtna_window_rows <- function(X, window_size, mode) {
  if (mode == "non-overlapping") return(X)
  .onehot_window_presence(X, window_size, mode) * 1
}


#' Compute co-occurrence counts
#'
#' For window_size <= 1: standard crossprod across all rows.
#' For window_size > 1: within-window pairwise co-occurrence matching tna's
#' windowed algorithm - every position in a window paired with every other
#' position in the same window.
#' @noRd
.wtna_cooccurrence <- function(X, window_size = 3L, mode = "non-overlapping") {
  n <- nrow(X)
  k <- ncol(X)

  if (window_size <= 1L) return(crossprod(X))

  weights <- matrix(0, k, k)
  X <- .wtna_window_rows(X, window_size, mode)
  n <- nrow(X)

  # Match tna: pair every position with every other position within
  # the same window (including self-pairs at same position)
  n_windows <- ceiling(n / window_size)
  # Loop kept: a running sum over windows (see .wtna_transitions()).
  for (i in seq_len(n_windows)) {
    idx <- seq((i - 1L) * window_size + 1L, min(n, i * window_size))
    s <- colSums(X[idx, , drop = FALSE])
    weights <- weights + tcrossprod(s)
  }
  weights
}


#' Dispatch weight computation by method
#'
#' @param X_raw Raw one-hot matrix (not collapsed). Transitions use this
#'   directly with pairwise between-window counting. Co-occurrence collapses
#'   windows first via \code{.wtna_to_matrix}.
#' @noRd
.wtna_compute_weights <- function(X_raw, method, window_size = 3L,
                                   mode = "non-overlapping") {
  switch(method,
    transition = .wtna_transitions(X_raw, window_size, mode),
    cooccurrence = .wtna_cooccurrence(X_raw, window_size, mode),
    both = list(
      transition = .wtna_transitions(X_raw, window_size, mode),
      cooccurrence = .wtna_cooccurrence(X_raw, window_size, mode)
    )
  )
}


#' Compute weights per actor and sum
#' @noRd
.wtna_compute_by_actor <- function(df, codes, window_size, mode, actor, method) {
  n_codes <- length(codes)
  init_mat <- matrix(0, n_codes, n_codes)
  # Raises on a missing identifier; split() below would drop those rows unseen.
  .observed_group_id(df[actor], context = "actor/session")
  if (length(actor) == 1L) {
    groups <- split(df, df[[actor]])
  } else {
    grp_key <- interaction(df[, actor, drop = FALSE], drop = TRUE)
    groups <- split(df, grp_key)
  }

  matrices <- lapply(groups, function(g) {
    X_raw <- .wtna_code_matrix(g, codes)
    .wtna_compute_weights(X_raw, method, window_size, mode)
  })

  if (method == "both") {
    list(
      transition = Reduce(`+`, lapply(matrices, `[[`, "transition"), init_mat),
      cooccurrence = Reduce(`+`, lapply(matrices, `[[`, "cooccurrence"), init_mat)
    )
  } else {
    Reduce(`+`, matrices, init_mat)
  }
}


#' Resolve codes specification to column names
#'
#' Supports multiple selection styles:
#' \itemize{
#'   \item \code{NULL} - auto-detect binary 0/1 columns (excluding \code{exclude})
#'   \item Character vector - column names (e.g. \code{c("A", "B", "C")})
#'   \item Numeric vector - column indices (e.g. \code{2:9})
#'   \item Single string with \code{:} - column name range
#'     (e.g. \code{"Planning:Evaluating"})
#' }
#' @param df Data frame.
#' @param codes Column specification (NULL, character, numeric, or range string).
#' @param exclude Character vector of column names to exclude from auto-detection.
#' @return Character vector of resolved column names.
#' @noRd
.resolve_codes <- function(df, codes = NULL, exclude = NULL) {
  if (is.null(codes)) {
    # Auto-detect binary columns, excluding specified columns
    check_df <- df[, setdiff(names(df), exclude), drop = FALSE]
    return(.wtna_auto_detect_codes(check_df))
  }

  if (is.numeric(codes)) {
    # Column indices
    .assert("`codes` positions must lie between 1 and the number of columns" =
              all(codes >= 1L & codes <= ncol(df)))
    return(names(df)[codes])
  }

  if (is.character(codes) && length(codes) == 1L && grepl(":", codes)) {
    # Name range "start:end"
    parts <- strsplit(codes, ":")[[1]]
    .assert("A `codes` range must be written \"first:last\"" =
              length(parts) == 2L)
    start <- match(parts[1], names(df))
    end <- match(parts[2], names(df))
    if (is.na(start)) {
      stop(errorCondition(paste0("Column '", parts[1], "' not found."),
                          class = "nestimate_bad_input", call = NULL))
    }
    if (is.na(end)) {
      stop(errorCondition(paste0("Column '", parts[2], "' not found."),
                          class = "nestimate_bad_input", call = NULL))
    }
    return(names(df)[start:end])
  }

  # Character vector of names
  .assert("`codes` must be column names, positions or a \"first:last\" range" =
            is.character(codes))
  missing <- setdiff(codes, names(df))
  if (length(missing) > 0) {
    stop(errorCondition(paste0("Columns not found: ",
                               paste(missing, collapse = ", ")),
                        class = "nestimate_bad_input", call = NULL))
  }
  codes
}


#' Auto-detect one-hot binary columns
#' @noRd
.wtna_auto_detect_codes <- function(df) {
  is_onehot <- vapply(df, function(x) {
    if (!is.numeric(x) && !is.logical(x)) return(FALSE)
    vals <- unique(x[!is.na(x)])
    length(vals) > 0 && all(vals %in% c(0, 1))
  }, logical(1L))

  codes <- names(df)[is_onehot]
  if (length(codes) == 0L) {
    stop(errorCondition(
      "No one-hot columns found.",
      class = "nestimate_bad_input", call = NULL)) # nocov
  }
  codes
}


#' Compute initial state probabilities for wtna transition networks
#'
#' For each actor (or the whole dataset if no actor), finds the first active
#' row and distributes probability uniformly across all simultaneously active
#' states. Results are averaged across actors to produce a distribution summing
#' to 1. This handles binary data where multiple states can be active at once.
#' @noRd
.wtna_initial_probs <- function(df, codes, actor) {
  X <- .wtna_code_matrix(df, codes)

  .first_row_probs <- function(mat) {
    active_rows <- which(rowSums(mat) > 0L)
    if (length(active_rows) == 0L) return(NULL)
    active_cols <- which(mat[active_rows[1L], ] > 0L)
    if (length(active_cols) == 0L) return(NULL)
    init <- setNames(numeric(length(codes)), codes)
    init[active_cols] <- 1.0 / length(active_cols)
    init
  }

  if (is.null(actor)) {
    return(.first_row_probs(X))
  }

  grp       <- df[[actor[1L]]]
  actor_ids <- unique(grp)

  prob_accum <- setNames(numeric(length(codes)), codes)
  n_valid    <- 0L

  # Loop kept: a running mean over actors that skips actors with no
  # observed first row; the order of addition is part of the result.
  for (id in actor_ids) {
    p <- .first_row_probs(X[grp == id, , drop = FALSE])
    if (is.null(p)) next
    prob_accum <- prob_accum + p
    n_valid    <- n_valid + 1L
  }

  if (n_valid == 0L) return(NULL)
  prob <- prob_accum / n_valid
  prob / sum(prob)   # ensure exact sum = 1 despite floating point
}

#' Finalize: row-normalize and build netobject
#' @noRd
.wtna_finalize <- function(weights, type, codes, data, method, initial = NULL,
                            window_size = 3L, mode = "non-overlapping",
                            actor = NULL) {
  if (type == "relative") {
    rs <- rowSums(weights)
    rs[rs == 0] <- 1
    weights <- weights / rs
  }

  dimnames(weights) <- list(codes, codes)
  directed <- method == "transition"

  # Extract edges
  edges <- .extract_edges_from_matrix(weights, directed = directed)

  nodes_df <- data.frame(
    id = seq_along(codes),
    label = codes,
    name = codes,
    x = NA_real_,
    y = NA_real_,
    stringsAsFactors = FALSE
  )
  wtna_method <- paste0("wtna_", method)

  structure(
    list(
      data = data,
      weights = weights,
      nodes = nodes_df,
      edges = edges,
      directed = directed,
      method = wtna_method,
      params = list(type = type, window_size = window_size, mode = mode,
                    codes = codes, actor = actor),
      scaling = NULL,
      threshold = 0,
      n_nodes = length(codes),
      n_edges = nrow(edges),
      level = NULL,
      initial = initial,
      meta = list(source = "nestimate", layout = NULL,
                  tna = list(method = wtna_method)),
      node_groups = NULL
    ),
    class = c("netobject", "cograph_network")
  )
}


# ---- Registry estimator wrappers ----
# These allow wtna to be used via build_network() and gain
# bootstrap_network() / permutation() support automatically.

#' Core wtna estimator for the registry
#'
#' @param data Data frame with one-hot columns.
#' @param codes Character vector or NULL. One-hot column names.
#' @param window_size Integer (>= 1). Window size. Default: 3.
#' @param mode Character. "non-overlapping" or "overlapping".
#' @param actor Character or NULL. Actor grouping column.
#' @param wtna_method Character. "transition" or "cooccurrence".
#' @param ... Ignored.
#' @return Standard estimator list (matrix, nodes, directed, cleaned_data).
#' @noRd
.estimator_wtna_core <- function(data, codes = NULL, window_size = 3L,
                                  mode = "non-overlapping", actor = NULL,
                                  wtna_method = "transition",
                                  type = "frequency", ...) {
  df <- as.data.frame(data)
  codes <- .resolve_codes(df, codes, exclude = actor)

  .assert("At least two one-hot code columns are needed" = length(codes) >= 2L)
  mode <- match.arg(mode, c("non-overlapping", "overlapping"))
  window_size <- as.integer(window_size)

  if (is.null(actor)) {
    X_raw <- .wtna_code_matrix(df, codes)
    weights <- .wtna_compute_weights(X_raw, wtna_method, window_size, mode)
  } else {
    .assert("`actor` must name columns of `data`" = all(actor %in% names(df)))
    weights <- .wtna_compute_by_actor(df, codes, window_size, mode,
                                       actor, wtna_method)
  }

  if (type == "relative") {
    rs <- rowSums(weights)
    rs[rs == 0] <- 1
    weights <- weights / rs
  }

  dimnames(weights) <- list(codes, codes)
  directed <- wtna_method == "transition"

  list(
    matrix = weights,
    nodes = codes,
    directed = directed,
    cleaned_data = data
  )
}


#' WTNA transition estimator (directed)
#'
#' Accepts both one-hot binary data and sequence data (wide/long).
#' One-hot: windowing + crossprod. Sequence: existing transition counting.
#' @noRd
.estimator_wtna <- function(data,
                             format = "auto",
                             action = NULL,
                             id = NULL,
                             time = NULL,
                             cols = NULL,
                             codes = NULL,
                             window_size = 3L,
                             mode = "non-overlapping",
                             actor = NULL,
                             ...) {
  df <- as.data.frame(data)

  # Before the format: a one-hot column may be called "Action".
  if (identical(format, "onehot") || !is.null(codes)) {
    return(.estimator_wtna_core(
      data, codes = codes, window_size = window_size,
      mode = mode, actor = actor, wtna_method = "transition", ...
    ))
  }

  # Detect format: one-hot vs sequence
  if (format == "auto") {
    format <- .format_from_action(action, df)
  }

  if (format == "wide") {
    state_cols <- .select_state_cols(df, c(id, actor), cols)
    is_binary <- length(state_cols) >= 2L && all(vapply(
      df[, state_cols, drop = FALSE],
      function(x) (is.numeric(x) || is.logical(x)) &&
          all(x[!is.na(x)] %in% c(0, 1)),
      logical(1)
    ))

    if (is_binary || !is.null(codes)) {
      # One-hot path: windowing + crossprod
      return(.estimator_wtna_core(
        data, codes = codes, window_size = window_size,
        mode = mode, actor = actor, wtna_method = "transition", ...
      ))
    }
  }

  # Sequence path: use existing transition counting directly
  .estimator_frequency(
    data, format = format, action = action, id = id,
    time = time, cols = cols, ...
  )
}




#' @rdname wtna
#' @param x For the `print()` method: a `wtna_mixed` object.
#' @param ... For the `print()` method: ignored.
#' @export
print.wtna_mixed <- function(x, ...) {
  cat("Mixed Window TNA (transition + co-occurrence)\n")
  cat("-- Transition (directed) --\n")
  t <- x$transition
  cat(sprintf("  Nodes: %d  |  Edges: %d\n", t$n_nodes, t$n_edges))
  cat("-- Co-occurrence (undirected) --\n")
  co <- x$cooccurrence
  cat(sprintf("  Nodes: %d  |  Edges: %d\n", co$n_nodes, co$n_edges))
  invisible(x)
}
