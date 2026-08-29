# ---- tna-parity sequence estimators: ngram, gap, reverse ----
#
# These mirror `tna::build_model()` types "n-gram", "gap" and "reverse".
# All three operate on a wide sequence matrix (rows = sequences, columns =
# time points) and count state pairs at chosen column offsets. The shared
# kernel `.count_column_pairs()` does the tabulation; each estimator only
# decides which column pairs to count and with what weight.

#' Coerce sequence input (wide or long) to a wide character matrix
#'
#' @return List with `matrix` (character, NA-padded) and `states`.
#' @noRd
.sequence_matrix <- function(data, format = "auto", action = "Action",
                             id = NULL, time = "Time", cols = NULL,
                             alphabet = NULL, concat = 1L,
                             begin_state = NULL, end_state = NULL) {
  stopifnot("`data` must be a data frame" = is.data.frame(data))
  if (identical(format, "auto")) {
    format <- if (action %in% names(data)) "long" else "wide"
  }
  if (identical(format, "long")) {
    id_col <- if (is.null(id)) ".sequence" else id[1L]
    long <- as.data.frame(data)
    if (is.null(id)) long$.sequence <- 1L
    if (length(id) > 1L) {
      long[[".sequence"]] <- .group_label(long[, id, drop = FALSE])
      id_col <- ".sequence"
    }
    time_col <- if (time %in% names(long)) time else NULL
    wide <- long_to_wide(long, id_col = id_col, time_col = time_col,
                         action_col = action, fill_na = TRUE)
    data <- wide
    id <- id_col
    cols <- NULL
  }
  .prepare_sequence_wide(data, id = id, cols = cols, alphabet = alphabet,
                         concat = concat, begin_state = begin_state,
                         end_state = end_state)
}

#' Count state pairs across a set of column offsets
#'
#' @param mat Character matrix (rows = sequences).
#' @param states Character vector of states (matrix dimnames).
#' @param from_cols,to_cols Integer vectors of equal length: column pairs.
#' @param pair_weight Numeric vector, one weight per column pair.
#' @param row_weight Optional numeric vector, one weight per row of `mat`.
#' @return Numeric `states x states` matrix of weighted pair counts.
#' @noRd
.count_column_pairs <- function(mat, states, from_cols, to_cols, pair_weight,
                                row_weight = NULL) {
  n_states <- length(states)
  nbins <- n_states * n_states
  empty <- matrix(0, nrow = n_states, ncol = n_states,
                  dimnames = list(states, states))
  if (n_states == 0L || length(from_cols) == 0L || nrow(mat) == 0L) {
    return(empty)
  }
  row_weight <- row_weight %||% rep(1, nrow(mat))
  contributions <- Map(function(fc, tc, w) {
    from <- match(mat[, fc], states)
    to <- match(mat[, tc], states)
    ok <- !is.na(from) & !is.na(to)
    if (!any(ok)) return(numeric(nbins))
    idx <- (from[ok] - 1L) * n_states + to[ok]
    sums <- rowsum(w * row_weight[ok], idx, reorder = FALSE)
    out <- numeric(nbins)
    out[as.integer(rownames(sums))] <- sums[, 1L]
    out
  }, from_cols, to_cols, pair_weight)
  counts <- Reduce(`+`, contributions, numeric(nbins))
  matrix(counts, nrow = n_states, ncol = n_states, byrow = TRUE,
         dimnames = list(states, states))
}

#' Assemble the estimator return list shared by the three estimators
#' @noRd
.sequence_estimator_result <- function(mat_counts, data, states, format,
                                       action, id, time, cols, concat,
                                       begin_state, end_state) {
  resolved_format <- if (identical(format, "auto")) {
    if (action %in% names(data)) "long" else "wide"
  } else format
  initial <- .compute_initial_probs(
    data, states, format = resolved_format,
    action = action, id = id, time = time, cols = cols, concat = concat,
    begin_state = begin_state, end_state = end_state
  )
  list(
    matrix = mat_counts,
    nodes = states,
    directed = TRUE,
    cleaned_data = data,
    frequency_matrix = mat_counts,
    initial = initial
  )
}

#' n-gram estimator (tna type "n-gram")
#'
#' Every window of `n_gram` consecutive positions contributes each of its
#' adjacent pairs once, so interior pairs are counted up to `n_gram - 1`
#' times. Matches `tna::build_model(type = "n-gram")`.
#' @noRd
.estimator_ngram <- function(data, format = "auto", action = "Action",
                             id = NULL, time = "Time", cols = NULL,
                             alphabet = NULL, concat = 1L,
                             begin_state = NULL, end_state = NULL,
                             n_gram = 2L, ...) {
  stopifnot(
    "`n_gram` must be a single integer >= 2" =
      is.numeric(n_gram) && length(n_gram) == 1L && !is.na(n_gram) &&
        n_gram >= 2
  )
  n_gram <- as.integer(n_gram)
  prep <- .sequence_matrix(data, format, action, id, time, cols, alphabet,
                           concat, begin_state, end_state)
  mat <- prep$matrix
  p <- ncol(mat)
  n_windows <- p - n_gram + 1L
  if (n_windows >= 1L) {
    # Pair (j, j+1) is counted once for each window that contains it.
    window_starts <- seq_len(n_windows)
    pair_starts <- unlist(lapply(window_starts, \(i) seq(i, i + n_gram - 2L)))
    multiplicity <- tabulate(pair_starts, nbins = p - 1L)
    from_cols <- which(multiplicity > 0L)
    counts <- .count_column_pairs(mat, prep$states, from_cols, from_cols + 1L,
                                  multiplicity[from_cols])
  } else {
    counts <- .count_column_pairs(mat, prep$states, integer(0), integer(0),
                                  numeric(0))
  }
  .sequence_estimator_result(counts, data, prep$states, format, action, id,
                             time, cols, concat, begin_state, end_state)
}

#' Gap estimator (tna type "gap")
#'
#' Pairs `(i, j)` with `1 <= j - i <= max_gap + 1` are counted with weight
#' `1 / (j - i)`. Matches `tna::build_model(type = "gap")`, where
#' `max_gap = 1` (default) allows one intervening state.
#' @noRd
.estimator_gap <- function(data, format = "auto", action = "Action",
                           id = NULL, time = "Time", cols = NULL,
                           alphabet = NULL, concat = 1L,
                           begin_state = NULL, end_state = NULL,
                           max_gap = 1L, ...) {
  stopifnot(
    "`max_gap` must be a single non-negative integer" =
      is.numeric(max_gap) && length(max_gap) == 1L && !is.na(max_gap) &&
        max_gap >= 0
  )
  max_gap <- as.integer(max_gap)
  prep <- .sequence_matrix(data, format, action, id, time, cols, alphabet,
                           concat, begin_state, end_state)
  mat <- prep$matrix
  p <- ncol(mat)
  distances <- seq_len(max_gap + 1L)
  pairs <- expand.grid(from = seq_len(max(p - 1L, 0L)), dist = distances)
  pairs$to <- pairs$from + pairs$dist
  pairs <- pairs[pairs$to <= p, , drop = FALSE]
  counts <- .count_column_pairs(mat, prep$states, pairs$from, pairs$to,
                                1 / pairs$dist)
  .sequence_estimator_result(counts, data, prep$states, format, action, id,
                             time, cols, concat, begin_state, end_state)
}

#' Reverse (reply-network) estimator (tna type "reverse")
#'
#' Each adjacent pair is counted from the later state to the earlier one:
#' the transpose of the frequency matrix. With `weighted = TRUE` each pair is
#' weighted by the inverse of its sequence length, as in tna.
#' @noRd
.estimator_reverse <- function(data, format = "auto", action = "Action",
                               id = NULL, time = "Time", cols = NULL,
                               alphabet = NULL, concat = 1L,
                               begin_state = NULL, end_state = NULL,
                               weighted = FALSE, ...) {
  stopifnot("`weighted` must be TRUE or FALSE" = isTRUE(weighted) ||
              isFALSE(weighted))
  prep <- .sequence_matrix(data, format, action, id, time, cols, alphabet,
                           concat, begin_state, end_state)
  mat <- prep$matrix
  p <- ncol(mat)
  from_cols <- seq_len(max(p - 1L, 0L)) + 1L
  row_weight <- if (weighted) 1 / rowSums(!is.na(mat)) else NULL
  counts <- .count_column_pairs(mat, prep$states, from_cols, from_cols - 1L,
                                rep(1, length(from_cols)), row_weight)
  .sequence_estimator_result(counts, data, prep$states, format, action, id,
                             time, cols, concat, begin_state, end_state)
}
