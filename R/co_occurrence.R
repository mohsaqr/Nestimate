# ---- Co-occurrence Network Construction ----

#' Build a Co-occurrence Network
#'
#' @description
#' `cooccurrence()` builds an undirected network of items that appear
#' together in the same transaction (a document, a session, a record): an
#' edge joins two items that occur together, weighted by how often they do,
#' or by a similarity measure that adjusts that count for how often each item
#' occurs overall. With `similarity = "relative"` the network is directed. It reads delimited columns (keywords separated by `";"`),
#' a long table of items per transaction, a 0/1 incidence table, wide
#' sequences or a list of transactions; the *Details* say how each is
#' recognised.
#'
#' Raw counts grow with how frequent the two items are, so the most frequent
#' items dominate a count network. The similarity measures correct for that
#' in different ways (van Eck & Waltman, 2009); compare a few before reading
#' which pairs are strongest. The network is descriptive: it keeps no source
#' data, so [bootstrap_network()] and [permutation()] cannot resample it.
#'
#' @param data The transactions, in one of these forms:
#'   \itemize{
#'     \item a data frame with a delimited column (`field` and `sep`);
#'     \item a data frame in long (bipartite) form, one row per item in a
#'       transaction (`field` and `by`);
#'     \item a 0/1 or `TRUE`/`FALSE` data frame or matrix, one row per
#'       transaction and one column per item, declared with
#'       `format = "binary"`;
#'     \item a wide sequence data frame or matrix with other values, whose
#'       rows are the transactions;
#'     \item a plain list of character vectors, one per transaction.
#'   }
#'   A network object (`netobject`, `netobject_group`, `cograph_network`,
#'   `tna`) raises an error of class `nestimate_bad_input`; pass the data it
#'   was built from.
#' @param field The column holding the items, which become the nodes. With
#'   `sep`, its values are split by `sep`; several column names pool the split
#'   values of each row. With `by`, it is the item column of a long table.
#'   Default `NULL`. Setting `by` or `sep` without `field`, or `field`
#'   without either, raises an error of class `nestimate_bad_input`, as does
#'   naming a column that `data` does not have, or setting any of the three
#'   for a list.
#' @param by The column that defines the transactions of a long table: each
#'   of its values is one transaction, holding the distinct `field` values of
#'   its rows (for example `"paper_id"` or `"session_id"`). Default `NULL`,
#'   where items co-occur within the same row. It needs a single `field`
#'   column (several raise an error of class `nestimate_bad_input`). Ignored,
#'   with a warning of class `nestimate_ignored_argument`, when `sep` is
#'   set, which reads `field` as delimited.
#' @param sep The separator that splits delimited `field` values, such as
#'   `";"` or `","`, matched literally; the split values are trimmed. Default
#'   `NULL`. See `field` for the combinations that raise an error.
#' @param similarity How the co-occurrence count \eqn{C_{ij}}{C_ij} of items
#'   \eqn{i} and \eqn{j}, with item frequencies \eqn{f_i} and \eqn{f_j}, is
#'   turned into an edge weight. Default `"none"`.
#'   \describe{
#'     \item{`"none"`}{Raw co-occurrence counts.}
#'     \item{`"jaccard"`}{\eqn{C_{ij} / (f_i + f_j - C_{ij})}{C_ij /
#'       (f_i + f_j - C_ij)}.}
#'     \item{`"cosine"`}{\eqn{C_{ij} / \sqrt{f_i \cdot f_j}}{C_ij /
#'       sqrt(f_i * f_j)}
#'       (Salton's cosine).}
#'     \item{`"inclusion"`}{\eqn{C_{ij} / \min(f_i, f_j)}{C_ij / min(f_i, f_j)}
#'       (Simpson coefficient).}
#'     \item{`"association"`}{\eqn{C_{ij} / (f_i \cdot f_j)}{C_ij /
#'       (f_i * f_j)}
#'       (association strength, or probabilistic affinity index;
#'       van Eck & Waltman, 2009).}
#'     \item{`"dice"`}{\eqn{2 C_{ij} / (f_i + f_j)}{2 C_ij / (f_i + f_j)}.}
#'     \item{`"equivalence"`}{\eqn{C_{ij}^2 / (f_i \cdot f_j)}{C_ij^2 /
#'       (f_i * f_j)}
#'       (Salton's cosine squared).}
#'     \item{`"relative"`}{Each row divided by its sum (the diagonal
#'       included when `diagonal = TRUE`), so each row sums to 1: the weight
#'       from `i` to `j` is the share of `i`'s co-occurrences that involve
#'       `j`. The matrix is no longer symmetric, so the network is directed
#'       and both directions of a pair are edges.}
#'   }
#'   A zero denominator is replaced by 1, so an item without co-occurrences
#'   gets weight 0 rather than `NaN`.
#' @param threshold The smallest weight kept: weights below it, on the scale
#'   of `similarity`, are set to zero. Default `0` (none).
#' @param min_occur The fewest transactions an item must appear in to be
#'   kept. Rarer items are dropped before the co-occurrences are counted, and
#'   transactions left empty are dropped with them. Default `1` (keep all).
#' @param diagonal Whether the diagonal is kept. With `TRUE` (default) it
#'   holds each item's frequency (the number of transactions it appears in),
#'   transformed by `similarity`; with `FALSE` it is set to zero. Item
#'   frequencies for the similarity measures are the same either way.
#' @param top_n How many of the strongest off-diagonal edges to keep: pairs
#'   in an undirected network, directed cells with `similarity =
#'   "relative"`. Edges tied with the `top_n`-th weight are all kept, so a
#'   few more may remain, and the diagonal is not affected. Default `NULL`
#'   (all edges).
#' @param format How 0/1 columns are read. Columns of only 0 and 1 (or
#'   `TRUE` and `FALSE`) could be items or recoded positions, so with
#'   `"auto"` (default) they raise an error of class `nestimate_bad_input`.
#'   `"binary"` declares them items, one row per transaction, where a missing
#'   cell counts as no occurrence; it raises the same error on data with any
#'   other value.
#' @param ... Not used; an argument given here is ignored with a warning of
#'   class `nestimate_ignored_argument`.
#'
#' @return A `netobject` (class `c("netobject", "cograph_network")`),
#'   undirected (directed with `similarity = "relative"`), with
#'   `method = "co_occurrence_fn"`. Print it for an
#'   overview; [extract_edges()] gives one row per edge,
#'   [extract_transition_matrix()] the weight matrix (one row and column per
#'   item kept), [net_centrality()] the node centralities, and
#'   `cograph::splot()` draws it. It is built without [build_network()] and
#'   holds no source data (`data` is `NULL`). Its `params` record
#'   `similarity`, `threshold`, `min_occur`, `diagonal`, `top_n`,
#'   `n_transactions` and `n_items`.
#'
#' @details
#' Six input forms are recognised from `data` and from which of `field`,
#' `by` and `sep` are set:
#'
#' \enumerate{
#'   \item \strong{Delimited}: `field` (one column) and `sep`. Each cell is
#'     split by `sep`, trimmed, and its distinct values form a transaction.
#'   \item \strong{Multi-column delimited}: `field` (several columns) and
#'     `sep`. The split values of a row's columns are pooled into one
#'     transaction.
#'   \item \strong{Long bipartite}: `field` and `by`. The distinct `field`
#'     values within each `by` group form a transaction.
#'   \item \strong{Binary}: none of the three, all values 0/1 or logical,
#'     and `format = "binary"`. Columns are items, rows are transactions;
#'     a missing cell is no occurrence. Without `format = "binary"` such
#'     data is refused, since 0/1 columns could also be recoded positions.
#'   \item \strong{Wide sequence}: none of the three, other values. The
#'     distinct states of each row form a transaction, leaving out missing
#'     values and the void markers (`"%"`, `"*"`, `""`, `"NA"`, `"NaN"`). A
#'     data frame must be a single column or a numbered run of position
#'     columns (`T1`, `T2`, ...) with nothing beside it; any other frame
#'     raises an error of class `nestimate_bad_input`, and a column called
#'     `action` and repeated column names are refused the same way.
#'   \item \strong{List}: a plain list of character vectors.
#' }
#'
#' To count co-occurrence of states within sequences through
#' [build_network()], with its resampling verbs, use
#' `build_network(method = "co_occurrence")` or [build_cna()] instead.
#'
#' All forms are turned into a list of transactions, which `min_occur`
#' filters. The raw counts are `crossprod(B)` of the binary transaction
#' matrix `B`; `similarity` is applied to them, then `threshold`, then
#' `top_n`.
#'
#' @references
#' van Eck, N. J., & Waltman, L. (2009). How to normalize co-occurrence
#' data? An analysis of some well-known similarity measures. \emph{Journal of
#' the American Society for Information Science and Technology}, 60(8),
#' 1635--1651.
#'
#' @seealso \code{\link{build_cna}} for sequence-positional co-occurrence via
#'   \code{build_network()}.
#'
#' @examples
#' # Delimited field (e.g., keyword co-occurrence)
#' df <- data.frame(
#'   id = 1:4,
#'   keywords = c("network; graph", "graph; matrix; network",
#'                "matrix; algebra", "network; algebra; graph")
#' )
#' net <- cooccurrence(df, field = "keywords", sep = ";")
#'
#' # Long/bipartite
#' long_df <- data.frame(
#'   paper = c(1, 1, 1, 2, 2, 3, 3),
#'   keyword = c("network", "graph", "matrix", "graph", "algebra",
#'               "network", "algebra")
#' )
#' net <- cooccurrence(long_df, field = "keyword", by = "paper")
#'
#' # List of transactions
#' transactions <- list(c("A", "B"), c("B", "C"), c("A", "B", "C"))
#' net <- cooccurrence(transactions, similarity = "jaccard")
#'
#' # Binary matrix
#' bin <- matrix(c(1,0,1, 1,1,0, 0,1,1), nrow = 3, byrow = TRUE,
#'               dimnames = list(NULL, c("X", "Y", "Z")))
#' net <- cooccurrence(bin, format = "binary")
#'
#' @export
cooccurrence <- function(data, field = NULL, by = NULL, sep = NULL,
                          similarity = c("none", "jaccard", "cosine",
                                         "inclusion", "association",
                                         "dice", "equivalence", "relative"),
                          threshold = 0, min_occur = 1L,
                          diagonal = TRUE, top_n = NULL, ...,
                          format = c("auto", "binary")) {
  similarity <- match.arg(similarity)
  format <- match.arg(format)
  threshold <- as.numeric(threshold)
  min_occur <- as.integer(min_occur)
  .assert(
    "`threshold` must be a single number of 0 or more" =
      length(threshold) == 1L && threshold >= 0,
    "`min_occur` must be a single whole number of 1 or more" =
      length(min_occur) == 1L && min_occur >= 1L,
    "`diagonal` must be TRUE or FALSE" =
      is.logical(diagonal) && length(diagonal) == 1L)
  .warn_ignored_dots("cooccurrence()", ...)

  # A network is a list, so its fields would otherwise be read as transactions
  if (inherits(data, c("netobject", "netobject_group", "cograph_network", "tna"))) {
    stop(errorCondition(paste0(
      "cooccurrence(): `data` is a network, not data. Pass the data it was ",
      "built from, or use build_network(method = \"co_occurrence\") for ",
      "positional co-occurrence of sequences."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (is.list(data) && !is.data.frame(data) &&
      (!is.null(oldClass(data)) || !all(vapply(data, is.atomic, logical(1L))))) {
    stop(errorCondition(
      "cooccurrence(): a list of transactions must be a plain list holding one vector of items each.",
      class = "nestimate_bad_input", call = NULL))
  }
  if (is.data.frame(data)) {
    .check_declared_columns(data, "cooccurrence()", field = field, by = by)
  }
  # A role that the format ignores is a mistake, not a request to read the input
  # some other way: `by` and `sep` both need `field`, and `field` needs one of them.
  # A list of transactions has no columns to name.
  roles_set <- !is.null(field) || !is.null(by) || !is.null(sep)
  if (is.list(data) && !is.data.frame(data) && roles_set) {
    stop(errorCondition(
      "cooccurrence(): `field`, `by` and `sep` name columns; a list has none.",
      class = "nestimate_bad_input", call = NULL))
  }
  if (roles_set && (is.null(field) || (is.null(by) && is.null(sep)))) {
    stop(errorCondition(
      paste0("cooccurrence(): set `field` and `by` for a long table, or `field` ",
             "and `sep` for delimited columns."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (!is.null(by) && is.null(sep) && length(field) != 1L) {
    stop(errorCondition(
      paste0("cooccurrence(): a long table has one item column; `field` names ",
             length(field), ". Name one, or set `sep` to pool delimited columns."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (!is.null(by) && !is.null(sep)) {
    warning(warningCondition(
      paste0("cooccurrence(): `by` is ignored when `sep` is set: `field` is ",
             "read as delimited, one transaction per row."),
      class = "nestimate_ignored_argument", call = NULL))
  }
  # Convert input to list of character vectors (transactions)
  fmt <- .co_detect_format(data, field, by, sep)
  # Columns of 0 and 1 could be items or recoded positions: no reading is guessed.
  if (identical(fmt, "binary") && identical(format, "auto")) {
    stop(errorCondition(
      paste0("cooccurrence(): these columns hold only 0 and 1, so they could be ",
             "items or positions. Set `format = \"binary\"` if they are items, where ",
             "a missing cell counts as no occurrence."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (identical(format, "binary") && !identical(fmt, "binary")) {
    stop(errorCondition(
      "cooccurrence(): `format = \"binary\"` needs a frame or matrix of only 0 and 1 (or TRUE and FALSE).",
      class = "nestimate_bad_input", call = NULL))
  }
  transactions <- switch(fmt,
    delimited       = .co_parse_delimited(data, field, sep),
    multi_delimited = .co_parse_multi_delimited(data, field, sep),
    long            = .co_parse_long(data, field, by),
    binary          = .co_parse_binary(data),
    wide            = {
      if (is.data.frame(data)) {
        .refuse_undeclared_roles(
          data, "cooccurrence()", whole_frame = TRUE,
          advice = paste0("Set `field` and `by` for a long table, or `field` ",
                          "and `sep` for delimited columns."))
      }
      .co_parse_wide(data)
    },
    list            = .co_parse_list(data)
  )

  # Drop empty transactions
  transactions <- transactions[vapply(transactions, length, integer(1)) > 0L]
  if (length(transactions) == 0L)
    stop(errorCondition("No non-empty transactions found in the input data.",
                        class = "nestimate_bad_input", call = NULL))

  # Apply min_occur filter (drop infrequent entities)
  if (min_occur > 1L) {
    all_items <- unlist(transactions)
    freq_table <- table(all_items)
    keep <- names(freq_table[freq_table >= min_occur])
    transactions <- lapply(transactions, function(t) t[t %in% keep])
    transactions <- transactions[vapply(transactions, length, integer(1)) > 0L]
    if (length(transactions) == 0L)
      stop(errorCondition("No transactions remain after min_occur filtering.",
                          class = "nestimate_bad_input", call = NULL))
  }

  # Build binary transaction matrix and compute co-occurrence
  B <- .co_transactions_to_matrix(transactions)
  C <- .co_compute_matrix(B)
  n_trans <- nrow(B)

  # Diagonal handling
  if (!diagonal) diag(C) <- 0

  # Item frequencies (from diagonal of raw co-occurrence = column sums of B)
  freq <- diag(C)
  if (!diagonal) freq <- colSums(B)

  # Normalize
  W <- .co_normalize(C, freq, n_trans, similarity)

  # Apply threshold (on final weight, after normalization)
  if (threshold > 0) W[W < threshold] <- 0

  # Row-normalised ("relative") weights are directional: W[i, j] is the
  # share of i's co-occurrences that involve j, so W[i, j] != W[j, i]. The
  # network is directed then, so every cell reaches the edge list.
  directed <- identical(similarity, "relative")

  # Apply top_n (keep only top N edges by weight). An edge is one cell of
  # the upper triangle when undirected, one off-diagonal cell when directed;
  # the cutoff is taken from those same cells. Ties with the cutoff stay.
  if (!is.null(top_n)) {
    .assert("`top_n` must be NULL or a positive number" =
              is.numeric(top_n) && all(top_n > 0))
    top_n <- as.integer(top_n)
    candidates <- if (directed) row(W) != col(W) else upper.tri(W)
    values <- W[candidates]
    if (length(values) > top_n) {
      cutoff <- sort(values, decreasing = TRUE)[top_n]
      drop <- candidates & W < cutoff
      if (!directed) drop <- drop | t(drop)
      W[drop] <- 0
    }
  }

  # Wrap as netobject
  net <- .wrap_netobject(
    W,
    data = NULL,
    method = "co_occurrence_fn",
    directed = directed
  )

  net$params <- list(
    similarity = similarity,
    threshold  = threshold,
    min_occur  = min_occur,
    diagonal   = diagonal,
    top_n      = top_n,
    n_transactions = n_trans,
    n_items = ncol(B)
  )

  net
}


# ---- Format detection ----

#' Detect input format from argument combination
#' @noRd
.co_detect_format <- function(data, field, by, sep) {
  if (is.list(data) && !is.data.frame(data) && !is.matrix(data))
    return("list")

  if (!is.null(sep) && !is.null(field) && length(field) > 1L)
    return("multi_delimited")

  if (!is.null(sep) && !is.null(field))
    return("delimited")

  if (!is.null(field) && !is.null(by))
    return("long")

  # Matrix or data.frame without field/by/sep
  if (is.matrix(data) || is.data.frame(data)) {
    mat <- if (is.data.frame(data)) as.matrix(data) else data
    if ((is.numeric(mat) || is.logical(mat)) && all(mat[!is.na(mat)] %in% c(0, 1)))
      return("binary")
    return("wide")
  }

  stop(errorCondition(paste0(
    "Cannot detect input format. Provide field/by/sep arguments or a ",
    "recognized data structure (data.frame, matrix, list)."),
    class = "nestimate_bad_input", call = NULL))
}


# ---- Parsers: each returns list of character vectors ----

#' Parse delimited field (single column)
#' @noRd
.co_parse_delimited <- function(data, field, sep) {
  .assert(is.data.frame(data), length(field) == 1L, field %in% names(data))
  vals <- as.character(data[[field]])
  lapply(strsplit(vals, sep, fixed = TRUE), function(items) {
    items <- trimws(items)
    items <- items[nzchar(items) & !is.na(items)]
    unique(items)
  })
}

#' Parse multi-column delimited (pool across columns)
#' @noRd
.co_parse_multi_delimited <- function(data, field, sep) {
  .assert(is.data.frame(data), all(field %in% names(data)))
  n <- nrow(data)
  lapply(seq_len(n), function(i) {
    items <- unlist(lapply(field, function(f) {
      strsplit(as.character(data[[f]][i]), sep, fixed = TRUE)[[1L]]
    }))
    items <- trimws(items)
    items <- items[nzchar(items) & !is.na(items)]
    unique(items)
  })
}

#' Parse long/bipartite format
#' @noRd
.co_parse_long <- function(data, field, by) {
  .assert(is.data.frame(data), field %in% names(data), by %in% names(data))
  groups <- split(as.character(data[[field]]), data[[by]])
  lapply(groups, function(items) {
    items <- items[nzchar(items) & !is.na(items)]
    unique(items)
  })
}

#' Parse binary matrix (columns = items, rows = transactions)
#' @noRd
.co_parse_binary <- function(data) {
  mat <- if (is.data.frame(data)) as.matrix(data) else data
  if (is.null(colnames(mat)))
    colnames(mat) <- paste0("V", seq_len(ncol(mat)))
  cn <- colnames(mat)
  # a missing cell is not an occurrence; a bare `==` would keep it as an NA item
  lapply(seq_len(nrow(mat)), function(i) cn[which(mat[i, ] == 1)])
}

#' Parse wide sequence data (unique values per row = transaction)
#' @noRd
.co_parse_wide <- function(data) {
  mat <- if (is.data.frame(data)) as.matrix(data) else data
  lapply(seq_len(nrow(mat)), function(i) {
    vals <- as.character(mat[i, ])
    vals <- vals[!is.na(vals) & nzchar(vals) & !(vals %in% .void_markers)]
    unique(vals)
  })
}

#' Parse list of character vectors
#' @noRd
.co_parse_list <- function(data) {
  lapply(data, function(items) {
    items <- as.character(items)
    items <- items[!is.na(items) & nzchar(items)]
    unique(items)
  })
}


# ---- Core computation ----

#' Build binary transaction matrix from list of character vectors
#' @return A logical matrix (rows = transactions, cols = items)
#' @noRd
.co_transactions_to_matrix <- function(transactions) {
  all_items <- sort(unique(unlist(transactions)))
  n <- length(transactions)
  k <- length(all_items)
  B <- matrix(FALSE, nrow = n, ncol = k,
              dimnames = list(NULL, all_items))
  B[cbind(rep(seq_len(n), lengths(transactions)),
          match(unlist(transactions, use.names = FALSE), all_items))] <- TRUE
  B
}

#' Compute raw co-occurrence matrix from binary transaction matrix
#' @return Numeric matrix (symmetric, items x items)
#' @noRd
.co_compute_matrix <- function(B) {
  C <- as.matrix(crossprod(B * 1L))
  storage.mode(C) <- "double"
  C
}

#' Normalize a co-occurrence matrix
#' @param C Raw co-occurrence matrix (symmetric).
#' @param freq Item frequency vector (length = ncol(C)).
#' @param N Number of transactions.
#' @param method Similarity method name.
#' @return Normalized matrix.
#' @noRd
.co_normalize <- function(C, freq, N, method) {
  if (method == "none") return(C)

  W <- C

  if (method == "jaccard") {
    denom <- outer(freq, freq, "+") - C
    denom[denom == 0] <- 1
    W <- C / denom

  } else if (method == "cosine") {
    denom <- outer(sqrt(freq), sqrt(freq), "*")
    denom[denom == 0] <- 1
    W <- C / denom

  } else if (method == "inclusion") {
    denom <- outer(freq, freq, pmin)
    denom[denom == 0] <- 1
    W <- C / denom

  } else if (method == "association") {
    denom <- outer(freq, freq, "*")
    denom[denom == 0] <- 1
    W <- C / denom

  } else if (method == "dice") {
    denom <- outer(freq, freq, "+")
    denom[denom == 0] <- 1
    W <- 2 * C / denom

  } else if (method == "equivalence") {
    denom <- outer(freq, freq, "*")
    denom[denom == 0] <- 1
    W <- C^2 / denom

  } else if (method == "relative") {
    rs <- rowSums(C)
    rs[rs == 0] <- 1
    W <- C / rs
  }

  W
}
