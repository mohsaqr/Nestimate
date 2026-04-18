# ---- Co-occurrence Network Construction ----

#' Build a Co-occurrence Network
#'
#' Constructs an undirected co-occurrence network from various input formats.
#' Entities that appear together in the same transaction, document, or record
#' are connected, with edge weights reflecting raw counts or a similarity
#' measure. Argument names follow the citenets convention.
#'
#' @param data Input data. Accepts:
#'   \itemize{
#'     \item A \code{data.frame} with a delimited column (\code{field} + \code{sep}).
#'     \item A \code{data.frame} in long/bipartite format (\code{field} + \code{by}).
#'     \item A binary (0/1) \code{data.frame} or \code{matrix} (auto-detected).
#'     \item A wide sequence \code{data.frame} or \code{matrix} (non-binary).
#'     \item A \code{list} of character vectors (each element is a transaction).
#'   }
#' @param field Character. The entity column — determines what the nodes are.
#'   For delimited format, a single column whose values are split by \code{sep}.
#'   For long/bipartite format, the item column. For multi-column delimited,
#'   a vector of column names whose split values are pooled per row.
#' @param by Character or \code{NULL}. What links the nodes. For
#'   long/bipartite format, the grouping column (e.g., \code{"paper_id"},
#'   \code{"session_id"}). Each unique value of \code{by} defines one
#'   transaction. If \code{NULL} (default), entities co-occur within the
#'   same row/document.
#' @param sep Character or \code{NULL}. Separator for splitting delimited
#'   fields (e.g., \code{";"}, \code{","}). Default \code{NULL}.
#' @param similarity Character. Similarity measure applied to the raw
#'   co-occurrence counts. One of:
#'   \describe{
#'     \item{\code{"none"}}{Raw co-occurrence counts.}
#'     \item{\code{"jaccard"}}{\eqn{C_{ij} / (f_i + f_j - C_{ij})}.}
#'     \item{\code{"cosine"}}{\eqn{C_{ij} / \sqrt{f_i \cdot f_j}}
#'       (Salton's cosine).}
#'     \item{\code{"inclusion"}}{\eqn{C_{ij} / \min(f_i, f_j)}
#'       (Simpson coefficient).}
#'     \item{\code{"association"}}{\eqn{C_{ij} / (f_i \cdot f_j)}
#'       (association strength / probabilistic affinity index;
#'       van Eck & Waltman, 2009).}
#'     \item{\code{"dice"}}{\eqn{2 C_{ij} / (f_i + f_j)}.}
#'     \item{\code{"equivalence"}}{\eqn{C_{ij}^2 / (f_i \cdot f_j)}
#'       (Salton's cosine squared).}
#'     \item{\code{"relative"}}{Row-normalized: each row sums to 1.}
#'   }
#' @param threshold Numeric. Minimum edge weight to retain. Edges below this
#'   value are set to zero. Applied \emph{after} similarity normalization.
#'   Default 0.
#' @param min_occur Integer. Minimum entity frequency (number of transactions
#'   an entity must appear in). Entities below this threshold are dropped
#'   before computing co-occurrence. Default 1 (keep all).
#' @param diagonal Logical. If \code{TRUE} (default), the diagonal of the
#'   co-occurrence matrix is kept (item self-co-occurrence = item frequency).
#'   If \code{FALSE}, the diagonal is zeroed.
#' @param top_n Integer or \code{NULL}. If specified, return only the top
#'   \code{top_n} edges by weight. Default \code{NULL} (all edges).
#' @param ... Currently unused.
#'
#' @return A \code{netobject} (undirected) with \code{method = "co_occurrence_fn"}.
#'   The \code{$weights} matrix contains similarity (or raw) co-occurrence values.
#'   The \code{$params} list stores the similarity method, threshold, and
#'   the number of transactions.
#'
#' @details
#' Six input formats are supported, auto-detected from the combination of
#' \code{field}, \code{by}, and \code{sep}:
#'
#' \enumerate{
#'   \item \strong{Delimited}: \code{field} + \code{sep} (single column).
#'     Each cell is split by \code{sep}, trimmed, and de-duplicated per row.
#'   \item \strong{Multi-column delimited}: \code{field} (vector) + \code{sep}.
#'     Values from multiple columns are split, pooled, and de-duplicated per row.
#'   \item \strong{Long bipartite}: \code{field} + \code{by}.
#'     Groups by \code{by}; unique values of \code{field} within each group
#'     form a transaction.
#'   \item \strong{Binary matrix}: No \code{field}/\code{by}/\code{sep}, all
#'     values 0/1. Columns are items, rows are transactions.
#'   \item \strong{Wide sequence}: No \code{field}/\code{by}/\code{sep},
#'     non-binary. Unique values across each row form a transaction.
#'   \item \strong{List}: A plain list of character vectors.
#' }
#'
#' The pipeline converts all formats into a list of character vectors
#' (transactions), optionally filters by \code{min_occur}, builds a binary
#' transaction matrix, computes \code{crossprod(B)} for the raw co-occurrence
#' counts, normalizes via the chosen \code{similarity}, then applies
#' \code{threshold} and \code{top_n} filtering.
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
#' net <- co_occurrence(df, field = "keywords", sep = ";")
#'
#' # Long/bipartite
#' long_df <- data.frame(
#'   paper = c(1, 1, 1, 2, 2, 3, 3),
#'   keyword = c("network", "graph", "matrix", "graph", "algebra",
#'               "network", "algebra")
#' )
#' net <- co_occurrence(long_df, field = "keyword", by = "paper")
#'
#' # List of transactions
#' transactions <- list(c("A", "B"), c("B", "C"), c("A", "B", "C"))
#' net <- co_occurrence(transactions, similarity = "jaccard")
#'
#' # Binary matrix
#' bin <- matrix(c(1,0,1, 1,1,0, 0,1,1), nrow = 3, byrow = TRUE,
#'               dimnames = list(NULL, c("X", "Y", "Z")))
#' net <- co_occurrence(bin)
#'
#' @export
co_occurrence <- function(data, field = NULL, by = NULL, sep = NULL,
                          similarity = c("none", "jaccard", "cosine",
                                         "inclusion", "association",
                                         "dice", "equivalence", "relative"),
                          threshold = 0, min_occur = 1L,
                          diagonal = TRUE, top_n = NULL, ...) {
  similarity <- match.arg(similarity)
  threshold <- as.numeric(threshold)
  min_occur <- as.integer(min_occur)
  stopifnot(threshold >= 0, min_occur >= 1L,
            is.logical(diagonal), length(diagonal) == 1L)

  # Convert input to list of character vectors (transactions)
  fmt <- .co_detect_format(data, field, by, sep)
  transactions <- switch(fmt,
    delimited       = .co_parse_delimited(data, field, sep),
    multi_delimited = .co_parse_multi_delimited(data, field, sep),
    long            = .co_parse_long(data, field, by),
    binary          = .co_parse_binary(data),
    wide            = .co_parse_wide(data),
    list            = .co_parse_list(data)
  )

  # Drop empty transactions
  transactions <- transactions[vapply(transactions, length, integer(1)) > 0L]
  if (length(transactions) == 0L)
    stop("No non-empty transactions found in the input data.", call. = FALSE)

  # Apply min_occur filter (drop infrequent entities)
  if (min_occur > 1L) {
    all_items <- unlist(transactions)
    freq_table <- table(all_items)
    keep <- names(freq_table[freq_table >= min_occur])
    transactions <- lapply(transactions, function(t) t[t %in% keep])
    transactions <- transactions[vapply(transactions, length, integer(1)) > 0L]
    if (length(transactions) == 0L)
      stop("No transactions remain after min_occur filtering.", call. = FALSE)
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

  # Apply top_n (keep only top N edges by weight)
  if (!is.null(top_n)) {
    stopifnot(is.numeric(top_n), top_n > 0)
    top_n <- as.integer(top_n)
    # Extract upper triangle values, find cutoff
    ut <- W[upper.tri(W)]
    if (length(ut) > top_n) {
      cutoff <- sort(ut, decreasing = TRUE)[top_n]
      W[W < cutoff & upper.tri(W)] <- 0
      W[W < cutoff & lower.tri(W)] <- 0
    }
  }

  # Wrap as netobject
  net <- .wrap_netobject(
    W,
    data = NULL,
    method = "co_occurrence_fn",
    directed = FALSE
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
    if (is.numeric(mat) && all(mat[!is.na(mat)] %in% c(0, 1)))
      return("binary")
    return("wide")
  }

  stop("Cannot detect input format. Provide field/by/sep arguments or a ",
       "recognized data structure (data.frame, matrix, list).", call. = FALSE)
}


# ---- Parsers: each returns list of character vectors ----

#' Parse delimited field (single column)
#' @noRd
.co_parse_delimited <- function(data, field, sep) {
  stopifnot(is.data.frame(data), length(field) == 1L, field %in% names(data))
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
  stopifnot(is.data.frame(data), all(field %in% names(data)))
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
  stopifnot(is.data.frame(data), field %in% names(data), by %in% names(data))
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
  lapply(seq_len(nrow(mat)), function(i) cn[mat[i, ] == 1])
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
  for (i in seq_len(n)) {
    B[i, transactions[[i]]] <- TRUE
  }
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
