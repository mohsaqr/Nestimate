# ---- Association Rule Mining ----

#' Discover Association Rules from Sequential or Transaction Data
#'
#' @description
#' `association_rules()` finds association rules of the form
#' \{A\} -> \{B\} with the Apriori algorithm (Agrawal & Srikant, 1994): sets
#' of items that occur together in at least `min_support` of the
#' transactions, and the rules among them that pass `min_confidence` and
#' `min_lift`. It reads a list of transactions, a wide data frame or matrix of
#' states, a 0/1 incidence table, or a network: the sequences behind a
#' `netobject` or `tna` model become transactions of their distinct states,
#' and a one-hot network's windows become transactions of their active codes.
#'
#' Support is how common the whole rule is, confidence how often B occurs
#' when A does, and lift how much more often A and B occur together than they
#' would if they were independent (see *Metrics*). A sequence read as a
#' transaction loses its order, so a rule says that states occur in the same
#' sequence, not that one follows the other; [build_network()] models order.
#' The rules are not tested for significance, and lowering the thresholds
#' quickly multiplies them, so check whether a rule of interest stays when
#' `min_support` changes.
#'
#' @param x The transactions, in one of these forms:
#'   \describe{
#'     \item{`netobject`}{Each sequence in its data becomes a transaction of
#'       its distinct states. In a one-hot network the codes are the items
#'       and each window is one transaction. A network without data, and an
#'       association network (`cor`, `pcor`, `glasso`, `ising`, `mgm`),
#'       whose rows are variables, raise an error; a `netobject_group` raises
#'       one of class `nestimate_bad_input` (pass one group's network). A
#'       bare `cograph_network` is read as a `netobject`, a `tna` model by
#'       its sequence data, and an `mcml` by its between-cluster network, so
#'       its items are the clusters.}
#'     \item{list}{A plain list holding one character vector of items per
#'       transaction.}
#'     \item{data frame}{Each row is a transaction of its distinct
#'       non-missing states. The frame must be a single column or a numbered
#'       run of position columns (`T1`, `T2`, ...) with nothing beside it;
#'       otherwise it raises an error of class `nestimate_bad_input`;
#'       declare the columns with `build_network(state_cols = )` and pass the
#'       network instead. A column called `action` and repeated column names
#'       are refused the same way. A frame of 0/1 or logical columns is a
#'       binary table, one column per item, read with `format = "binary"`.}
#'     \item{matrix}{A character matrix is read like a character data frame.
#'       A numeric or logical matrix is a binary transaction table (rows are
#'       transactions, columns are items, named `I1`, `I2`, ... when the
#'       matrix has no column names).}
#'   }
#'   Binary input must hold only 0/1 or `TRUE`/`FALSE`: a missing cell or any
#'   other value raises an error of class `nestimate_bad_input`. For the
#'   `print()` and `plot()` methods: a `net_association_rules` object.
#' @param min_support The smallest share of transactions an itemset must
#'   occur in to be frequent, above 0 and at most 1. Default `0.1`. It is
#'   applied as a count, `ceiling(min_support * n_transactions)`.
#' @param min_confidence The smallest confidence a rule must have, between
#'   0 and 1. Default `0.5`.
#' @param min_lift The smallest lift a rule must have, at least 0. Default
#'   `1`, which keeps rules whose items occur together at least as often as
#'   independence predicts.
#' @param max_length The largest itemset size searched, at least 2. Default
#'   `5`.
#' @param format How 0/1 columns are read. Columns of only 0 and 1 could be
#'   items or recoded positions, so with `"auto"` (default) they raise an
#'   error of class `nestimate_bad_input`. `"binary"` declares them items,
#'   one row per transaction; it raises the same error on a data frame with
#'   any other value. A list of transactions or a network has no 0/1
#'   columns, so `"binary"` is ignored for them with a warning of class
#'   `nestimate_ignored_argument`.
#'
#' @return A `net_association_rules` object. `print()` shows the counts, the
#'   thresholds and the ten rules with the highest lift; `summary()` returns
#'   the rules table, one row per rule; `plot()` draws the rules by support
#'   and confidence. The object's fields:
#' \describe{
#'   \item{rules}{Data frame, one row per rule, ordered by decreasing lift
#'     and then confidence, with columns `antecedent` and `consequent` (the
#'     itemsets as comma-separated strings), `support`, `confidence`, `lift`,
#'     `conviction`, `count` (transactions holding the whole rule) and
#'     `n_transactions`.}
#'   \item{frequent}{Data frame, one row per frequent itemset, with columns
#'     `itemset`, `size`, `support` and `count`.}
#'   \item{frequent_itemsets}{List of the frequent itemsets per size k, each
#'     entry a list of `items`, `count` and `support`.}
#'   \item{items}{Character vector of the frequent single items the mining
#'     ran on (all items when none reaches `min_support`).}
#'   \item{n_transactions}{Number of transactions.}
#'   \item{n_rules}{Number of rules.}
#'   \item{params}{List of `min_support`, `min_confidence`, `min_lift` and
#'     `max_length`.}
#' }
#'
#' @details
#' ## Algorithm
#'
#' Level-wise Apriori (Agrawal & Srikant, 1994) with the full pruning step:
#' after the join step generates the candidates of size k, every subset of
#' size k - 1 is checked to be frequent before support is counted, which
#' keeps the search fast at k >= 4. Support of pairs is counted with
#' `crossprod()`, of larger itemsets by logical indexing. Every split of a
#' frequent itemset into a non-empty antecedent and consequent is a
#' candidate rule.
#'
#' ## Metrics
#'
#' \describe{
#'   \item{support}{P(A and B): the share of transactions holding both the
#'     antecedent and the consequent.}
#'   \item{confidence}{P(B | A): the share of transactions holding the
#'     antecedent that also hold the consequent.}
#'   \item{lift}{P(A and B) / (P(A) P(B)). Values above 1 indicate positive
#'     association, below 1 negative association.}
#'   \item{conviction}{(1 - P(B)) / (1 - confidence), the departure from
#'     independence in the direction of the rule (Brin et al., 1997); higher
#'     means a stronger implication, and a rule with confidence 1 has
#'     conviction `Inf`.}
#' }
#'
#' @references
#' Agrawal, R. & Srikant, R. (1994). Fast algorithms for mining association
#' rules. In \emph{Proc. 20th VLDB Conference}, 487--499.
#'
#' Brin, S., Motwani, R., Ullman, J. D. & Tsur, S. (1997). Dynamic itemset
#' counting and implication rules for market basket data. In \emph{Proc. ACM
#' SIGMOD}, 255--264. (lift and conviction)
#'
#' @examples
#' # From a list of transactions
#' trans <- list(
#'   c("plan", "discuss", "execute"),
#'   c("plan", "research", "analyze"),
#'   c("discuss", "execute", "reflect"),
#'   c("plan", "discuss", "execute", "reflect"),
#'   c("research", "analyze", "reflect")
#' )
#' rules <- association_rules(trans, min_support = 0.3, min_confidence = 0.5)
#' print(rules)
#'
#' # From a netobject (sequences as transactions)
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:5], 50, TRUE),
#'   V2 = sample(LETTERS[1:5], 50, TRUE),
#'   V3 = sample(LETTERS[1:5], 50, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' rules <- association_rules(net, min_support = 0.1)
#'
#' @seealso \code{\link{build_network}}, \code{\link{predict_links}}
#'
#' @export
association_rules <- function(x,
                              min_support = 0.1,
                              min_confidence = 0.5,
                              min_lift = 1.0,
                              max_length = 5L,
                              format = c("auto", "binary")) {
  format <- match.arg(format)

  # ---- Input validation ----
  .assert(
    "`min_support` must be a single number above 0 and at most 1" =
      is.numeric(min_support) && length(min_support) == 1 &&
      min_support > 0 && min_support <= 1,
    "`min_confidence` must be a single number from 0 to 1" =
      is.numeric(min_confidence) && length(min_confidence) == 1 &&
      min_confidence >= 0 && min_confidence <= 1,
    "`min_lift` must be a single number of 0 or more" =
      is.numeric(min_lift) && length(min_lift) == 1 && min_lift >= 0,
    "`max_length` must be a single number of 2 or more" =
      is.numeric(max_length) && length(max_length) == 1 && max_length >= 2
  )
  max_length <- as.integer(max_length)
  if (identical(format, "binary") && is.list(x) && !is.data.frame(x)) {
    warning(warningCondition(paste0(
      "association_rules(): `format = \"binary\"` is ignored: it declares ",
      "the 0/1 columns of a data frame or matrix, and `x` is a list of ",
      "transactions or a network."),
      class = "nestimate_ignored_argument", call = NULL))
  }

  # ---- Parse input to binary transaction matrix ----
  parsed <- .ar_parse_input(x, format)
  trans_mat <- parsed$matrix  # logical matrix: rows=transactions, cols=items
  items <- parsed$items
  n_trans <- nrow(trans_mat)
  n_items <- length(items)
  min_count <- ceiling(min_support * n_trans)

  # ---- Level 1: frequent 1-itemsets ----
  item_counts <- colSums(trans_mat)
  freq1_mask <- item_counts >= min_count
  freq1_items <- items[freq1_mask]
  if (length(freq1_items) == 0) {
    return(.ar_empty_result(items, n_trans, min_support, min_confidence,
                            min_lift, max_length))
  }
  # Prune transaction matrix to frequent items only
  trans_mat <- trans_mat[, freq1_items, drop = FALSE]
  items <- freq1_items
  n_items <- length(items)
  item_counts <- item_counts[freq1_mask]

  frequent_itemsets <- list()
  frequent_itemsets[[1]] <- lapply(seq_along(items), function(i) {
    list(items = items[i], count = item_counts[i],
         support = item_counts[i] / n_trans)
  })

  # ---- Level 2: frequent 2-itemsets via crossprod ----
  if (max_length >= 2 && n_items >= 2) {
    co_counts <- as.matrix(crossprod(trans_mat * 1L))
    # Upper triangle only (avoid duplicates)
    pairs <- unname(which(upper.tri(co_counts) & co_counts >= min_count,
                          arr.ind = TRUE))
    # Row by row, as the pairs are listed: (1, 2), (1, 3), ..., (2, 3), ...
    pairs <- pairs[order(pairs[, 1L], pairs[, 2L]), , drop = FALSE]
    freq2 <- lapply(seq_len(nrow(pairs)), function(r) {
      i <- pairs[r, 1L]
      j <- pairs[r, 2L]
      list(items = c(items[i], items[j]),
           count = co_counts[i, j],
           support = co_counts[i, j] / n_trans)
    })
    if (length(freq2) > 0) frequent_itemsets[[2]] <- freq2
  }

  # ---- Level k >= 3: join + prune + count ----
  if (max_length >= 3) {
    k <- 3L
    while (k <= max_length) {
      if (length(frequent_itemsets) < k - 1L) break
      prev <- frequent_itemsets[[k - 1L]]
      if (is.null(prev) || length(prev) < 2L) break

      prev_items <- lapply(prev, `[[`, "items")
      # Join step: merge (k-1)-itemsets sharing first k-2 items
      candidates <- .ar_generate_candidates(prev_items, k)
      if (length(candidates) == 0) break

      # Prune step: verify all (k-1)-subsets are frequent
      prev_set <- vapply(prev_items, paste, character(1), collapse = "\t")
      candidates <- Filter(function(cand) {
        subsets <- combn(cand, k - 1L, simplify = FALSE)
        all(vapply(subsets, function(s) paste(s, collapse = "\t"), character(1))
            %in% prev_set)
      }, candidates)
      if (length(candidates) == 0) break

      # Count support via logical AND
      freq_k <- Filter(Negate(is.null), lapply(candidates, function(cand) {
        mask <- Reduce(`&`, lapply(cand, function(it) trans_mat[, it]))
        cnt <- sum(mask)
        if (cnt >= min_count) {
          list(items = cand, count = cnt, support = cnt / n_trans)
        }
      }))
      if (length(freq_k) == 0) break
      frequent_itemsets[[k]] <- freq_k
      k <- k + 1L
    }
  }

  # ---- Generate rules from frequent itemsets ----
  rules <- .ar_generate_rules(frequent_itemsets, item_counts, items,
                               n_trans, min_confidence, min_lift)

  # ---- Build tidy rules data frame ----
  if (nrow(rules) > 0) {
    tidy_rules <- data.frame(
      antecedent = vapply(rules$antecedent, paste, character(1), collapse = ", "),
      consequent = vapply(rules$consequent, paste, character(1), collapse = ", "),
      support    = rules$support,
      confidence = rules$confidence,
      lift       = rules$lift,
      conviction = rules$conviction,
      count      = rules$count,
      n_transactions = rules$n_transactions,
      stringsAsFactors = FALSE
    )
  } else {
    tidy_rules <- data.frame(
      antecedent = character(0), consequent = character(0),
      support = numeric(0), confidence = numeric(0),
      lift = numeric(0), conviction = numeric(0),
      count = integer(0), n_transactions = integer(0),
      stringsAsFactors = FALSE
    )
  }

  # ---- Build tidy frequent itemsets data frame ----
  fi_rows <- lapply(seq_along(frequent_itemsets), function(k) {
    level <- frequent_itemsets[[k]]
    data.frame(
      itemset = vapply(level, function(fi) paste(fi$items, collapse = ", "), character(1)),
      size    = k,
      support = vapply(level, `[[`, numeric(1), "support"),
      count   = vapply(level, `[[`, numeric(1), "count"),
      stringsAsFactors = FALSE
    )
  })
  frequent <- do.call(rbind, fi_rows)
  if (is.null(frequent)) {
    frequent <- data.frame(itemset = character(0), size = integer(0),
                           support = numeric(0), count = numeric(0),
                           stringsAsFactors = FALSE)
  }

  structure(list(
    rules              = tidy_rules,
    frequent           = frequent,
    frequent_itemsets   = frequent_itemsets,
    items              = items,
    n_transactions     = n_trans,
    n_rules            = nrow(tidy_rules),
    params             = list(min_support = min_support,
                              min_confidence = min_confidence,
                              min_lift = min_lift,
                              max_length = max_length)
  ), class = "net_association_rules")
}


# ---- Input Parsing ----

#' @noRd
# Columns of 0 and 1 could be items or recoded positions: no reading is guessed.
#' @noRd
.ar_stop_undeclared_binary <- function() {
  stop(errorCondition(
    paste0("association_rules(): these columns hold only 0 and 1, so they could be ",
           "items or positions. Set `format = \"binary\"` if they are items."),
    class = "nestimate_bad_input", call = NULL))
}

.ar_parse_input <- function(x, format = "auto") {
  # mcml -> netobject_group -> use first element
  if (inherits(x, "mcml")) {
    x <- as_tna(x)
    if (inherits(x, "netobject_group")) x <- x[[1]]
  }
  if (inherits(x, "netobject_group")) {
    stop(errorCondition(paste(
      "association_rules(): a netobject_group holds one network per group;",
      "pass one group's network."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (inherits(x, "cograph_network") && !inherits(x, "netobject")) {
    x <- .as_netobject(x)
  }

  # tna object: extract $data, decode integer labels if needed
  if (inherits(x, "tna")) {
    if (!is.null(x$data)) {
      df <- as.data.frame(x$data, stringsAsFactors = FALSE)
      labels <- rownames(x$weights)
      # Decode integer-encoded data to state names
      if (!is.null(labels) && length(labels) > 0 &&
          (is.integer(df[[1]]) || is.numeric(df[[1]]))) {
        df[] <- lapply(df, function(col) {
          idx <- as.integer(col)
          ifelse(is.na(idx) | idx < 1L | idx > length(labels),
                 NA_character_, labels[idx])
        })
      }
      transactions <- lapply(seq_len(nrow(df)), function(i) {
        vals <- as.character(unlist(df[i, ], use.names = FALSE))
        unique(vals[!is.na(vals) & vals != ""])
      })
      return(.ar_transactions_to_matrix(transactions))
    }
    stop(errorCondition("tna object has no $data. Build from sequence data.",
                        class = "nestimate_bad_input", call = NULL))
  }

  # netobject: sequences -> transactions (unique states per row)
  if (inherits(x, "netobject")) {
    if (is.null(x$data)) {
      stop(errorCondition("netobject has no $data. Build from sequence data.",
                          class = "nestimate_bad_input", call = NULL))
    }
    if (.is_onehot_params(x$params)) return(.ar_onehot_input(x))
    if (isTRUE(x$method %in% .association_methods)) {
      stop(errorCondition(sprintf(paste(
        "association_rules(): this \"%s\" network holds variables, not",
        "transactions. Pass its data with `format = \"binary\"` if the",
        "variables are 0/1 items."), x$method),
        class = "nestimate_bad_input", call = NULL))
    }
    df <- as.data.frame(x$data, stringsAsFactors = FALSE)
    # .sequence_rows() reads factor columns by their labels; unlist() of a
    # mixed row would turn a factor into its integer codes
    transactions <- lapply(.sequence_rows(df), function(vals) {
      unique(vals[!is.na(vals) & vals != ""])
    })
    return(.ar_transactions_to_matrix(transactions))
  }

  # List of character vectors
  if (is.list(x) && !is.data.frame(x)) {
    # A classed list (packageDescription, htest, ...) is a result, not
    # transactions, even when every field is character.
    if (!is.null(oldClass(x)) || !all(vapply(x, is.character, logical(1L)))) {
      stop(errorCondition(paste(
        "association_rules(): a list of transactions must be a plain list",
        "holding one character vector of items each."),
        class = "nestimate_bad_input", call = NULL))
    }
    return(.ar_transactions_to_matrix(x))
  }

  # Matrix: a character matrix holds states (wide sequences), so it takes the
  # data-frame path; a numeric or logical matrix is a binary item matrix
  if (is.matrix(x)) {
    if (is.character(x)) {
      return(.ar_parse_input(as.data.frame(x, stringsAsFactors = FALSE), format))
    }
    if (identical(format, "auto") && all(x[!is.na(x)] %in% c(0, 1))) {
      .ar_stop_undeclared_binary()
    }
    return(.ar_binary_input(x, colnames(x) %||% paste0("I", seq_len(ncol(x)))))
  }

  # Data frame
  if (is.data.frame(x)) {
    # Binary (0/1 numeric or logical) columns are items
    all_binary_type <- all(vapply(x, \(col) is.numeric(col) || is.logical(col),
                                  logical(1)))
    if (all_binary_type && all(unlist(x) %in% c(0, 1, NA))) {
      if (identical(format, "auto")) .ar_stop_undeclared_binary()
      return(.ar_binary_input(as.matrix(x), colnames(x)))
    }
    if (identical(format, "binary")) {
      stop(errorCondition(
        "association_rules(): `format = \"binary\"` needs a frame or matrix of only 0 and 1.",
        class = "nestimate_bad_input", call = NULL))
    }
    .refuse_undeclared_roles(x, "association_rules()", via_network = TRUE,
                             whole_frame = TRUE)
    # Wide character data: each row is a transaction
    # .sequence_rows() reads factor columns by their labels; unlist() of a
    # mixed row would turn a factor into its integer codes
    transactions <- lapply(.sequence_rows(x), function(vals) {
      unique(vals[!is.na(vals) & vals != ""])
    })
    return(.ar_transactions_to_matrix(transactions))
  }

  stop(errorCondition("x must be a netobject, list, data.frame, or matrix.",
                      class = "nestimate_bad_input", call = NULL))
}


# One-hot network: the codes are the items and each window is a transaction,
# the unit the network itself counts co-occurrence in.
#' @noRd
.ar_onehot_input <- function(x) {
  data <- .resampling_transition_data(x$data, x$metadata, x$params)
  codes <- .param_get(x$params, "codes")
  actor <- .param_get(x$params, "actor")
  size <- as.integer(.param_get(x$params, "window_size", 1L))
  mode <- .param_get(x$params, "mode", "non-overlapping")
  X <- .wtna_code_matrix(data, codes)
  units <- if (is.null(actor)) list(seq_len(nrow(X))) else
    .split_first_seen(.observed_group_id(data[actor], context = "actor/session"))
  windows <- lapply(units, function(rows) {
    .onehot_window_presence(X[rows, , drop = FALSE], size, mode)
  })
  .ar_binary_input(do.call(rbind, windows), codes)
}

# One logical row per window: which codes are active anywhere in it. The windows
# are those of .wtna_cooccurrence().
#' @noRd
.onehot_window_presence <- function(X, size, mode) {
  n <- nrow(X)
  if (size <= 1L) return(X > 0L)
  if (mode == "non-overlapping") {
    return(unname(rowsum(X, (seq_len(n) - 1L) %/% size, reorder = FALSE) > 0L))
  }
  starts <- seq_len(max(n - size + 1L, 0L))
  matrix(vapply(starts, function(s) {
    colSums(X[s:(s + size - 1L), , drop = FALSE]) > 0L
  }, logical(ncol(X))), nrow = length(starts), ncol = ncol(X), byrow = TRUE)
}

#' @noRd
.ar_binary_input <- function(mat, items) {
  # A missing cell is neither present nor absent; reading it as either would
  # change every support count, so the caller must decide
  n_missing <- sum(is.na(mat))
  if (n_missing > 0L) {
    stop(errorCondition(
      sprintf(paste0("binary transaction input `x` has %d missing cell(s); ",
                     "each cell must be 1/TRUE (item present) or 0/FALSE ",
                     "(absent)."), n_missing),
      class = "nestimate_bad_input", call = NULL))
  }
  if (!all(mat %in% c(0, 1))) {
    bad <- utils::head(setdiff(unique(as.vector(mat)), c(0, 1)), 3L)
    stop(errorCondition(
      sprintf(paste0("a numeric transaction matrix `x` must be binary (0/1); ",
                     "found %s. Pass `x > 0` with `format = \"binary\"` to read ",
                     "counts as presence, or a character matrix of states."),
              paste(bad, collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }
  colnames(mat) <- items
  ord <- order(items)
  list(matrix = (mat > 0)[, ord, drop = FALSE], items = items[ord])
}


#' @noRd
.ar_transactions_to_matrix <- function(transactions) {
  items <- sort(unique(unlist(transactions, use.names = FALSE)))
  n_trans <- length(transactions)
  n_items <- length(items)

  # Vectorized construction
  row_idx <- rep(seq_along(transactions), lengths(transactions))
  col_idx <- match(unlist(transactions, use.names = FALSE), items)

  mat <- matrix(FALSE, n_trans, n_items, dimnames = list(NULL, items))
  mat[cbind(row_idx, col_idx)] <- TRUE

  list(matrix = mat, items = items)
}


# ---- Candidate Generation ----

#' @noRd
.ar_generate_candidates <- function(prev_items, k) {
  # Standard Apriori join: two (k-1)-itemsets sharing first k-2 items
  n <- length(prev_items)
  candidates <- lapply(seq_len(n - 1L), function(i) {
    a <- prev_items[[i]]
    prefix_a <- a[seq_len(k - 2L)]
    joined <- lapply(seq(i + 1L, n), function(j) {
      b <- prev_items[[j]]
      prefix_b <- b[seq_len(k - 2L)]
      if (identical(prefix_a, prefix_b) && a[k - 1L] < b[k - 1L]) {
        c(a, b[k - 1L])
      }
    })
    Filter(Negate(is.null), joined)
  })
  unlist(candidates, recursive = FALSE) %||% list()
}


# ---- Rule Generation ----

#' @noRd
.ar_generate_rules <- function(frequent_itemsets, item_counts, items,
                                n_trans, min_confidence, min_lift) {
  # Pre-build support lookup for all frequent itemsets
  all_sets <- unlist(frequent_itemsets, recursive = FALSE)
  support_lookup <- list2env(
    stats::setNames(
      lapply(all_sets, function(fi) fi$support),
      vapply(all_sets, function(fi) paste(fi$items, collapse = "\t"),
             character(1L))),
    hash = TRUE, parent = emptyenv())

  # Item-level support for lift/conviction
  item_support <- setNames(item_counts / n_trans, items)

  rules_list <- list()

  # Generate rules from itemsets of size >= 2
  if (length(frequent_itemsets) < 2L) {
    return(data.frame(
      antecedent = I(list()), consequent = I(list()),
      support = numeric(0), confidence = numeric(0),
      lift = numeric(0), conviction = numeric(0),
      count = integer(0), n_transactions = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  # Loops on purpose: every itemset is split into each antecedent and
  # consequent in turn, with four early exits; the nesting is the algorithm.
  for (k in seq(2L, length(frequent_itemsets))) {
    level <- frequent_itemsets[[k]]
    if (is.null(level)) next
    for (fi in level) {
      itemset <- fi$items
      sup_ab <- fi$support
      cnt_ab <- fi$count
      n_items_k <- length(itemset)

      # Generate all non-empty proper subsets as antecedents
      for (size in seq_len(n_items_k - 1L)) {
        subsets <- combn(n_items_k, size, simplify = FALSE)
        for (idx in subsets) {
          ante <- itemset[idx]
          cons <- itemset[-idx]

          # Antecedent support
          ante_key <- paste(ante, collapse = "\t")
          sup_a <- support_lookup[[ante_key]]
          if (is.null(sup_a)) next

          # Confidence = P(A,B) / P(A)
          confidence <- sup_ab / sup_a
          if (confidence < min_confidence) next

          # Consequent support (for lift/conviction)
          cons_key <- paste(cons, collapse = "\t")
          sup_b <- support_lookup[[cons_key]]
          if (is.null(sup_b)) {
            # Single-item consequent
            if (length(cons) == 1L) {
              sup_b <- item_support[cons]
            } else {
              next
            }
          }

          # Lift = P(A,B) / (P(A) * P(B))
          lift <- sup_ab / (sup_a * sup_b)
          if (lift < min_lift) next

          # Conviction = (1 - P(B)) / (1 - confidence)
          conviction <- if (confidence >= 1) Inf
                        else (1 - sup_b) / (1 - confidence)

          rules_list <- c(rules_list, list(data.frame(
            antecedent     = I(list(ante)),
            consequent     = I(list(cons)),
            support        = sup_ab,
            confidence     = confidence,
            lift           = lift,
            conviction     = conviction,
            count          = cnt_ab,
            n_transactions = n_trans,
            stringsAsFactors = FALSE
          )))
        }
      }
    }
  }

  if (length(rules_list) == 0) {
    return(data.frame(
      antecedent = I(list()), consequent = I(list()),
      support = numeric(0), confidence = numeric(0),
      lift = numeric(0), conviction = numeric(0),
      count = integer(0), n_transactions = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  rules <- do.call(rbind, rules_list)
  rules <- rules[order(-rules$lift, -rules$confidence), , drop = FALSE]
  rownames(rules) <- NULL
  rules
}


#' @noRd
.ar_empty_result <- function(items, n_trans, min_support, min_confidence,
                              min_lift, max_length) {
  empty_rules <- data.frame(
    antecedent = character(0), consequent = character(0),
    support = numeric(0), confidence = numeric(0),
    lift = numeric(0), conviction = numeric(0),
    count = integer(0), n_transactions = integer(0),
    stringsAsFactors = FALSE
  )
  empty_freq <- data.frame(
    itemset = character(0), size = integer(0),
    support = numeric(0), count = numeric(0),
    stringsAsFactors = FALSE
  )
  structure(list(
    rules = empty_rules,
    frequent = empty_freq,
    frequent_itemsets = list(),
    items = items,
    n_transactions = n_trans,
    n_rules = 0L,
    params = list(min_support = min_support, min_confidence = min_confidence,
                  min_lift = min_lift, max_length = max_length)
  ), class = "net_association_rules")
}


# ---- S3 Methods ----

#' @rdname association_rules
#' @param ... For the `print()`, `summary()` and `plot()` methods: ignored.
#' @param object For the `summary()` method: a `net_association_rules`
#'   object.
#' @return The `print()` method returns its input invisibly. `summary()`
#'   returns the rules data frame described above. `plot()` prints a scatter
#'   plot of support against confidence, with point size and colour showing
#'   lift, and returns the `ggplot` object invisibly; with no rules it says so
#'   in a message and returns `NULL` invisibly.
#' @export
print.net_association_rules <- function(x, ...) {
  cat(sprintf("Association Rules  [%d rules | %d items | %d transactions]\n",
              x$n_rules, length(x$items), x$n_transactions))
  cat(sprintf("  Support >= %.2f  |  Confidence >= %.2f  |  Lift >= %.2f\n",
              x$params$min_support, x$params$min_confidence, x$params$min_lift))

  if (x$n_rules > 0) {
    top <- utils::head(x$rules, 10)
    cat("\n  Top rules (by lift):\n")
    cat(sprintf("    %d. %s -> %s  (sup=%.3f conf=%.3f lift=%.2f)\n",
                seq_len(nrow(top)), top$antecedent, top$consequent,
                top$support, top$confidence, top$lift), sep = "")
    if (x$n_rules > 10) {
      cat(sprintf("    ... and %d more rules\n", x$n_rules - 10))
    }
  }
  invisible(x)
}


#' @rdname association_rules
#' @export
summary.net_association_rules <- function(object, ...) {
  r <- object$rules
  row.names(r) <- NULL
  r
}


#' @rdname association_rules
#' @import ggplot2
#' @export
plot.net_association_rules <- function(x, ...) {
  r <- x$rules
  if (nrow(r) == 0) {
    message("No rules to plot.")
    return(invisible(NULL))
  }

  df <- data.frame(
    support    = r$support,
    confidence = r$confidence,
    lift       = r$lift,
    rule       = paste0("{", r$antecedent, "} => {", r$consequent, "}"),
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = support, y = confidence,
                                         size = lift, color = lift)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::scale_color_gradient(low = "#2E7D32", high = "#C62828") +
    ggplot2::scale_size_continuous(range = c(2, 8)) +
    ggplot2::labs(x = "Support", y = "Confidence",
                  title = sprintf("Association Rules (%d rules)", nrow(r)),
                  size = "Lift", color = "Lift") +
    ggplot2::theme_minimal(base_size = 12)

  print(p)
  invisible(p)
}
