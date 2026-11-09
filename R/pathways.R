# ---- Pathway Extraction for Higher-Order Visualization ----

#' Extract Pathways from Higher-Order Network Objects
#'
#' @description
#' `pathways()` turns a higher-order result into pathway strings in arrow
#' notation, `"A B -> C"`: the states before the arrow lead, together, to the
#' state after it. The strings are what `cograph::plot_simplicial()` draws,
#' one simplex per pathway. Which pathways are returned depends on the input:
#' \itemize{
#'   \item a `net_hon` gives its higher-order edges, those whose source node
#'     has order above 1, i.e. the transitions that differ from a first-order
#'     Markov chain, ranked by count;
#'   \item a `net_hypa` gives its anomalous paths, over- or under-represented
#'     relative to the hypergeometric null model, the over-represented ones
#'     first;
#'   \item a `net_mogen` gives every transition at the optimal order (or the
#'     order given), ranked by count;
#'   \item a `netobject` built from sequences is first turned into a HON or a
#'     HYPA, and that result's pathways are returned;
#'   \item a `net_association_rules` gives its rules `{A, B} => {C}` as
#'     `"A B -> C"`: antecedent items become the sources, consequent items the
#'     target. Rules whose antecedent and consequent together hold the same
#'     items describe the same simplex, so only the highest-lift rule per item
#'     set is kept; the rules are ranked by lift, then confidence;
#'   \item a `net_link_prediction` gives its predicted links, ranked by score;
#'     with `evidence = TRUE` each predicted edge `A -> B` is extended with the
#'     common neighbours that support it, `"A cn1 cn2 -> B"`.
#' }
#'
#' The strings keep only which pathways passed the filters, not the counts,
#' probabilities, scores or tests behind them; those stay in the source
#' object (for a MOGen, [mogen_transitions()] gives them per order).
#'
#' @param x The object to read: a `net_hon`, `net_hypa` or `net_mogen` (from
#'   [build_hon()], [build_hypa()] or [build_mogen()]), a `netobject` built
#'   from sequences, a `net_association_rules` (from [association_rules()])
#'   or a `net_link_prediction` (from [predict_links()]). A `netobject` of an
#'   association method (whose rows are variables) or a one-hot `netobject`
#'   raises an error of class `nestimate_bad_input`; a network without its
#'   sequences stops `pathways()`.
#' @param ... For the `netobject` method: passed to the `pathways()` method of
#'   the HON or HYPA built from it, so `min_count`, `min_prob`, `top` and
#'   `order` apply to a HON and `type` to a HYPA. Ignored by the other
#'   methods.
#'
#' @return A character vector of pathway strings in arrow notation (e.g.
#'   `"A B -> C"`), one per pathway and in the ranking order described above,
#'   ready for `cograph::plot_simplicial()`. Every method returns this shape,
#'   and `character(0)` when nothing passes its filters.
#'
#' @examples
#' \donttest{
#' seqs <- list(c("A","B","C","D"), c("A","B","C","A"))
#' hon <- build_hon(seqs, max_order = 3)
#' pw <- pathways(hon)
#' }
#'
#' @export
pathways <- function(x, ...) {
  UseMethod("pathways")
}


#' @describeIn pathways Higher-order edges of a HON.
#'
#' @param min_count For a HON or a MOGen: the smallest transition count
#'   kept, which filters out rare observations. Default `1` (all).
#' @param min_prob For a HON or a MOGen: the smallest transition probability
#'   kept. Default `0` (all).
#' @param top How many pathways to keep after ranking: by count for a HON or
#'   a MOGen, by lift then confidence for association rules (after the item
#'   sets are de-duplicated), and by score for link predictions. Default
#'   `NULL` (all), except for link predictions, whose default is `10`. Not
#'   used for a HYPA.
#' @param order For a HON: keep only pathways whose source has this order.
#'   Default `NULL`, every order above 1; order 1 returns nothing. For a
#'   MOGen: the Markov order to extract. Default `NULL`, the optimal order
#'   from model selection; an order below 1 returns `character(0)`.
#'
#' @export
pathways.net_hon <- function(x, min_count = 1L, min_prob = 0,
                             top = NULL, order = NULL, ...) {
  edges <- x$ho_edges
  if (is.null(edges) || nrow(edges) == 0L) return(character(0)) # nocov

  # Higher-order edges: from_order > 1
  ho <- edges[edges$from_order > 1L, , drop = FALSE]
  if (!is.null(order)) {
    ho <- ho[ho$from_order == order, , drop = FALSE]
  }
  if (min_count > 1L) {
    ho <- ho[ho$count >= min_count, , drop = FALSE]
  }
  if (min_prob > 0) {
    ho <- ho[ho$probability >= min_prob, , drop = FALSE]
  }
  # Rank by count descending
  ho <- ho[order(-ho$count), , drop = FALSE]
  if (!is.null(top) && nrow(ho) > top) {
    ho <- ho[seq_len(top), , drop = FALSE]
  }
  if (nrow(ho) == 0L) return(character(0))

  # Convert "A -> B -> C" path format to "A B -> C" for plot_simplicial
  # path column is already "A -> B -> C" - split into states
  vapply(ho$path, function(p) {
    parts <- trimws(strsplit(p, "->", fixed = TRUE)[[1]])
    if (length(parts) < 2L) return(p) # nocov
    sources <- paste(parts[-length(parts)], collapse = " ")
    paste(sources, "->", parts[length(parts)])
  }, character(1), USE.NAMES = FALSE)
}


#' @describeIn pathways Anomalous paths of a HYPA.
#'
#' @param type For a HYPA: which anomalies to include, `"all"` (default),
#'   `"over"` (over-represented paths) or `"under"` (under-represented ones).
#'
#' @export
pathways.net_hypa <- function(x, type = "all", ...) {
  type <- .ho_match_arg(type, "type", c("all", "over", "under"))
  scores <- x$scores
  if (is.null(scores) || nrow(scores) == 0L) return(character(0)) # nocov

  if (type == "all") {
    anom <- rbind(x$over, x$under)
  } else if (type == "over") {
    anom <- x$over
  } else {
    anom <- x$under
  }
  if (nrow(anom) == 0L) return(character(0))

  vapply(anom$path, function(p) {
    parts <- trimws(strsplit(p, "->", fixed = TRUE)[[1]])
    if (length(parts) < 2L) return(p) # nocov
    sources <- paste(parts[-length(parts)], collapse = " ")
    paste(sources, "->", parts[length(parts)])
  }, character(1), USE.NAMES = FALSE)
}


#' @describeIn pathways Higher-order pathways of a network's sequences.
#'
#' @param ho_method For a `netobject`: the higher-order model built from its
#'   sequences, `"hon"` (default, [build_hon()]) or `"hypa"`
#'   ([build_hypa()]), each with its default settings. The network must carry
#'   the sequences it was built from.
#'
#' @export
pathways.netobject <- function(x, ho_method = c("hon", "hypa"), ...) {
  ho_method <- .ho_match_arg(ho_method, "ho_method")
  if (ho_method == "hon") {
    pathways(build_hon(x), ...)
  } else {
    pathways(build_hypa(x), ...)
  }
}


#' @describeIn pathways Association rules as pathways.
#'
#' @param min_lift For association rules: the smallest lift kept, on top of
#'   the threshold the rules were mined with. Default `NULL` (no further
#'   filter).
#' @param min_confidence For association rules: the smallest confidence
#'   kept, on top of the mining threshold. Default `NULL` (no further
#'   filter).
#'
#' @examples
#' trans <- list(c("A","B","C"), c("A","B"), c("B","C","D"), c("A","C","D"))
#' rules <- association_rules(trans, min_support = 0.3, min_confidence = 0.3,
#'                            min_lift = 0)
#' pathways(rules)
#'
#' @export
pathways.net_association_rules <- function(x, top = NULL, min_lift = NULL,
                                           min_confidence = NULL, ...) {
  rules <- x$rules
  if (nrow(rules) == 0L) return(character(0))

  if (!is.null(min_lift)) {
    rules <- rules[rules$lift >= min_lift, , drop = FALSE]
  }
  if (!is.null(min_confidence)) {
    rules <- rules[rules$confidence >= min_confidence, , drop = FALSE]
  }
  if (nrow(rules) == 0L) return(character(0))

  rules <- rules[order(-rules$lift, -rules$confidence), , drop = FALSE]

  # Deduplicate: each rule's antecedent + consequent is a simplex (unordered).
  # {A,B} => {C} and {A,C} => {B} are the same simplex {A,B,C}.
  # Keep the highest-lift version of each unique itemset.
  itemsets <- vapply(seq_len(nrow(rules)), function(i) {
    items <- sort(unique(c(
      strsplit(rules$antecedent[i], ", ", fixed = TRUE)[[1]],
      strsplit(rules$consequent[i], ", ", fixed = TRUE)[[1]]
    )))
    paste(items, collapse = "\t")
  }, character(1))
  keep <- !duplicated(itemsets)
  rules <- rules[keep, , drop = FALSE]

  if (!is.null(top) && nrow(rules) > top) {
    rules <- rules[seq_len(top), , drop = FALSE]
  }

  vapply(seq_len(nrow(rules)), function(i) {
    ante <- gsub(", ", " ", rules$antecedent[i])
    cons <- gsub(", ", " ", rules$consequent[i])
    paste(ante, "->", cons)
  }, character(1), USE.NAMES = FALSE)
}


#' @describeIn pathways Predicted links as pathways.
#'
#' @param method For link predictions: the prediction method whose links are
#'   returned. Default `NULL`, the first method in the object. A method the
#'   object does not hold returns `character(0)`.
#' @param evidence For link predictions: whether to extend each predicted
#'   link `A -> B` with its common neighbours, the nodes with an edge from
#'   `A` and an edge to `B`, giving `"A cn1 cn2 -> B"`. Default `TRUE`.
#'   Without a stored adjacency matrix the links are returned as `"A -> B"`.
#' @param max_evidence For link predictions with `evidence = TRUE`: the most
#'   common neighbours added to one link; when there are more, those with
#'   the largest summed weight of the edges from `A` and to `B` are kept.
#'   Default `3`.
#'
#' @examples
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:5], 50, TRUE),
#'   V2 = sample(LETTERS[1:5], 50, TRUE),
#'   V3 = sample(LETTERS[1:5], 50, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' pred <- predict_links(net, methods = "common_neighbors")
#' pathways(pred)
#'
#' @export
pathways.net_link_prediction <- function(x, method = NULL, top = 10L,
                                          evidence = TRUE, max_evidence = 3L,
                                          ...) {
  if (is.null(method)) method <- x$methods[1]
  df <- x$predictions[x$predictions$method == method, , drop = FALSE]
  df <- df[order(-df$score), , drop = FALSE]

  if (!is.null(top) && nrow(df) > top) {
    df <- df[seq_len(top), , drop = FALSE]
  }
  if (nrow(df) == 0L) return(character(0))

  # Simple mode: just "from -> to"
  if (!isTRUE(evidence) || is.null(x$adjacency)) {
    return(paste(df$from, "->", df$to))
  }

  # Evidence mode: include common neighbors as structural bridge
  A <- x$adjacency
  nodes <- x$nodes

  vapply(seq_len(nrow(df)), function(i) {
    from_idx <- match(df$from[i], nodes)
    to_idx <- match(df$to[i], nodes)

    # Common neighbors: reachable from source AND reaching target
    from_out <- A[from_idx, ] > 0
    to_in <- A[, to_idx] > 0
    cn_mask <- from_out & to_in
    cn_mask[from_idx] <- FALSE
    cn_mask[to_idx] <- FALSE
    cn_indices <- which(cn_mask)
    cn_nodes <- nodes[cn_indices]

    if (length(cn_nodes) > max_evidence) {
      # Rank by combined weight to source and target
      cn_weights <- A[from_idx, cn_indices] + A[cn_indices, to_idx]
      cn_nodes <- cn_nodes[order(-cn_weights)][seq_len(max_evidence)]
    }

    sources <- c(df$from[i], cn_nodes)
    paste(paste(sources, collapse = " "), "->", df$to[i])
  }, character(1), USE.NAMES = FALSE)
}


#' @describeIn pathways Transitions of a MOGen at one order.
#'
#' @export
pathways.net_mogen <- function(x, order = NULL, min_count = 1L,
                               min_prob = 0, top = NULL, ...) {
  if (is.null(order)) order <- x$optimal_order
  if (order < 1L) return(character(0))

  trans <- mogen_transitions(x, order = order)
  if (is.null(trans) || nrow(trans) == 0L) return(character(0)) # nocov

  if (min_count > 1L) {
    trans <- trans[trans$count >= min_count, , drop = FALSE] # nocov
  }
  if (min_prob > 0) {
    trans <- trans[trans$probability >= min_prob, , drop = FALSE]
  }
  trans <- trans[order(-trans$count), , drop = FALSE]
  if (!is.null(top) && nrow(trans) > top) {
    trans <- trans[seq_len(top), , drop = FALSE] # nocov
  }
  if (nrow(trans) == 0L) return(character(0)) # nocov

  vapply(trans$path, function(p) {
    parts <- trimws(strsplit(p, "->", fixed = TRUE)[[1]])
    if (length(parts) < 2L) return(p) # nocov
    sources <- paste(parts[-length(parts)], collapse = " ")
    paste(sources, "->", parts[length(parts)])
  }, character(1), USE.NAMES = FALSE)
}
