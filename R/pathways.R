# ---- Pathway Extraction for Higher-Order Visualization ----

#' Extract Pathways from Higher-Order Network Objects
#'
#' @description
#' Extracts higher-order pathway strings suitable for
#' \code{cograph::plot_simplicial()}. Each pathway represents a
#' multi-step dependency: source states lead to a target state.
#'
#' For \code{net_hon}: extracts edges where the source node is
#' higher-order (order > 1), i.e., the transitions that differ from
#' first-order Markov.
#'
#' For \code{net_hypa}: extracts anomalous paths (over- or
#' under-represented relative to the hypergeometric null model).
#'
#' For \code{net_mogen}: extracts all transitions at the optimal order
#' (or a specified order).
#'
#' @param x A higher-order network object (\code{net_hon},
#'   \code{net_hypa}, or \code{net_mogen}).
#' @param ... Additional arguments.
#'
#' @return A character vector of pathway strings in arrow notation
#'   (e.g. \code{"A B -> C"}), suitable for
#'   \code{cograph::plot_simplicial()}.
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


#' @describeIn pathways Extract higher-order pathways from HON
#'
#' @param min_count Integer. Minimum transition count to include
#'   (default: 1). Filters noise from rare observations.
#' @param min_prob Numeric. Minimum transition probability to include
#'   (default: 0). Useful for filtering weak transitions.
#' @param top Integer or NULL. Return only the top N pathways ranked
#'   by count (default: NULL = all).
#' @param order Integer or NULL. If specified, only include pathways
#'   at this source order. Default: all orders > 1.
#'
#' @return A character vector of pathway strings.
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
  # path column is already "A -> B -> C" — split into states
  vapply(ho$path, function(p) {
    parts <- trimws(strsplit(p, "->", fixed = TRUE)[[1]])
    if (length(parts) < 2L) return(p) # nocov
    sources <- paste(parts[-length(parts)], collapse = " ")
    paste(sources, "->", parts[length(parts)])
  }, character(1), USE.NAMES = FALSE)
}


#' @describeIn pathways Extract anomalous pathways from HYPA
#'
#' @param type Character. Which anomalies to include: \code{"all"}
#'   (default), \code{"over"}, or \code{"under"}.
#'
#' @return A character vector of pathway strings.
#'
#' @export
pathways.net_hypa <- function(x, type = "all", ...) {
  type <- match.arg(type, c("all", "over", "under"))
  scores <- x$scores
  if (is.null(scores) || nrow(scores) == 0L) return(character(0)) # nocov

  if (type == "all") {
    anom <- scores[scores$anomaly != "normal", , drop = FALSE]
  } else {
    anom <- scores[scores$anomaly == type, , drop = FALSE]
  }
  if (nrow(anom) == 0L) return(character(0))

  vapply(anom$path, function(p) {
    parts <- trimws(strsplit(p, "->", fixed = TRUE)[[1]])
    if (length(parts) < 2L) return(p) # nocov
    sources <- paste(parts[-length(parts)], collapse = " ")
    paste(sources, "->", parts[length(parts)])
  }, character(1), USE.NAMES = FALSE)
}


#' @describeIn pathways Extract pathways from a netobject
#'
#' Builds a Higher-Order Network (HON) from the netobject's sequence data
#' and returns the higher-order pathways. Requires that the netobject was
#' built from sequence data (has \code{$data}).
#'
#' @param ho_method Character. Higher-order method: \code{"hon"} (default) or
#'   \code{"hypa"}.
#'
#' @return A character vector of pathway strings.
#'
#' @export
pathways.netobject <- function(x, ho_method = c("hon", "hypa"), ...) {
  ho_method <- match.arg(ho_method)
  if (ho_method == "hon") {
    pathways(build_hon(x), ...)
  } else {
    pathways(build_hypa(x), ...)
  }
}


#' @describeIn pathways Extract pathways from association rules
#'
#' Converts association rules \code{{A, B} => {C}} into pathway strings
#' \code{"A B -> C"} suitable for \code{cograph::plot_simplicial()}.
#' Antecedent items become source nodes; consequent items become the target.
#'
#' @param top Integer or NULL. Return only the top N rules ranked by lift
#'   (default: NULL = all).
#' @param min_lift Numeric or NULL. Additional lift filter applied on top of
#'   the object's original threshold (default: NULL).
#' @param min_confidence Numeric or NULL. Additional confidence filter
#'   (default: NULL).
#'
#' @return A character vector of pathway strings.
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

  if (!is.null(top) && nrow(rules) > top) {
    rules <- rules[seq_len(top), , drop = FALSE]
  }

  vapply(seq_len(nrow(rules)), function(i) {
    ante <- paste(rules$antecedent[[i]], collapse = " ")
    cons <- paste(rules$consequent[[i]], collapse = " ")
    paste(ante, "->", cons)
  }, character(1), USE.NAMES = FALSE)
}


#' @describeIn pathways Extract pathways from link predictions
#'
#' Converts predicted links into pathway strings for
#' \code{cograph::plot_simplicial()}. When \code{evidence = TRUE}
#' (default), each predicted edge \code{A -> B} is enriched with common
#' neighbors that structurally support the prediction, producing
#' \code{"A cn1 cn2 -> B"}.
#'
#' @param method Character or NULL. Which prediction method to use.
#'   Default: first method in the object.
#' @param top Integer. Number of top predictions to include (default: 10).
#' @param evidence Logical. If TRUE, include common neighbor evidence
#'   nodes in each pathway. Default: TRUE.
#' @param max_evidence Integer. Maximum number of evidence nodes per
#'   pathway (default: 3).
#'
#' @return A character vector of pathway strings.
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


#' @describeIn pathways Extract transition pathways from MOGen
#'
#' @param order Integer or NULL. Markov order to extract. Default:
#'   optimal order from model selection.
#' @param min_count Integer. Minimum transition count to include
#'   (default: 1).
#' @param min_prob Numeric. Minimum transition probability to include
#'   (default: 0).
#' @param top Integer or NULL. Return only the top N pathways ranked
#'   by count (default: NULL = all).
#'
#' @return A character vector of pathway strings.
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
