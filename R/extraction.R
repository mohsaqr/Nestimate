#' Model Extraction Functions
#'
#' @description
#' Functions for extracting components from TNA model objects.
#'
#' @name extraction
#' @keywords internal
NULL

#' Extract Transition Matrix from Model
#'
#' @description
#' Extract the transition probability matrix from a TNA model object.
#'
#' @param model A TNA model object or a list containing a 'weights' element.
#' @param type Character. Type of matrix to return:
#' \describe{
#'   \item{"raw"}{The raw weight matrix as stored in the model.}
#'   \item{"scaled"}{Row-normalized to ensure rows sum to 1.}
#' }
#' Default: "raw".
#'
#' @return A square numeric matrix with row and column names as state names.
#'
#' @details
#' TNA models store transition weights in different locations depending on the
#' model type. This function handles the extraction automatically.
#'
#' For "scaled" type, each row is divided by its sum to create valid transition
#' probabilities. This is useful when the original weights don't sum to 1.
#'
#' @examples
#' \donttest{
#' seqs <- data.frame(V1 = c("A","B","A"), V2 = c("B","A","C"), V3 = c("A","C","B"))
#' net <- build_network(seqs, method = "relative")
#' trans_mat <- extract_transition_matrix(net)
#' print(trans_mat)
#' }
#'
#' @seealso \code{\link{extract_initial_probs}} for extracting initial probabilities,
#'   \code{\link{extract_edges}} for extracting an edge list.
#'
#' @export
extract_transition_matrix <- function(model, type = c("raw", "scaled")) {
  type <- match.arg(type)

  # Extract weights matrix
  weights <- NULL

  # Try different possible locations
  if (inherits(model, "tna") || inherits(model, "ftna") ||
      inherits(model, "ctna") || inherits(model, "atna")) {
    # Standard TNA model objects
    if (!is.null(model$weights)) {
      weights <- model$weights
    } else if (!is.null(model$transition)) {
      weights <- model$transition
    }
  } else if (is.list(model)) {
    # Generic list
    if (!is.null(model$weights)) {
      weights <- model$weights
    } else if (!is.null(model$transition_matrix)) {
      weights <- model$transition_matrix
    } else if (!is.null(model$transition)) {
      weights <- model$transition
    }
  } else if (is.matrix(model)) {
    # Direct matrix input
    weights <- model
  }

  if (is.null(weights)) {
    stop("Could not extract transition matrix from model. ",
         "Expected a TNA model object, a list with 'weights', or a matrix.")
  }

  # Ensure it's a matrix
  weights <- as.matrix(weights)

  # Scale if requested
  if (type == "scaled") {
    row_sums <- rowSums(weights)
    # Avoid division by zero
    row_sums[row_sums == 0] <- 1
    weights <- weights / row_sums
  }

  return(weights)
}


#' Extract Initial Probabilities from Model
#'
#' @description
#' Extract the initial state probability vector from a TNA model object.
#'
#' @param model A TNA model object or a list containing an 'initial' element.
#'
#' @return A named numeric vector of initial state probabilities.
#'
#' @details
#' Initial probabilities represent the probability of starting a sequence in
#' each state. If the model doesn't have explicit initial probabilities,
#' this function attempts to estimate them from the data or use uniform
#' probabilities.
#'
#' @examples
#' \donttest{
#' seqs <- data.frame(V1 = c("A","B","A"), V2 = c("B","A","C"), V3 = c("A","C","B"))
#' net <- build_network(seqs, method = "relative")
#' init_probs <- extract_initial_probs(net)
#' print(init_probs)
#' }
#'
#' @seealso \code{\link{extract_transition_matrix}} for extracting the transition matrix,
#'   \code{\link{extract_edges}} for extracting an edge list.
#'
#' @export
extract_initial_probs <- function(model) {
  initial <- NULL

  # Try different possible locations
  if (inherits(model, "tna") || inherits(model, "ftna") ||
      inherits(model, "ctna") || inherits(model, "atna")) {
    if (!is.null(model$initial)) {
      initial <- model$initial
    } else if (!is.null(model$initial_probs)) { # nocov start
      initial <- model$initial_probs # nocov end
    }
  } else if (is.list(model)) {
    if (!is.null(model$initial)) {
      initial <- model$initial
    } else if (!is.null(model$initial_probs)) {
      initial <- model$initial_probs # nocov
    } else if (!is.null(model$initial_probabilities)) {
      initial <- model$initial_probabilities # nocov
    }
  }

  if (is.null(initial)) {
    # Try to infer from transition matrix
    weights <- tryCatch(
      extract_transition_matrix(model),
      error = function(e) NULL
    )
    if (!is.null(weights)) {
      n_states <- nrow(weights)
      states <- rownames(weights)
      if (is.null(states)) states <- paste0("S", seq_len(n_states))
      # Use uniform distribution
      initial <- rep(1 / n_states, n_states)
      names(initial) <- states
      warning("Initial probabilities not found in model. Using uniform distribution.")
    } else {
      stop("Could not extract initial probabilities from model.")
    }
  }

  # Ensure it's a named vector
  if (is.null(names(initial))) {
    names(initial) <- paste0("S", seq_along(initial))
  }

  # Normalize to sum to 1
  if (abs(sum(initial) - 1) > 1e-6) {
    initial <- initial / sum(initial)
  }

  return(initial)
}


#' Extract Edge List with Weights
#'
#' @description
#' Extract an edge list from a TNA model, representing the network as a
#' data frame of from-to-weight tuples.
#'
#' @param model A TNA model object or a matrix of weights.
#' @param threshold Numeric. Minimum weight to include an edge. Default: 0.
#' @param include_self Logical. Whether to include self-loops. Default: FALSE.
#' @param sort_by Character. Column to sort by: "weight" (descending),
#'   "from", "to", or NULL for no sorting. Default: "weight".
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{from}{Source state name.}
#'   \item{to}{Target state name.}
#'   \item{weight}{Edge weight (transition probability).}
#' }
#'
#' @details
#' This function converts the transition matrix into an edge list format,
#' which is useful for visualization, analysis with igraph, or export to
#' other network tools.
#'
#' @examples
#' \donttest{
#' seqs <- data.frame(V1 = c("A","B","A"), V2 = c("B","A","C"), V3 = c("A","C","B"))
#' net <- build_network(seqs, method = "relative")
#' edges <- extract_edges(net, threshold = 0.05)
#' head(edges)
#' }
#'
#' @seealso \code{\link{extract_transition_matrix}} for the full matrix,
#'   \code{\link{build_network}} for network estimation.
#'
#' @export
extract_edges <- function(model,
                          threshold = 0,
                          include_self = FALSE,
                          sort_by = "weight") {
  # Extract transition matrix
  weights <- extract_transition_matrix(model)

  # Get state names
  states <- rownames(weights)
  if (is.null(states)) {
    states <- paste0("S", seq_len(nrow(weights)))
  }

  # Create edge list
  edges <- expand.grid(
    from = states,
    to = states,
    stringsAsFactors = FALSE
  )
  edges$weight <- as.vector(weights)

  # Filter by threshold
  edges <- edges[edges$weight >= threshold, ]

  # Remove self-loops if requested
  if (!include_self) {
    edges <- edges[edges$from != edges$to, ]
  }

  # Sort if requested
  if (!is.null(sort_by)) {
    if (sort_by == "weight") {
      edges <- edges[order(-edges$weight), ]
    } else if (sort_by == "from") {
      edges <- edges[order(edges$from, edges$to), ]
    } else if (sort_by == "to") {
      edges <- edges[order(edges$to, edges$from), ]
    }
  }

  rownames(edges) <- NULL
  return(edges)
}
