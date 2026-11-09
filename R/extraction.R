#' Extract Transition Matrix from Model
#'
#' @description
#' `extract_transition_matrix()` returns the weight matrix of a network: for
#' a transition network, the transition probabilities or counts between
#' states, with rows as the states moved from and columns as the states
#' moved to. It reads a `netobject` (or any list with a `weights` element), a
#' `tna`, `ftna`, `ctna` or `atna` model, an `mcml` (one matrix per layer) or
#' a bare matrix, and can rescale the rows to sum to 1.
#'
#' @param model The network: a `netobject`, a `tna`-family model, an `mcml`,
#'   a list holding the matrix as `weights` (or, failing that,
#'   `transition_matrix` or `transition`), or a weight matrix. Anything else
#'   raises an error.
#' @param type Which matrix to return. `"raw"` (default) returns the weights
#'   as stored. `"scaled"` divides each row by its sum, so the rows of a
#'   count matrix become transition probabilities; a row that sums to zero
#'   is left as it is.
#'
#' @return For a single network: a square numeric matrix with the model's
#'   state names (if any) as row and column names, of class
#'   `c("nest_transition_matrix", "matrix", "array")`. It behaves as an
#'   ordinary matrix in arithmetic and indexing; `summary()`
#'   returns a data frame with one row per non-zero entry (self-transitions
#'   included) and the columns `from`, `to` and `weight`, sorted by
#'   decreasing absolute weight.
#'
#'   For an `mcml`: a named list of such matrices, `macro` first and then
#'   one per cluster.
#'
#' @details
#' Models store their weights in different places; the function looks in
#' each in turn. `"scaled"` is useful when the stored weights are counts or
#' otherwise do not sum to 1 per row.
#'
#' @examples
#' seqs <- data.frame(V1 = c("A","B","A"), V2 = c("B","A","C"), V3 = c("A","C","B"))
#' net <- build_network(seqs, method = "relative")
#' trans_mat <- extract_transition_matrix(net)
#' print(trans_mat)
#'
#' @seealso \code{\link{extract_initial_probs}} for extracting initial probabilities,
#'   \code{\link{extract_edges}} for extracting an edge list.
#'
#' @export
extract_transition_matrix <- function(model, type = c("raw", "scaled")) {
  type <- match.arg(type)

  # mcml dispatch: return named list of matrices per layer
  if (inherits(model, "mcml")) {
    result <- list(macro = extract_transition_matrix(model$macro, type = type))
    if (!is.null(model$clusters)) {
      # Loop kept: `[[<-` by cluster name, so a cluster called "macro"
      # replaces the macro entry exactly as before.
      for (nm in names(model$clusters)) {
        result[[nm]] <- extract_transition_matrix(model$clusters[[nm]], type = type)
      }
    }
    return(result)
  }

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
    stop(errorCondition(
      paste0(
        "Could not extract transition matrix from model. ",
        "Expected a TNA model object, a list with 'weights', or a matrix."),
      class = "nestimate_bad_input", call = NULL))
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

  class(weights) <- c("nest_transition_matrix", class(weights))
  weights
}


#' @rdname extract_transition_matrix
#' @param x For the `print()` method: a `nest_transition_matrix`.
#' @export
print.nest_transition_matrix <- function(x, ...) {
  # the plain matrix, without a trailing attr(,"class") block
  print(unclass(x), ...)
  invisible(x)
}


#' @rdname extract_transition_matrix
#' @param object For the `summary()` method: a `nest_transition_matrix`.
#' @param ... For the `summary()` method: ignored. For the `print()` method:
#'   passed on to the matrix print.
#' @export
summary.nest_transition_matrix <- function(object, ...) {
  .matrix_to_long_df(unclass(object), value_col = "weight")
}


#' Extract Initial Probabilities from Model
#'
#' @description
#' `extract_initial_probs()` returns the initial-state probabilities of a
#' network: the share of sequences that start in each state. It reads a
#' `netobject` (or any list) that stores them as `initial`, `initial_probs`,
#' `inits` or `initial_probabilities`, a `tna`-family object that stores
#' them as `initial`, `initial_probs` or `inits` (models from the tna
#' package store `inits`), or an `mcml` (one vector per layer).
#'
#' Of the `build_network()` methods, the transition methods (`"relative"`,
#' `"frequency"`, `"attention"`, `"ngram"`, `"gap"` and `"reverse"`) store
#' initial probabilities, as do [wtna()] transition networks. For
#' co-occurrence and association networks, and for a model that keeps them
#' under another name, the function returns a uniform distribution over the
#' states of the weight matrix, with a warning: those values are a
#' placeholder, not an estimate.
#'
#' @param model The network: a `netobject`, a `tna`-family model, an `mcml`,
#'   or a list holding an `initial` element. Input from which neither
#'   initial probabilities nor a weight matrix can be read raises an error.
#'
#' @return For a single network: a named numeric vector of initial-state
#'   probabilities that sums to 1 (stored values are rescaled when their sum
#'   differs from 1 by more than 1e-6), of class `c("nest_initial_probs",
#'   "numeric")`. `summary()` returns a data frame with one row per state and
#'   the columns `state` and `prob`, sorted by decreasing probability.
#'
#'   For an `mcml`: a named list of plain numeric vectors, `macro` first and
#'   then one per cluster, taken as stored in each layer (`inits`).
#'
#' @examples
#' seqs <- data.frame(V1 = c("A","B","A"), V2 = c("B","A","C"), V3 = c("A","C","B"))
#' net <- build_network(seqs, method = "relative")
#' init_probs <- extract_initial_probs(net)
#' print(init_probs)
#'
#' @seealso \code{\link{extract_transition_matrix}} for extracting the transition matrix,
#'   \code{\link{extract_edges}} for extracting an edge list.
#'
#' @export
extract_initial_probs <- function(model) {
  # mcml dispatch: return named list of init vectors per layer
  if (inherits(model, "mcml")) {
    result <- list(macro = model$macro$inits)
    if (!is.null(model$clusters)) {
      # Loop kept: `[[<-` skips a cluster without `$inits`; lapply() would
      # add a NULL entry for it.
      for (nm in names(model$clusters)) {
        result[[nm]] <- model$clusters[[nm]]$inits
      }
    }
    return(result)
  }

  initial <- NULL

  # Try different possible locations
  if (inherits(model, "tna") || inherits(model, "ftna") ||
      inherits(model, "ctna") || inherits(model, "atna")) {
    # Models from the tna package store the initial probabilities as $inits.
    initial <- model$initial %||% model$initial_probs %||% model$inits
  } else if (is.list(model)) {
    initial <- model$initial %||% model$initial_probs %||% model$inits %||%
      model$initial_probabilities
  }

  if (is.null(initial)) {
    # Try to infer from transition matrix
    # only "no weights in this object" falls through to the error below
    weights <- tryCatch(
      extract_transition_matrix(model),
      nestimate_bad_input = function(e) NULL
    )
    if (!is.null(weights)) {
      n_states <- nrow(weights)
      states <- rownames(weights)
      if (is.null(states)) states <- paste0("S", seq_len(n_states))
      # Use uniform distribution
      initial <- rep(1 / n_states, n_states)
      names(initial) <- states
      warning(warningCondition(
        "Initial probabilities not found in model. Using uniform distribution.",
        class = "nestimate_uniform_initial", call = NULL))
    } else {
      stop(errorCondition("Could not extract initial probabilities from model.",
                          class = "nestimate_bad_input", call = NULL))
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

  class(initial) <- c("nest_initial_probs", class(initial))
  initial
}


#' @rdname extract_initial_probs
#' @param x For the `print()` method: a `nest_initial_probs`.
#' @export
print.nest_initial_probs <- function(x, ...) {
  # the plain named vector, without a trailing attr(,"class") block
  print(unclass(x), ...)
  invisible(x)
}


#' @rdname extract_initial_probs
#' @param object For the `summary()` method: a `nest_initial_probs`.
#' @param ... For the `summary()` method: ignored. For the `print()` method:
#'   passed on to the vector print.
#' @export
summary.nest_initial_probs <- function(object, ...) {
  v <- unclass(object)
  states <- names(v)
  if (is.null(states)) states <- paste0("S", seq_along(v))
  df <- data.frame(
    state = states,
    prob  = as.numeric(v),
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
  df[order(-df$prob), ]
}


#' Extract Edge List with Weights
#'
#' @description
#' `extract_edges()` returns the edges of a network as a data frame, one row
#' per edge with its weight, read from the weight matrix (rows are the
#' `from` states, columns the `to` states). It reads every input
#' [extract_transition_matrix()] reads; an `mcml` gives one edge list per
#' layer.
#'
#' An edge is a cell whose absolute weight is above `threshold`, so with
#' the default `threshold = 0` every non-zero cell is an edge, negative
#' weights (in correlation, partial-correlation and mlVAR networks)
#' included, and zero cells are not. A directed network lists every
#' ordered pair; an undirected one lists each edge once, with `from` the
#' earlier state in the matrix order. A network's own `directed` field
#' decides which; a `tna` model is directed unless it is a co-occurrence
#' model, and a bare matrix is undirected when it is symmetric.
#'
#' @param model The network: a `netobject`, a `tna`-family model, an `mcml`
#'   or a weight matrix, as for [extract_transition_matrix()].
#' @param threshold Cells whose absolute weight is at most `threshold` are
#'   dropped. A single finite number of at least 0, otherwise an error is
#'   raised (class `nestimate_bad_input` for a negative value); default `0`.
#' @param include_self Whether to keep self-loops (a state followed by
#'   itself), `TRUE` or `FALSE`. Default `FALSE`.
#' @param sort_by The row order: `"weight"` (default) by decreasing weight,
#'   `"from"` by source then target state, `"to"` by target then source
#'   state, or `NULL` for the order of the weight matrix (column by column).
#'   Any other value raises an error.
#'
#' @return For a single network: a data frame with one row per edge and
#'   the columns:
#' \describe{
#'   \item{from}{The source state.}
#'   \item{to}{The target state.}
#'   \item{weight}{The edge weight (for a `"relative"` network, the
#'     transition probability).}
#' }
#'
#'   For an `mcml`: a named list of such data frames, `macro` first and then
#'   one per cluster.
#'
#' @details
#' An edge list is the form most other network tools (igraph, export to
#' Gephi or a spreadsheet) take.
#'
#' @examples
#' seqs <- data.frame(V1 = c("A","B","A"), V2 = c("B","A","C"), V3 = c("A","C","B"))
#' net <- build_network(seqs, method = "relative")
#' edges <- extract_edges(net, threshold = 0.05)
#' head(edges)
#'
#' @seealso \code{\link{extract_transition_matrix}} for the full matrix,
#'   \code{\link{build_network}} for network estimation.
#'
#' @note The \pkg{idiographic} package exports a function with the same
#'   name. With both packages attached, the one attached last masks the
#'   other; write `Nestimate::extract_edges()` to call this one.
#' @export
extract_edges <- function(model,
                          threshold = 0,
                          include_self = FALSE,
                          sort_by = "weight") {
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold)) {
    stop(errorCondition("'threshold' must be a single finite numeric value.",
                        class = "nestimate_bad_input", call = NULL))
  }
  if (threshold < 0) {
    stop(errorCondition(
      "'threshold' must be a single finite numeric value of at least 0.",
      class = "nestimate_bad_input", call = NULL))
  }
  if (!is.logical(include_self) || length(include_self) != 1L ||
      is.na(include_self)) {
    stop(errorCondition("'include_self' must be TRUE or FALSE.",
                        class = "nestimate_bad_input", call = NULL))
  }
  if (!is.null(sort_by) &&
      (!is.character(sort_by) || length(sort_by) != 1L ||
       is.na(sort_by) ||
       !sort_by %in% c("weight", "from", "to"))) {
    stop(errorCondition(
      "'sort_by' must be one of 'weight', 'from', 'to', or NULL.",
      class = "nestimate_bad_input", call = NULL))
  }

  # mcml dispatch: return named list of edge data frames per layer
  if (inherits(model, "mcml")) {
    result <- list(macro = extract_edges(model$macro, threshold = threshold,
                                         include_self = include_self,
                                         sort_by = sort_by))
    if (!is.null(model$clusters)) {
      # Loop kept: `[[<-` by cluster name (see extract_transition_matrix()).
      for (nm in names(model$clusters)) {
        result[[nm]] <- extract_edges(model$clusters[[nm]],
                                       threshold = threshold,
                                       include_self = include_self,
                                       sort_by = sort_by)
      }
    }
    return(result)
  }

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

  # An edge is a non-zero cell: |weight| > threshold keeps negative edges
  # (correlation, partial-correlation, mlVAR networks) and drops zero
  # cells. An undirected network lists each edge once, from the upper
  # triangle (row <= column), as .extract_edges_from_matrix() does.
  keep <- abs(edges$weight) > threshold
  if (!.edges_directed(model, weights)) {
    keep <- keep & as.vector(row(weights) <= col(weights))
  }
  edges <- edges[keep, , drop = FALSE]

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

#' Whether a model's edges are directed
#'
#' A declared `$directed` wins (netobject and its subclasses); a tna model
#' is directed unless its type is co-occurrence; a bare matrix is
#' undirected only when it is symmetric.
#' @noRd
.edges_directed <- function(model, weights) {
  declared <- if (is.list(model)) model$directed else NULL
  if (is.logical(declared) && length(declared) == 1L && !is.na(declared)) {
    return(declared)
  }
  if (inherits(model, c("tna", "ftna", "ctna", "atna"))) {
    return(!identical(attr(model, "type"), "co-occurrence"))
  }
  w <- unclass(weights)
  !isSymmetric(unname(w))
}
