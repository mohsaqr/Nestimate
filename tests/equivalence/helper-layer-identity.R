# Helpers of the composite-identity (Layer 9) and invariance (Layer 5) tests. A composite route (a wrapper, a
# grouped fit, the networks of a clustering, a resampling verb's observed network) must return exactly what the
# direct build_network() fit returns on the same rows with the same settings; a refit after a transformation of
# the input (rows permuted, states renamed, sequences duplicated, variables rescaled) must return the original
# carried through that transformation. Both sides are computed by Nestimate, so every reference is `internal`.
# These helpers only take objects apart for expect_reference(); they compute nothing of their own.

#' The reference name of an identity or invariance check
#'
#' @param what The route or transformation, in a few words.
#' @return A length-one character vector starting `internal: `.
li_ref <- function(what) paste0("internal: ", what)

#' The edge table of a network keyed by node labels
#'
#' An undirected edge is written with its two labels in alphabetical order, so that the key does not depend on
#' the order of the nodes.
#'
#' @param net A `netobject`.
#' @return A data frame `from`, `to` (labels) and `weight`, one row per edge.
li_edges <- function(net) {
  labels <- net$nodes$label %||% rownames(net$weights)
  edges <- net$edges
  from <- labels[edges$from]
  to <- labels[edges$to]
  if (!isTRUE(net$directed)) {
    first <- pmin(from, to)
    to <- pmax(from, to)
    from <- first
  }
  data.frame(from = from, to = to, weight = edges$weight, stringsAsFactors = FALSE)
}

#' The settings of a network as one text value each
#'
#' @param net A `netobject`.
#' @param drop Parameter names left out (bookkeeping a route adds, such as the shared group alphabet).
#' @return A named character vector: `method`, `directed`, `n_nodes`, `n_edges`, `threshold`, `scaling`, and
#'   one `params$<name>` entry per stored estimator parameter.
li_settings <- function(net, drop = character()) {
  params <- net$params %||% list()
  kept <- sort(setdiff(as.character(names(params)), drop))
  stored <- stats::setNames(vapply(params[kept], \(value) paste(deparse(value), collapse = " "), character(1L)),
                            paste0("params$", kept, recycle0 = TRUE))
  c(method = as.character(net$method), directed = as.character(isTRUE(net$directed)),
    n_nodes = as.character(net$n_nodes), n_edges = as.character(net$n_edges),
    threshold = format(net$threshold %||% 0, digits = 17L),
    scaling = paste(net$scaling %||% "none", collapse = "+"),
    stored)
}

#' The initial probabilities of a network, or a marker when it has none
#'
#' @param net A `netobject`.
#' @return A named numeric vector, or the text `"absent"` keyed `initial`.
li_initial <- function(net) {
  initial <- net$initial %||% net$inits
  if (is.null(initial)) return(c(initial = "absent"))
  stats::setNames(as.numeric(initial), names(initial))
}

#' Assert that two networks are the same network, part by part
#'
#' Weights (keyed by row and column label), edges (keyed by from and to), initial probabilities and settings
#' are each one expect_reference() with tolerance 0. A part named in `defect` is a declared defect.
#'
#' @param actual The network of the composite route.
#' @param expected The direct fit.
#' @param route What the composite route is, for the reference name.
#' @param covers The exported functions checked.
#' @param quantity What is compared, in the reader's words; the part's name is appended.
#' @param parts The parts compared.
#' @param drop_params Parameter names left out of the settings.
#' @param defect A named list: part -> one-line defect description.
#' @param tolerance The largest difference allowed; 0 unless a solver's precision is the reason.
#' @return `actual`, invisibly.
li_expect_same_network <- function(actual, expected, route, covers, quantity,
                                   parts = c("weights", "edges", "initial", "settings"),
                                   drop_params = character(), defect = list(), tolerance = 0) {
  reference <- li_ref(route)
  check <- function(part, a, e, key = NULL) {
    expect_reference(a, e, reference = reference, covers = covers, quantity = sprintf("%s: %s", quantity, part),
                     tolerance = tolerance, key = key, authority = "internal", defect = defect[[part]])
  }
  if ("weights" %in% parts) check("weights", actual$weights, expected$weights)
  if ("edges" %in% parts) {
    edges_actual <- li_edges(actual)
    edges_expected <- li_edges(expected)
    # a network without edges has an empty table, which holds no keyed cell to compare: compare the counts
    if (nrow(edges_actual) == 0L || nrow(edges_expected) == 0L) {
      check("edges", c(n_edges = nrow(edges_actual)), c(n_edges = nrow(edges_expected)))
    } else {
      check("edges", edges_actual, edges_expected, key = c("from", "to"))
    }
  }
  if ("initial" %in% parts) check("initial", li_initial(actual), li_initial(expected))
  if ("settings" %in% parts) {
    check("settings", li_settings(actual, drop_params), li_settings(expected, drop_params))
  }
  invisible(actual)
}

#' Relabel the rows and columns of a matrix, or the names of a vector, through a map
#'
#' @param x A matrix with dimnames, or a named vector.
#' @param map A named character vector: old name -> new name. Names not in `map` stay.
#' @return `x` with its labels replaced.
li_relabel <- function(x, map) {
  swap <- function(labels) ifelse(labels %in% names(map), unname(map[labels]), labels)
  if (is.matrix(x)) {
    dimnames(x) <- lapply(dimnames(x), swap)
    return(x)
  }
  stats::setNames(x, swap(names(x)))
}

#' Apply a state map to every cell of a wide sequence table
#'
#' @param wide A data frame of state columns.
#' @param map A named character vector: old state -> new state.
#' @return The table with every state replaced; missing cells stay missing.
li_rename_states <- function(wide, map) {
  out <- wide
  out[] <- lapply(wide, \(column) unname(map[as.character(column)]))
  out
}

#' The inverse of a map
#'
#' @param map A named character vector, one-to-one.
#' @return The map from its values back to its names.
li_invert <- function(map) stats::setNames(names(map), unname(map))

#' Rewrite the states inside pattern strings (`"A->B->C"`) through a map
#'
#' @param patterns A character vector of patterns.
#' @param map A named character vector: state -> new state.
#' @param sep The separator between states.
#' @return The patterns with each state replaced.
li_map_patterns <- function(patterns, map, sep = "->") {
  vapply(strsplit(patterns, sep, fixed = TRUE),
         \(states) paste(ifelse(states %in% names(map), unname(map[states]), states), collapse = sep),
         character(1L))
}

#' Two-regime sequences: half the rows drawn from each of two transition matrices
#'
#' Gives a mixture model a well-separated optimum, so that a refit after a transformation reaches the same one.
#'
#' @param n_per Rows per regime.
#' @param length Columns.
#' @param seed Seed of the draw.
#' @return A data frame of characters over A-D, columns T1..T<length>.
li_two_regimes <- function(n_per, length, seed) {
  set.seed(seed)
  states <- LETTERS[1:4]
  stay <- matrix(0.05, 4L, 4L) + diag(0.75, 4L)
  cycle <- matrix(0.1, 4L, 4L)
  cycle[cbind(1:4, c(2:4, 1L))] <- 0.7
  draw <- function(P) {
    Reduce(\(previous, step) sample(states, 1L, prob = P[match(previous, states), ]),
           seq_len(length - 1L), accumulate = TRUE, init = sample(states, 1L))
  }
  rows <- c(replicate(n_per, draw(stay), simplify = FALSE), replicate(n_per, draw(cycle), simplify = FALSE))
  wide <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(wide) <- paste0("T", seq_len(length))
  wide
}

#' The order of clusters of a refit that matches the original: for each original cluster, the refit cluster
#' sharing most of its members
#'
#' @param original,refit Assignment vectors of the same units in the same order.
#' @return An integer vector: position `j` is the refit cluster matching original cluster `j`.
li_match_clusters <- function(original, refit) {
  k <- max(original)
  vapply(seq_len(k), \(cluster) {
    hits <- tabulate(refit[original == cluster], nbins = max(refit))
    which.max(hits)
  }, integer(1L))
}
