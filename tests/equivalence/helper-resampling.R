# Reference-side helpers of the resampling, comparison, centrality and pruning tests (bootstrap_network(),
# permutation(), centrality_stability(), network_reliability(), net_centrality(), net_prune(), compare_model(),
# ...). Base R only; nothing here calls Nestimate.

#' Random sequences in wide format: one row per actor, columns T1..Tn, states drawn uniformly
#'
#' @param n_actors,n_states,length Sizes.
#' @param seed The seed of the draw (the input only; no package randomness is replayed here).
#' @return A data frame of character columns.
rs_random_sequences <- function(n_actors, n_states, length, seed) {
  set.seed(seed)
  states <- LETTERS[seq_len(n_states)]
  draws <- matrix(sample(states, n_actors * length, replace = TRUE), n_actors, length,
                  dimnames = list(NULL, paste0("T", seq_len(length))))
  as.data.frame(draws, stringsAsFactors = FALSE)
}

#' Matrices with the same dimnames as one long table: one row per cell, from = row, to = column
#'
#' @param ... Named matrices; each becomes a column.
#' @param drop_diagonal Leave the self-loops out.
#' @return A data frame with `from`, `to` and one column per matrix.
rs_edge_table <- function(..., drop_diagonal = FALSE) {
  matrices <- list(...)
  first <- matrices[[1L]]
  cells <- expand.grid(from = rownames(first), to = colnames(first), stringsAsFactors = FALSE)
  values <- lapply(matrices, \(m) m[cbind(cells$from, cells$to)])
  table <- cbind(cells, as.data.frame(values))
  if (drop_diagonal) table <- table[table$from != table$to, , drop = FALSE]
  rownames(table) <- NULL
  table
}

#' tna's "from -> to" edge names split into two key columns
#'
#' @param stats A data frame with a column `edge_name`.
#' @return The data frame with `from` and `to` and without `edge_name`.
rs_split_edge_names <- function(stats) {
  stats <- as.data.frame(stats)
  parts <- strsplit(stats$edge_name, " -> ", fixed = TRUE)
  stats$from <- vapply(parts, `[`, character(1), 1L)
  stats$to <- vapply(parts, `[`, character(1), 2L)
  stats$edge_name <- NULL
  stats
}
