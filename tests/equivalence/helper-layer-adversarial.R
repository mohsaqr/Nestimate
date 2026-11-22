# Layer 6 of the verification suite: adversarial data. The inputs on which a sequence or network verb is most
# likely to fail (a missing first state, a state seen once, a state never left, one actor, two states, twenty
# states with sparse counts, numeric-looking and odd labels, a group of one, duplicated time stamps, a constant or
# near-collinear variable, more variables than rows) and the base-R references the layer-6 tests compare with.
# Nothing here calls Nestimate: each function is one side of a comparison, the package is the other.

#' A wide table of sequences from a list of character vectors, padded with NA to the longest
#'
#' @param rows A list of character vectors (NA allowed anywhere).
#' @return A data frame of characters, columns T1..T<longest>.
adv_wide <- function(rows) {
  width <- max(lengths(rows))
  padded <- lapply(rows, \(r) c(r, rep(NA_character_, width - length(r))))
  wide <- as.data.frame(do.call(rbind, padded), stringsAsFactors = FALSE)
  names(wide) <- paste0("T", seq_len(width))
  wide
}

#' The adversarial sequence sets, each a wide data frame of characters (columns T1..Tn)
#'
#' Every set is fixed by its seed, so the comparisons are reproducible.
#' @return A named list of data frames.
adv_sequence_sets <- function() {
  set.seed(6001L)
  draw <- \(n, states = c("A", "B", "C")) sample(states, n, replace = TRUE)
  base_rows <- lapply(seq_len(12L), \(i) draw(8L))
  list(
    # the first one to three positions of four rows are missing: the first observed state is not at T1
    leading_na = adv_wide(Map(\(r, gap) c(rep(NA_character_, gap), r[seq_len(8L - gap)]),
                              base_rows, c(1L, 2L, 3L, 1L, integer(8L)))),
    # runs of NA at the end of half of the rows (drop-out)
    trailing_na = adv_wide(Map(\(r, cut) c(r[seq_len(8L - cut)], rep(NA_character_, cut)),
                               base_rows, rep(c(0L, 3L, 5L), 4L))),
    # gaps inside the rows: a transition across a gap is not a transition
    internal_na = adv_wide(Map(\(r, at) replace(r, at, NA_character_), base_rows,
                               rep(list(3L, c(2L, 5L), 4:5, integer()), 3L))),
    # Q occurs exactly once, in the middle of a row
    singleton_state = adv_wide(c(list(c("A", "B", "Q", "C", "A", "B", "C", "A")), base_rows[-1L])),
    # Z is entered but never left (it ends every row it occurs in); S is left but never entered (only at T1):
    # a zero row and a zero column in the transition table
    absorbing = adv_wide(c(lapply(base_rows[1:4], \(r) c(r[1:6], "Z", NA_character_)),
                           lapply(base_rows[5:8], \(r) c("S", r[2:8])), base_rows[9:12])),
    # one row of length 1 and one of length 2 among ordinary rows
    short = adv_wide(c(list("B", c("C", "A")), base_rows[3:12])),
    # a single actor
    one_actor = adv_wide(list(c("A", "B", "B", "C", "A", "C", "B", "A", "A", "C"))),
    # every row the same
    identical = adv_wide(rep(list(c("A", "B", "C", "A", "C", "B")), 8L)),
    # an alphabet of two states
    two_states = adv_wide(lapply(seq_len(10L), \(i) draw(7L, c("A", "B")))),
    # 22 states over 10 short rows: most cells of the table are zero
    many_states = adv_wide(lapply(seq_len(10L), \(i) draw(6L, sprintf("s%02d", seq_len(22L))))),
    # states that look like numbers
    numeric_labels = adv_wide(lapply(base_rows, \(r) c(A = "1", B = "2", C = "10")[r])),
    # states with spaces and punctuation
    odd_labels = adv_wide(lapply(base_rows, \(r) c(A = "state A", B = "b-c.d", C = "(e)/f")[r])),
    # a deterministic cycle: every state has the same in- and out-strength, so the centralities tie
    ties = adv_wide(rep(list(c("A", "B", "C", "A", "B", "C", "A")), 6L)),
    # one row with no observed state at all among rows with drop-out
    empty_row = adv_wide(c(Map(\(r, cut) c(r[seq_len(8L - cut)], rep(NA_character_, cut)), base_rows[1:4],
                               c(0L, 3L, 5L, 0L)), list(rep(NA_character_, 8L)), base_rows[6:12]))
  )
}

#' Two groups of very different sizes, and a group of one, beside the sequences of `wide`
#'
#' @param wide A wide table of sequences.
#' @return A named list of data frames, each `wide` with a character column `grp`.
adv_group_sets <- function(wide) {
  n <- nrow(wide)
  list(imbalanced = cbind(grp = c(rep("big", n - 2L), "small", "small"), wide, stringsAsFactors = FALSE),
       singleton_group = cbind(grp = c(rep("many", n - 1L), "one"), wide, stringsAsFactors = FALSE))
}

#' A long event log with time stamps duplicated within actors: rows are in log order, ties are kept in that order
#'
#' @return A data frame `actor`, `time`, `action`.
adv_long_duplicated_times <- function() {
  data.frame(actor = rep(c("u1", "u2", "u3"), c(6L, 5L, 4L)),
             time = c(1, 2, 2, 3, 3, 4, 1, 1, 2, 3, 4, 5, 5, 5, 6),
             action = c("A", "B", "C", "A", "B", "C", "C", "A", "A", "B", "C", "B", "A", "C", "A"),
             stringsAsFactors = FALSE)
}

#' Adversarial association data: a constant column, two near-collinear columns, more variables than rows, and
#' binary columns (one of them nearly constant) for the Ising model
#'
#' @return A named list of numeric data frames.
adv_association_sets <- function() {
  set.seed(6002L)
  n <- 60L
  base <- matrix(stats::rnorm(n * 4L), n, 4L, dimnames = list(NULL, paste0("V", 1:4)))
  near <- base[, "V1"] + stats::rnorm(n, sd = 1e-4)
  wide_p <- matrix(stats::rnorm(8L * 12L), 8L, 12L, dimnames = list(NULL, paste0("V", 1:12)))
  latent <- stats::rnorm(200L)
  binary <- vapply(1:5, \(j) as.integer(latent + stats::rnorm(200L) > 0), integer(200L))
  colnames(binary) <- paste0("V", 1:5)
  rare <- binary
  rare[, "V5"] <- c(1L, 1L, 1L, integer(197L))
  list(constant = data.frame(base, V5 = rep(3, n)),
       near_collinear = data.frame(base, V5 = near),
       n_less_p = as.data.frame(wide_p),
       binary = as.data.frame(binary),
       binary_rare = as.data.frame(rare))
}

# ---- references ---------------------------------------------------------------------------------------------

#' The alphabet of a wide table of sequences, sorted as Nestimate and tna sort it
#'
#' @param wide A wide table of sequences.
#' @return A character vector.
adv_states <- function(wide) sort(unique(stats::na.omit(unlist(wide, use.names = FALSE))))

#' Consecutive pairs of a wide table: both positions observed, adjacent columns (a gap breaks the pair)
#'
#' @param wide A wide table of sequences.
#' @return A data frame `row`, `from`, `to`, one row per pair.
adv_pairs <- function(wide) {
  m <- as.matrix(wide)
  if (ncol(m) < 2L) return(data.frame(row = integer(), from = character(), to = character()))
  from <- m[, -ncol(m), drop = FALSE]
  to <- m[, -1L, drop = FALSE]
  keep <- !is.na(from) & !is.na(to)
  data.frame(row = row(from)[keep], from = from[keep], to = to[keep], stringsAsFactors = FALSE)
}

#' Transition counts of a wide table, rows = from, columns = to, over `states`
#'
#' @param wide A wide table of sequences.
#' @param states The alphabet.
#' @return A numeric matrix.
adv_counts <- function(wide, states = adv_states(wide)) {
  pairs <- adv_pairs(wide)
  counts <- table(factor(pairs$from, states), factor(pairs$to, states))
  matrix(as.numeric(counts), length(states), dimnames = list(states, states))
}

#' Row-normalised counts; a row without transitions stays 0
#'
#' @param counts A count matrix.
#' @return A matrix of the same shape.
adv_relative <- function(counts) {
  totals <- rowSums(counts)
  counts / ifelse(totals > 0, totals, 1)
}

#' Share of rows whose first OBSERVED state is each state (a row with no observed state is left out)
#'
#' @param wide A wide table of sequences.
#' @param states The alphabet.
#' @return A named numeric vector.
adv_initial <- function(wide, states = adv_states(wide)) {
  first <- apply(as.matrix(wide), 1L, \(r) r[!is.na(r)][1L])
  first <- first[!is.na(first)]
  counts <- table(factor(first, states))
  stats::setNames(as.numeric(counts) / sum(counts), states)
}

#' The observed states of each row, NA dropped (positions are lost)
#'
#' @param wide A wide table of sequences.
#' @return A list of character vectors.
adv_rows <- function(wide) lapply(seq_len(nrow(wide)), \(i) {
  r <- unlist(wide[i, ], use.names = FALSE)
  r[!is.na(r)]
})

#' A matrix as one row per cell: from (row name), to (column name), weight
#'
#' @param m A matrix with dimnames.
#' @return A data frame.
adv_cells <- function(m) {
  data.frame(from = rownames(m)[as.vector(row(m))], to = colnames(m)[as.vector(col(m))], weight = as.vector(m),
             stringsAsFactors = FALSE)
}

#' Every ordered pair of observed positions i < j within a row, with their distance j - i
#'
#' @param wide A wide table of sequences.
#' @return A data frame `row`, `from`, `to`, `distance`.
adv_position_pairs <- function(wide) {
  m <- as.matrix(wide)
  p <- ncol(m)
  if (p < 2L) return(data.frame(row = integer(), from = character(), to = character(), distance = integer()))
  index <- utils::combn(p, 2L)
  parts <- lapply(seq_len(ncol(index)), \(k) {
    i <- index[1L, k]
    j <- index[2L, k]
    keep <- !is.na(m[, i]) & !is.na(m[, j])
    data.frame(row = which(keep), from = m[keep, i], to = m[keep, j], distance = rep(j - i, sum(keep)),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, parts)
}

#' A weighted table of state pairs: sum of `w` over the rows of `pairs`, rows = from, columns = to
#'
#' @param pairs A data frame with `from`, `to`.
#' @param w Weights, one per row of `pairs`.
#' @param states The alphabet.
#' @return A numeric matrix.
adv_weighted_table <- function(pairs, w, states) {
  sums <- tapply(w, list(factor(pairs$from, states), factor(pairs$to, states)), sum)
  sums[is.na(sums)] <- 0
  matrix(as.numeric(sums), length(states), dimnames = list(states, states))
}

#' Co-occurrence of a wide table: every unordered pair of observed positions of a row adds one to [a, b] and
#' [b, a] (once to the diagonal when a == b), tna's ctna() definition
#'
#' @param wide A wide table of sequences.
#' @param states The alphabet.
#' @return A symmetric numeric matrix.
adv_cooccurrence <- function(wide, states = adv_states(wide)) {
  pairs <- adv_position_pairs(wide)
  counts <- adv_weighted_table(pairs, rep(1, nrow(pairs)), states)
  counts + t(counts) - diag(diag(counts), length(states))
}

#' Attention weights: every ordered pair of observed positions i < j adds exp(-(j - i) / lambda) to [s_i, s_j]
#'
#' @param wide A wide table of sequences.
#' @param lambda The decay.
#' @param states The alphabet.
#' @return A numeric matrix.
adv_attention <- function(wide, lambda = 1, states = adv_states(wide)) {
  pairs <- adv_position_pairs(wide)
  adv_weighted_table(pairs, exp(-pairs$distance / lambda), states)
}

#' Gap-allowed transitions: an ordered pair of observed positions at distance d <= max_gap + 1 adds 1 / d
#'
#' @param wide A wide table of sequences.
#' @param max_gap The largest gap.
#' @param states The alphabet.
#' @return A numeric matrix.
adv_gap <- function(wide, max_gap = 1L, states = adv_states(wide)) {
  pairs <- adv_position_pairs(wide)
  near <- pairs$distance <= max_gap + 1L
  adv_weighted_table(pairs[near, , drop = FALSE], 1 / pairs$distance[near], states)
}

#' n-gram transitions: each window of n consecutive positions adds one for every adjacent observed pair in it
#'
#' @param wide A wide table of sequences.
#' @param n_gram The window length.
#' @param states The alphabet.
#' @return A numeric matrix.
adv_ngram <- function(wide, n_gram = 2L, states = adv_states(wide)) {
  m <- as.matrix(wide)
  p <- ncol(m)
  starts <- seq_len(max(p - n_gram + 1L, 0L))
  # the adjacent pairs (j, j + 1) of every window, one entry per window that holds the pair
  pair_starts <- unlist(lapply(starts, \(s) seq(s, s + n_gram - 2L)))
  parts <- lapply(pair_starts, \(j) {
    keep <- !is.na(m[, j]) & !is.na(m[, j + 1L])
    data.frame(from = m[keep, j], to = m[keep, j + 1L], stringsAsFactors = FALSE)
  })
  pairs <- do.call(rbind, c(list(data.frame(from = character(), to = character())), parts))
  adv_weighted_table(pairs, rep(1, nrow(pairs)), states)
}

#' igraph's metrics of a weighted directed network, in the rows summary() of a netobject prints
#'
#' Self-loops count as edges, strengths and degrees; mean distance and centralization are igraph's without loops,
#' the weights read as edge lengths; density is the edge count (loops included) over n (n - 1), capped at 1
#' (Nestimate's documented definition, which igraph::edge_density() does not cap).
#' @param weights A weight matrix with dimnames.
#' @return A data frame `metric`, `value`.
adv_igraph_summary <- function(weights) {
  with_loops <- igraph::graph_from_adjacency_matrix(weights, mode = "directed", weighted = TRUE, diag = TRUE)
  without_loops <- igraph::graph_from_adjacency_matrix(weights, mode = "directed", weighted = TRUE, diag = FALSE)
  n <- nrow(weights)
  out_strength <- igraph::strength(with_loops, mode = "out")
  in_strength <- igraph::strength(with_loops, mode = "in")
  out_degree <- igraph::degree(with_loops, mode = "out")
  data.frame(
    metric = c("Node Count", "Edge Count", "Network Density", "Mean Distance", "Mean Out-Strength", "SD Out-Strength",
               "Mean In-Strength", "SD In-Strength", "Mean Out-Degree", "SD Out-Degree", "Centralization (Out-Degree)",
               "Centralization (In-Degree)", "Reciprocity"),
    value = c(n, igraph::ecount(with_loops), min(1, igraph::ecount(with_loops) / (n * (n - 1))),
              igraph::mean_distance(without_loops), mean(out_strength), stats::sd(out_strength), mean(in_strength),
              stats::sd(in_strength), mean(out_degree), stats::sd(out_degree),
              igraph::centr_degree(without_loops, mode = "out", loops = FALSE)$centralization,
              igraph::centr_degree(without_loops, mode = "in", loops = FALSE)$centralization,
              igraph::reciprocity(with_loops)),
    stringsAsFactors = FALSE)
}

#' The non-zero cells of a weight matrix as an edge list, self-loops left out unless asked for
#'
#' @param weights A weight matrix with dimnames.
#' @param include_self Keep the diagonal.
#' @param threshold Keep cells whose weight is above it.
#' @return A data frame `from`, `to`, `weight`.
adv_edge_list <- function(weights, include_self = FALSE, threshold = 0) {
  cells <- adv_cells(weights)
  keep <- cells$weight != 0 & abs(cells$weight) > threshold & (include_self | cells$from != cells$to)
  cells[keep, , drop = FALSE]
}

#' Transition counts with a start state before the first OBSERVED state of each row and an end state after its
#' last observed state (the documented boundary markers); a row with no observed state adds nothing
#'
#' @param wide A wide table of sequences.
#' @param start,end The labels of the boundary states.
#' @return A numeric matrix over c(start, alphabet, end).
adv_boundary_counts <- function(wide, start = "Start", end = "End") {
  states <- c(start, adv_states(wide), end)
  rows <- adv_rows(wide)
  inner <- adv_pairs(wide)
  seen <- Filter(length, rows)
  pairs <- rbind(inner[c("from", "to")],
                 data.frame(from = rep(start, length(seen)), to = vapply(seen, `[`, character(1), 1L)),
                 data.frame(from = vapply(seen, \(r) r[length(r)], character(1)), to = rep(end, length(seen))))
  adv_weighted_table(pairs, rep(1, nrow(pairs)), states)
}
