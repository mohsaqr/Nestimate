# Inputs and reference computations shared by the sequence-analysis, clustering and mixture-model equivalence
# tests, written in base R. Nothing here calls Nestimate.

#' Random wide sequences: `n` rows of `length` states drawn from the first `n_states` letters
#'
#' @param n,length Rows and columns.
#' @param n_states Size of the alphabet.
#' @param seed Seed of the draw.
#' @return A data frame of characters, columns T1..T<length>.
sq_random <- function(n, length, n_states, seed) {
  set.seed(seed)
  data <- as.data.frame(matrix(sample(LETTERS[seq_len(n_states)], n * length, replace = TRUE), n, length),
                        stringsAsFactors = FALSE)
  names(data) <- paste0("T", seq_len(length))
  data
}

#' A distance object as one row per unordered pair of sequences (i < j)
#'
#' @param distance A `dist` or a square matrix.
#' @return A data frame `i`, `j`, `distance`.
sq_pairs <- function(distance) {
  m <- as.matrix(distance)
  index <- which(upper.tri(m), arr.ind = TRUE)
  data.frame(i = index[, "row"], j = index[, "col"], distance = m[index])
}

#' Transition probabilities of wide sequences: consecutive observed pairs counted and row-normalised
#'
#' A pair is counted when both positions hold a state of `nodes`; rows without a pair stay 0.
#' @param wide A data frame of states (NA for missing).
#' @param nodes The states kept, in order.
#' @return A matrix `nodes` by `nodes`.
sq_transitions <- function(wide, nodes) {
  from <- unlist(wide[-ncol(wide)], use.names = FALSE)
  to <- unlist(wide[-1L], use.names = FALSE)
  keep <- from %in% nodes & to %in% nodes
  counts <- table(factor(from[keep], nodes), factor(to[keep], nodes))
  m <- matrix(as.numeric(counts), length(nodes), dimnames = list(nodes, nodes))
  totals <- rowSums(m)
  m / ifelse(totals == 0, 1, totals)
}

#' Wide sequences with every state recoded to its cluster
#'
#' @param wide A data frame of states.
#' @param clusters A named list of states per cluster.
#' @return A data frame of cluster names.
sq_recode <- function(wide, clusters) {
  lookup <- stats::setNames(rep(names(clusters), lengths(clusters)), unlist(clusters, use.names = FALSE))
  out <- as.data.frame(lapply(wide, \(column) unname(lookup[column])), stringsAsFactors = FALSE)
  names(out) <- names(wide)
  out
}

#' Wide sequences with the states outside `nodes` set to NA
#'
#' @param wide A data frame of states.
#' @param nodes The states kept.
#' @return A data frame.
sq_mask <- function(wide, nodes) {
  as.data.frame(lapply(wide, \(column) ifelse(column %in% nodes, column, NA_character_)), stringsAsFactors = FALSE)
}

#' A matrix as one row per cell: `from` (row name), `to` (column name), `value`
#'
#' @param m A matrix with dimnames.
#' @param name The name of the value column.
#' @return A data frame.
sq_cells <- function(m, name = "weight") {
  out <- data.frame(from = rownames(m)[as.vector(row(m))], to = colnames(m)[as.vector(col(m))], value = as.vector(m))
  names(out)[[3L]] <- name
  out
}

#' seqHMM's mixture Markov model at the parameters of a fitted mixture (no fitting): one hidden state per
#' observed state, the identity emission, and the mixing proportions as multinomial log-odds against cluster 1
#'
#' @param fit A fitted mixture: `states`, `mixing`, and per cluster `models[[k]]$weights` and `$initial`.
#' @param data The wide sequences it was fitted to.
#' @return A seqHMM `mmm` object.
sq_seqhmm <- function(fit, data) {
  states <- fit$states
  sequences <- suppressMessages(TraMineR::seqdef(data, alphabet = states))
  transitions <- lapply(fit$models, \(model) model$weights[states, states])
  initials <- lapply(fit$models, \(model) model$initial[states])
  seqHMM::build_mmm(sequences, transition_probs = transitions, initial_probs = initials,
                    coefficients = matrix(log(fit$mixing / fit$mixing[[1L]]), nrow = 1L))
}

#' The posterior cluster probabilities of a seqHMM mixture as a matrix, sequences by clusters
#'
#' @param model A seqHMM `mmm` object.
#' @return A numeric matrix.
sq_posterior <- function(model) {
  long <- as.data.frame(seqHMM::posterior_cluster_probabilities(model))
  matrix(long$probability, ncol = model$n_clusters, byrow = TRUE)
}

#' Quality statistics of a mixture from its log-likelihood and posterior (definitions as in Latent GOLD and
#' Biernacki, Celeux and Govaert 2000)
#'
#' @param loglik Log-likelihood.
#' @param posterior Sequences by clusters.
#' @param n_states Number of states.
#' @return A list: n_params, AIC, BIC, ICL (BIC + 2 sum -log of each sequence's modal posterior), avepp per cluster,
#'   AvePP (mean modal posterior), Entropy (sum -p log p over N log k), assignments.
sq_mixture_quality <- function(loglik, posterior, n_states) {
  k <- ncol(posterior)
  n <- nrow(posterior)
  n_params <- (k - 1) + k * (n_states - 1) + k * n_states * (n_states - 1)
  assignments <- max.col(posterior, ties.method = "first")
  modal <- posterior[cbind(seq_len(n), assignments)]
  bic <- -2 * loglik + n_params * log(n)
  list(n_params = n_params, AIC = -2 * loglik + 2 * n_params, BIC = bic, ICL = bic - 2 * sum(log(modal)),
       avepp = vapply(seq_len(k), \(cluster) mean(posterior[assignments == cluster, cluster]), numeric(1)),
       AvePP = mean(modal), Entropy = -sum(posterior * log(posterior)) / (n * log(k)), assignments = assignments)
}

#' Sequences from two groups with different state frequencies (an identifiable two-cluster mixture)
#'
#' @param n_per_group Rows per group.
#' @param length Columns.
#' @param seed Seed of the draw.
#' @return A data frame of characters.
sq_two_groups <- function(n_per_group, length, seed) {
  set.seed(seed)
  draw <- \(prob) matrix(sample(c("A", "B", "C"), n_per_group * length, replace = TRUE, prob = prob), n_per_group)
  data <- as.data.frame(rbind(draw(c(0.7, 0.2, 0.1)), draw(c(0.1, 0.2, 0.7))), stringsAsFactors = FALSE)
  names(data) <- paste0("T", seq_len(length))
  data
}

#' Contiguous n-grams of wide sequences (trailing and internal NA ends a run of observed states)
#'
#' @param wide A data frame of states, one sequence per row.
#' @param len The n-gram length.
#' @return A data frame: `sequence` (row number), `pattern` ("a->b->c"), one row per occurrence.
sq_ngrams <- function(wide, len) {
  rows <- lapply(seq_len(nrow(wide)), \(i) {
    states <- unlist(wide[i, ], use.names = FALSE)
    starts <- seq_len(max(0L, length(states) - len + 1L))
    grams <- vapply(starts, \(s) {
      window <- states[s:(s + len - 1L)]
      if (anyNA(window)) NA_character_ else paste(window, collapse = "->")
    }, character(1))
    grams <- grams[!is.na(grams)]
    data.frame(sequence = rep(i, length(grams)), pattern = grams)
  })
  do.call(rbind, rows)
}

#' Counts of each n-gram: occurrences, sequences containing it, support, share of all n-grams of that length
#'
#' @param grams The result of sq_ngrams().
#' @param n_sequences Number of sequences.
#' @return A data frame `pattern`, `occurrences`, `n_sequences`, `support`, `proportion`.
sq_pattern_counts <- function(grams, n_sequences) {
  occurrences <- table(grams$pattern)
  containing <- table(unique(grams)$pattern)
  data.frame(pattern = names(occurrences), occurrences = as.vector(occurrences),
             n_sequences = as.vector(containing[names(occurrences)]),
             support = as.vector(containing[names(occurrences)]) / n_sequences,
             proportion = as.vector(occurrences) / nrow(grams))
}
