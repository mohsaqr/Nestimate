# Inputs and reference computations shared by the higher-order, path, Markov-chain and entropy equivalence tests,
# written in base R (markovchain where named). Nothing here calls Nestimate.

#' A row-stochastic matrix with strictly positive entries, states s1..sn
#'
#' @param n_states Number of states.
#' @param seed Seed of the draw.
#' @return A named square matrix whose rows sum to 1.
ho_stochastic <- function(n_states, seed) {
  set.seed(seed)
  weights <- matrix(stats::runif(n_states^2, 0.05, 1), n_states, n_states)
  states <- paste0("s", seq_len(n_states))
  matrix(weights / rowSums(weights), n_states, n_states, dimnames = list(states, states))
}

#' markovchain's view of a transition matrix
#'
#' @param matrix_p A named row-stochastic matrix.
#' @return A markovchain object.
ho_chain <- function(matrix_p) {
  methods::new(methods::getClassDef("markovchain", package = "markovchain"),
               states = rownames(matrix_p), byrow = TRUE, transitionMatrix = matrix_p)
}

#' The stationary distribution of a chain by markovchain::steadyStates(), named by state
#'
#' @param matrix_p A named row-stochastic matrix.
#' @return A named numeric vector.
ho_stationary <- function(matrix_p) {
  stationary <- markovchain::steadyStates(ho_chain(matrix_p))
  stats::setNames(as.vector(stationary), colnames(stationary))
}

#' Random sequences of a first-order chain, one per row of a wide data frame
#'
#' @param n_seq,length Number and length of the sequences.
#' @param matrix_p A named row-stochastic matrix.
#' @param seed Seed of the draw.
#' @return A data frame with columns T1..T<length>.
ho_sequences <- function(n_seq, length, matrix_p, seed) {
  set.seed(seed)
  states <- rownames(matrix_p)
  rows <- lapply(seq_len(n_seq), \(i) {
    Reduce(\(previous, step) sample(states, 1L, prob = matrix_p[previous, ]), seq_len(length - 1L),
           accumulate = TRUE, init = sample(states, 1L))
  })
  wide <- as.data.frame(do.call(rbind, lapply(rows, unlist)), stringsAsFactors = FALSE)
  names(wide) <- paste0("T", seq_len(length))
  wide
}
