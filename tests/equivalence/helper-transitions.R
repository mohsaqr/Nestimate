# Reference computations for the transition-network tests, written from the definitions in base R. Nothing here
# calls Nestimate: each function is one side of a comparison, the package is the other.

# Random sequences in wide format: one actor per row, states drawn uniformly from the first n_states letters,
# columns T1..Tn. With ragged = TRUE each row keeps a random length (at least 2) and is padded with NA.
tn_simulate <- function(n_actors, n_states, seq_length, seed, ragged = FALSE) {
  set.seed(seed)
  states <- LETTERS[seq_len(n_states)]
  rows <- lapply(seq_len(n_actors), \(i) {
    row <- sample(states, seq_length, replace = TRUE)
    keep <- if (ragged) sample(2:seq_length, 1L) else seq_length
    c(row[seq_len(keep)], rep(NA_character_, seq_length - keep))
  })
  wide <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(wide) <- paste0("T", seq_len(seq_length))
  wide
}

# The sequences of a wide table: a list of character vectors with missing values dropped.
tn_sequences <- function(wide) {
  lapply(seq_len(nrow(wide)), \(row) as.character(stats::na.omit(unlist(wide[row, ], use.names = FALSE))))
}

# Counts of consecutive pairs (from row, to column) over all sequences, on the sorted alphabet.
tn_pair_counts <- function(sequences, states = sort(unique(unlist(sequences)))) {
  pairs <- do.call(rbind, lapply(sequences, \(s) {
    n <- length(s)
    if (n < 2L) return(NULL)
    data.frame(from = s[-n], to = s[-1L])
  }))
  unclass(table(factor(pairs$from, levels = states), factor(pairs$to, levels = states), dnn = NULL))
}

# Row-normalised counts: the maximum likelihood transition probabilities; a row without transitions stays 0.
tn_row_normalise <- function(counts) {
  totals <- rowSums(counts)
  counts / ifelse(totals > 0, totals, 1)
}

# Proportion of sequences that start in each state.
tn_initial <- function(sequences, states = sort(unique(unlist(sequences)))) {
  first <- vapply(Filter(length, sequences), \(s) s[1L], character(1L))
  counts <- table(factor(first, levels = states), dnn = NULL)
  stats::setNames(as.numeric(counts) / sum(counts), states)
}

# Symmetric co-occurrence counts: every unordered pair of positions (i < j) within a sequence adds one to
# [a, b] and [b, a] (once to the diagonal when a == b), the definition of tna's ctna() without windows.
tn_cooccurrence <- function(sequences, states = sort(unique(unlist(sequences)))) {
  pairs <- do.call(rbind, lapply(sequences, \(s) {
    n <- length(s)
    if (n < 2L) return(NULL)
    index <- utils::combn(n, 2L)
    data.frame(a = s[index[1L, ]], b = s[index[2L, ]])
  }))
  counts <- unclass(table(factor(pairs$a, levels = states), factor(pairs$b, levels = states), dnn = NULL))
  counts + t(counts) - diag(diag(counts), nrow(counts))
}

# A long event log of wide sequences: actor, time (position), action.
tn_long <- function(wide) {
  sequences <- tn_sequences(wide)
  do.call(rbind, lapply(seq_along(sequences), \(i) {
    data.frame(actor = i, time = seq_along(sequences[[i]]), action = sequences[[i]])
  }))
}

# Attention-weighted transitions: every ordered pair of positions i < j within a sequence adds
# exp(-(j - i) / lambda) to [s_i, s_j], the decay weighting of tna's atna().
tn_attention <- function(sequences, lambda = 1, states = sort(unique(unlist(sequences)))) {
  pairs <- do.call(rbind, lapply(sequences, \(s) {
    n <- length(s)
    if (n < 2L) return(NULL)
    index <- utils::combn(n, 2L)
    data.frame(a = s[index[1L, ]], b = s[index[2L, ]], w = exp(-(index[2L, ] - index[1L, ]) / lambda))
  }))
  sums <- tapply(pairs$w, list(factor(pairs$a, levels = states), factor(pairs$b, levels = states)), sum)
  sums[is.na(sums)] <- 0
  unclass(sums)
}
