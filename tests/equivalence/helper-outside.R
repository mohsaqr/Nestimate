# Readers of outside reference software used by several files: TraMineR state sequences and their state
# frequencies, tna's networks as edge tables.

# A TraMineR state-sequence object of a wide table; trailing missing cells are the end of the sequence.
out_seqdef <- function(wide, states = sort(unique(stats::na.omit(unlist(wide, use.names = FALSE))))) {
  suppressMessages(TraMineR::seqdef(as.data.frame(wide), alphabet = states, right = "DEL", cpal = NULL))
}

# TraMineR::seqstatf() of a wide table as one row per state: count and proportion of all events.
out_state_table <- function(wide, states = sort(unique(stats::na.omit(unlist(wide, use.names = FALSE))))) {
  table <- TraMineR::seqstatf(out_seqdef(wide, states))
  data.frame(state = states, count = as.integer(table$Freq), proportion = table$Percent / 100)
}

# A weight matrix as one row per cell (from, to, weight), every cell or the non-zero ones only.
out_edges <- function(weights, nonzero = FALSE) {
  weights <- unclass(weights)
  edges <- data.frame(from = rep(rownames(weights), times = ncol(weights)),
                      to = rep(colnames(weights), each = nrow(weights)), weight = as.vector(weights))
  if (nonzero) subset(edges, weight != 0) else edges
}

# Draws of the row-wise Dirichlet posterior of a matrix of concentration parameters by gtools::rdirichlet(): an
# array draws x rows x columns.
out_dirichlet_draws <- function(alpha, draws, seed) {
  set.seed(seed)
  rows <- lapply(seq_len(nrow(alpha)), \(i) gtools::rdirichlet(draws, alpha[i, ]))
  aperm(simplify2array(rows), c(1L, 3L, 2L))
}

# A list of sequences as a wide table, NA after each one's end.
out_pad <- function(sequences) {
  width <- max(lengths(sequences))
  as.data.frame(do.call(rbind, lapply(sequences, \(s) c(s, rep(NA_character_, width - length(s))))))
}
