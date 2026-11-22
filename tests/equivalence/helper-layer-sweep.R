# Reference-side helpers of the argument sweep (layer 8, test-layer8-sweep-*.R). Base R only: nothing here calls
# Nestimate. Each function computes one side of a comparison from the definition the documentation states; the
# package is the other side. The sweep fits every option of every argument, so the references take the options
# as plain arguments.

#' Skip a comparison whose tna definition changed in tna 1.5
#'
#' Windows on sequence data, the n-gram end column and range-normalised Diffusion are tna >= 1.5 behaviour; the
#' default library of the verify job may hold an older tna, whose definitions these comparisons do not test.
#' @return Invisibly `NULL`; skips the test when tna is older than 1.5.7.
sw_skip_old_tna <- function() {
  testthat::skip_if(utils::packageVersion("tna") < "1.5.7", "needs tna >= 1.5.7 (the definition compared here)")
}

#' Positions of wide sequences with the boundary states the documentation defines
#'
#' `start` is put before every sequence; `end` in the cell after each sequence's last observed state (a column
#' is added for it).
#' @param wide A data frame of character position columns.
#' @param start,end A label, or `NULL` for none.
#' @return A character matrix without dimnames.
sw_positions <- function(wide, start = NULL, end = NULL) {
  mat <- as.matrix(wide)
  mat[!is.na(mat) & mat == ""] <- NA
  if (!is.null(start)) mat <- cbind(start, mat)
  if (!is.null(end)) {
    last <- apply(!is.na(mat), 1L, \(observed) max(c(0L, which(observed))))
    mat <- cbind(mat, NA_character_)
    mat[cbind(which(last > 0L), last[last > 0L] + 1L)] <- end
  }
  unname(mat)
}

#' The alphabet of a position matrix: the start label, the sorted observed states, the end label
#'
#' @param mat A character matrix.
#' @param start,end The boundary labels, or `NULL`.
#' @return A character vector.
sw_states <- function(mat, start = NULL, end = NULL) {
  observed <- sort(setdiff(unique(as.vector(mat[!is.na(mat)])), c(start, end)))
  c(start, observed, end)
}

#' Weighted counts of the state pairs found at the given column pairs, both cells observed
#'
#' @param mat A character matrix.
#' @param states The alphabet (row and column order of the result).
#' @param from,to Column indices, one pair per element.
#' @param weight One weight per column pair.
#' @param row_weight One weight per row of `mat`.
#' @return A numeric `states x states` matrix with dimnames.
sw_pair_table <- function(mat, states, from, to, weight, row_weight = rep(1, nrow(mat))) {
  empty <- matrix(0, length(states), length(states), dimnames = list(states, states))
  if (!length(from)) return(empty)
  cells <- do.call(rbind, Map(\(f, t, w) data.frame(a = mat[, f], b = mat[, t], w = w * row_weight,
                                                    stringsAsFactors = FALSE), from, to, weight))
  cells <- cells[!is.na(cells$a) & !is.na(cells$b), , drop = FALSE]
  if (!nrow(cells)) return(empty)
  sums <- tapply(cells$w, list(factor(cells$a, levels = states), factor(cells$b, levels = states)), sum)
  sums[is.na(sums)] <- 0
  matrix(as.vector(sums), length(states), dimnames = list(states, states))
}

#' Row-normalised weights: each row with a positive sum divided by it
#'
#' @param m A numeric matrix.
#' @return The matrix.
sw_row_normalise <- function(m) {
  totals <- rowSums(abs(m))
  m / ifelse(totals > 0, totals, 1)
}

#' The weight matrix of a sequence method, written from the definition in ?build_network and in tna
#'
#' frequency: adjacent observed pairs; relative: frequency row-normalised; reverse: frequency transposed;
#' co_occurrence: every pair of positions i < j, both directions (the diagonal once), `type = "relative"`
#' row-normalised; attention: every ordered pair i < j weighted exp(-(j - i) / lambda) (`direction` "backward" is
#' the transpose, "both" the sum); gap: pairs 1 <= j - i <= max_gap + 1 weighted 1 / (j - i); ngram: each adjacent
#' pair once per window of n_gram consecutive columns that contains it. `weighted = TRUE` weights every pair by
#' the inverse of the number of observed cells of its row.
#' @param mat A character position matrix (from sw_positions()).
#' @param states The alphabet.
#' @param method One of the method names above.
#' @param params The method's settings.
#' @return A numeric matrix with dimnames.
sw_reference <- function(mat, states, method, params = list()) {
  p <- ncol(mat)
  row_weight <- if (isTRUE(params$weighted)) 1 / rowSums(!is.na(mat)) else rep(1, nrow(mat))
  adjacent <- seq_len(max(p - 1L, 0L))
  every <- if (p >= 2L) utils::combn(p, 2L) else matrix(integer(0), 2L, 0L)
  frequency <- \() sw_pair_table(mat, states, adjacent, adjacent + 1L, rep(1, length(adjacent)), row_weight)
  switch(method,
    frequency = frequency(),
    relative = sw_row_normalise(frequency()),
    reverse = t(frequency()),
    co_occurrence = {
      counts <- sw_pair_table(mat, states, every[1L, ], every[2L, ], rep(1, ncol(every)), row_weight)
      counts <- counts + t(counts) - diag(diag(counts), nrow(counts))
      if (identical(params$type, "relative")) sw_row_normalise(counts) else counts
    },
    attention = {
      forward <- sw_pair_table(mat, states, every[1L, ], every[2L, ],
                               exp(-(every[2L, ] - every[1L, ]) / (params$lambda %||% 1)))
      switch(params$direction %||% "forward", forward = forward, backward = t(forward),
             both = forward + t(forward))
    },
    gap = {
      grid <- expand.grid(from = adjacent, distance = seq_len((params$max_gap %||% 1L) + 1L))
      grid <- grid[grid$from + grid$distance <= p, , drop = FALSE]
      sw_pair_table(mat, states, grid$from, grid$from + grid$distance, 1 / grid$distance)
    },
    ngram = {
      n <- params$n_gram %||% 2L
      windows <- seq_len(max(p - n + 1L, 0L))
      covered <- unlist(lapply(windows, \(i) seq(i, i + n - 2L)))
      multiplicity <- tabulate(covered, nbins = max(p - 1L, 0L))
      used <- which(multiplicity > 0L)
      sw_pair_table(mat, states, used, used + 1L, multiplicity[used])
    },
    stop(sprintf("sw_reference(): no reference for method %s", method), call. = FALSE)
  )
}

#' Share of sequences whose first observed cell is each state
#'
#' @param mat A character position matrix.
#' @param states The alphabet.
#' @return A named numeric vector.
sw_initial <- function(mat, states) {
  first <- apply(mat, 1L, \(row) row[!is.na(row)][1L])
  first <- first[!is.na(first)]
  counts <- table(factor(first, levels = states), dnn = NULL)
  stats::setNames(as.numeric(counts) / sum(counts), states)
}

#' The scalings of ?build_network, applied in the order given
#'
#' minmax on a transition method spans every cell, zeros included (tna's rule); on any other method it maps the
#' non-zero cells into (eps, 1] with eps = sqrt(.Machine$double.eps); max divides by the largest absolute weight;
#' rank replaces the non-zero weights by their average ranks; normalize divides each row by its absolute sum.
#' @param m A numeric matrix.
#' @param scaling A character vector of scalings.
#' @param include_zeros Whether min-max spans the zero cells.
#' @return The scaled matrix.
sw_scale <- function(m, scaling, include_zeros = TRUE) {
  one <- function(m, s) {
    switch(s,
      minmax = {
        values <- if (include_zeros) as.vector(m) else m[m != 0]
        if (!length(values) || isTRUE(all.equal(min(values), max(values), tolerance = 0))) return(m)
        if (include_zeros) {
          (m - min(values)) / (max(values) - min(values))
        } else {
          eps <- sqrt(.Machine$double.eps)
          nonzero <- m != 0
          m[nonzero] <- eps + (1 - eps) * (m[nonzero] - min(values)) / (max(values) - min(values))
          m
        }
      },
      max = if (max(abs(m)) > 0) m / max(abs(m)) else m,
      rank = {
        nonzero <- m != 0
        m[nonzero] <- rank(m[nonzero], ties.method = "average")
        m
      },
      normalize = sw_row_normalise(m))
  }
  Reduce(one, scaling, m)
}

#' Weights below a threshold in absolute value set to zero
#'
#' @param m A numeric matrix.
#' @param threshold A number.
#' @return The matrix.
sw_threshold <- function(m, threshold) {
  m[abs(m) < threshold] <- 0
  m
}

#' A long event log of wide sequences: actor, time (seconds from a fixed origin, one per position), action
#'
#' @param wide A data frame of position columns (trailing missing values only).
#' @return A data frame.
sw_long <- function(wide) {
  rows <- lapply(seq_len(nrow(wide)), \(i) {
    states <- as.character(stats::na.omit(unlist(wide[i, ], use.names = FALSE)))
    data.frame(actor = sprintf("a%02d", i), time = 1700000000 + 60 * seq_along(states), action = states,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' tna's bootstrap summary with its bootstrap means and SDs added, one row per edge of its summary
#'
#' The same table as test-bootstrap_network.R's reference, kept here so the sweep files run alone.
#' @param boot A `tna::bootstrap()` result.
#' @return A data frame keyed by `from` and `to`.
sw_tna_bootstrap <- function(boot) {
  moments <- rs_edge_table(mean = boot$weights_mean, sd = boot$weights_sd)
  merge(as.data.frame(boot$summary), moments, by = c("from", "to"))
}

#' summary.net_stability()'s table written from tna's correlations (draws x drop proportions)
#'
#' The mean, SD and share at or above `threshold` of each measure's correlations per drop proportion; the same
#' reference as test-centrality_stability.R, kept here so the sweep files run alone.
#' @param estimate A `tna::estimate_cs()` result.
#' @param threshold The correlation threshold.
#' @return A data frame keyed by `measure` and `drop_prop`.
sw_stability_summary <- function(estimate, threshold = 0.7) {
  drop_prop <- attr(estimate, "drop_prop")
  do.call(rbind, lapply(names(estimate), \(m) {
    correlations <- estimate[[m]]$correlations
    data.frame(measure = m, drop_prop = drop_prop,
               mean_cor = colMeans(correlations, na.rm = TRUE),
               sd_cor = apply(correlations, 2L, stats::sd, na.rm = TRUE),
               prop_above = colMeans(!is.na(correlations) & correlations >= threshold))
  }))
}
