# Layer 4 of the verification: recovery of known generating processes. The simulators below draw data from a
# process whose parameters are set in the test (a first-order Markov chain, a second-order chain with one planted
# rule, a mixture of first-order chains, a multivariate normal with a known precision matrix, an Ising model), and
# the criteria judge what Nestimate estimates from them against those parameters: bias within Monte Carlo error
# (Morris, White and Crowther 2019, Statistics in Medicine 38:2074-2102), rejection rates and interval coverage
# within the central binomial range of their nominal value. Everything here is base R and nothing calls Nestimate.

# ---- Randomness and replication ------------------------------------------------------------------------------

#' Evaluate an expression under a seed, leaving the session's random stream as it was
#'
#' @param seed The seed.
#' @param code The expression, evaluated after `set.seed(seed)`.
#' @return The value of `code`.
rc_with_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
          else rm(".Random.seed", envir = globalenv()), add = TRUE, after = FALSE)
  set.seed(seed)
  force(code)
}

#' The number of worker processes for the replicates: 4 on a Unix machine, 1 elsewhere or under a core limit
#'
#' @return A length-one integer.
rc_cores <- function() {
  limited <- isTRUE(as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_", "FALSE")))
  if (.Platform$OS.type == "unix" && !limited) 4L else 1L
}

#' Run the replicates, in parallel where the platform allows, and stop on the first failed replicate
#'
#' Every replicate draws its own data under its own seed (`rc_with_seed()` inside `fun`), so the results do not
#' depend on the number of workers or on the order the workers run in. A forked worker's warnings would be lost,
#' so each replicate records its own and they are raised again here, once per distinct message with its count.
#' @param replicates The replicate numbers, `seq_len(R)`.
#' @param fun A function of one replicate number.
#' @param cores Worker processes.
#' @return A list, one element per replicate.
rc_replicate <- function(replicates, fun, cores = rc_cores()) {
  run <- function(r) {
    warned <- character()
    value <- withCallingHandlers(fun(r), warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    list(value = value, warned = warned)
  }
  results <- if (cores > 1L) {
    parallel::mclapply(replicates, run, mc.cores = cores, mc.preschedule = TRUE)
  } else {
    lapply(replicates, run)
  }
  failed <- vapply(results, \(r) inherits(r, "try-error") || is.null(r), logical(1L))
  if (any(failed)) {
    first <- results[[which(failed)[[1L]]]]
    stop(errorCondition(sprintf("%d of %d replicates failed; the first: %s", sum(failed), length(results),
                                if (is.null(first)) "no result (worker died)" else as.character(first)),
                        class = "rc_replicate_failed", call = NULL))
  }
  warned <- table(unlist(lapply(results, \(r) unique(r$warned))))
  invisible(lapply(names(warned), \(message) {
    warning(warningCondition(sprintf("%d of %d replicates warned: %s", warned[[message]], length(results), message),
                             class = "rc_replicate_warning", call = NULL))
  }))
  lapply(results, \(r) r$value)
}

# ---- Criteria --------------------------------------------------------------------------------------------------

#' The central binomial range of a rate observed over R independent replicates
#'
#' @param replicates R.
#' @param rate The nominal rate (a rejection rate alpha, or a coverage 1 - alpha).
#' @param level The probability content of the range.
#' @return A named numeric vector: `lower`, `upper` (as rates) and `tolerance`, the larger of the two
#'   half-widths, the single number the harness judges `|observed - rate|` against.
rc_binomial_band <- function(replicates, rate, level = 0.999) {
  tail <- (1 - level) / 2
  bounds <- stats::qbinom(c(tail, 1 - tail), replicates, rate) / replicates
  c(lower = bounds[[1L]], upper = bounds[[2L]], tolerance = max(rate - bounds[[1L]], bounds[[2L]] - rate))
}

#' Bias of replicated estimates in units of its Monte Carlo standard error
#'
#' @param estimates A list of numeric vectors or matrices of one shape, one per replicate.
#' @param truth The generating values, the same shape.
#' @return A list: `mean` (the mean estimate), `mcse` (the Monte Carlo SE of that mean, sd / sqrt(R)) and `z`
#'   ((mean - truth) / mcse; 0 where the estimate never varied and equals the truth).
rc_standardised_bias <- function(estimates, truth) {
  stack <- simplify2array(estimates)
  margin <- seq_along(dim(as.array(truth)))
  replicates <- length(estimates)
  mean_estimate <- apply(stack, margin, mean)
  mcse <- apply(stack, margin, stats::sd) / sqrt(replicates)
  difference <- mean_estimate - truth
  z <- ifelse(mcse > 0, difference / pmax(mcse, .Machine$double.xmin), ifelse(abs(difference) < 1e-12, 0, Inf))
  if (is.matrix(truth)) {
    mean_estimate <- matrix(mean_estimate, nrow(truth), dimnames = dimnames(truth))
    mcse <- matrix(mcse, nrow(truth), dimnames = dimnames(truth))
    z <- matrix(z, nrow(truth), dimnames = dimnames(truth))
  } else {
    names(mean_estimate) <- names(mcse) <- names(z) <- names(truth)
  }
  list(mean = mean_estimate, mcse = mcse, z = z)
}

#' A matrix as one row per cell, keyed `from`, `to`
#'
#' @param ... Named matrices with one set of dimnames; each becomes a column.
#' @return A data frame with `from`, `to` and one column per matrix.
rc_cells <- function(...) {
  matrices <- list(...)
  first <- matrices[[1L]]
  cells <- expand.grid(from = rownames(first), to = colnames(first), stringsAsFactors = FALSE)
  cells[names(matrices)] <- lapply(matrices, \(m) as.vector(m))
  cells
}

# ---- Markov chains ---------------------------------------------------------------------------------------------

#' Draw the next state of every sequence at once
#'
#' @param probabilities A matrix with one row per sequence: the distribution of its next state.
#' @return Integer state indices.
rc_draw_rows <- function(probabilities) {
  cumulative <- t(apply(probabilities, 1L, cumsum))
  if (ncol(probabilities) == 1L) cumulative <- t(cumulative)
  u <- stats::runif(nrow(probabilities))
  pmin(1L + rowSums(u > cumulative[, -ncol(cumulative), drop = FALSE]), ncol(probabilities))
}

#' Sequences of a first-order Markov chain, or of a mixture of chains
#'
#' @param n_seq,len Number and length of the sequences.
#' @param transition A named row-stochastic matrix, or a list of them (one per mixture component).
#' @param initial The initial distribution, or a list of them.
#' @param component Integer component of each sequence (mixtures); `NULL` for one chain.
#' @return A data frame with columns T1..T<len>, one sequence per row, states as character.
rc_markov_sequences <- function(n_seq, len, transition, initial, component = NULL) {
  transitions <- if (is.list(transition)) transition else list(transition)
  initials <- if (is.list(initial)) initial else list(initial)
  component <- component %||% rep(1L, n_seq)
  states <- rownames(transitions[[1L]])
  n_states <- length(states)
  stacked <- do.call(rbind, transitions)
  first <- rc_draw_rows(do.call(rbind, initials)[component, , drop = FALSE])
  path <- Reduce(\(current, step) rc_draw_rows(stacked[(component - 1L) * n_states + current, , drop = FALSE]),
                 seq_len(len - 1L), accumulate = TRUE, init = first)
  wide <- as.data.frame(matrix(states[unlist(path)], n_seq, len), stringsAsFactors = FALSE)
  names(wide) <- paste0("T", seq_len(len))
  wide
}

#' Sequences of a first-order chain with one planted second-order rule
#'
#' After the two states `context`, the next state is `target` with probability `strength` (the other states share
#' the rest in proportion to the first-order row); everywhere else the first-order matrix applies.
#' @inheritParams rc_markov_sequences
#' @param context Two states.
#' @param target The state that follows them.
#' @param strength Its probability after `context`.
#' @return A data frame with columns T1..T<len>.
rc_planted_sequences <- function(n_seq, len, transition, initial, context, target, strength) {
  states <- rownames(transition)
  rule_row <- transition[context[[2L]], ]
  rule_row[target] <- 0
  rule_row <- rule_row / sum(rule_row) * (1 - strength)
  rule_row[target] <- strength
  first <- rc_draw_rows(matrix(initial, n_seq, length(states), byrow = TRUE))
  second <- rc_draw_rows(transition[first, , drop = FALSE])
  context_index <- match(context, states)
  path <- Reduce(\(window, step) {
    planted <- window$previous == context_index[[1L]] & window$current == context_index[[2L]]
    probabilities <- transition[window$current, , drop = FALSE]
    probabilities[planted, ] <- matrix(rule_row, sum(planted), length(states), byrow = TRUE)
    list(previous = window$current, current = rc_draw_rows(probabilities))
  }, seq_len(len - 2L), accumulate = TRUE, init = list(previous = first, current = second))
  columns <- c(list(first), lapply(path, \(window) window$current))
  wide <- as.data.frame(matrix(states[unlist(columns)], n_seq, len), stringsAsFactors = FALSE)
  names(wide) <- paste0("T", seq_len(len))
  wide
}

#' Expected transition counts of n_seq chains of length len: n_seq * sum_t diag(pi_t) P, pi_1 = initial
#'
#' @param n_seq,len Number and length of the sequences.
#' @param transition The transition matrix.
#' @param initial The initial distribution.
#' @return The matrix of expected counts.
rc_expected_counts <- function(n_seq, len, transition, initial) {
  marginals <- Reduce(\(pi_t, step) as.vector(pi_t %*% transition), seq_len(len - 2L), accumulate = TRUE,
                      init = initial)
  visits <- Reduce(`+`, marginals)
  n_seq * visits * transition
}

# ---- Mixtures ----------------------------------------------------------------------------------------------------

#' All permutations of 1..k
#'
#' @param k A small whole number.
#' @return A matrix with one permutation per row.
rc_permutations <- function(k) {
  if (k == 1L) return(matrix(1L, 1L, 1L))
  smaller <- rc_permutations(k - 1L)
  do.call(rbind, lapply(seq_len(k), \(first) cbind(first, matrix(setdiff(seq_len(k), first)[smaller], ncol = k - 1L))))
}

#' The relabelling of estimated components that brings them closest to the generating ones
#'
#' @param estimated,truth Lists of k transition matrices with one set of dimnames.
#' @return An integer vector: element j is the estimated component matched to generating component j.
rc_match_components <- function(estimated, truth) {
  orders <- rc_permutations(length(truth))
  distance <- apply(orders, 1L, \(o) sum(vapply(seq_along(truth), \(j) sum((estimated[[o[[j]]]] - truth[[j]])^2),
                                                numeric(1L))))
  orders[which.min(distance), ]
}

#' Adjusted Rand index of two partitions (Hubert and Arabie 1985)
#'
#' @param a,b Two labelings of the same units.
#' @return A number.
rc_adjusted_rand <- function(a, b) {
  table_ab <- table(a, b)
  pairs <- \(x) sum(x * (x - 1) / 2)
  index <- pairs(table_ab)
  rows <- pairs(rowSums(table_ab))
  cols <- pairs(colSums(table_ab))
  expected <- rows * cols / pairs(sum(table_ab))
  (index - expected) / ((rows + cols) / 2 - expected)
}

#' The Bayes classification of sequences under the generating mixture (the oracle the fit is measured against)
#'
#' @param wide Sequences, one per row.
#' @param transition,initial Lists of the generating matrices and initial distributions.
#' @param mixing The generating mixing proportions.
#' @return Integer component of each sequence with the highest posterior probability.
rc_oracle_classes <- function(wide, transition, initial, mixing) {
  states <- rownames(transition[[1L]])
  codes <- matrix(match(as.matrix(wide), states), nrow(wide))
  from <- codes[, -ncol(codes), drop = FALSE]
  to <- codes[, -1L, drop = FALSE]
  log_joint <- vapply(seq_along(transition), \(j) {
    log(mixing[[j]]) + log(initial[[j]][codes[, 1L]]) +
      rowSums(matrix(log(transition[[j]][cbind(as.vector(from), as.vector(to))]), nrow(wide)))
  }, numeric(nrow(wide)))
  max.col(log_joint, ties.method = "first")
}

# ---- Psychometric networks ---------------------------------------------------------------------------------------

#' A precision matrix of a ring of p variables with given partial correlations of neighbours
#'
#' Unit diagonal, -pcor on the ring's edges, zero elsewhere; positive definite while every row's absolute
#' off-diagonal sum is below 1.
#' @param p Number of variables.
#' @param pcor The partial correlation of each ring edge (i, i + 1), length p (the last closes the ring).
#' @param drop Ring edges (by their first variable) set to zero.
#' @return A named precision matrix.
rc_ring_precision <- function(p, pcor, drop = integer()) {
  pcor[drop] <- 0
  precision <- diag(p)
  edges <- cbind(seq_len(p), c(seq_len(p - 1L) + 1L, 1L))
  precision[edges] <- -pcor
  precision[edges[, 2:1]] <- -pcor
  stopifnot("the precision matrix must be positive definite" =
              min(eigen(precision, symmetric = TRUE, only.values = TRUE)$values) > 0)
  dimnames(precision) <- list(paste0("V", seq_len(p)), paste0("V", seq_len(p)))
  precision
}

#' The partial correlations of a precision matrix, -K_ij / sqrt(K_ii K_jj), zero diagonal
#'
#' @param precision A precision matrix.
#' @return The matrix of partial correlations.
rc_partial <- function(precision) {
  scale <- sqrt(diag(precision))
  partial <- -precision / outer(scale, scale)
  diag(partial) <- 0
  partial
}

#' Multivariate normal draws with a given precision matrix
#'
#' @param n Rows.
#' @param precision The precision matrix (named).
#' @return A data frame with the precision matrix's names as columns.
rc_normal <- function(n, precision) {
  covariance <- chol2inv(chol(precision))
  draws <- matrix(stats::rnorm(n * nrow(precision)), n) %*% chol(covariance)
  colnames(draws) <- colnames(precision)
  as.data.frame(draws)
}

#' Exact draws from an Ising model on {0, 1}^p, P(x) proportional to exp(tau'x + sum_{i<j} omega_ij x_i x_j)
#'
#' The 2^p states are enumerated, so p is small.
#' @param n Rows.
#' @param omega Symmetric interaction matrix with zero diagonal (named).
#' @param tau Thresholds.
#' @return A data frame of 0/1 integers.
rc_ising <- function(n, omega, tau) {
  p <- nrow(omega)
  states <- as.matrix(expand.grid(rep(list(0:1), p)))
  potential <- as.vector(states %*% tau) + rowSums((states %*% omega) * states) / 2
  probability <- exp(potential - max(potential))
  rows <- sample.int(nrow(states), n, replace = TRUE, prob = probability / sum(probability))
  draws <- states[rows, , drop = FALSE]
  dimnames(draws) <- list(NULL, colnames(omega))
  storage.mode(draws) <- "integer"
  as.data.frame(draws)
}

#' Sensitivity and specificity of an estimated edge set against the generating one (upper triangle)
#'
#' @param estimate,truth Symmetric weight matrices.
#' @return A named vector: `sensitivity` (share of true edges estimated non-zero) and `specificity` (share of
#'   absent edges estimated zero).
rc_edge_recovery <- function(estimate, truth) {
  upper <- upper.tri(truth)
  present <- truth[upper] != 0
  found <- abs(estimate[upper]) > 0
  c(sensitivity = mean(found[present]), specificity = mean(!found[!present]))
}
