#' Per-Sequence Indices
#'
#' Computes descriptive indices for every sequence: length and validity,
#' composition and diversity, spells, transitions, the fate of the initial
#' state, the dominant and emergent states, cyclic recurrence, and a
#' composite complexity index. The indices are those of
#' \code{codyna::sequence_indices()}, with three corrections listed under
#' Details. One row per sequence, so the result can be described, plotted,
#' or used as predictors.
#'
#' @details
#' Positions after a sequence's last observed state are not part of it;
#' \code{p} below is the position of the last observed state.
#' \describe{
#'   \item{valid_n, valid_proportion}{Observed states, and their share of
#'     the first \code{p} positions.}
#'   \item{unique_states}{Number of distinct states.}
#'   \item{mean_spell_duration, max_spell_duration}{Mean and longest run of
#'     one state.}
#'   \item{longitudinal_entropy}{Shannon entropy of the state proportions,
#'     divided by the log of the number of states in the data (0 to 1).}
#'   \item{simpson_diversity}{One minus the sum of squared proportions.}
#'   \item{self_loop_tendency, transition_rate}{Share of transitions (pairs
#'     of adjacent observed states) that repeat the state, and share that
#'     change it; the two sum to 1.}
#'   \item{transition_complexity}{Distinct state-changing transitions, as a
#'     share of the \eqn{a(a - 1)} possible ones.}
#'   \item{initial_state_persistence}{Length of the opening spell as a share
#'     of \code{p}.}
#'   \item{initial_state_proportion}{Share of the sequence in its first
#'     state.}
#'   \item{initial_state_influence_decay}{Share of the first third spent in
#'     the first state minus the share of the last third.}
#'   \item{cyclic_feedback_strength}{Largest share, over lags of 2 or more,
#'     of positions whose state recurs that many steps later.}
#'   \item{first_state, last_state, dominant_state, dominant_proportion,
#'     dominant_max_spell}{Opening and closing states; the most frequent
#'     state, its share, and its longest spell.}
#'   \item{emergent_state, emergent_state_persistence,
#'     emergent_state_proportion}{A state that settles in later (a spell of
#'     at least 3 that outlasts the dominant spell, a dominant state that is
#'     not the first, or a return to the first state in a spell of at least 3
#'     and more than half the opening one), its spell and its share;
#'     \code{NA} when none.}
#'   \item{integrative_potential}{With \code{favorable}: position-weighted
#'     share of favourable states over the observed positions, weights
#'     \code{position^omega} (Brzinsky-Fay 2018).}
#'   \item{complexity_index}{\code{0.4 * longitudinal_entropy + 0.4 *
#'     transition_rate + 0.2 * min(CV of spell durations, 1)}; 0 for a
#'     sequence of one spell (one state, no change). This is codyna's
#'     composite, not the complexity index of Gabadinho et al. (2011) that
#'     \code{TraMineR::seqici()} computes.}
#' }
#'
#' Missing values: a missing position is not an observation. The first
#' state is the first observed one, a missing position inside a sequence
#' ends a spell, and a pair of states across it is not a transition.
#' Indices that are undefined for a sequence (for example a rate with no
#' transition) are \code{NA}. On sequences without internal gaps the
#' indices shared with TraMineR (\code{seqindic()}, \code{seqistatd()},
#' \code{seqdur()}, \code{seqintegr()}) agree exactly; at an internal gap
#' TraMineR joins the states on either side into one sequence, so spells
#' and transitions differ there.
#'
#' Corrections relative to \code{codyna::sequence_indices()} 1.0.0, which
#' otherwise agrees value for value: \code{dominant_max_spell} is computed
#' per sequence (codyna reports the last sequence's value for every row);
#' \code{longitudinal_entropy} is normalised by the log of the number of
#' states (codyna counts \code{NA} as one when the data have missing
#' values); and \code{complexity_index} uses the normalised entropy once
#' (codyna divides it by the log of the number of states a second time) and
#' is 0, not \code{NA}, for a single spell.
#' Leading and internal missing values follow the rule above, where codyna
#' reads position 1 as the first state and counts missing positions as
#' unfavourable and as possible transitions.
#'
#' @param x Sequence data, in any form accepted by
#'   \code{\link{extract_sequences}()}.
#' @param favorable Optional character vector of favourable states; adds
#'   \code{integrative_potential}. States not in the data raise
#'   \code{nestimate_bad_input} when none of them occur.
#' @param omega Exponent of the position weights in
#'   \code{integrative_potential}. Default \code{1}.
#' @param meta Logical. Put each sequence's covariates (actor, session,
#'   group, metadata) before the indices. Default \code{FALSE}.
#' @param cols State columns of a wide data frame; see
#'   \code{\link{extract_sequences}()}.
#' @param ... Passed to \code{\link{extract_sequences}()}, e.g.
#'   \code{actor}, \code{action} and \code{time} for a long event log.
#'
#' @return A \code{data.frame} with one row per sequence, in the order of
#'   \code{\link{extract_sequences}()}, and the columns described in Details
#'   (\code{integrative_potential} only with \code{favorable}).
#'
#' @references
#' Saqr, M., Tikka, S., & López-Pernas, S. (2026). codyna: Complex Dynamic
#' Systems. R package version 1.0.0.
#'
#' Brzinsky-Fay, C. (2018). Unused resources: Sequence and trajectory
#' indicators. International Symposium on Sequence Analysis and Related
#' Methods, Monte Verità.
#'
#' Gabadinho, A., Ritschard, G., Müller, N. S., & Studer, M. (2011).
#' Analyzing and visualizing state sequences in R with TraMineR. Journal of
#' Statistical Software, 40(4), 1-37. \doi{10.18637/jss.v040.i04}
#'
#' @examples
#' head(sequence_metrics(trajectories))
#'
#' head(sequence_metrics(trajectories, favorable = "Active"))
#'
#' @seealso \code{\link{sequence_patterns}}, \code{\link{extract_sequences}}
#' @export
sequence_metrics <- function(x, favorable = NULL, omega = 1, meta = FALSE,
                             cols = NULL, ...) {
  stopifnot(
    "`favorable` must be NULL or character" = is.null(favorable) ||
      is.character(favorable),
    "`omega` must be a single finite number" = is.numeric(omega) &&
      length(omega) == 1L && is.finite(omega),
    "`meta` must be TRUE or FALSE" = isTRUE(meta) || isFALSE(meta)
  )
  seqs <- .as_sequences(x, cols = cols, ...)
  alphabet <- seqs$alphabet
  a <- length(alphabet)
  codes <- matrix(match(seqs$states, alphabet), nrow = nrow(seqs$states))
  fav <- if (is.null(favorable)) integer(0) else which(alphabet %in% favorable)
  if (!is.null(favorable) && length(fav) == 0L) {
    stop(errorCondition(
      sprintf("none of `favorable` (%s) occur among the states (%s).",
              paste(favorable, collapse = ", "), paste(alphabet, collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }

  observed <- !is.na(codes)
  last_obs <- max.col(observed, ties.method = "last")
  rows <- lapply(seq_len(nrow(codes)), \(i) {
    .sm_row(codes[i, ], last_obs[i], a, fav, omega)
  })
  num <- \(field) vapply(rows, \(r) as.numeric(r[[field]]), numeric(1))
  state <- \(field) alphabet[vapply(rows, \(r) r[[field]], integer(1))]

  out <- data.frame(
    valid_n = vapply(rows, \(r) r$valid, integer(1)),
    valid_proportion = num("valid") / last_obs,
    unique_states = vapply(rows, \(r) r$unique_states, integer(1)),
    mean_spell_duration = num("mean_spell"),
    max_spell_duration = num("max_spell"),
    longitudinal_entropy = num("entropy"),
    simpson_diversity = num("simpson"),
    self_loop_tendency = num("loops"),
    transition_rate = num("rate"),
    transition_complexity = num("trans_complexity"),
    initial_state_persistence = num("init_persistence"),
    initial_state_proportion = num("init_proportion"),
    initial_state_influence_decay = num("init_decay"),
    cyclic_feedback_strength = .sm_cyclic(codes, last_obs),
    first_state = state("first"),
    last_state = state("last"),
    dominant_state = state("dominant"),
    dominant_proportion = num("dominant_proportion"),
    dominant_max_spell = num("dominant_spell"),
    emergent_state = state("emergent"),
    emergent_state_persistence = num("emergent_spell"),
    emergent_state_proportion = num("emergent_proportion"),
    integrative_potential = num("integrative"),
    complexity_index = num("complexity"),
    stringsAsFactors = FALSE
  )
  # a row with no observed state: max.col() still returns the last column,
  # so the two indices computed outside .sm_row() must be set to NA here
  empty <- rowSums(observed) == 0L
  out$valid_proportion[empty] <- NA_real_
  out$cyclic_feedback_strength[empty] <- NA_real_
  if (length(fav) == 0L) out$integrative_potential <- NULL
  if (meta && ncol(seqs$units) > 0L) out <- cbind(seqs$units, out)
  rownames(out) <- NULL
  out
}


#' Indices of one integer-coded sequence
#'
#' @param r Integer state codes, NA after the end or at gaps.
#' @param p Position of the last observed state.
#' @param a Alphabet size.
#' @param fav Codes of favourable states (possibly empty).
#' @param omega Position-weight exponent.
#' @return Named list of scalars (state fields are integer codes).
#' @noRd
.sm_row <- function(r, p, a, fav, omega) {
  valid <- sum(!is.na(r))
  if (valid == 0L) return(.sm_empty_row())
  first <- r[!is.na(r)][1L]
  freq <- tabulate(r, nbins = a)
  prop <- freq / valid
  runs <- rle(r)
  values <- runs$values[!is.na(runs$values)]
  spells <- runs$lengths[!is.na(runs$values)]

  from <- r[-length(r)]
  to <- r[-1L]
  ok <- !is.na(from) & !is.na(to)
  n_trans <- sum(ok)
  n_self <- sum(from[ok] == to[ok])
  # distinct state-changing transitions, each pair encoded as one integer
  pairs <- (from[ok] - 1L) * a + to[ok]
  n_change_types <- length(unique(pairs[from[ok] != to[ok]]))
  shares <- prop[freq > 0]

  # opening spell: from the first observed state until the state changes or
  # a position is missing
  opening <- r[seq(which(!is.na(r))[1L], p)]
  first_change <- which(is.na(opening) | opening != first)[1L]
  persistence <- (if (is.na(first_change)) length(opening) else
    first_change - 1L) / p
  first_third <- r[seq_len(ceiling(p / 3))]
  last_third <- r[seq(ceiling(2 * p / 3), p)]
  decay <- sum(first_third == first, na.rm = TRUE) / length(first_third) -
    sum(last_third == first, na.rm = TRUE) / length(last_third)

  dominant <- which.max(freq)
  dominant_spell <- max(spells[values == dominant])
  emergent <- .sm_emergent(values, spells, dominant, dominant_spell, first,
                           prop)

  entropy <- if (a > 1L) -sum(shares * log(shares)) / log(a) else 0
  # transitions are pairs of adjacent observed states, so loops + rate = 1
  rate <- if (n_trans > 0L) (n_trans - n_self) / n_trans else NA_real_
  # one spell: one state, no change, no spread of durations -> 0, as in
  # TraMineR's seqici() (codyna: NA, the CV of a single duration)
  complexity <- if (length(spells) == 1L) {
    0
  } else {
    0.4 * entropy + 0.4 * rate +
      0.2 * min(stats::sd(spells) / mean(spells), 1)
  }
  # position weights over observed positions only (Brzinsky-Fay; TraMineR's
  # seqintegr): a missing position is not an unfavourable one
  integrative <- if (length(fav) > 0L) {
    seen <- which(!is.na(r))
    w <- seen^omega
    sum((r[seen] %in% fav) * w) / sum(w)
  } else {
    NA_real_
  }

  list(
    valid = valid,
    unique_states = sum(freq > 0),
    mean_spell = mean(spells),
    max_spell = max(spells),
    entropy = entropy,
    simpson = 1 - sum(prop^2),
    loops = if (n_trans > 0L) n_self / n_trans else NA_real_,
    rate = rate,
    trans_complexity = if (a > 1L) n_change_types / (a * (a - 1)) else NA_real_,
    init_persistence = persistence,
    init_proportion = prop[first],
    init_decay = decay,
    first = first,
    last = r[p],
    dominant = dominant,
    dominant_proportion = prop[dominant],
    dominant_spell = dominant_spell,
    emergent = emergent$state,
    emergent_spell = emergent$spell,
    emergent_proportion = emergent$proportion,
    integrative = integrative,
    complexity = complexity
  )
}


#' Indices of a sequence with no observed state: all undefined
#' @noRd
.sm_empty_row <- function() {
  num <- c("mean_spell", "max_spell", "entropy", "simpson", "loops", "rate",
           "trans_complexity", "init_persistence", "init_proportion",
           "init_decay", "dominant_proportion", "dominant_spell",
           "emergent_spell", "emergent_proportion", "integrative",
           "complexity")
  state <- c("first", "last", "dominant", "emergent")
  c(list(valid = 0L, unique_states = 0L),
    stats::setNames(rep(list(NA_real_), length(num)), num),
    stats::setNames(rep(list(NA_integer_), length(state)), state))
}


#' Emergent state of one sequence (codyna's rule)
#'
#' Candidates, in order: a spell of at least 3 that outlasts the dominant
#' spell; the dominant state when it is not the first and has a spell of at
#' least 3; the first state when it returns in a spell of at least 3 and
#' more than half the opening spell (codyna's \code{2 * return > opening}).
#' The candidate with the longest spell wins.
#' @noRd
.sm_emergent <- function(values, spells, dominant, dominant_spell, first,
                         prop) {
  none <- list(state = NA_integer_, spell = NA_real_, proportion = NA_real_)
  persisting <- spells >= 3L
  if (!any(persisting)) return(none)
  outlasting <- persisting & spells > dominant_spell
  true_state <- 0L
  true_spell <- 0
  if (any(outlasting)) {
    longest <- which.max(spells[outlasting])
    true_state <- values[outlasting][longest]
    true_spell <- spells[outlasting][longest]
  }
  dominant_emergent <- dominant_spell >= 3L && dominant != first
  returns <- spells[-1L][values[-1L] == first]
  max_return <- if (length(returns) > 0L) max(returns) else 0
  return_emergent <- max_return * 2 > spells[1L] && max_return >= 3L
  candidates <- c(true_state, dominant, first)
  candidate_spells <- c(true_spell, dominant_spell * dominant_emergent,
                        max_return * return_emergent)
  if (!any(candidate_spells > 0)) return(none)
  best <- which.max(candidate_spells)
  list(state = as.integer(candidates[best]), spell = candidate_spells[best],
       proportion = prop[candidates[best]])
}


#' Cyclic feedback strength for all sequences at once
#'
#' For each lag of 2 to k - 1: positions whose state recurs that many steps
#' later, divided by the number of comparable positions (p - lag). The
#' largest value over lags, 0 when no lag applies.
#' @noRd
.sm_cyclic <- function(codes, last_obs) {
  n <- nrow(codes)
  k <- ncol(codes)
  if (k < 3L) return(numeric(n))
  lags <- seq(2L, k - 1L)
  strength <- matrix(vapply(lags, \(lag) {
    span <- seq_len(k - lag)
    recurs <- rowSums(codes[, span, drop = FALSE] ==
                        codes[, span + lag, drop = FALSE], na.rm = TRUE)
    ifelse(last_obs > lag, recurs / (last_obs - lag), recurs)
  }, numeric(n)), nrow = n)
  pmax(0, apply(strength, 1L, max))
}
