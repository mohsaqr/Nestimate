#' Compute descriptive indices for each sequence
#'
#' @description
#' `sequence_metrics()` computes descriptive indices for every sequence and
#' returns one row per sequence: length and validity, composition and
#' diversity, spells, transitions, the fate of the initial state, the
#' dominant and emergent states, cyclic recurrence and a composite
#' complexity index, plus the integrative potential of a set of favourable
#' states when `favorable` is given. It reads any input
#' [extract_sequences()] reads: a wide table of positions, an event log, the
#' output of [prepare()], a sequence network or group of networks, or a
#' `tna` model. The indices are those of `codyna::sequence_indices()`, with
#' three corrections listed under *Details*.
#'
#' The indices describe each sequence as it was observed, and several of
#' them (spells, transitions, recurrence) rest on few positions in a short
#' sequence; read them beside `valid_n`. Because there is one row per
#' sequence, the result can be summarised, plotted or used as predictors;
#' [sequence_patterns()] and [pattern_outcome()] look at which subsequences
#' recur and how they relate to an outcome.
#'
#' @details
#' Positions after a sequence's last observed state are not part of it;
#' \code{p} below is the position of the last observed state.
#' \describe{
#'   \item{valid_n, valid_proportion}{The number of observed states, and
#'     their share of the first \code{p} positions.}
#'   \item{unique_states}{The number of distinct states.}
#'   \item{mean_spell_duration, max_spell_duration}{The mean and the longest
#'     run of one state.}
#'   \item{longitudinal_entropy}{The Shannon entropy of the state
#'     proportions, divided by the log of the number of states in the data
#'     (0 to 1; 0 when the data hold one state).}
#'   \item{simpson_diversity}{One minus the sum of squared proportions.}
#'   \item{self_loop_tendency, transition_rate}{The share of transitions
#'     (pairs of adjacent observed states) that repeat the state, and the
#'     share that change it; the two sum to 1.}
#'   \item{transition_complexity}{The number of distinct state-changing
#'     transitions, as a share of the \eqn{a(a - 1)} possible ones, where
#'     \eqn{a} is the number of states in the data; \code{NA} when the data
#'     hold one state.}
#'   \item{initial_state_persistence}{The length of the opening spell as a
#'     share of \code{p}.}
#'   \item{initial_state_proportion}{The share of the observed states that
#'     are the first state.}
#'   \item{initial_state_influence_decay}{The share of the first third spent
#'     in the first state minus the share of the last third.}
#'   \item{cyclic_feedback_strength}{The largest share, over lags of 2 or
#'     more, of positions whose state recurs that many steps later.}
#'   \item{first_state, last_state, dominant_state, dominant_proportion,
#'     dominant_max_spell}{The opening and closing states; the most frequent
#'     state (ties go to the state that sorts first), its share and its
#'     longest spell.}
#'   \item{emergent_state, emergent_state_persistence,
#'     emergent_state_proportion}{A state that settles in later, its spell
#'     and its share; \code{NA} when there is none. The candidates are a
#'     spell of at least 3 that outlasts the dominant spell, the dominant
#'     state when it is not the first and has a spell of at least 3, and a
#'     return to the first state in a spell of at least 3 and more than half
#'     the opening one; the longest of them is taken.}
#'   \item{integrative_potential}{With \code{favorable}: the
#'     position-weighted share of favourable states over the observed
#'     positions, with weights \code{position^omega} (Brzinsky-Fay 2018).}
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
#' transition) are \code{NA}; a sequence with no observed state has
#' \code{valid_n} and \code{unique_states} 0 and every other index
#' \code{NA}. On sequences without
#' internal gaps the indices shared with TraMineR (\code{seqindic()},
#' \code{seqistatd()}, \code{seqdur()}, \code{seqintegr()}) agree exactly;
#' at an internal gap TraMineR joins the states on either side into one
#' sequence, so spells and transitions differ there.
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
#' @param x The sequence data, in any form [extract_sequences()] accepts.
#' @param favorable The favourable states; giving them adds
#'   `integrative_potential`. States that do not occur in the data are
#'   ignored, and an error of class `nestimate_bad_input` is raised when
#'   none of them occur. Default `NULL`.
#' @param omega The exponent of the position weights in
#'   `integrative_potential`; a positive value gives later positions more
#'   weight, `0` weights every position equally and a negative value gives
#'   earlier positions more weight; any single finite number is accepted.
#'   Ignored without `favorable`. Default `1`.
#' @param meta Whether to put each sequence's covariates (actor, session,
#'   group, metadata) before the indices. Default `FALSE`.
#' @param cols The state columns of a wide data frame, as in
#'   [extract_sequences()]. Default `NULL`.
#' @param actor,action,time The actor, action and timestamp columns of a
#'   long event log; setting `action` marks `x` as one. Default `NULL`.
#'   With `group` and sequences that are not an event log (a network, a
#'   wide table), `actor` instead names the covariate that nests the
#'   sequences (the session of a half-session, the student of a course), and
#'   the permutation test shuffles the group labels within it (see
#'   *Comparing groups*). For a network built with `actor`, the default
#'   nests the sequences in that actor (its `.session_label`), so the two
#'   halves of one session are compared within the session.
#' @param group A covariate of the sequences holding each sequence's group:
#'   for a grouped network its group column (for example `"phase"` after
#'   [mark_phases()]), for a wide table one of its non-state columns. With
#'   `group`, the result compares the indices between the groups instead of
#'   listing them per sequence. A grouped network is compared by its own
#'   groups by default; `meta = TRUE` lists its sequences instead. A name that is not a covariate of the
#'   sequences, or fewer than two groups, raises an error of class
#'   `nestimate_bad_input`. Default `NULL`.
#' @param indices The indices to return, by name (for example
#'   `c("valid_n", "transition_rate")`). Default `NULL` returns them all. A
#'   name that is not an index raises an error of class
#'   `nestimate_bad_input`.
#' @param iter The number of permutations of the group comparison, a whole
#'   number of at least 1. Default `1000`.
#' @param adjust The method of [stats::p.adjust()] applied to the p-values of
#'   the group comparison. Default `"BH"`.
#' @param format The layout of a group comparison: `"long"` (default), one
#'   row per index and group, or `"wide"`, one row per index with a column
#'   per group. Ignored without `group`.
#' @param ... Further arguments passed to [prepare()] for a long event log,
#'   such as `time_threshold` or `order`.
#'
#' @return A data frame with one row per sequence, in the order of
#'   [extract_sequences()], and the columns described in *Details*
#'   (`integrative_potential` only with `favorable`). With `meta = TRUE`,
#'   the input's covariate columns, if it has any, come first.
#'
#'   With `group`, a data frame of class `c("net_sequence_metrics",
#'   "data.frame")` with one row per numeric index and group: `index`,
#'   `group`, `n` (the group's sequences with a value of the index), `mean`
#'   and `sd`; then the index's comparison, which repeats on each of its
#'   rows: with two groups `difference` (the second group's mean minus the
#'   first's, in sorted order), `hedges_g` and its 95% interval `g_lower`,
#'   `g_upper`; for any number of groups `eta_squared`, the one-way
#'   `statistic` (F), the permutation `p_value` and the adjusted `p_adj`.
#'   With `format = "wide"`, one row per index with `n_<group>`,
#'   `mean_<group>` and `sd_<group>` columns instead. The settings are
#'   attributes (`groups`, `group`, `actor`, `iter`, `adjust`, `format`,
#'   `n_indices`). `print()` shows the first `n` indices, rounded;
#'   `summary()` returns the whole table as a plain data frame, and a subset
#'   is a plain data frame.
#'
#' @section Comparing groups:
#' For each numeric index, the groups are compared with the one-way F
#' statistic, which for two groups is the square of the t statistic. Its
#' p-value is `(b + 1) / (iter + 1)`, where `b` counts the permutations of
#' the group labels whose F is at least the observed one (Good, 2005); call
#' [base::set.seed()] first for reproducible p-values. A sequence with a
#' missing value of an index (a one-state sequence has no transition rate)
#' takes no part in that index. With `actor`, the labels are shuffled within
#' each actor whose sequences are in several groups and between the actors
#' whose sequences are all in one group, as in [sequence_compare()]. With
#' two groups, `hedges_g` is the standardised mean difference corrected for
#' small-sample bias, with an interval from its normal approximation
#' (Hedges, 1981). `eta_squared` is the share of the index's variance
#' between the groups.
#'
#' @references
#' Saqr, M., Tikka, S., & L\enc{ó}{o}pez-Pernas, S. (2026). codyna: Complex Dynamic
#' Systems. R package version 1.0.0.
#'
#' Brzinsky-Fay, C. (2018). Unused resources: Sequence and trajectory
#' indicators. International Symposium on Sequence Analysis and Related
#' Methods, Monte Verit\enc{à}{a}.
#'
#' Good, P. (2005). \emph{Permutation, parametric and bootstrap tests of
#' hypotheses}. Springer. \doi{10.1007/b138696}
#'
#' Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect
#' size and related estimators. \emph{Journal of Educational Statistics},
#' 6(2), 107-128. \doi{10.3102/10769986006002107}
#'
#' Gabadinho, A., Ritschard, G., M\enc{ü}{u}ller, N. S., & Studer, M. (2011).
#' Analyzing and visualizing state sequences in R with TraMineR. Journal of
#' Statistical Software, 40(4), 1-37. \doi{10.18637/jss.v040.i04}
#'
#' @examples
#' head(sequence_metrics(trajectories))
#'
#' head(sequence_metrics(trajectories, favorable = "Active"))
#'
#' # compare the indices of the first and second halves of each session
#' halves <- mark_phases(human_long, actor = "session_id",
#'                       time = "order_in_session",
#'                       labels = c("early", "late"))
#' net <- build_network(halves, actor = "session_id", action = "cluster",
#'                      group = "phase", method = "relative")
#' set.seed(1)
#' sequence_metrics(net, iter = 99)
#'
#' @seealso \code{\link{sequence_patterns}}, \code{\link{extract_sequences}}
#' @export
sequence_metrics <- function(x, favorable = NULL, omega = 1, meta = FALSE,
                             cols = NULL, actor = NULL, action = NULL,
                             time = NULL, group = NULL, indices = NULL,
                             iter = 1000L, adjust = "BH",
                             format = c("long", "wide"), ...) {
  format <- match.arg(format)
  .assert(
    "`favorable` must be NULL or character" = is.null(favorable) ||
      is.character(favorable),
    "`omega` must be a single finite number" = is.numeric(omega) &&
      length(omega) == 1L && is.finite(omega),
    "`meta` must be TRUE or FALSE" = isTRUE(meta) || isFALSE(meta),
    "`group` must be NULL or one column name" = is.null(group) ||
      (is.character(group) && length(group) == 1L && !is.na(group)),
    "`iter` must be a single whole number of at least 1" = is.numeric(iter) &&
      length(iter) == 1L && is.finite(iter) && iter >= 1 && iter == round(iter),
    "`adjust` must be one of p.adjust.methods" = is.character(adjust) &&
      length(adjust) == 1L && adjust %in% stats::p.adjust.methods,
    "`indices` must be NULL or index names" = is.null(indices) ||
      (is.character(indices) && length(indices) > 0L && !anyNA(indices))
  )
  # an event log's actor defines its sequences; for other input it is a
  # covariate, so it is kept out of the states like `group`
  is_log <- !is.null(action)
  seqs <- .as_sequences(x, cols = cols, actor = if (is_log) actor,
                        action = action, time = time,
                        exclude = c(group, if (!is_log) actor), ...)
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
  if (!is.null(indices)) {
    unknown <- setdiff(indices, names(out))
    if (length(unknown) > 0L) {
      stop(errorCondition(sprintf(
        "unknown index(es) in `indices`: %s. Available: %s.",
        paste(unknown, collapse = ", "), paste(names(out), collapse = ", ")),
        class = "nestimate_bad_input", call = NULL))
    }
    out <- out[unique(indices)]
  }
  # a grouped network is compared by its own groups unless its sequences
  # are asked for (meta = TRUE)
  if (is.null(group) && !meta && inherits(x, "netobject_group")) {
    group <- .group_col_name(x)
  }
  if (!is.null(group)) {
    # a network built with an actor nests its sequences in that actor
    nest <- if (!is.null(actor) && !is_log) actor
            else if (".session_label" %in% names(seqs$units)) ".session_label"
    return(.sm_compare(out, seqs$units, group = group, actor = nest,
                       iter = as.integer(iter), adjust = adjust,
                       format = format))
  }
  if (meta && ncol(seqs$units) > 0L) out <- cbind(seqs$units, out)
  rownames(out) <- NULL
  out
}


#' Compare sequence indices between groups
#'
#' @param metrics The per-sequence index table of `sequence_metrics()`.
#' @param units The sequences' covariates, aligned with `metrics`.
#' @param group,actor Covariate names: the group, and `NULL` or the nesting.
#' @param iter,adjust Permutations and p-value adjustment.
#' @return One row per numeric index; see `sequence_metrics()`.
#' @noRd
.sm_compare <- function(metrics, units, group, actor, iter, adjust,
                        format = "long") {
  if (!group %in% names(units)) {
    stop(errorCondition(sprintf(paste0(
      "`group` column '%s' is not a covariate of the sequences; available: %s."),
      group, paste(names(units), collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }
  labels <- as.character(units[[group]])
  groups <- sort(unique(labels[!is.na(labels)]))
  if (length(groups) < 2L) {
    stop(errorCondition(
      sprintf("`group` column '%s' has fewer than two groups.", group),
      class = "nestimate_bad_input", call = NULL))
  }
  nest <- !is.null(actor) && actor %in% names(units)
  design <- if (nest) {
    .cs_actor_design(as.character(units[[actor]]), labels, actor)
  } else NULL
  index_cols <- names(metrics)[vapply(metrics, is.numeric, logical(1))]
  values <- as.matrix(metrics[index_cols])

  f_stats <- function(lab) {
    vapply(seq_along(index_cols), function(j) {
      .sm_f_statistic(values[, j], lab, groups)
    }, numeric(1))
  }
  observed <- f_stats(labels)
  permute_once <- function() {
    if (nest) return(.cs_actor_permute(design))
    shuffled <- labels
    used <- which(!is.na(labels))
    shuffled[used] <- labels[used[sample.int(length(used))]]
    shuffled
  }
  permuted <- vapply(seq_len(iter), function(i) f_stats(permute_once()),
                     numeric(length(index_cols)))
  permuted <- matrix(permuted, nrow = length(index_cols))
  # ties within floating-point noise count as at least as large
  tol <- sqrt(.Machine$double.eps) * pmax(abs(observed), 1)
  exceed <- rowSums(permuted >= observed - tol, na.rm = TRUE)
  p_value <- ifelse(is.na(observed), NA_real_, (exceed + 1) / (iter + 1))

  per_group <- lapply(groups, function(g) {
    in_g <- !is.na(labels) & labels == g
    stats <- vapply(seq_along(index_cols), function(j) {
      v <- values[in_g, j]
      v <- v[!is.na(v)]
      c(length(v), if (length(v)) mean(v) else NA_real_,
        if (length(v) > 1L) stats::sd(v) else NA_real_)
    }, numeric(3))
    stats::setNames(as.data.frame(t(stats)),
                    paste0(c("n_", "mean_", "sd_"), g))
  })
  out <- cbind(data.frame(index = index_cols, stringsAsFactors = FALSE),
               do.call(cbind, per_group))
  if (length(groups) == 2L) {
    g <- .sm_hedges_g(out[[paste0("n_", groups[1])]],
                      out[[paste0("mean_", groups[1])]],
                      out[[paste0("sd_", groups[1])]],
                      out[[paste0("n_", groups[2])]],
                      out[[paste0("mean_", groups[2])]],
                      out[[paste0("sd_", groups[2])]])
    out$difference <- out[[paste0("mean_", groups[2])]] -
      out[[paste0("mean_", groups[1])]]
    out <- cbind(out, g)
  }
  out$eta_squared <- vapply(seq_along(index_cols), function(j) {
    .sm_eta_squared(values[, j], labels, groups)
  }, numeric(1))
  out$statistic <- observed
  out$p_value <- p_value
  out$p_adj <- stats::p.adjust(p_value, method = adjust)
  count_cols <- paste0("n_", groups)
  out[count_cols] <- lapply(out[count_cols], as.integer)
  rownames(out) <- NULL
  table <- if (identical(format, "wide")) out else .sm_long(out, groups)
  .result_df(table, "net_sequence_metrics", groups = groups, group = group,
             actor = actor, iter = iter, adjust = adjust, format = format,
             n_indices = nrow(out))
}

#' One row per index and group
#'
#' @param wide The wide comparison: `index`, per group `n_`, `mean_`,
#'   `sd_<group>`, then the index-level columns.
#' @param groups The group names, in the order of the columns.
#' @return `index`, `group`, `n`, `mean`, `sd`, then the index-level
#'   columns, which repeat on each group's row; indices in the order of
#'   `wide`, groups in `groups` order within each.
#' @noRd
.sm_long <- function(wide, groups) {
  n_idx <- nrow(wide)
  n_grp <- length(groups)
  row_id <- rep(seq_len(n_idx), each = n_grp)
  group_id <- rep(seq_len(n_grp), times = n_idx)
  pick <- function(prefix) {
    columns <- lapply(paste0(prefix, groups), function(cn) wide[[cn]])
    vapply(seq_along(row_id), function(k) {
      columns[[group_id[k]]][row_id[k]]
    }, if (identical(prefix, "n_")) integer(1L) else numeric(1L))
  }
  long <- data.frame(index = wide$index[row_id], group = groups[group_id],
                     n = pick("n_"), mean = pick("mean_"), sd = pick("sd_"),
                     stringsAsFactors = FALSE)
  per_index <- setdiff(names(wide),
                       c("index", paste0(c("n_", "mean_", "sd_"),
                                         rep(groups, each = 3L))))
  long[per_index] <- lapply(wide[per_index], `[`, row_id)
  stopifnot(nrow(long) == n_idx * n_grp)
  long
}


#' @rdname sequence_metrics
#' @param n For `print()` of a group comparison: the number of indices
#'   shown. Default `10`.
#' @param digits For `print()`: the decimals of the means, SDs and effect
#'   sizes. Default `2`; p-values show 3 significant digits.
#' @export
print.net_sequence_metrics <- function(x, n = 10L, digits = 2L, ...) {
  .check_print_args(n, digits)
  groups <- attr(x, "groups", exact = TRUE)
  if (is.null(groups)) {
    # a subset made with `[` keeps the class but not the settings
    print(.plain_df(x), ...)
    return(invisible(x))
  }
  n_indices <- attr(x, "n_indices", exact = TRUE)
  actor <- attr(x, "actor", exact = TRUE)
  cat(sprintf("Sequence metrics by %s  [%d indices | %d groups: %s]\n",
              attr(x, "group", exact = TRUE), n_indices, length(groups),
              paste(groups, collapse = ", ")))
  cat(sprintf("  permutation of the one-way F: %d iter (%s)%s\n",
              attr(x, "iter", exact = TRUE), attr(x, "adjust", exact = TRUE),
              if (is.null(actor)) "" else
                sprintf(", within %s", sub("^\\.session_label$",
                                           "session", actor))))
  table <- summary(x)
  shown <- utils::head(unique(table$index), n)
  top <- table[table$index %in% shown, , drop = FALSE]
  # display only: short headers fit the console; summary() keeps the names
  short <- c(difference = "diff", hedges_g = "g", g_lower = "g_lo",
             g_upper = "g_hi", eta_squared = "eta_sq", statistic = "F",
             p_value = "p")
  hit <- names(top) %in% names(short)
  names(top)[hit] <- short[names(top)[hit]]
  cat("\n")
  .print_result_rows(top, digits, p_cols = c("p", "p_adj"))
  if (n_indices > n) {
    cat(sprintf("  ... and %d more indices (%d rows); summary() returns them all\n",
                n_indices - n,
                nrow(table) - sum(table$index %in% shown)))
  }
  invisible(x)
}


#' @rdname sequence_metrics
#' @param object For `summary()`: a group comparison from
#'   `sequence_metrics(group = )`.
#' @export
summary.net_sequence_metrics <- function(object, ...) {
  .plain_df(object)
}


#' @rdname sequence_metrics
#' @param i,j For `[`: rows and columns, as for a data frame. A subset of a
#'   group comparison is a plain data frame.
#' @export
`[.net_sequence_metrics` <- function(x, i, j, ...) {
  out <- NextMethod()
  if (is.data.frame(out)) .plain_df(out) else out
}

#' One-way F statistic of one index across the groups
#' @noRd
.sm_f_statistic <- function(y, lab, groups) {
  keep <- !is.na(y) & !is.na(lab)
  y <- y[keep]; lab <- factor(lab[keep], levels = groups)
  n_k <- tabulate(lab, nbins = length(groups))
  present <- n_k > 0L
  if (sum(present) < 2L || length(y) <= sum(present)) return(NA_real_)
  means <- rowsum(y, lab, reorder = TRUE)[, 1] / pmax(n_k, 1L)
  grand <- mean(y)
  ss_between <- sum(n_k[present] * (means[present] - grand)^2)
  ss_within <- sum((y - means[as.integer(lab)])^2)
  df_between <- sum(present) - 1L
  df_within <- length(y) - sum(present)
  if (ss_within <= 0) return(if (ss_between > 0) Inf else NA_real_)
  (ss_between / df_between) / (ss_within / df_within)
}

#' Share of one index's variance between the groups
#' @noRd
.sm_eta_squared <- function(y, lab, groups) {
  keep <- !is.na(y) & !is.na(lab)
  y <- y[keep]; lab <- factor(lab[keep], levels = groups)
  if (length(y) < 2L) return(NA_real_)
  total <- sum((y - mean(y))^2)
  if (total <= 0) return(NA_real_)
  means <- tapply(y, lab, mean)
  n_k <- tabulate(lab, nbins = length(groups))
  sum(n_k * (means - mean(y))^2, na.rm = TRUE) / total
}

#' Hedges' g of the second group against the first, with a 95% interval
#' @noRd
.sm_hedges_g <- function(n1, m1, s1, n2, m2, s2) {
  df <- n1 + n2 - 2
  pooled <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / df)
  d <- ifelse(pooled > 0, (m2 - m1) / pooled, NA_real_)
  correction <- 1 - 3 / (4 * (n1 + n2) - 9)
  g <- correction * d
  se <- sqrt((n1 + n2) / (n1 * n2) + g^2 / (2 * (n1 + n2)))
  z <- stats::qnorm(0.975)
  data.frame(hedges_g = g, g_lower = g - z * se, g_upper = g + z * se)
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
    # in log space, scaled by the largest weight: seen^omega overflows to
    # Inf for a large omega and Inf / Inf gave NaN
    log_w <- omega * log(seen)
    w <- exp(log_w - max(log_w))
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
