#' Discover Sequential Patterns
#'
#' @description
#' `sequence_patterns()` finds the ordered patterns of states that occur in
#' a set of sequences and returns one row per pattern with how often it
#' occurs, how many sequences contain it, and its lift. It reads sequences
#' in every form [extract_sequences()] accepts, and finds three kinds of
#' pattern:
#' \describe{
#'   \item{`"ngram"`}{Contiguous runs of `len` states, e.g. `A->B->C`.}
#'   \item{`"gapped"`}{Two states separated by `gap` positions of any
#'     state, e.g. `A->*->B` (gap 1) or `A->**->B` (gap 2).}
#'   \item{`"repeated"`}{Runs of `len` identical states, e.g. `A->A->A`.}
#' }
#' A `pattern` template searches for the instances of one pattern instead,
#' and `group` adds a test of whether the share of sequences containing a
#' pattern differs between groups (see *Group test*).
#'
#' Support and lift describe the sequences at hand. A lift above 1 says the
#' states of a pattern occur in this order in more sequences than their
#' separate presence predicts; it does not account for the transition
#' structure, so a pattern that the transition probabilities make likely
#' can have a lift above 1. Many patterns are searched at once, so read
#' the group test's `p_adj`, not `p_value`, and its effect sizes beside
#' them. To relate patterns to an outcome, use [pattern_outcome()]; to
#' compare how often each run of states occurs in each group,
#' [sequence_compare()].
#'
#' @details
#' `support` is the share of sequences that contain the pattern at least
#' once. `lift` compares it with the support expected if the states in the
#' pattern appeared in sequences independently:
#' `support / prod(state_support)`, where the support of a state is the
#' share of sequences in which it occurs and wildcard positions are left
#' out.
#'
#' A window that contains a missing value, at a fixed or a wildcard
#' position, is skipped, so a pattern never spans a gap or the end of a
#' sequence.
#'
#' `proportion` is computed over all patterns of the same length before
#' `min_freq`, `min_support`, `start`, `end` and `contain` are applied, so
#' it does not change with the filters. Counts, support and lift agree with
#' `codyna::discover_patterns()`, which computes `proportion` after the
#' filters instead.
#'
#' @param x The sequences, in any form [extract_sequences()] accepts:
#'   prepared data, a long event log (with `action`), a wide data frame or
#'   character matrix, a `netobject`, `netobject_group` or `tna` model.
#' @param type The kind of pattern to find: `"ngram"` (default),
#'   `"gapped"` or `"repeated"`. Ignored when `pattern` is given.
#' @param len The pattern lengths, in states, for `"ngram"` and
#'   `"repeated"`: positive whole numbers. Default `2:5`.
#' @param gap The gap sizes, in positions, for `"gapped"`: positive whole
#'   numbers. Default `1:3`.
#' @param min_freq The fewest times a pattern must occur to be kept, at
#'   least 1. Default `2`.
#' @param min_support The smallest share of sequences that must contain a
#'   pattern for it to be kept, between 0 and 1. Default `0.01`.
#' @param start,end States a pattern must begin (`start`) or end (`end`)
#'   with, one of them at least; wildcards do not count. Default `NULL`
#'   (any).
#' @param contain States a pattern must include, one of them at least.
#'   Default `NULL` (any).
#' @param pattern A template to search for instead of discovering all
#'   patterns, such as `"discuss->*->consensus"`: each `*` is one position
#'   holding any state (`**` is two). Every concrete instance is returned
#'   as its own row (`"discuss->plan->consensus"`, ...). States are joined
#'   by `"->"` (spaces around it are allowed), and a template may mix fixed
#'   states and wildcards or hold no wildcard at all. `type`, `len` and
#'   `gap` are then ignored. An empty template, or one with an empty step
#'   such as `"A->->B"` or a leading or trailing `"->"`, raises an error of
#'   class `nestimate_bad_input`. Default `NULL`.
#' @param group The grouping of the sequences for the group test: the name
#'   of a column stored with them (a covariate of a wide data frame or long
#'   event log, or of a network's metadata), or a vector with one value per
#'   sequence. A `netobject_group` is grouped by its own groups when
#'   `group` is `NULL`. `"last_obs"` takes each sequence's last observed
#'   state as its group and removes every state that ends some sequence
#'   from all sequences, for data whose final state codes the outcome
#'   (such as an appended `"pass"` or `"fail"`), as codyna's `last_obs`
#'   does. Sequences with a missing group are dropped with a warning of
#'   class `nestimate_missing_group`. A column not found, a vector of the
#'   wrong length or a single group raise an error of class
#'   `nestimate_bad_input`. Default `NULL` (no test).
#' @param reference For two groups: the reference group, so the odds ratio
#'   is that of the other group. Default `NULL`, the first group in sorted
#'   order; a value that is not a group raises an error of class
#'   `nestimate_bad_input`.
#' @param adjust The multiplicity correction of the group test's p-values,
#'   one of `p.adjust.methods` (see [stats::p.adjust()]). Default `"BH"`.
#' @param top The number of rows to return, a positive whole number: the
#'   first in the output order. Default `NULL` (all).
#' @param cols The state columns of a wide data frame, as in
#'   [extract_sequences()]. Default `NULL`.
#' @param actor,action,time The actor, action and timestamp columns of a
#'   long event log, as in [extract_sequences()]: setting `action` declares
#'   `x` an event log. Default `NULL`.
#' @param ... Further arguments passed to [prepare()] for a long event log,
#'   such as `time_threshold` or `order`.
#'
#' @section Group test:
#' With `group`, each pattern's presence (in a sequence or not) is
#' cross-tabulated against the group. For two groups, Fisher's exact test
#' gives the p-value, the conditional maximum-likelihood odds ratio of
#' containing the pattern in the non-reference group, and its 95%
#' confidence interval. A pattern that occurs in one group only has an
#' infinite (or zero) odds ratio, with a finite lower (or upper) limit; a
#' pattern present in every sequence has none (`NA`). For three or more
#' groups, the Pearson chi-square of the presence-by-group table gives the
#' p-value (G - 1 degrees of freedom), with Cram\enc{é}{e}r's V as the
#' effect size. Unlike a comparison of counts against an equal split, both
#' tests account for unequal group sizes. The p-values are adjusted with
#' `adjust` over all patterns that pass the filters, before `top` is
#' applied, and the rows are ordered by `p_value`.
#'
#' @return A data frame with one row per pattern, ordered by `support`
#'   (then `occurrences`, then `pattern`):
#'   \describe{
#'     \item{pattern}{The states joined by `"->"`; `*` marks a wildcard
#'       position.}
#'     \item{length}{The number of positions the pattern spans, gaps
#'       included.}
#'     \item{occurrences}{The number of times it occurs.}
#'     \item{proportion}{`occurrences` as a share of all patterns of the
#'       same length.}
#'     \item{n_sequences}{The number of sequences that contain it.}
#'     \item{support}{`n_sequences` divided by the number of sequences.}
#'     \item{lift}{Support relative to the independent presence of its
#'       states.}
#'   }
#'   With `group`, the columns `support_<group>` (the support within each
#'   group) follow; then `odds_ratio`, `or_lower` and `or_upper` for two
#'   groups, or `chisq` and `cramers_v` for more; then `p_value` and
#'   `p_adj`. The rows are then ordered by `p_value` (then `support`).
#'
#'   A state name that contains `"->"` or consists only of asterisks
#'   (reserved for wildcards) raises an error of class
#'   `nestimate_bad_input`.
#'
#' @references
#' Agrawal, R., & Srikant, R. (1995). Mining sequential patterns.
#' \emph{Proceedings of the Eleventh International Conference on Data
#' Engineering}, 3--14.
#'
#' Fisher, R. A. (1935). The logic of inductive inference. \emph{Journal of
#' the Royal Statistical Society}, 98(1), 39--82.
#'
#' @examples
#' sequence_patterns(trajectories, len = 2:3, top = 10)
#'
#' sequence_patterns(trajectories, type = "gapped", gap = 1)
#'
#' sequence_patterns(trajectories, type = "repeated", len = 3:5,
#'                   start = "Disengaged")
#'
#' sequence_patterns(group_regulation_long, actor = "Actor",
#'                   action = "Action", time = "Time", len = 2,
#'                   group = "Achiever", top = 10)
#'
#' sequence_patterns(trajectories, pattern = "Active->*->Disengaged")
#'
#' @seealso \code{\link{pattern_outcome}}, \code{\link{sequence_compare}}
#' @export
sequence_patterns <- function(x, type = c("ngram", "gapped", "repeated"),
                              len = 2:5, gap = 1:3, min_freq = 2,
                              min_support = 0.01, start = NULL, end = NULL,
                              contain = NULL, pattern = NULL, group = NULL,
                              reference = NULL, adjust = "BH", top = NULL,
                              cols = NULL, actor = NULL, action = NULL,
                              time = NULL, ...) {
  type <- match.arg(type)
  adjust <- match.arg(adjust, stats::p.adjust.methods)
  .sp_check_args(len, gap, min_freq, min_support, start, end, contain)
  .assert(
    "`top` must be NULL or a single positive whole number" = is.null(top) ||
      (is.numeric(top) && length(top) == 1L && top >= 1 && top == round(top)),
    "`pattern` must be NULL or a single string" = is.null(pattern) ||
      (is.character(pattern) && length(pattern) == 1L && !is.na(pattern))
  )
  if (is.null(group) && inherits(x, "netobject_group")) {
    group <- .group_col_name(x)
  }
  last_obs <- identical(group, "last_obs")
  group_is_column <- !last_obs && is.character(group) && length(group) == 1L
  seqs <- .as_sequences(x, cols = cols, actor = actor, action = action,
                        time = time,
                        exclude = if (group_is_column) group, ...)
  states <- seqs$states

  labels <- NULL
  if (last_obs) {
    split <- .sp_last_obs(states)
    states <- split$states
    labels <- split$labels
  } else if (!is.null(group)) {
    labels <- .sp_group_labels(group, group_is_column, seqs$units,
                               nrow(states))
  }
  if (!is.null(labels)) {
    keep <- !is.na(labels)
    if (!all(keep)) {
      warning(warningCondition(
        sprintf("%d of %d sequences have a missing group and are dropped.",
                sum(!keep), length(keep)),
        class = "nestimate_missing_group"))
    }
    states <- states[keep, , drop = FALSE]
    labels <- labels[keep]
  }

  found <- .sp_find(states, type = type, len = len, gap = gap,
                    template = pattern)
  out <- .sp_filter(found$table, min_freq = min_freq,
                    min_support = min_support, start = start, end = end,
                    contain = contain)
  if (!is.null(labels)) {
    out <- .sp_group_test(out, found, labels, reference = reference,
                          adjust = adjust)
  }
  out$states <- NULL
  if (!is.null(top)) out <- utils::head(out, top)
  rownames(out) <- NULL
  out
}


#' Resolve the group of each sequence as a character vector
#' @noRd
.sp_group_labels <- function(group, group_is_column, units, n) {
  if (group_is_column) {
    if (!group %in% names(units)) {
      stop(errorCondition(
        sprintf("group column `%s` not found with the sequences. %s",
                group,
                if (ncol(units) == 0L) {
                  "This input stores no covariates; pass `group` as a vector."
                } else {
                  sprintf("Available: %s.", paste(names(units), collapse = ", "))
                }),
        class = "nestimate_bad_input", call = NULL))
    }
    return(as.character(units[[group]]))
  }
  if (length(group) != n) {
    stop(errorCondition(
      sprintf("`group` has %d values; the data has %d sequences.",
              length(group), n),
      class = "nestimate_bad_input", call = NULL))
  }
  as.character(group)
}


#' Sequences containing each pattern, per group (patterns x groups)
#' @noRd
.sp_group_counts <- function(found, patterns, labels, levels) {
  occ <- found$occurrences
  ids <- match(patterns, found$table$pattern)
  matrix(vapply(levels, \(lv) {
    in_level <- which(labels == lv)
    vapply(ids, \(pid) sum(occ$seq[occ$id == pid] %in% in_level),
           numeric(1))
  }, numeric(length(ids))), nrow = length(ids), ncol = length(levels),
  dimnames = list(NULL, levels))
}


#' Add per-group support and the presence-by-group test
#' @noRd
.sp_group_test <- function(tab, found, labels, reference = NULL,
                           adjust = "BH") {
  values <- sort(unique(labels))
  if (length(values) < 2L) {
    stop(errorCondition(
      sprintf("`group` must have at least two values; it has %d.",
              length(values)),
      class = "nestimate_bad_input", call = NULL))
  }
  if (!is.null(reference) && !reference %in% values) {
    stop(errorCondition(
      sprintf("`reference` must be one of: %s.",
              paste(values, collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }
  reference <- reference %||% values[1L]
  levels <- c(reference, setdiff(values, reference))
  sizes <- as.numeric(table(factor(labels, levels = levels)))
  present <- .sp_group_counts(found, tab$pattern, labels, levels)

  # rep(each =) rather than matrix(byrow =): no warning at zero patterns
  size_mat <- matrix(rep(sizes, each = nrow(present)), nrow(present),
                     length(sizes))
  support_by <- as.data.frame(present / size_mat, check.names = FALSE)
  names(support_by) <- paste0("support_", levels)
  out <- cbind(tab, support_by)

  if (length(levels) == 2L) {
    tests <- lapply(seq_len(nrow(present)), \(i) {
      # rows: other group, reference; columns: contains, does not
      ft <- stats::fisher.test(matrix(
        c(present[i, 2L], present[i, 1L],
          sizes[2L] - present[i, 2L], sizes[1L] - present[i, 1L]), 2L))
      c(unname(ft$estimate), ft$conf.int, ft$p.value)
    })
    tests <- do.call(rbind, c(list(matrix(numeric(0), 0L, 4L)), tests))
    # A pattern in every sequence has no contrast; fisher.test() would
    # report an odds ratio of 0 with an unbounded interval
    no_contrast <- out$n_sequences >= sum(sizes)
    tests[no_contrast, 1:3] <- NA_real_
    out$odds_ratio <- tests[, 1L]
    out$or_lower <- tests[, 2L]
    out$or_upper <- tests[, 3L]
    out$p_value <- tests[, 4L]
  } else {
    n_total <- sum(sizes)
    expected <- outer(out$n_sequences, sizes) / n_total
    absent_expected <- outer(n_total - out$n_sequences, sizes) / n_total
    absent <- size_mat - present
    # A pattern in every (or no) sequence has zero expected cells, whose
    # observed counts are zero too: they contribute nothing, not 0 / 0
    cell <- \(obs, exp) ifelse(exp > 0, (obs - exp)^2 / exp, 0)
    out$chisq <- rowSums(cell(present, expected) +
                           cell(absent, absent_expected))
    out$cramers_v <- sqrt(out$chisq / n_total)
    out$p_value <- stats::pchisq(out$chisq, length(levels) - 1L,
                                 lower.tail = FALSE)
  }
  out$p_adj <- stats::p.adjust(out$p_value, method = adjust)
  out[order(out$p_value, -out$support, out$pattern), , drop = FALSE]
}


#' Parse a search template such as "A->*->B" or "A->**->B->C"
#'
#' Each "*" is one position of any state.
#' @return list(span, positions of the fixed states, the fixed states)
#' @noRd
.sp_parse_template <- function(template) {
  # "A -> * -> B" (the HON notation) is read like "A->*->B"
  template <- trimws(template)
  tokens <- trimws(strsplit(template, "->", fixed = TRUE)[[1L]])
  # strsplit() drops a trailing empty token, so test the end explicitly
  if (length(tokens) == 0L || any(!nzchar(tokens)) ||
      endsWith(template, "->")) {
    stop(errorCondition(
      sprintf("`pattern` \"%s\" is not of the form \"A->*->B\".", template),
      class = "nestimate_bad_input", call = NULL))
  }
  wild <- grepl("^\\*+$", tokens)
  widths <- ifelse(wild, nchar(tokens), 1L)
  starts <- cumsum(widths) - widths + 1L
  list(span = sum(widths), positions = starts[!wild], states = tokens[!wild])
}


#' Split each sequence into its last observed state and the rest
#'
#' codyna's `last_obs` outcome: the last observed state of each sequence is
#' its outcome, and every state that ends some sequence is removed from all
#' sequences. Meant for sequences whose final state codes the outcome (an
#' appended "pass"/"fail").
#' @return list(states = matrix without those states, labels = last states)
#' @noRd
.sp_last_obs <- function(states) {
  observed <- !is.na(states)
  last_pos <- max.col(observed, ties.method = "last")
  labels <- states[cbind(seq_len(nrow(states)), last_pos)]
  labels[rowSums(observed) == 0L] <- NA_character_
  states[states %in% stats::na.omit(unique(labels))] <- NA_character_
  list(states = states, labels = labels)
}


#' Validate the pattern arguments shared by the pattern verbs
#' @noRd
.sp_check_args <- function(len, gap, min_freq, min_support,
                           start = NULL, end = NULL, contain = NULL) {
  is_counts <- function(v) is.numeric(v) && length(v) >= 1L &&
    all(is.finite(v)) && all(v >= 1) && all(v == round(v))
  .assert(
    "`len` must be positive whole numbers" = is_counts(len),
    "`gap` must be positive whole numbers" = is_counts(gap),
    "`min_freq` must be a single number >= 1" =
      is.numeric(min_freq) && length(min_freq) == 1L && min_freq >= 1,
    "`min_support` must be a single number in [0, 1]" =
      is.numeric(min_support) && length(min_support) == 1L &&
      min_support >= 0 && min_support <= 1,
    "`start` must be NULL or character" = is.null(start) || is.character(start),
    "`end` must be NULL or character" = is.null(end) || is.character(end),
    "`contain` must be NULL or character" =
      is.null(contain) || is.character(contain)
  )
  invisible(TRUE)
}


#' Find every pattern occurrence and tabulate it
#'
#' @return A list with \code{table} (one row per distinct pattern, with a
#'   list column \code{states}) and \code{occurrences} (data.frame of
#'   \code{seq}, integer pattern \code{id}; one row per occurrence).
#' @noRd
.sp_find <- function(states, type, len, gap, template = NULL) {
  alphabet <- sort(unique(states[!is.na(states)]))
  if (any(grepl("->", alphabet, fixed = TRUE))) {
    stop(errorCondition(
      "state names must not contain \"->\", the pattern separator.",
      class = "nestimate_bad_input", call = NULL))
  }
  if (any(grepl("^\\*+$", alphabet))) {
    stop(errorCondition(
      "state names must not consist only of asterisks, which mark wildcards.",
      class = "nestimate_bad_input", call = NULL))
  }
  n <- nrow(states)
  m <- ncol(states)
  search <- if (is.null(template)) NULL else .sp_parse_template(template)
  gapped <- is.null(search) && identical(type, "gapped")
  # a repeated length or gap would count every occurrence twice
  sizes <- unique(if (!is.null(search)) {
    search$span
  } else if (gapped) {
    as.integer(gap)
  } else {
    as.integer(len)
  })

  # One window spec per (size, start position); each yields the occurrences
  # of every sequence at that position
  windows <- unlist(lapply(sizes, \(s) {
    span <- if (gapped) s + 2L else s
    n_pos <- m - span + 1L
    if (n_pos < 1L) return(list())
    lapply(seq_len(n_pos), \(i) list(size = s, span = span, start = i))
  }), recursive = FALSE)

  hits <- lapply(windows, \(w) {
    idx <- if (gapped) {
      c(w$start, w$start + w$size + 1L)
    } else {
      w$start + seq_len(w$size) - 1L
    }
    block <- states[, idx, drop = FALSE]
    # the whole window must be observed, the gap positions included
    window <- states[, w$start + seq_len(w$span) - 1L, drop = FALSE]
    valid <- rowSums(is.na(window)) == 0L
    if (!is.null(search)) {
      fixed <- matrix(search$states, nrow(block), length(search$states),
                      byrow = TRUE)
      valid <- valid & rowSums(block[, search$positions, drop = FALSE] ==
                                 fixed) == length(search$states)
    } else if (identical(type, "repeated")) {
      valid <- valid & rowSums(block == block[, 1L]) == length(idx)
    }
    valid <- which(valid)
    if (length(valid) == 0L) return(NULL)
    block <- block[valid, , drop = FALSE]
    sep <- if (gapped) {
      paste0("->", strrep("*", w$size), "->")
    } else {
      "->"
    }
    list(seq = valid,
         pattern = do.call(paste, c(as.data.frame(block), sep = sep)),
         span = rep(w$span, length(valid)))
  })
  hits <- Filter(Negate(is.null), hits)

  empty <- data.frame(pattern = character(0), length = integer(0),
                      occurrences = integer(0), proportion = numeric(0),
                      n_sequences = integer(0), support = numeric(0),
                      lift = numeric(0))
  if (length(hits) == 0L) {
    empty$states <- list()
    return(list(table = empty,
                occurrences = data.frame(seq = integer(0), id = integer(0),
                                         times = integer(0))))
  }

  seq_id <- unlist(lapply(hits, `[[`, "seq"), use.names = FALSE)
  pattern <- unlist(lapply(hits, `[[`, "pattern"), use.names = FALSE)
  span <- unlist(lapply(hits, `[[`, "span"), use.names = FALSE)

  patterns <- unique(pattern)
  id <- match(pattern, patterns)
  k <- length(patterns)
  occurrences <- tabulate(id, nbins = k)
  pair <- (as.numeric(id) - 1) * n + seq_id
  first_in_seq <- !duplicated(pair)
  times_in_seq <- tabulate(match(pair, pair[first_in_seq]),
                           nbins = sum(first_in_seq))
  n_sequences <- tabulate(id[first_in_seq], nbins = k)
  lengths <- span[match(seq_len(k), id)]

  pattern_states <- lapply(strsplit(patterns, "->", fixed = TRUE),
                           \(s) s[!grepl("^\\*+$", s)])
  state_support <- vapply(alphabet, \(a) {
    mean(rowSums(states == a, na.rm = TRUE) > 0)
  }, numeric(1))
  support <- n_sequences / n
  expected <- vapply(pattern_states, \(s) prod(state_support[s]), numeric(1))
  length_total <- tapply(occurrences, lengths, sum)

  table <- data.frame(
    pattern = patterns,
    length = as.integer(lengths),
    occurrences = as.integer(occurrences),
    proportion = occurrences / as.numeric(length_total[as.character(lengths)]),
    n_sequences = as.integer(n_sequences),
    support = support,
    lift = support / expected,
    stringsAsFactors = FALSE
  )
  table$states <- pattern_states
  list(table = table,
       occurrences = data.frame(seq = seq_id[first_in_seq],
                                id = id[first_in_seq],
                                times = times_in_seq))
}


#' Apply the pattern filters and the deterministic ordering
#' @noRd
.sp_filter <- function(table, min_freq, min_support, start = NULL,
                       end = NULL, contain = NULL, sort = TRUE) {
  keep <- table$occurrences >= min_freq & table$support >= min_support
  first <- vapply(table$states, \(s) s[1L], character(1))
  last <- vapply(table$states, \(s) s[length(s)], character(1))
  if (!is.null(start)) keep <- keep & first %in% start
  if (!is.null(end)) keep <- keep & last %in% end
  if (!is.null(contain)) {
    keep <- keep & vapply(table$states, \(s) any(s %in% contain), logical(1))
  }
  out <- table[keep, , drop = FALSE]
  # sort = FALSE keeps first-appearance order (by length, position, row),
  # which pattern_outcome() needs to break ties the way codyna does
  if (sort) {
    out <- out[order(-out$support, -out$occurrences, out$pattern), ,
               drop = FALSE]
  }
  rownames(out) <- NULL
  out
}
