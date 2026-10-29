#' Discover Sequential Patterns
#'
#' Finds the ordered patterns that occur in a set of sequences and reports,
#' for each one, how often it occurs, how many sequences contain it, and its
#' lift. Three pattern types are supported:
#' \describe{
#'   \item{\code{"ngram"}}{Contiguous runs of \code{len} states,
#'     e.g. \code{A->B->C}.}
#'   \item{\code{"gapped"}}{Two states separated by \code{gap} positions of
#'     any state, e.g. \code{A->*->B} (gap 1) or \code{A->**->B} (gap 2).}
#'   \item{\code{"repeated"}}{Runs of \code{len} identical states,
#'     e.g. \code{A->A->A}.}
#' }
#' A window that contains a missing value is skipped, so a pattern never
#' spans a gap or the end of a sequence.
#'
#' @details
#' \code{support} is the share of sequences that contain the pattern at
#' least once. \code{lift} compares it with the support expected if the
#' states in the pattern appeared in sequences independently:
#' \code{support / prod(state_support)}, where the support of a state is the
#' share of sequences in which it occurs and wildcard positions are left out.
#' A lift above 1 means the states co-occur in this order more often than
#' independent presence predicts. It does not account for the transition
#' structure: a pattern that the transition probabilities make likely can
#' have a lift above 1.
#'
#' \code{proportion} is computed over all patterns of the same length before
#' \code{min_freq}, \code{min_support}, \code{start}, \code{end} and
#' \code{contain} are applied, so it does not change with the filters.
#'
#' @param x Sequence data, in any form accepted by
#'   \code{\link{extract_sequences}()}: prepared data, a long event log
#'   (with \code{action}), a wide data frame or character matrix, a
#'   \code{netobject}, \code{netobject_group} or \code{tna} model.
#' @param type Pattern type: \code{"ngram"} (default), \code{"gapped"} or
#'   \code{"repeated"}.
#' @param len Integer vector of pattern lengths for \code{"ngram"} and
#'   \code{"repeated"}. Default \code{2:5}.
#' @param gap Integer vector of gap sizes for \code{"gapped"}. Default
#'   \code{1:3}.
#' @param min_freq Minimum number of occurrences. Default \code{2}.
#' @param min_support Minimum share of sequences containing the pattern.
#'   Default \code{0.01}.
#' @param start,end Optional character vectors of states. Keep only patterns
#'   that begin (\code{start}) or end (\code{end}) with one of them.
#' @param contain Optional character vector of states. Keep only patterns
#'   that include at least one of them.
#' @param pattern Optional template to search for instead of discovering
#'   all patterns, e.g. \code{"discuss->*->consensus"}: each \code{*} is one
#'   position holding any state (\code{**} is two). Every concrete instance
#'   is returned as its own row (\code{"discuss->plan->consensus"}, ...).
#'   \code{type}, \code{len} and \code{gap} are then ignored.
#' @param group Optional grouping of the sequences: the name of a column
#'   stored with them (a covariate of a wide data frame or long event log,
#'   or of a network's metadata), or a vector with one value per sequence.
#'   For a \code{netobject_group} its own grouping is used when \code{group}
#'   is \code{NULL}. Adds a test of pattern presence by group; see Details.
#'   \code{group = "last_obs"} uses each sequence's last observed state as
#'   its group and removes the states that end some sequence from all
#'   sequences (for data whose final state codes the outcome, e.g. an
#'   appended \code{"pass"}/\code{"fail"}; codyna's \code{last_obs}).
#' @param reference For two groups, the reference group; the odds ratio is
#'   for the other group. Default: the first group in sorted order.
#' @param adjust Multiplicity correction for the group test, passed to
#'   \code{\link[stats]{p.adjust}}. Default \code{"BH"}.
#' @param top Optional number of rows to return (the first \code{top} in the
#'   output order). Default \code{NULL}, all rows.
#' @param cols State columns of a wide data frame; see
#'   \code{\link{extract_sequences}()}.
#' @param ... Passed to \code{\link{extract_sequences}()}, e.g.
#'   \code{actor}, \code{action} and \code{time} for a long event log.
#'
#' @section Group test:
#' With \code{group}, each pattern's presence (in a sequence or not) is
#' cross-tabulated against the group. For two groups, Fisher's exact test
#' gives the p-value, the conditional maximum-likelihood odds ratio of
#' containing the pattern in the non-reference group, and its 95\%
#' confidence interval. A pattern that occurs in one group only has an
#' infinite (or zero) odds ratio, with a finite lower (or upper) limit; a
#' pattern present in every sequence has none (\code{NA}). For three or more
#' groups, the Pearson chi-square of
#' the presence-by-group table gives the p-value (G - 1 degrees of freedom),
#' with Cramér's V as the effect size. Unlike a comparison of counts against
#' an equal split, both tests account for unequal group sizes. p-values are
#' adjusted across the returned patterns' family (all patterns passing the
#' filters) with \code{adjust}. The rows are then ordered by \code{p_value}.
#'
#' @return A \code{data.frame} with one row per pattern, ordered by
#'   \code{support} (then \code{occurrences}, then \code{pattern}):
#'   \describe{
#'     \item{pattern}{States joined by \code{"->"}; \code{*} marks a
#'       wildcard position.}
#'     \item{length}{Number of positions the pattern spans.}
#'     \item{occurrences}{Total number of times it occurs.}
#'     \item{proportion}{\code{occurrences} as a share of all patterns of
#'       the same length.}
#'     \item{n_sequences}{Number of sequences that contain it.}
#'     \item{support}{\code{n_sequences} divided by the number of
#'       sequences.}
#'     \item{lift}{Support relative to independent presence of its states.}
#'   }
#'   With \code{group}, also \code{support_<group>} (the support within each
#'   group); for two groups \code{odds_ratio}, \code{or_lower},
#'   \code{or_upper}; for more, \code{chisq} and \code{cramers_v}; then
#'   \code{p_value} and \code{p_adj}. Rows are ordered by \code{p_value}
#'   (then \code{support}).
#'
#'   Conditions: \code{nestimate_bad_input} for a state name containing
#'   \code{"->"}, a missing group column, a group vector of the wrong length,
#'   or a group with fewer than two values; warning
#'   \code{nestimate_missing_group} when sequences with a missing group are
#'   dropped.
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
                              cols = NULL, ...) {
  type <- match.arg(type)
  adjust <- match.arg(adjust, stats::p.adjust.methods)
  .sp_check_args(len, gap, min_freq, min_support, start, end, contain)
  stopifnot(
    "`top` must be NULL or a single positive whole number" = is.null(top) ||
      (is.numeric(top) && length(top) == 1L && top >= 1 && top == round(top)),
    "`pattern` must be NULL or a single string" = is.null(pattern) ||
      (is.character(pattern) && length(pattern) == 1L && !is.na(pattern))
  )
  if (is.null(group) && inherits(x, "netobject_group")) {
    group <- attr(x, "group_col") %||% "group"
  }
  last_obs <- identical(group, "last_obs")
  group_is_column <- !last_obs && is.character(group) && length(group) == 1L
  seqs <- .as_sequences(x, cols = cols,
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

  support_by <- as.data.frame(
    present / matrix(sizes, nrow(present), length(sizes), byrow = TRUE),
    check.names = FALSE)
  names(support_by) <- paste0("support_", levels)
  out <- cbind(tab, support_by)
  if (nrow(out) == 0L) return(out)

  if (length(levels) == 2L) {
    tests <- lapply(seq_len(nrow(present)), \(i) {
      # rows: other group, reference; columns: contains, does not
      ft <- stats::fisher.test(matrix(
        c(present[i, 2L], present[i, 1L],
          sizes[2L] - present[i, 2L], sizes[1L] - present[i, 1L]), 2L))
      c(unname(ft$estimate), ft$conf.int, ft$p.value)
    })
    tests <- do.call(rbind, tests)
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
    absent <- matrix(sizes, nrow(present), length(sizes), byrow = TRUE) -
      present
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
  tokens <- strsplit(template, "->", fixed = TRUE)[[1L]]
  if (length(tokens) == 0L || any(!nzchar(tokens))) {
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
  stopifnot(
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
  n <- nrow(states)
  m <- ncol(states)
  search <- if (is.null(template)) NULL else .sp_parse_template(template)
  gapped <- is.null(search) && identical(type, "gapped")
  sizes <- if (!is.null(search)) {
    search$span
  } else if (gapped) {
    as.integer(gap)
  } else {
    as.integer(len)
  }

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
    valid <- rowSums(is.na(block)) == 0L
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
