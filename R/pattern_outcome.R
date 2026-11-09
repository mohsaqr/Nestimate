#' Model a Binary Outcome From Sequential Patterns
#'
#' @description
#' `pattern_outcome()` selects the `n` sequential patterns most associated
#' with a two-valued outcome and fits one logistic regression of the
#' outcome on their presence (or their counts), returning one row per
#' selected pattern with its log odds ratio, odds ratio and Wald interval.
#' It reads any sequence data [extract_sequences()] reads, takes the
#' outcome from a stored covariate, a vector or each sequence's last state,
#' and can adjust for covariates and add a random intercept per actor.
#'
#' Each estimate is the association of one pattern with the outcome,
#' adjusted for the other selected patterns and the covariates; it is not a
#' causal effect. The patterns are chosen by their association with the
#' outcome in the same data the model is fitted to, so the p-values do not
#' account for that selection and overstate the evidence: treat them as a
#' ranking, or confirm the selected patterns on new data. A pattern and a
#' longer pattern containing it often occur together, and each coefficient
#' is then adjusted for the other; [sequence_patterns()] shows how often
#' each pattern occurs, overall and by group.
#'
#' @details
#' `pattern_outcome()` reproduces `codyna::analyze_outcome()`: the same
#' pattern discovery and filters, the same ranking and the same model, so
#' the coefficients agree. The steps:
#' \enumerate{
#'   \item Patterns are found as in [sequence_patterns()] and filtered by
#'     `min_freq`, `min_support`, `start`, `end` and `contain`. Patterns
#'     absent from all sequences of one outcome value are dropped.
#'   \item Each pattern gets the chi-square statistic of the number of
#'     sequences containing it under each outcome value, against an equal
#'     split across the values (codyna's statistic; it does not weight by
#'     the size of each outcome group).
#'   \item The top `n` patterns by `priority` enter one model,
#'     `outcome ~ adjust + pattern_1 + ... + pattern_n`, as presence
#'     indicators (or occurrence counts with `occurrences = TRUE`). Ties
#'     keep codyna's order: more occurrences first, then first appearance.
#'     With `group`, a random intercept per group is added
#'     (`lme4::glmer()`); otherwise the model is a plain `glm()`.
#' }
#'
#' @param x The sequence data, in any form [extract_sequences()] accepts:
#'   prepared data, a long event log, a wide data frame, a `netobject`,
#'   `netobject_group` or `tna` model.
#' @param outcome The outcome, with exactly two distinct values: the name of
#'   a column stored with the sequences (a covariate of a wide data frame or
#'   long event log, or of a network's metadata), or a vector with one value
#'   per sequence; use a vector for a model that stores no covariates, such
#'   as a `tna` model. The default, `"last_obs"`, uses each sequence's last
#'   observed state as the outcome and removes the states that end some
#'   sequence from all sequences, as codyna does; it suits data whose final
#'   state codes the outcome, such as an appended `"pass"` or `"fail"`. An
#'   outcome without exactly two values, a vector of the wrong length or a
#'   column that is not stored with the sequences raises an error of class
#'   `nestimate_bad_input`.
#' @param type,len,gap,min_freq,min_support,start,end,contain The pattern
#'   definition and filters; see [sequence_patterns()]. The defaults are
#'   codyna's: `type = "ngram"`, `len = 1:2`, `gap = 1`, `min_freq = 5` and
#'   `min_support = 0.01`. When no pattern passes the filters under both
#'   outcome values, an error of class `nestimate_no_patterns` is raised.
#' @param n The number of patterns entered into the model, a positive whole
#'   number. Default `10`.
#' @param priority The columns ranking the patterns for selection, applied
#'   in order: `"chisq"` (default), `"p_value"`, `"support"`,
#'   `"n_sequences"`, `"occurrences"`, `"proportion"`, `"lift"` and
#'   `"length"`.
#' @param desc Whether each `priority` column ranks in descending order,
#'   recycled to the length of `priority`. Default `TRUE`, so `"p_value"`
#'   needs `desc = FALSE` to select the smallest p-values.
#' @param occurrences Whether each pattern enters the model as its number of
#'   occurrences in the sequence instead of its presence (0/1). Default
#'   `FALSE`.
#' @param adjust Covariate columns added to the model. Default `NULL`. A
#'   name that is not a column stored with the sequences raises an error of
#'   class `nestimate_bad_input`.
#' @param adjust_length Whether to add the number of observed states per
#'   sequence as a covariate. Default `FALSE`, as in codyna.
#' @param group The column identifying the actor (or other cluster) each
#'   sequence belongs to; it adds a random intercept per group. It needs
#'   \pkg{lme4}, otherwise an error of class `nestimate_no_lme4` is raised.
#'   Default `NULL`. A name that is not a column stored with the sequences
#'   raises an error of class `nestimate_bad_input`.
#' @param reference The reference outcome value; the odds ratios are for
#'   the other value. Default `NULL`, the first value in sorted order; a
#'   value that does not occur raises an error of class
#'   `nestimate_bad_input`.
#' @param ci_level The confidence level of the Wald intervals, a number
#'   between 0 and 1. Default `0.95`.
#' @param cols The state columns of a wide data frame; see
#'   [extract_sequences()]. Default `NULL`.
#' @param actor,action,time The actor, action and timestamp columns of a
#'   long event log, as in [extract_sequences()]: setting `action` declares
#'   `x` an event log. Default `NULL`.
#' @param ... Further arguments passed to [prepare()] for a long event log,
#'   such as `time_threshold` or `order`.
#'
#' @return A data frame with one row per selected pattern, ordered by
#'   pattern name, with the columns `pattern`, `length`, `occurrences`,
#'   `n_sequences`, `support` and `lift` (as in [sequence_patterns()]); the
#'   share of the sequences of each outcome value that contain the pattern
#'   (`support_<value>`); `chisq` (the selection statistic); `estimate` (the
#'   log odds ratio, adjusted for the other terms), `std_error`, `ci_lower`,
#'   `ci_upper` and `p_value` (Wald z); and `odds_ratio`, `or_lower` and
#'   `or_upper`.
#'
#'   Sequences with a missing outcome, covariate or group are dropped with a
#'   warning of class `nestimate_missing_outcome`. A selected pattern
#'   collinear with the other terms gets `NA` estimates (and `NA` standard
#'   error, intervals, p-value and odds ratios), its descriptive columns
#'   kept, and a warning of class
#'   `nestimate_aliased_pattern`. Warnings from the model fit (such as
#'   separation) are passed on unchanged.
#'
#' @references
#' Saqr, M., Tikka, S., & L\enc{ó}{o}pez-Pernas, S. (2026). codyna: Complex Dynamic
#' Systems. R package version 1.0.0.
#'
#' @examples
#' pattern_outcome(group_regulation_long, outcome = "Achiever",
#'                 actor = "Actor", action = "Action", time = "Time")
#'
#' @seealso \code{\link{sequence_patterns}}
#' @export
pattern_outcome <- function(x, outcome = "last_obs",
                            type = c("ngram", "gapped", "repeated"),
                            len = 1:2, gap = 1, min_freq = 5,
                            min_support = 0.01, start = NULL, end = NULL,
                            contain = NULL, n = 10, priority = "chisq",
                            desc = TRUE, occurrences = FALSE, adjust = NULL,
                            adjust_length = FALSE, group = NULL,
                            reference = NULL, ci_level = 0.95, cols = NULL,
                            actor = NULL, action = NULL, time = NULL, ...) {
  type <- match.arg(type)
  priority <- match.arg(priority, c("chisq", "p_value", "support",
                                    "n_sequences", "occurrences",
                                    "proportion", "lift", "length"),
                        several.ok = TRUE)
  .assert(
    "`outcome` must be a column name or a vector of outcome values" =
      is.atomic(outcome) && length(outcome) >= 1L,
    "`n` must be a single positive whole number" = is.numeric(n) &&
      length(n) == 1L && n >= 1 && n == round(n),
    "`desc` must be logical" = is.logical(desc) && !anyNA(desc),
    "`occurrences` must be TRUE or FALSE" =
      isTRUE(occurrences) || isFALSE(occurrences),
    "`adjust` must be NULL or character" = is.null(adjust) ||
      is.character(adjust),
    "`adjust_length` must be TRUE or FALSE" =
      isTRUE(adjust_length) || isFALSE(adjust_length),
    "`group` must be NULL or a single column name" = is.null(group) ||
      (is.character(group) && length(group) == 1L),
    "`ci_level` must be a single number in (0, 1)" = is.numeric(ci_level) &&
      length(ci_level) == 1L && ci_level > 0 && ci_level < 1
  )
  .sp_check_args(len, gap, min_freq, min_support, start, end, contain)
  if (!is.null(group) && !requireNamespace("lme4", quietly = TRUE)) {
    stop(errorCondition("`group` needs the lme4 package; install it.",
                        class = "nestimate_no_lme4", call = NULL))
  }

  last_obs <- identical(outcome, "last_obs")
  outcome_is_column <- !last_obs && is.character(outcome) &&
    length(outcome) == 1L
  seqs <- .as_sequences(x, cols = cols, actor = actor, action = action,
                        time = time,
                        exclude = c(if (outcome_is_column) outcome, adjust,
                                    group), ...)
  units <- seqs$units
  if (last_obs) {
    split <- .sp_last_obs(seqs$states)
    seqs$states <- split$states
    outcome_values <- split$labels
  } else if (!outcome_is_column) {
    if (length(outcome) != nrow(seqs$states)) {
      stop(errorCondition(
        sprintf("`outcome` has %d values; the data has %d sequences.",
                length(outcome), nrow(seqs$states)),
        class = "nestimate_bad_input", call = NULL))
    }
    outcome_values <- outcome
  }
  needed <- c(if (outcome_is_column) outcome, adjust, group)
  .po_check_columns(x, setdiff(needed, names(units)), names(units))
  if (outcome_is_column) outcome_values <- units[[outcome]]
  # Keep a generated outcome separate from the user's covariates: even a
  # column named .outcome can be an adjustment or random-intercept group.
  complete <- stats::complete.cases(outcome_values)
  if (length(needed) > 0L) {
    complete <- complete & stats::complete.cases(units[needed])
  }
  if (!all(complete)) {
    warning(warningCondition(
      sprintf(paste("%d of %d sequences have a missing outcome, covariate",
                    "or group and are dropped."),
              sum(!complete), length(complete)),
      class = "nestimate_missing_outcome"))
  }
  states <- seqs$states[complete, , drop = FALSE]
  units <- units[complete, , drop = FALSE]
  y <- .po_outcome(outcome_values[complete], reference)

  # ---- Discover, filter, rank and select (codyna's order of operations) ----
  found <- .sp_find(states, type = type, len = len, gap = gap)
  tab <- .sp_filter(found$table, min_freq = min_freq,
                    min_support = min_support, start = start, end = end,
                    contain = contain, sort = FALSE)
  occ <- found$occurrences
  # Sequences containing each pattern, per outcome value (patterns x values)
  by_level <- .sp_group_counts(found, tab$pattern, y$label, y$levels)
  in_all <- rowSums(by_level > 0) == ncol(by_level)
  tab <- tab[in_all, , drop = FALSE]
  by_level <- by_level[in_all, , drop = FALSE]
  if (nrow(tab) == 0L) {
    stop(errorCondition(
      "no pattern passes the filters under both outcome values.",
      class = "nestimate_no_patterns", call = NULL))
  }
  expected <- tab$n_sequences / ncol(by_level)
  tab$chisq <- rowSums((by_level - expected)^2 / expected)
  tab$p_value <- stats::pchisq(tab$chisq, ncol(by_level) - 1L,
                               lower.tail = FALSE)

  # Ties keep codyna's order: occurrences (desc), then first appearance
  base_order <- order(-tab$occurrences, method = "radix")
  tab <- tab[base_order, , drop = FALSE]
  by_level <- by_level[base_order, , drop = FALSE]
  desc <- rep_len(desc, length(priority))
  keys <- Map(\(col, d) if (d) -xtfrm(tab[[col]]) else xtfrm(tab[[col]]),
              priority, desc)
  ranked <- do.call(order, c(unname(keys), list(method = "radix")))
  chosen <- utils::head(ranked, n)
  chosen <- chosen[order(tab$pattern[chosen], method = "radix")]
  tab <- tab[chosen, , drop = FALSE]
  by_level <- by_level[chosen, , drop = FALSE]

  # ---- One model with every selected pattern ----
  sel_ids <- match(tab$pattern, found$table$pattern)
  # Internal names for every model column: a data column can then neither
  # collide with a pattern term nor break the formula (spaces, symbols)
  terms <- sprintf(".p%d", seq_along(sel_ids))
  adjust_terms <- sprintf(".a%d", seq_along(adjust))
  model_data <- units[c(adjust, group)]
  names(model_data) <- c(adjust_terms, if (!is.null(group)) ".g")
  model_data$.y <- y$value
  if (adjust_length) model_data$.length <- rowSums(!is.na(states))
  model_data[terms] <- lapply(sel_ids, \(pid) {
    rows <- occ$id == pid
    value <- integer(nrow(states))
    value[occ$seq[rows]] <- if (occurrences) occ$times[rows] else 1L
    value
  })
  rhs <- c(adjust_terms, if (adjust_length) ".length", terms,
           if (!is.null(group)) "(1 | .g)")
  form <- stats::reformulate(rhs, response = ".y")
  fit <- if (is.null(group)) {
    stats::glm(form, data = model_data, family = stats::binomial())
  } else {
    lme4::glmer(form, data = model_data, family = stats::binomial())
  }
  # An aliased pattern (e.g. present in every sequence) has no estimate:
  # glm returns NA, glmer drops the column; both become NA here
  all_coef <- if (is.null(group)) stats::coef(fit) else lme4::fixef(fit)
  coef_tab <- stats::coef(summary(fit))
  estimate <- unname(all_coef[terms])
  std_error <- unname(coef_tab[match(terms, rownames(coef_tab)), "Std. Error"])
  if (anyNA(estimate)) {
    warning(warningCondition(
      sprintf(paste("pattern(s) collinear with the other terms have no",
                    "estimate: %s."),
              paste(tab$pattern[is.na(estimate)], collapse = ", ")),
      class = "nestimate_aliased_pattern"))
  }

  z <- stats::qnorm(1 - (1 - ci_level) / 2)
  support_by <- as.data.frame(
    by_level / matrix(table(factor(y$label, levels = y$levels)),
                      nrow(by_level), ncol(by_level), byrow = TRUE),
    check.names = FALSE)
  names(support_by) <- paste0("support_", y$levels)

  out <- data.frame(
    tab[c("pattern", "length", "occurrences", "n_sequences", "support",
          "lift")],
    support_by,
    chisq = tab$chisq,
    estimate = estimate,
    std_error = std_error,
    ci_lower = estimate - z * std_error,
    ci_upper = estimate + z * std_error,
    p_value = 2 * stats::pnorm(abs(estimate / std_error), lower.tail = FALSE),
    check.names = FALSE
  )
  out$odds_ratio <- exp(out$estimate)
  out$or_lower <- exp(out$ci_lower)
  out$or_upper <- exp(out$ci_upper)
  rownames(out) <- NULL
  out
}


#' Refuse missing covariate columns, saying why they are missing
#' @noRd
.po_check_columns <- function(x, absent, available) {
  if (length(absent) == 0L) return(invisible(TRUE))
  as_states <- inherits(x, "netobject") && is.data.frame(x$data) &&
    any(absent %in% names(x$data))
  hint <- if (as_states) {
    sprintf(paste("%s was read as a state column when the network was",
                  "built; rebuild it with `metadata_cols = %s`."),
            paste(intersect(absent, names(x$data)), collapse = ", "),
            deparse(intersect(absent, names(x$data))))
  } else if (length(available) == 0L) {
    paste("This input stores no covariates; pass `outcome` as a vector",
          "with one value per sequence.")
  } else {
    sprintf("Available: %s.", paste(available, collapse = ", "))
  }
  stop(errorCondition(
    sprintf("column(s) not found with the sequences: %s. %s",
            paste(absent, collapse = ", "), hint),
    class = "nestimate_bad_input", call = NULL))
}


#' Encode a two-valued outcome (1 = the non-reference value)
#' @noRd
.po_outcome <- function(v, reference = NULL) {
  values <- sort(unique(as.character(v)))
  if (length(values) != 2L) {
    stop(errorCondition(
      sprintf("`outcome` must have exactly two distinct values; it has %d.",
              length(values)),
      class = "nestimate_bad_input", call = NULL))
  }
  reference <- reference %||% values[1L]
  if (!reference %in% values) {
    stop(errorCondition(
      sprintf("`reference` must be one of: %s.",
              paste(values, collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }
  label <- as.character(v)
  list(value = as.integer(label != reference), label = label,
       levels = c(reference, setdiff(values, reference)))
}
